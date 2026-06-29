"""
SYBR API — Job runner service.

Manages the lifecycle of pipeline jobs: starting sybr.sh as a subprocess,
monitoring progress, and handling cancellation.
"""

import json
import logging
import os
import signal
import subprocess
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from api import database as db
from api.config import (
    MAX_CONCURRENT_JOBS,
    MAX_CORES_PER_JOB,
    PIPELINE_PATHS_YAML,
    SYBR_JOBS_DIR,
    SYBR_PIPELINE_DIR,
    SYBR_RUNNER,
)
from api.services.file_manager import get_job_log_path

logger = logging.getLogger("sybr.runner")

# Thread pool for running pipeline jobs
_executor: ThreadPoolExecutor | None = None
_lock = threading.Lock()


def init_runner():
    """Initialize the thread pool executor."""
    global _executor
    _executor = ThreadPoolExecutor(
        max_workers=MAX_CONCURRENT_JOBS,
        thread_name_prefix="sybr-job",
    )
    logger.info(f"Job runner initialized (max_workers={MAX_CONCURRENT_JOBS})")


def shutdown_runner():
    """Shut down the thread pool (waits for running jobs)."""
    global _executor
    if _executor:
        logger.info("Shutting down job runner...")
        _executor.shutdown(wait=False, cancel_futures=True)
        _executor = None


def submit_job(job_id: str):
    """
    Submit a job for execution in the background thread pool.

    The job must already exist in the database with status 'queued' or 'validating'.
    """
    if _executor is None:
        raise RuntimeError("Job runner not initialized. Call init_runner() first.")

    _executor.submit(_run_pipeline, job_id)
    logger.info(f"Job {job_id} submitted to thread pool")


def _run_pipeline(job_id: str):
    """
    Execute the SYBR pipeline for a single job.

    This runs in a worker thread. It:
    1. Updates job status to 'running'
    2. Invokes sybr.sh via subprocess
    3. Monitors stdout for progress
    4. Updates status to 'completed' or 'failed'
    """
    job = db.get_job(job_id)
    if job is None:
        logger.error(f"Job {job_id} not found in database")
        return

    job_dir = Path(job["input_dir"]).parent  # job_dir = SYBR_JOBS_DIR/job_id
    config_path = job_dir / "run_sybr_config.yaml"
    log_path = get_job_log_path(job_id)

    # Validate prerequisites
    if not config_path.exists():
        db.update_job(
            job_id,
            status="failed",
            error_message=f"Config file not found: {config_path}",
        )
        db.add_job_log(job_id, "ERROR", f"Config file not found: {config_path}")
        return

    if not SYBR_RUNNER.exists():
        db.update_job(
            job_id,
            status="failed",
            error_message=f"Pipeline runner not found: {SYBR_RUNNER}",
        )
        db.add_job_log(job_id, "ERROR", f"Pipeline runner not found: {SYBR_RUNNER}")
        return

    # Build the sybr.sh command
    # Parse window sizes and step size from the config
    try:
        config_data = json.loads(job["config_json"])
    except (json.JSONDecodeError, TypeError):
        config_data = {}

    requested_cores = int(job.get("cores", 4) or 4)
    cores = max(1, min(requested_cores, MAX_CORES_PER_JOB))

    cmd = [
        str(SYBR_RUNNER),
        "-c", str(config_path),
        "-P", str(PIPELINE_PATHS_YAML),
        "-j", str(cores),
        "-s",  # skip validation (API does its own)
    ]

    # Add custom window sizes if provided
    window_sizes = config_data.get("window_sizes")
    if window_sizes and isinstance(window_sizes, list):
        cmd.extend(["-w", ",".join(str(w) for w in window_sizes)])

    step_size = config_data.get("step_size")
    if step_size:
        cmd.extend(["-p", str(step_size)])

    logger.info(f"Starting job {job_id}: {' '.join(cmd)}")
    db.update_job(job_id, status="running", started_at=db.now_utc())
    db.add_job_log(job_id, "INFO", f"Pipeline started: {' '.join(cmd)}")

    try:
        with open(log_path, "w") as log_file:
            process = subprocess.Popen(
                cmd,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                cwd=str(SYBR_PIPELINE_DIR),
                env={
                    **os.environ,
                    "PAGER": "cat",  # Prevent interactive pagers
                },
                start_new_session=True,  # Own process group for clean kill
            )

        # Store PID for cancellation
        db.update_job(job_id, pid=process.pid)
        db.add_job_log(job_id, "INFO", f"Pipeline PID: {process.pid}")

        # Wait for completion
        return_code = process.wait()

        if return_code == 0:
            db.update_job(
                job_id,
                status="completed",
                completed_at=db.now_utc(),
                progress_pct=100.0,
                pid=None,
            )
            db.add_job_log(job_id, "INFO", "Pipeline completed successfully")
            logger.info(f"Job {job_id} completed successfully")
        else:
            # Read last 20 lines of log for error context
            error_tail = _tail_file(log_path, 20)
            db.update_job(
                job_id,
                status="failed",
                completed_at=db.now_utc(),
                error_message=f"Pipeline exited with code {return_code}",
                pid=None,
            )
            db.add_job_log(
                job_id,
                "ERROR",
                f"Pipeline failed (exit code {return_code}):\n{error_tail}",
            )
            logger.error(f"Job {job_id} failed with exit code {return_code}")

    except Exception as exc:
        db.update_job(
            job_id,
            status="failed",
            completed_at=db.now_utc(),
            error_message=str(exc),
            pid=None,
        )
        db.add_job_log(job_id, "ERROR", f"Runner exception: {exc}")
        logger.exception(f"Job {job_id} runner exception")


def cancel_job(job_id: str) -> bool:
    """
    Cancel a running job by killing its process group.

    Returns True if the job was successfully cancelled.
    """
    job = db.get_job(job_id)
    if job is None:
        return False

    if job["status"] not in ("running", "queued", "validating"):
        return False

    pid = job.get("pid")
    if pid:
        try:
            # Kill the entire process group
            os.killpg(os.getpgid(pid), signal.SIGTERM)
            logger.info(f"Sent SIGTERM to process group of PID {pid}")
        except (ProcessLookupError, PermissionError):
            logger.warning(f"Could not kill PID {pid} — process may have already exited")

    db.update_job(
        job_id,
        status="cancelled",
        completed_at=db.now_utc(),
        pid=None,
    )
    db.add_job_log(job_id, "INFO", "Job cancelled by user")
    return True


def get_active_job_count() -> int:
    """Return the number of currently running jobs."""
    jobs = db.list_jobs(status="running")
    return len(jobs)


def _tail_file(path: Path, n: int = 20) -> str:
    """Read the last n lines of a file."""
    try:
        with open(path) as f:
            lines = f.readlines()
        return "".join(lines[-n:])
    except Exception:
        return "(could not read log file)"
