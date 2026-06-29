"""
SYBR API — Jobs router.

Handles job creation, listing, status polling, starting, and cancellation.
"""

import json
import uuid
from datetime import datetime, timezone

from fastapi import APIRouter, Depends, HTTPException, status

from api import database as db
from api.auth import get_current_key
from api.config import SYBR_JOBS_DIR
from api.models import (
    AuthVerifyResponse,
    JobCreateResponse,
    JobListResponse,
    JobResponse,
    JobStatus,
    JobSubmitRequest,
)
from api.services.config_generator import write_job_config
from api.services.file_manager import create_job_directory, cleanup_job
from api.services import job_runner

router = APIRouter(tags=["jobs"])


def _generate_job_id() -> str:
    """Generate a unique, human-readable job ID."""
    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    short_uuid = uuid.uuid4().hex[:8]
    return f"sybr_{ts}_{short_uuid}"


# ── Auth verification ────────────────────────────────────────────────────────

@router.get("/auth/verify", response_model=AuthVerifyResponse)
async def verify_api_key(key: dict = Depends(get_current_key)):
    """Verify that the provided API key is valid."""
    return AuthVerifyResponse(
        valid=True,
        key_name=key["name"],
        permissions=key.get("permissions", "read,write"),
    )


# ── Job CRUD ─────────────────────────────────────────────────────────────────

@router.post(
    "/jobs",
    response_model=JobCreateResponse,
    status_code=status.HTTP_202_ACCEPTED,
)
async def create_job(
    request: JobSubmitRequest,
    key: dict = Depends(get_current_key),
):
    """
    Create a new pipeline job.

    The job is created in 'queued' status. Upload input files, then
    call POST /jobs/{job_id}/start to begin execution.
    """
    job_id = _generate_job_id()

    # Create job directory tree
    job_dir = create_job_directory(job_id)

    # Generate run_sybr_config.yaml
    write_job_config(request, job_dir)

    # Store job in database
    config_json = json.dumps(request.model_dump(), default=str)
    job = db.create_job(
        job_id=job_id,
        api_key_id=key["id"],
        name=request.job_name,
        config_json=config_json,
        input_dir=str(job_dir / "inputs"),
        output_dir=str(job_dir / "outputs"),
        cores=request.cores,
    )

    db.add_job_log(job_id, "INFO", f"Job created by key '{key['name']}'")

    return JobCreateResponse(
        job_id=job_id,
        status=JobStatus.queued,
        created_at=job["created_at"],
        message="Job created. Upload input files, then POST /jobs/{job_id}/start to run.",
    )


@router.get("/jobs", response_model=JobListResponse)
async def list_jobs(
    status_filter: str | None = None,
    limit: int = 50,
    offset: int = 0,
    key: dict = Depends(get_current_key),
):
    """List jobs belonging to the authenticated user."""
    jobs = db.list_jobs(
        api_key_id=key["id"],
        status=status_filter,
        limit=limit,
        offset=offset,
    )
    return JobListResponse(
        jobs=[JobResponse.from_db_row(j) for j in jobs],
        total=len(jobs),
    )


@router.get("/jobs/{job_id}", response_model=JobResponse)
async def get_job(
    job_id: str,
    key: dict = Depends(get_current_key),
):
    """
    Get detailed status of a specific job by ID.

    Any authenticated user can look up a job by its ID.
    This enables the 'Check Job' workflow where users share Job IDs.
    """
    job = db.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    return JobResponse.from_db_row(job)


@router.post("/jobs/{job_id}/start", response_model=JobResponse)
async def start_job(
    job_id: str,
    key: dict = Depends(get_current_key),
):
    """
    Start pipeline execution for a job.

    The job must be in 'queued' or 'uploading' status with input files uploaded.
    """
    job = db.get_job(job_id)
    if job is None or job["api_key_id"] != key["id"]:
        raise HTTPException(status_code=404, detail="Job not found")

    if job["status"] not in ("queued", "uploading"):
        raise HTTPException(
            status_code=400,
            detail=f"Cannot start job in '{job['status']}' status. "
                   f"Job must be in 'queued' or 'uploading' status.",
        )

    # Update status and submit to runner
    db.update_job(job_id, status="validating")
    db.add_job_log(job_id, "INFO", "Job submitted for execution")

    job_runner.submit_job(job_id)

    return JobResponse.from_db_row(db.get_job(job_id))


@router.get("/jobs/{job_id}/logs")
async def get_job_logs(
    job_id: str,
    tail: int = 100,
    key: dict = Depends(get_current_key),
):
    """
    Get the pipeline log output for a job.

    Returns the last `tail` lines from the sybr.log file and DB log entries.
    """
    job = db.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")

    # Read the sybr.log file
    from api.services.file_manager import get_job_log_path

    log_path = get_job_log_path(job_id)
    pipeline_log = ""
    if log_path.exists():
        try:
            with open(log_path) as f:
                lines = f.readlines()
            pipeline_log = "".join(lines[-tail:])
        except Exception:
            pipeline_log = "(could not read log file)"

    # Also get structured DB logs
    db_logs = db.get_job_logs(job_id, limit=tail)

    return {
        "job_id": job_id,
        "pipeline_log": pipeline_log,
        "events": db_logs,
    }


@router.delete("/jobs/{job_id}")
async def delete_job(
    job_id: str,
    key: dict = Depends(get_current_key),
):
    """
    Cancel a running job, or delete a completed/failed job and its files.
    """
    job = db.get_job(job_id)
    if job is None or job["api_key_id"] != key["id"]:
        raise HTTPException(status_code=404, detail="Job not found")

    # If running, cancel first
    if job["status"] in ("running", "validating"):
        job_runner.cancel_job(job_id)

    # Clean up files
    cleanup_job(job_id)

    # Delete from database
    db.delete_job(job_id)

    return {"message": f"Job {job_id} deleted", "job_id": job_id}
