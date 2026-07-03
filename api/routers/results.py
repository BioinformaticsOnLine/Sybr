"""
SYBR API — Results router.

Download pipeline outputs for completed jobs.
"""

from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import FileResponse

from api import database as db
from api.auth import get_current_key
from api.models import ResultEntry, ResultListResponse
from api.services.file_manager import (
    FileManagerError,
    get_output_file_path,
    list_job_outputs,
)

router = APIRouter(tags=["results"])


@router.get("/jobs/{job_id}/results", response_model=ResultListResponse)
async def list_results(
    job_id: str,
    path: str = "",
    key: dict = Depends(get_current_key),
):
    """
    List output files and directories for a completed job.

    Use the `path` query parameter to browse subdirectories.
    """
    job = db.get_job(job_id)
    if job is None or job["api_key_id"] != key["id"]:
        raise HTTPException(status_code=404, detail="Job not found")

    try:
        entries = list_job_outputs(job_id, subpath=path)
    except FileManagerError as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    return ResultListResponse(
        job_id=job_id,
        entries=[ResultEntry(**e) for e in entries],
    )


@router.get("/jobs/{job_id}/results/{file_path:path}")
async def download_result(
    job_id: str,
    file_path: str,
    key: dict = Depends(get_current_key),
):
    """
    Download a specific output file.

    The `file_path` is relative to the job's outputs directory.
    Example: `eba_analysis/EBRs/EBRs_stats.html`
    """
    job = db.get_job(job_id)
    if job is None or job["api_key_id"] != key["id"]:
        raise HTTPException(status_code=404, detail="Job not found")

    try:
        full_path = get_output_file_path(job_id, file_path)
    except FileManagerError as exc:
        raise HTTPException(status_code=404, detail=str(exc))

    return FileResponse(
        path=str(full_path),
        filename=full_path.name,
        media_type="application/octet-stream",
    )


@router.get("/jobs/{job_id}/results_archive")
async def download_results_archive(
    job_id: str,
    key: dict = Depends(get_current_key),
):
    """
    Download the entire outputs directory as a ZIP archive.
    The archive is cached in the job directory after creation.
    """
    import os
    import shutil
    from api.config import SYBR_JOBS_DIR

    job = db.get_job(job_id)
    if job is None or job["api_key_id"] != key["id"]:
        raise HTTPException(status_code=404, detail="Job not found")

    job_dir = SYBR_JOBS_DIR / job_id
    output_dir = job_dir / "outputs"
    archive_path = job_dir / f"{job_id}_results.zip"

    if not output_dir.exists():
        raise HTTPException(status_code=404, detail="Outputs directory not found")

    # Only recreate the zip if it doesn't exist (acting as a cache)
    if not archive_path.exists():
        # shutil.make_archive adds the .zip extension automatically
        base_name = str(job_dir / f"{job_id}_results")
        shutil.make_archive(base_name, "zip", root_dir=str(output_dir))

    return FileResponse(
        path=str(archive_path),
        filename=archive_path.name,
        media_type="application/zip",
    )
