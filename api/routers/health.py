"""
SYBR API — Health check router.
"""

from fastapi import APIRouter

from api import database as db
from api.config import SYBR_JOBS_DIR, SYBR_PIPELINE_DIR
from api.models import HealthResponse
from api.services import job_runner

router = APIRouter(tags=["health"])


@router.get("/health", response_model=HealthResponse)
async def health_check():
    """Health check — no authentication required."""
    return HealthResponse(
        status="ok",
        version="0.1.0",
        pipeline_dir=str(SYBR_PIPELINE_DIR),
        jobs_dir=str(SYBR_JOBS_DIR),
        active_jobs=job_runner.get_active_job_count(),
    )
