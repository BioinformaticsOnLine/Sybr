"""
SYBR API — Pydantic request / response models.
"""

from __future__ import annotations

from datetime import datetime
from enum import Enum
from typing import Any, Optional

from pydantic import BaseModel, Field

from api.config import DEFAULT_CORES, MAX_CORES_PER_JOB


# ── Enums ────────────────────────────────────────────────────────────────────

class JobStatus(str, Enum):
    queued = "queued"
    uploading = "uploading"
    validating = "validating"
    running = "running"
    completed = "completed"
    failed = "failed"
    cancelled = "cancelled"


class FileCategory(str, Enum):
    fasta = "fasta"
    satsuma_alignments = "satsuma_alignments"
    sequence_lengths = "sequence_lengths"
    scaffolds = "scaffolds"
    lastz_alignments = "lastz_alignments"
    species_info = "species_info"
    tree = "tree"
    classification = "classification"
    annotation = "annotation"
    kegg = "kegg"
    hgt = "hgt"


# ── Request models ───────────────────────────────────────────────────────────

class RunStages(BaseModel):
    synteny_processing: bool = True
    eba_analysis: bool = True
    enrichment_analysis: bool = True
    chainNet_generation: bool = True
    Ancestor_seq_recunstruction: bool = True
    hgt_overlap_analysis: bool = False


class EBAParams(BaseModel):
    n: int = Field(default=5, ge=1, description="Number of EBA iterations")
    r: str = Field(description="Reference species name (e.g. Adineta_vaga)")
    p: int = Field(default=300, ge=1, description="Resolution parameter")


class EnrichParams(BaseModel):
    r: str = Field(default="ko", description="KEGG organism code or 'ko'")


class JobSubmitRequest(BaseModel):
    """Body for POST /jobs."""
    job_name: str = Field(
        default="sybr_job",
        min_length=1,
        max_length=128,
        description="Human-readable job name",
    )
    email: Optional[str] = Field(
        default=None,
        description="Submitter email address for notifications and admin tracking",
    )
    run_stages: RunStages = Field(default_factory=RunStages)
    reference_name: str = Field(
        description="Full Genus_species reference name matching FASTA filename",
        examples=["Adineta_vaga"],
    )
    reference_species: str = Field(
        description="Short species name used internally",
        examples=["vaga"],
    )
    eba: EBAParams
    getenrich: EnrichParams = Field(default_factory=EnrichParams)
    window_sizes: list[int] = Field(
        default=[100000, 300000, 500000],
        description="Window sizes in bp for synteny assignment",
    )
    step_size: int = Field(
        default=30000,
        ge=1,
        description="Step size in bp for synteny assignment",
    )
    cores: int = Field(
        default=DEFAULT_CORES,
        ge=1,
        le=MAX_CORES_PER_JOB,
        description=f"Number of CPU cores for this job (max {MAX_CORES_PER_JOB})",
    )


# ── Response models ──────────────────────────────────────────────────────────

class JobProgress(BaseModel):
    current_stage: Optional[str] = None
    percent: float = 0.0


class JobResponse(BaseModel):
    job_id: str
    job_name: Optional[str] = None
    email: Optional[str] = None
    status: JobStatus
    progress: JobProgress = Field(default_factory=JobProgress)
    config: Optional[dict[str, Any]] = None
    cores: int = 4
    created_at: Optional[str] = None
    started_at: Optional[str] = None
    completed_at: Optional[str] = None
    error: Optional[str] = None

    @classmethod
    def from_db_row(cls, row: dict) -> "JobResponse":
        """Construct from a database row dict."""
        import json

        config = None
        if row.get("config_json"):
            try:
                config = json.loads(row["config_json"])
            except (json.JSONDecodeError, TypeError):
                config = None

        return cls(
            job_id=row["id"],
            job_name=row.get("name"),
            email=row.get("email"),
            status=JobStatus(row.get("status", "queued")),
            progress=JobProgress(
                current_stage=row.get("current_stage"),
                percent=row.get("progress_pct", 0.0),
            ),
            config=config,
            cores=row.get("cores", 4),
            created_at=row.get("created_at"),
            started_at=row.get("started_at"),
            completed_at=row.get("completed_at"),
            error=row.get("error_message"),
        )


class JobCreateResponse(BaseModel):
    job_id: str
    status: JobStatus
    created_at: str
    message: str


class JobListResponse(BaseModel):
    jobs: list[JobResponse]
    total: int


class FileInfo(BaseModel):
    name: str
    category: str
    size_bytes: int
    uploaded_at: Optional[str] = None


class FileListResponse(BaseModel):
    job_id: str
    files: list[FileInfo]
    total: int


class FileUploadResponse(BaseModel):
    job_id: str
    filename: str
    category: str
    size_bytes: int
    destination: str
    message: str


class ResultEntry(BaseModel):
    path: str
    name: str
    is_dir: bool
    size_bytes: Optional[int] = None


class ResultListResponse(BaseModel):
    job_id: str
    entries: list[ResultEntry]


class HealthResponse(BaseModel):
    status: str = "ok"
    version: str
    pipeline_dir: str
    jobs_dir: str
    active_jobs: int


class AuthVerifyResponse(BaseModel):
    valid: bool
    key_name: str
    permissions: str


class ErrorResponse(BaseModel):
    detail: str
