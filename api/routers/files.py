"""
SYBR API — File upload router.
"""

from fastapi import APIRouter, Depends, HTTPException, UploadFile, File, Form

from api import database as db
from api.auth import get_current_key
from api.models import FileCategory, FileListResponse, FileUploadResponse, FileInfo
from api.services.file_manager import (
    FileManagerError,
    list_job_files,
    save_upload,
)

router = APIRouter(tags=["files"])


@router.post(
    "/jobs/{job_id}/upload",
    response_model=FileUploadResponse,
)
async def upload_file(
    job_id: str,
    category: FileCategory = Form(..., description="File category"),
    file: UploadFile = File(..., description="File to upload"),
    key: dict = Depends(get_current_key),
):
    """
    Upload an input file to a job.

    The `category` determines where the file is placed in the input directory:
    - **fasta**: `inputs/fasta/` (*.fa, *.fasta, *.fna)
    - **satsuma_alignments**: `inputs/synteny_processing/Satsuma_alignments/` (*.txt)
    - **sequence_lengths**: `inputs/synteny_processing/` (all_sequence_lengths.txt)
    - **scaffolds**: `inputs/synteny_processing/` (Scaffolds.txt)
    - **lastz_alignments**: `inputs/Ancestor_seq_recunstruction/LastZ_alignments/` (*.axt)
    - **species_info**: `inputs/Ancestor_seq_recunstruction/` (species_info.txt)
    - **tree**: `inputs/Ancestor_seq_recunstruction/` (tree.txt)
    - **classification**: `inputs/eba_analysis/` (classification.eba)
    - **annotation**: `inputs/enrichment_analysis/` (protein_annotation.tsv)
    - **kegg**: `inputs/enrichment_analysis/` (3kegg_annotationTOgenes.txt)
    - **hgt**: `inputs/HGTs/` (hgt.txt)
    """
    # Verify job exists and belongs to this key
    job = db.get_job(job_id)
    if job is None or job["api_key_id"] != key["id"]:
        raise HTTPException(status_code=404, detail="Job not found")

    # Only allow uploads for jobs not yet running
    if job["status"] not in ("queued", "uploading"):
        raise HTTPException(
            status_code=400,
            detail=f"Cannot upload files to a job in '{job['status']}' status",
        )

    try:
        result = await save_upload(
            job_id=job_id,
            category=category.value,
            filename=file.filename or "unnamed",
            file_obj=file,
        )
    except FileManagerError as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    # Update job status to 'uploading' if it was 'queued'
    if job["status"] == "queued":
        db.update_job(job_id, status="uploading")

    db.add_job_log(
        job_id,
        "INFO",
        f"Uploaded {result['filename']} ({result['size_bytes']} bytes) "
        f"to {result['destination']}",
    )

    return FileUploadResponse(
        job_id=job_id,
        filename=result["filename"],
        category=category.value,
        size_bytes=result["size_bytes"],
        destination=result["destination"],
        message=f"File uploaded successfully to {result['destination']}",
    )


@router.get("/jobs/{job_id}/files", response_model=FileListResponse)
async def list_files(
    job_id: str,
    key: dict = Depends(get_current_key),
):
    """List all uploaded input files for a job."""
    job = db.get_job(job_id)
    if job is None or job["api_key_id"] != key["id"]:
        raise HTTPException(status_code=404, detail="Job not found")

    files = list_job_files(job_id)

    return FileListResponse(
        job_id=job_id,
        files=[FileInfo(**f) for f in files],
        total=len(files),
    )
