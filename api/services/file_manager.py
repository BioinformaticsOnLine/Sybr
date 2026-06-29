"""
SYBR API — File manager service.

Handles creating job directory trees, placing uploaded files into the
correct subdirectory, listing files, and cleanup.
"""

import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

from api.config import (
    ALLOWED_EXTENSIONS,
    CATEGORY_PATHS,
    MAX_UPLOAD_SIZE_BYTES,
    SYBR_JOBS_DIR,
)


class FileManagerError(Exception):
    """Raised when file operations fail."""


def create_job_directory(job_id: str) -> Path:
    """
    Create the full input directory tree for a new job.

    Returns the job root directory.
    """
    job_dir = SYBR_JOBS_DIR / job_id
    input_dir = job_dir / "inputs"

    # Create the standard input subdirectories
    subdirs = [
        "fasta",
        "synteny_processing/Satsuma_alignments",
        "Ancestor_seq_recunstruction/LastZ_alignments",
        "eba_analysis",
        "enrichment_analysis",
        "HGTs",
    ]

    for subdir in subdirs:
        (input_dir / subdir).mkdir(parents=True, exist_ok=True)

    # Create outputs directory
    (job_dir / "outputs").mkdir(parents=True, exist_ok=True)

    return job_dir


def get_upload_destination(job_id: str, category: str, filename: str) -> Path:
    """
    Determine the full path where an uploaded file should be stored.

    Validates category and file extension.
    """
    # Validate category
    if category not in CATEGORY_PATHS:
        raise FileManagerError(
            f"Unknown file category: '{category}'. "
            f"Valid categories: {list(CATEGORY_PATHS.keys())}"
        )

    # Sanitize filename (prevent path traversal)
    safe_name = Path(filename).name
    if not safe_name or safe_name.startswith("."):
        raise FileManagerError(f"Invalid filename: '{filename}'")

    # Validate extension
    ext = Path(safe_name).suffix.lower()
    allowed = ALLOWED_EXTENSIONS.get(category, set())
    if allowed and ext not in allowed:
        raise FileManagerError(
            f"Extension '{ext}' not allowed for category '{category}'. "
            f"Allowed: {allowed}"
        )

    # Build destination path
    job_dir = SYBR_JOBS_DIR / job_id
    subdir = CATEGORY_PATHS[category]
    dest = job_dir / "inputs" / subdir / safe_name

    return dest


async def save_upload(job_id: str, category: str, filename: str, file_obj) -> dict:
    """
    Save an uploaded file to the correct location.

    Args:
        job_id: The job identifier.
        category: File category (fasta, satsuma_alignments, etc.).
        filename: Original filename from the upload.
        file_obj: FastAPI UploadFile or file-like object.

    Returns:
        Dict with filename, category, size_bytes, destination.
    """
    dest = get_upload_destination(job_id, category, filename)
    dest.parent.mkdir(parents=True, exist_ok=True)

    # Stream the file to disk in chunks to handle large files
    total_bytes = 0
    chunk_size = 1024 * 1024  # 1 MB chunks

    with open(dest, "wb") as f:
        while True:
            chunk = await file_obj.read(chunk_size)
            if not chunk:
                break
            total_bytes += len(chunk)
            if total_bytes > MAX_UPLOAD_SIZE_BYTES:
                # Clean up partial file
                f.close()
                dest.unlink(missing_ok=True)
                raise FileManagerError(
                    f"File exceeds maximum upload size "
                    f"({MAX_UPLOAD_SIZE_BYTES // (1024*1024)} MB)"
                )
            f.write(chunk)

    return {
        "filename": dest.name,
        "category": category,
        "size_bytes": total_bytes,
        "destination": str(dest.relative_to(SYBR_JOBS_DIR / job_id)),
    }


def list_job_files(job_id: str) -> list[dict]:
    """List all uploaded files in a job's input directory."""
    job_dir = SYBR_JOBS_DIR / job_id
    input_dir = job_dir / "inputs"

    if not input_dir.exists():
        return []

    files = []
    # Build a reverse map: subdir → category
    reverse_map = {}
    for cat, subdir in CATEGORY_PATHS.items():
        reverse_map[subdir] = cat

    for path in sorted(input_dir.rglob("*")):
        if path.is_file():
            rel = path.relative_to(input_dir)
            # Determine category from parent directory
            category = "unknown"
            for subdir, cat in reverse_map.items():
                if str(rel).startswith(subdir):
                    category = cat
                    break

            stat = path.stat()
            files.append({
                "name": path.name,
                "category": category,
                "size_bytes": stat.st_size,
                "uploaded_at": datetime.fromtimestamp(
                    stat.st_mtime, tz=timezone.utc
                ).strftime("%Y-%m-%dT%H:%M:%SZ"),
            })

    return files


def list_job_outputs(job_id: str, subpath: str = "") -> list[dict]:
    """List files/dirs in a job's output directory."""
    job_dir = SYBR_JOBS_DIR / job_id
    output_dir = job_dir / "outputs"

    if subpath:
        # Prevent path traversal
        target = (output_dir / subpath).resolve()
        if not str(target).startswith(str(output_dir.resolve())):
            raise FileManagerError("Invalid path: directory traversal detected")
        output_dir = target

    if not output_dir.exists():
        return []

    entries = []
    for item in sorted(output_dir.iterdir()):
        entry = {
            "path": str(item.relative_to(SYBR_JOBS_DIR / job_id / "outputs")),
            "name": item.name,
            "is_dir": item.is_dir(),
        }
        if item.is_file():
            entry["size_bytes"] = item.stat().st_size
        entries.append(entry)

    return entries


def get_output_file_path(job_id: str, file_path: str) -> Path:
    """
    Resolve and validate a file path within a job's output directory.

    Raises FileManagerError if the file doesn't exist or path traversal is detected.
    """
    job_dir = SYBR_JOBS_DIR / job_id
    output_dir = job_dir / "outputs"
    target = (output_dir / file_path).resolve()

    # Security: ensure we're still within the outputs directory
    if not str(target).startswith(str(output_dir.resolve())):
        raise FileManagerError("Invalid path: directory traversal detected")

    if not target.exists():
        raise FileManagerError(f"File not found: {file_path}")

    if not target.is_file():
        raise FileManagerError(f"Not a file: {file_path}")

    return target


def cleanup_job(job_id: str):
    """Remove all files for a job."""
    job_dir = SYBR_JOBS_DIR / job_id
    if job_dir.exists():
        shutil.rmtree(job_dir, ignore_errors=True)


def get_job_log_path(job_id: str) -> Path:
    """Return the path to a job's pipeline log file."""
    return SYBR_JOBS_DIR / job_id / "sybr.log"
