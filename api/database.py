"""
SYBR API — SQLite database layer.

Handles connection management, schema creation, and common queries.
"""

import sqlite3
import threading
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

from api.config import SYBR_DB_PATH

# Thread-local storage for connections (SQLite is not thread-safe by default)
_local = threading.local()

SCHEMA_SQL = """
-- API keys table
CREATE TABLE IF NOT EXISTS api_keys (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    key_hash    TEXT    NOT NULL UNIQUE,
    name        TEXT    NOT NULL,
    created_at  TEXT    DEFAULT (strftime('%Y-%m-%dT%H:%M:%SZ', 'now')),
    last_used   TEXT,
    is_active   INTEGER DEFAULT 1,
    permissions TEXT    DEFAULT 'read,write'
);

-- Jobs table
CREATE TABLE IF NOT EXISTS jobs (
    id              TEXT    PRIMARY KEY,
    api_key_id      INTEGER REFERENCES api_keys(id),
    name            TEXT,
    email           TEXT,
    status          TEXT    DEFAULT 'queued',
    config_json     TEXT    NOT NULL,
    input_dir       TEXT    NOT NULL,
    output_dir      TEXT    NOT NULL,
    pid             INTEGER,
    cores           INTEGER DEFAULT 4,
    progress_pct    REAL    DEFAULT 0.0,
    current_stage   TEXT,
    error_message   TEXT,
    created_at      TEXT    DEFAULT (strftime('%Y-%m-%dT%H:%M:%SZ', 'now')),
    started_at      TEXT,
    completed_at    TEXT,
    updated_at      TEXT    DEFAULT (strftime('%Y-%m-%dT%H:%M:%SZ', 'now'))
);

-- Job logs table
CREATE TABLE IF NOT EXISTS job_logs (
    id         INTEGER PRIMARY KEY AUTOINCREMENT,
    job_id     TEXT    REFERENCES jobs(id) ON DELETE CASCADE,
    timestamp  TEXT    DEFAULT (strftime('%Y-%m-%dT%H:%M:%SZ', 'now')),
    level      TEXT,
    message    TEXT
);

-- Index for fast job lookup by status
CREATE INDEX IF NOT EXISTS idx_jobs_status ON jobs(status);
CREATE INDEX IF NOT EXISTS idx_jobs_api_key ON jobs(api_key_id);
CREATE INDEX IF NOT EXISTS idx_job_logs_job_id ON job_logs(job_id);
"""


def _get_connection() -> sqlite3.Connection:
    """Get or create a thread-local SQLite connection."""
    if not hasattr(_local, "conn") or _local.conn is None:
        SYBR_DB_PATH.parent.mkdir(parents=True, exist_ok=True)
        _local.conn = sqlite3.connect(
            str(SYBR_DB_PATH),
            check_same_thread=False,
            timeout=30,
        )
        _local.conn.row_factory = sqlite3.Row
        _local.conn.execute("PRAGMA journal_mode=WAL")
        _local.conn.execute("PRAGMA foreign_keys=ON")
    return _local.conn


@contextmanager
def get_db():
    """Context manager that yields a database connection and auto-commits."""
    conn = _get_connection()
    try:
        yield conn
        conn.commit()
    except Exception:
        conn.rollback()
        raise


def init_db():
    """Create tables if they don't exist, and run incremental migrations."""
    with get_db() as conn:
        conn.executescript(SCHEMA_SQL)
        # Migration: add email column to existing databases that predate this field
        existing_cols = {
            row[1]
            for row in conn.execute("PRAGMA table_info(jobs)").fetchall()
        }
        if "email" not in existing_cols:
            conn.execute("ALTER TABLE jobs ADD COLUMN email TEXT")


def now_utc() -> str:
    """Return current UTC time as ISO-8601 string."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ── Job queries ──────────────────────────────────────────────────────────────

def create_job(
    job_id: str,
    api_key_id: int,
    name: str,
    config_json: str,
    input_dir: str,
    output_dir: str,
    cores: int,
    email: str | None = None,
) -> dict:
    """Insert a new job row and return it as a dict."""
    with get_db() as conn:
        conn.execute(
            """INSERT INTO jobs
               (id, api_key_id, name, email, config_json, input_dir, output_dir, cores)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
            (job_id, api_key_id, name, email, config_json, input_dir, output_dir, cores),
        )
    return get_job(job_id)


def get_job(job_id: str) -> dict | None:
    """Fetch a single job by ID."""
    with get_db() as conn:
        row = conn.execute("SELECT * FROM jobs WHERE id = ?", (job_id,)).fetchone()
    return dict(row) if row else None


def list_jobs(
    api_key_id: int | None = None,
    status: str | None = None,
    limit: int = 50,
    offset: int = 0,
) -> list[dict]:
    """List jobs with optional filters."""
    query = "SELECT * FROM jobs WHERE 1=1"
    params: list = []
    if api_key_id is not None:
        query += " AND api_key_id = ?"
        params.append(api_key_id)
    if status:
        query += " AND status = ?"
        params.append(status)
    query += " ORDER BY created_at DESC LIMIT ? OFFSET ?"
    params.extend([limit, offset])
    with get_db() as conn:
        rows = conn.execute(query, params).fetchall()
    return [dict(r) for r in rows]


def update_job(job_id: str, **fields) -> dict | None:
    """Update specific fields on a job."""
    if not fields:
        return get_job(job_id)
    fields["updated_at"] = now_utc()
    set_clause = ", ".join(f"{k} = ?" for k in fields)
    values = list(fields.values()) + [job_id]
    with get_db() as conn:
        conn.execute(f"UPDATE jobs SET {set_clause} WHERE id = ?", values)
    return get_job(job_id)


def delete_job(job_id: str):
    """Delete a job row (cascades to job_logs)."""
    with get_db() as conn:
        conn.execute("DELETE FROM jobs WHERE id = ?", (job_id,))


# ── API key queries ──────────────────────────────────────────────────────────

def create_api_key(key_hash: str, name: str) -> dict:
    """Insert a new API key."""
    with get_db() as conn:
        cursor = conn.execute(
            "INSERT INTO api_keys (key_hash, name) VALUES (?, ?)",
            (key_hash, name),
        )
        row = conn.execute(
            "SELECT * FROM api_keys WHERE id = ?", (cursor.lastrowid,)
        ).fetchone()
    return dict(row)


def get_api_key_by_hash(key_hash: str) -> dict | None:
    """Look up an active API key by its hash."""
    with get_db() as conn:
        row = conn.execute(
            "SELECT * FROM api_keys WHERE key_hash = ? AND is_active = 1",
            (key_hash,),
        ).fetchone()
    if row:
        # Update last_used timestamp
        conn.execute(
            "UPDATE api_keys SET last_used = ? WHERE id = ?",
            (now_utc(), row["id"]),
        )
        conn.commit()
    return dict(row) if row else None


def list_api_keys() -> list[dict]:
    """List all API keys (without hashes)."""
    with get_db() as conn:
        rows = conn.execute(
            "SELECT id, name, created_at, last_used, is_active, permissions FROM api_keys"
        ).fetchall()
    return [dict(r) for r in rows]


def revoke_api_key(key_id: int):
    """Deactivate an API key."""
    with get_db() as conn:
        conn.execute("UPDATE api_keys SET is_active = 0 WHERE id = ?", (key_id,))


# ── Log queries ──────────────────────────────────────────────────────────────

def add_job_log(job_id: str, level: str, message: str):
    """Append a log entry for a job."""
    with get_db() as conn:
        conn.execute(
            "INSERT INTO job_logs (job_id, level, message) VALUES (?, ?, ?)",
            (job_id, level, message),
        )


def get_job_logs(job_id: str, limit: int = 200, offset: int = 0) -> list[dict]:
    """Get log entries for a job, newest first."""
    with get_db() as conn:
        rows = conn.execute(
            "SELECT * FROM job_logs WHERE job_id = ? ORDER BY id DESC LIMIT ? OFFSET ?",
            (job_id, limit, offset),
        ).fetchall()
    return [dict(r) for r in rows]
