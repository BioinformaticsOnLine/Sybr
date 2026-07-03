"""
SYBR API — FastAPI application assembly.

This is the main FastAPI app that wires together all routers,
middleware, and lifecycle hooks.
"""

import logging
from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from api import database as db
from api.config import API_PREFIX, LOG_LEVEL, SYBR_JOBS_DIR, SYBR_PIPELINE_DIR
from api.routers import files, health, jobs, results
from api.services import job_runner

# ── Logging setup ────────────────────────────────────────────────────────────

logging.basicConfig(
    level=getattr(logging, LOG_LEVEL.upper(), logging.INFO),
    format="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("sybr.api")


# ── Application lifespan ─────────────────────────────────────────────────────

@asynccontextmanager
async def lifespan(app: FastAPI):
    """Startup and shutdown hooks."""
    # ── Startup ──────────────────────────────────────────────────────────
    logger.info("=" * 60)
    logger.info("  SYBR API starting up")
    logger.info(f"  Pipeline dir : {SYBR_PIPELINE_DIR}")
    logger.info(f"  Jobs dir     : {SYBR_JOBS_DIR}")
    logger.info("=" * 60)

    # Ensure directories exist
    SYBR_JOBS_DIR.mkdir(parents=True, exist_ok=True)

    # Initialize database
    db.init_db()
    logger.info("Database initialized")

    # Initialize job runner thread pool
    job_runner.init_runner()
    logger.info("Job runner initialized")

    yield

    # ── Shutdown ─────────────────────────────────────────────────────────
    logger.info("SYBR API shutting down...")
    job_runner.shutdown_runner()
    logger.info("Job runner shut down")


# ── FastAPI app ──────────────────────────────────────────────────────────────

app = FastAPI(
    title="SYBR API",
    description=(
        "REST API for the SYBR Snakemake pipeline — "
        "Synteny Block discovery, evolutionary Breakpoint identification, "
        "and ancestral genome Reconstruction."
    ),
    version="0.1.0",
    lifespan=lifespan,
    docs_url="/docs",
    redoc_url="/redoc",
)

# ── CORS middleware ──────────────────────────────────────────────────────────

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Tighten in production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ── Mount routers ────────────────────────────────────────────────────────────

app.include_router(health.router, prefix=API_PREFIX)
app.include_router(jobs.router, prefix=API_PREFIX)
app.include_router(files.router, prefix=API_PREFIX)
app.include_router(results.router, prefix=API_PREFIX)


# ── Root redirect ────────────────────────────────────────────────────────────

@app.get("/", include_in_schema=False)
async def root():
    """Redirect root to API docs."""
    from fastapi.responses import RedirectResponse
    return RedirectResponse(url="/docs")
