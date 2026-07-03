"""
SYBR API — API key authentication.

Validates the X-API-Key header on every protected request.
"""

import hashlib
import secrets

from fastapi import Depends, HTTPException, Security, status, Query
from fastapi.security import APIKeyHeader

from api import database as db

# Header name used by clients
API_KEY_HEADER = APIKeyHeader(name="X-API-Key", auto_error=False)


def hash_key(raw_key: str) -> str:
    """SHA-256 hash of a raw API key."""
    return hashlib.sha256(raw_key.encode()).hexdigest()


def generate_api_key() -> str:
    """Generate a cryptographically secure API key (48 hex chars)."""
    return secrets.token_hex(24)


async def get_current_key(
    api_key_header: str | None = Security(API_KEY_HEADER),
    api_key_query: str | None = Query(None, alias="api_key"),
) -> dict:
    """
    FastAPI dependency that validates the API key.
    Checks the X-API-Key header first, then the api_key query parameter.

    Returns the api_keys row dict if valid; raises 401/403 otherwise.
    """
    api_key = api_key_header or api_key_query

    if api_key is None:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Missing API key (header or query parameter)",
        )

    key_row = db.get_api_key_by_hash(hash_key(api_key))

    if key_row is None:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Invalid or revoked API key",
        )

    return key_row
