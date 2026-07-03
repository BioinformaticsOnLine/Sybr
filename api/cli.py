"""
SYBR API — CLI tool for managing the API server.

Usage:
    python -m api.cli create-key --name "my-frontend"
    python -m api.cli list-keys
    python -m api.cli revoke-key --id 1
    python -m api.cli serve --port 8000
"""

import argparse
import sys

from api.auth import generate_api_key, hash_key
from api import database as db
from api.config import API_HOST, API_PORT


def cmd_create_key(args):
    """Create a new API key."""
    db.init_db()

    raw_key = generate_api_key()
    key_hash = hash_key(raw_key)

    row = db.create_api_key(key_hash=key_hash, name=args.name)

    print()
    print("=" * 60)
    print("  NEW API KEY CREATED")
    print("=" * 60)
    print(f"  Name    : {row['name']}")
    print(f"  ID      : {row['id']}")
    print(f"  Key     : {raw_key}")
    print()
    print("  ⚠  Save this key now — it cannot be recovered!")
    print("=" * 60)
    print()


def cmd_list_keys(args):
    """List all API keys."""
    db.init_db()

    keys = db.list_api_keys()
    if not keys:
        print("No API keys found. Create one with: python -m api.cli create-key --name 'my-key'")
        return

    print()
    print(f"{'ID':<5} {'Name':<25} {'Active':<8} {'Created':<22} {'Last Used':<22}")
    print("-" * 82)
    for k in keys:
        active = "✓" if k["is_active"] else "✗"
        last_used = k.get("last_used") or "never"
        print(f"{k['id']:<5} {k['name']:<25} {active:<8} {k['created_at']:<22} {last_used:<22}")
    print()


def cmd_revoke_key(args):
    """Revoke an API key."""
    db.init_db()

    db.revoke_api_key(args.id)
    print(f"API key {args.id} revoked.")


def cmd_serve(args):
    """Start the API server."""
    import uvicorn

    print()
    print("┏━┓╻ ╻┏┓ ┏━┓    ╻ ╻┏━┓╻")
    print("┗━┓┗┳┛┣┻┓┣┳┛    ┣━┫┣━┛┃")
    print("┗━┛ ╹ ┗━┛╹┗╸    ╹ ╹╹  ╹")
    print()
    print(f"Starting SYBR API server on {args.host}:{args.port}")
    print(f"Docs: http://{args.host}:{args.port}/docs")
    print()

    uvicorn.run(
        "api.main:app",
        host=args.host,
        port=args.port,
        reload=args.reload,
        log_level="info",
    )


def main():
    parser = argparse.ArgumentParser(
        prog="sybr-api",
        description="SYBR API management CLI",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # create-key
    p_create = subparsers.add_parser("create-key", help="Create a new API key")
    p_create.add_argument("--name", required=True, help="Human-readable key name")
    p_create.set_defaults(func=cmd_create_key)

    # list-keys
    p_list = subparsers.add_parser("list-keys", help="List all API keys")
    p_list.set_defaults(func=cmd_list_keys)

    # revoke-key
    p_revoke = subparsers.add_parser("revoke-key", help="Revoke an API key")
    p_revoke.add_argument("--id", required=True, type=int, help="Key ID to revoke")
    p_revoke.set_defaults(func=cmd_revoke_key)

    # serve
    p_serve = subparsers.add_parser("serve", help="Start the API server")
    p_serve.add_argument("--host", default=API_HOST, help=f"Bind host (default: {API_HOST})")
    p_serve.add_argument("--port", type=int, default=API_PORT, help=f"Bind port (default: {API_PORT})")
    p_serve.add_argument("--reload", action="store_true", help="Enable auto-reload for development")
    p_serve.set_defaults(func=cmd_serve)

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    args.func(args)


if __name__ == "__main__":
    main()
