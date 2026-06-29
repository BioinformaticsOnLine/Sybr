#!/usr/bin/env python3
"""
SYBR API — Quick-start entry point.

Usage:
    python sybr_api.py                    # Start server on :8000
    python sybr_api.py --port 9000        # Custom port
    python sybr_api.py create-key         # Create API key
    python sybr_api.py list-keys          # List keys
"""

import sys

if __name__ == "__main__":
    # If no args or first arg looks like a flag, start the server
    if len(sys.argv) == 1 or sys.argv[1].startswith("-"):
        # Quick serve mode
        sys.argv = [sys.argv[0], "serve"] + sys.argv[1:]

    from api.cli import main
    main()
