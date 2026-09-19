#!/usr/bin/env python3
"""Fail when source topology or immutable build-input authority regresses."""

from governance.repository import main


if __name__ == "__main__":
    raise SystemExit(main())
