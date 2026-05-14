#!/usr/bin/env python3
"""Thin wrapper around ``crypticip download-validation``."""
from crypticip.cli import main
import sys
if __name__ == "__main__":
    raise SystemExit(main(["download-validation", *sys.argv[1:]]))
