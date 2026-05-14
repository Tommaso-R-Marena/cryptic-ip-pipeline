#!/usr/bin/env python3
from crypticip.cli import main
import sys
if __name__ == "__main__":
    raise SystemExit(main(["validate", *sys.argv[1:]]))
