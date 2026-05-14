"""Lightweight, JSON-or-text logging used everywhere in the pipeline."""
from __future__ import annotations

import json
import logging
import os
import sys
import time
from pathlib import Path


_CONFIGURED = False


class _JsonFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        payload = {
            "ts": time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(record.created)),
            "level": record.levelname,
            "name": record.name,
            "msg": record.getMessage(),
        }
        if record.exc_info:
            payload["exc"] = self.formatException(record.exc_info)
        for k, v in record.__dict__.items():
            if k in ("args", "asctime", "created", "exc_info", "exc_text", "filename",
                     "funcName", "levelname", "levelno", "lineno", "message", "module",
                     "msecs", "msg", "name", "pathname", "process", "processName",
                     "relativeCreated", "stack_info", "thread", "threadName"):
                continue
            payload[k] = v
        return json.dumps(payload, default=str)


def configure_logging(level: str = "INFO", json_output: bool = False,
                      logfile: Path | None = None) -> None:
    """Idempotently configure the root logger for the CLI."""
    global _CONFIGURED
    root = logging.getLogger()
    root.setLevel(level.upper())
    if _CONFIGURED:
        return

    handler: logging.Handler = logging.StreamHandler(sys.stderr)
    if json_output or os.environ.get("CRYPTICIP_LOG_JSON"):
        handler.setFormatter(_JsonFormatter())
    else:
        handler.setFormatter(logging.Formatter(
            fmt="%(asctime)s %(levelname)-7s %(name)s | %(message)s",
            datefmt="%H:%M:%S",
        ))
    root.handlers.clear()
    root.addHandler(handler)

    if logfile:
        logfile.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(logfile, mode="a")
        fh.setFormatter(_JsonFormatter())
        root.addHandler(fh)

    _CONFIGURED = True


def get_logger(name: str) -> logging.Logger:
    if not _CONFIGURED:
        configure_logging()
    return logging.getLogger(name)
