"""Shared checkout paths and source-token parsing for repository governance."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = ROOT / "src" / "sirius"
OPERATING_MODEL = ROOT / "tests" / "operating_model.json"
CPP_NON_CODE = re.compile(
    r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|//[^\n]*|/\*.*?\*/',
    re.MULTILINE | re.DOTALL,
)


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()
