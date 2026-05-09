"""Internal utilities for shared CLI config loading."""

from __future__ import annotations

import json
from dataclasses import fields
from pathlib import Path

from .config import AstraConfig


def load_config_json(path: Path) -> AstraConfig:
    """Load AstraConfig values from a JSON file."""
    path = Path(path).expanduser()
    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError("Config JSON must contain an object.")
    allowed = {field.name for field in fields(AstraConfig)}
    unknown = sorted(set(data) - allowed)
    if unknown:
        raise ValueError(f"Unknown AstraConfig field(s): {', '.join(unknown)}")
    config = AstraConfig()
    for key, value in data.items():
        setattr(config, key, value)
    return config
