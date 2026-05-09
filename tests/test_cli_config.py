from __future__ import annotations

import argparse
import json

import pytest

from astra_xas._config_utils import load_config_json
from astra_xas.cli import _resolve_config
from astra_xas.config import AstraConfig


def test_load_config_json_applies_valid_fields(tmp_path):
    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "e0": 7050.0,
                "analysis_mode": "trans",
            }
        ),
        encoding="utf-8",
    )
    cfg = load_config_json(config_path)
    assert cfg.e0 == 7050.0
    assert cfg.analysis_mode == "trans"


def test_resolve_config_uses_file_values_when_no_overrides(tmp_path):
    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "e0": 7050.0,
                "analysis_mode": "trans",
                "foil_alignment_mode": "ref",
            }
        ),
        encoding="utf-8",
    )
    args = argparse.Namespace(
        config=config_path,
        mode=None,
        foil_mode=None,
        e0=None,
    )
    cfg = _resolve_config(args)
    assert cfg.e0 == 7050.0
    assert cfg.analysis_mode == "trans"
    assert cfg.foil_alignment_mode == "ref"


def test_resolve_config_cli_flags_override_file(tmp_path):
    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "e0": 7050.0,
                "analysis_mode": "trans",
                "foil_alignment_mode": "trans",
            }
        ),
        encoding="utf-8",
    )
    args = argparse.Namespace(
        config=config_path,
        mode="fluo",
        foil_mode="ref",
        e0=7100.0,
    )
    cfg = _resolve_config(args)
    assert cfg.e0 == 7100.0
    assert cfg.analysis_mode == "fluo"
    assert cfg.foil_alignment_mode == "ref"


def test_resolve_config_no_file_uses_defaults():
    args = argparse.Namespace(
        config=None,
        mode=None,
        foil_mode=None,
        e0=None,
    )
    cfg = _resolve_config(args)
    defaults = AstraConfig()
    assert cfg.e0 == defaults.e0
    assert cfg.analysis_mode == defaults.analysis_mode
    assert cfg.foil_alignment_mode == defaults.foil_alignment_mode


def test_load_config_json_unknown_key_raises(tmp_path):
    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "this_is_not_a_real_field": 42,
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError) as excinfo:
        load_config_json(config_path)
    assert "this_is_not_a_real_field" in str(excinfo.value)


def test_load_config_json_missing_file_raises(tmp_path):
    missing = tmp_path / "definitely_does_not_exist.json"
    with pytest.raises(FileNotFoundError):
        load_config_json(missing)


def test_load_config_json_non_object_raises(tmp_path):
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps([1, 2, 3]), encoding="utf-8")
    with pytest.raises(ValueError, match="Config JSON must contain an object"):
        load_config_json(config_path)
