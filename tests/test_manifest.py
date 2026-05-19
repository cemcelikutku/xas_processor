from __future__ import annotations

import json

import pytest

from astra_xas import manifest as manifest_module
from astra_xas.manifest import (
    SCHEMA_VERSION,
    Comment,
    ConfigRef,
    GroupEntry,
    InputSpec,
    Manifest,
    ManifestError,
    ManifestValidationError,
    QCSummary,
    ScanEntry,
    SessionMeta,
    load_manifest,
    save_manifest,
    validate_manifest,
)


def _minimal_manifest() -> Manifest:
    return Manifest(
        input=InputSpec(base_dir="data"),
        config=ConfigRef(source="path", path="cfg.json"),
        groups=[GroupEntry(id="g1", label="Group 1")],
        scans=[
            ScanEntry(
                filename="a_001.xasd",
                path="a_001.xasd",
                status="unreviewed",
                group="g1",
            ),
        ],
    )


# ---------------------------------------------------------------------------
# Round-trip
# ---------------------------------------------------------------------------

def test_save_then_load_roundtrip(tmp_path):
    manifest = _minimal_manifest()
    out = tmp_path / "session.json"
    save_manifest(manifest, out)

    loaded = load_manifest(out)
    assert loaded.schema_version == SCHEMA_VERSION
    assert loaded.input.base_dir == "data"
    assert loaded.config.source == "path"
    assert loaded.config.path == "cfg.json"
    assert loaded.config.inline is None
    assert [g.id for g in loaded.groups] == ["g1"]
    assert [s.filename for s in loaded.scans] == ["a_001.xasd"]
    assert loaded.scans[0].status == "unreviewed"
    assert loaded.scans[0].group == "g1"
    assert loaded.scans[0].is_foil is None
    assert loaded.scans[0].assigned_foil is None
    assert loaded.scans[0].comments == []
    assert loaded.scans[0].qc is None


def test_saved_file_uses_indented_json(tmp_path):
    out = tmp_path / "session.json"
    save_manifest(_minimal_manifest(), out)
    text = out.read_text(encoding="utf-8")
    # Indented, sorted keys, trailing newline.
    assert text.endswith("\n")
    assert "\n  " in text
    parsed = json.loads(text)
    assert parsed["schema_version"] == SCHEMA_VERSION


# ---------------------------------------------------------------------------
# Timestamps
# ---------------------------------------------------------------------------

def test_first_save_sets_both_timestamps(tmp_path, monkeypatch):
    monkeypatch.setattr(manifest_module, "_now_iso", lambda: "2026-01-01T10:00:00")
    manifest = _minimal_manifest()
    save_manifest(manifest, tmp_path / "session.json")
    assert manifest.created_iso == "2026-01-01T10:00:00"
    assert manifest.updated_iso == "2026-01-01T10:00:00"


def test_resave_preserves_created_updates_updated(tmp_path, monkeypatch):
    times = iter(["2026-01-01T10:00:00", "2026-01-01T11:00:00"])
    monkeypatch.setattr(manifest_module, "_now_iso", lambda: next(times))
    manifest = _minimal_manifest()
    out = tmp_path / "session.json"
    save_manifest(manifest, out)
    save_manifest(manifest, out)
    assert manifest.created_iso == "2026-01-01T10:00:00"
    assert manifest.updated_iso == "2026-01-01T11:00:00"


def test_load_preserves_created_iso_from_disk(tmp_path, monkeypatch):
    times = iter(["2026-01-01T10:00:00", "2026-02-02T20:00:00"])
    monkeypatch.setattr(manifest_module, "_now_iso", lambda: next(times))
    out = tmp_path / "session.json"
    save_manifest(_minimal_manifest(), out)
    loaded = load_manifest(out)
    assert loaded.created_iso == "2026-01-01T10:00:00"
    save_manifest(loaded)
    assert loaded.created_iso == "2026-01-01T10:00:00"
    assert loaded.updated_iso == "2026-02-02T20:00:00"


# ---------------------------------------------------------------------------
# Validation: schema version, statuses, required fields
# ---------------------------------------------------------------------------

def test_invalid_schema_version_rejected(tmp_path):
    manifest = _minimal_manifest()
    manifest.schema_version = 99
    with pytest.raises(ManifestValidationError) as exc:
        save_manifest(manifest, tmp_path / "session.json")
    assert "schema_version" in str(exc.value)


@pytest.mark.parametrize("status", ["unreviewed", "good", "bad", "excluded"])
def test_all_valid_statuses_accepted(tmp_path, status):
    manifest = _minimal_manifest()
    manifest.scans[0].status = status
    save_manifest(manifest, tmp_path / "session.json")


def test_invalid_status_rejected(tmp_path):
    manifest = _minimal_manifest()
    manifest.scans[0].status = "maybe"
    with pytest.raises(ManifestValidationError) as exc:
        save_manifest(manifest, tmp_path / "session.json")
    assert "status" in str(exc.value)


def test_missing_base_dir_rejected(tmp_path):
    manifest = _minimal_manifest()
    manifest.input.base_dir = ""
    with pytest.raises(ManifestValidationError):
        save_manifest(manifest, tmp_path / "session.json")


def test_missing_scan_filename_rejected(tmp_path):
    manifest = _minimal_manifest()
    manifest.scans[0].filename = ""
    with pytest.raises(ManifestValidationError):
        save_manifest(manifest, tmp_path / "session.json")


def test_missing_scan_path_rejected(tmp_path):
    manifest = _minimal_manifest()
    manifest.scans[0].path = ""
    with pytest.raises(ManifestValidationError):
        save_manifest(manifest, tmp_path / "session.json")


def test_duplicate_scan_filename_rejected(tmp_path):
    manifest = _minimal_manifest()
    manifest.scans.append(
        ScanEntry(filename="a_001.xasd", path="other/a_001.xasd", group="g1")
    )
    with pytest.raises(ManifestValidationError) as exc:
        save_manifest(manifest, tmp_path / "session.json")
    assert "duplicates" in str(exc.value)


def test_duplicate_group_id_rejected(tmp_path):
    manifest = _minimal_manifest()
    manifest.groups.append(GroupEntry(id="g1", label="Conflicting"))
    with pytest.raises(ManifestValidationError) as exc:
        save_manifest(manifest, tmp_path / "session.json")
    assert "g1" in str(exc.value)


def test_empty_group_id_rejected(tmp_path):
    manifest = _minimal_manifest()
    manifest.groups.append(GroupEntry(id="", label="Nameless"))
    with pytest.raises(ManifestValidationError):
        save_manifest(manifest, tmp_path / "session.json")


# ---------------------------------------------------------------------------
# Validation: group cross-references
# ---------------------------------------------------------------------------

def test_scan_referencing_unknown_group_rejected(tmp_path):
    manifest = _minimal_manifest()
    manifest.scans[0].group = "ghost"
    with pytest.raises(ManifestValidationError) as exc:
        save_manifest(manifest, tmp_path / "session.json")
    assert "ghost" in str(exc.value)


def test_scan_with_null_group_is_allowed(tmp_path):
    manifest = _minimal_manifest()
    manifest.scans[0].group = None
    save_manifest(manifest, tmp_path / "session.json")


# ---------------------------------------------------------------------------
# Validation: config source variants
# ---------------------------------------------------------------------------

def test_config_source_path_requires_path(tmp_path):
    manifest = _minimal_manifest()
    manifest.config = ConfigRef(source="path", path=None)
    with pytest.raises(ManifestValidationError):
        save_manifest(manifest, tmp_path / "session.json")


def test_config_source_path_rejects_inline(tmp_path):
    manifest = _minimal_manifest()
    manifest.config = ConfigRef(source="path", path="cfg.json", inline={"foo": 1})
    with pytest.raises(ManifestValidationError):
        save_manifest(manifest, tmp_path / "session.json")


def test_config_source_inline_requires_inline(tmp_path):
    manifest = _minimal_manifest()
    manifest.config = ConfigRef(source="inline", inline=None)
    with pytest.raises(ManifestValidationError):
        save_manifest(manifest, tmp_path / "session.json")


def test_config_source_inline_rejects_path(tmp_path):
    manifest = _minimal_manifest()
    manifest.config = ConfigRef(source="inline", path="cfg.json", inline={"e0": 7000})
    with pytest.raises(ManifestValidationError):
        save_manifest(manifest, tmp_path / "session.json")


def test_config_source_inline_accepted_and_preserved(tmp_path):
    """v1 validates the shape of inline configs and round-trips them.

    Phase 2 will actually apply the inline config; for now we just keep it
    intact so a hand-edited or watcher-written inline manifest is durable.
    """
    manifest = _minimal_manifest()
    manifest.config = ConfigRef(
        source="inline",
        inline={"analysis_mode": "fluo", "e0": 7121.03},
    )
    out = tmp_path / "session.json"
    save_manifest(manifest, out)
    loaded = load_manifest(out)
    assert loaded.config.source == "inline"
    assert loaded.config.path is None
    assert loaded.config.inline == {"analysis_mode": "fluo", "e0": 7121.03}


def test_config_source_unknown_rejected(tmp_path):
    manifest = _minimal_manifest()
    manifest.config = ConfigRef(source="magic")
    with pytest.raises(ManifestValidationError):
        save_manifest(manifest, tmp_path / "session.json")


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------

def test_relative_base_dir_resolves_to_manifest_parent(tmp_path):
    manifest = _minimal_manifest()
    manifest.input.base_dir = "incoming"
    out = tmp_path / "session.json"
    save_manifest(manifest, out)
    loaded = load_manifest(out)
    assert loaded.resolve_base_dir() == (tmp_path / "incoming").resolve()


def test_absolute_base_dir_used_as_is(tmp_path):
    elsewhere = (tmp_path / "elsewhere" / "data").resolve()
    elsewhere.mkdir(parents=True)
    manifest = _minimal_manifest()
    manifest.input.base_dir = str(elsewhere)
    out = tmp_path / "session.json"
    save_manifest(manifest, out)
    loaded = load_manifest(out)
    assert loaded.resolve_base_dir() == elsewhere


def test_relative_config_path_resolves_to_manifest_parent(tmp_path):
    manifest = _minimal_manifest()
    manifest.config.path = "configs/k.json"
    out = tmp_path / "session.json"
    save_manifest(manifest, out)
    loaded = load_manifest(out)
    assert loaded.resolve_config_path() == (tmp_path / "configs" / "k.json").resolve()


def test_absolute_config_path_used_as_is(tmp_path):
    cfg_path = (tmp_path / "shared" / "cfg.json").resolve()
    manifest = _minimal_manifest()
    manifest.config.path = str(cfg_path)
    out = tmp_path / "session.json"
    save_manifest(manifest, out)
    loaded = load_manifest(out)
    assert loaded.resolve_config_path() == cfg_path


def test_resolve_config_path_none_when_source_is_inline(tmp_path):
    manifest = _minimal_manifest()
    manifest.config = ConfigRef(source="inline", inline={"e0": 7000.0})
    out = tmp_path / "session.json"
    save_manifest(manifest, out)
    loaded = load_manifest(out)
    assert loaded.resolve_config_path() is None


def test_scan_path_resolves_relative_to_base_dir(tmp_path):
    manifest = _minimal_manifest()
    manifest.input.base_dir = "incoming"
    out = tmp_path / "session.json"
    save_manifest(manifest, out)
    loaded = load_manifest(out)
    expected = (tmp_path / "incoming" / "a_001.xasd").resolve()
    assert loaded.resolve_scan_path(loaded.scans[0]) == expected


def test_scan_path_absolute_used_as_is(tmp_path):
    abs_path = (tmp_path / "elsewhere" / "weird_name.xasd").resolve()
    manifest = _minimal_manifest()
    manifest.scans[0].path = str(abs_path)
    out = tmp_path / "session.json"
    save_manifest(manifest, out)
    loaded = load_manifest(out)
    assert loaded.resolve_scan_path(loaded.scans[0]) == abs_path


def test_manifest_dir_requires_source_path():
    manifest = _minimal_manifest()
    with pytest.raises(ManifestError):
        manifest.manifest_dir()


# ---------------------------------------------------------------------------
# Field preservation through round-trip
# ---------------------------------------------------------------------------

def test_assigned_foil_preserved_through_roundtrip(tmp_path):
    """Reserved-for-v2: the field is accepted and preserved even though
    the offline pipeline ignores it."""
    manifest = _minimal_manifest()
    manifest.scans[0].assigned_foil = "foil_001.xasd"
    out = tmp_path / "session.json"
    save_manifest(manifest, out)
    loaded = load_manifest(out)
    assert loaded.scans[0].assigned_foil == "foil_001.xasd"


def test_comments_preserved_through_roundtrip(tmp_path):
    manifest = _minimal_manifest()
    manifest.scans[0].comments = [
        Comment(iso="2026-05-13T10:00:00", author="Cem", text="Looks good"),
        Comment(iso="2026-05-13T11:00:00", author="BLS", text="Confirmed."),
    ]
    out = tmp_path / "session.json"
    save_manifest(manifest, out)
    loaded = load_manifest(out)
    assert len(loaded.scans[0].comments) == 2
    assert loaded.scans[0].comments[0].text == "Looks good"
    assert loaded.scans[0].comments[1].author == "BLS"


def test_qc_summary_preserved_through_roundtrip(tmp_path):
    manifest = _minimal_manifest()
    manifest.scans[0].qc = QCSummary(
        timestamp_iso="2026-05-13T10:00:00",
        auto_status="warn",
        n_warnings=2,
        n_jumps=1,
        notes="detector jump near edge",
    )
    out = tmp_path / "session.json"
    save_manifest(manifest, out)
    loaded = load_manifest(out)
    assert loaded.scans[0].qc is not None
    assert loaded.scans[0].qc.auto_status == "warn"
    assert loaded.scans[0].qc.n_warnings == 2
    assert loaded.scans[0].qc.n_jumps == 1
    assert loaded.scans[0].qc.notes == "detector jump near edge"


def test_qc_can_disagree_with_status(tmp_path):
    """Roadmap requirement: qc (automated) and status (curation) are
    independent. A scan can have qc.auto_status='warn' and status='good'."""
    manifest = _minimal_manifest()
    manifest.scans[0].status = "good"
    manifest.scans[0].qc = QCSummary(auto_status="warn", n_warnings=3)
    out = tmp_path / "session.json"
    save_manifest(manifest, out)
    loaded = load_manifest(out)
    assert loaded.scans[0].status == "good"
    assert loaded.scans[0].qc.auto_status == "warn"


def test_qc_can_be_null(tmp_path):
    manifest = _minimal_manifest()
    manifest.scans[0].qc = None
    out = tmp_path / "session.json"
    save_manifest(manifest, out)
    loaded = load_manifest(out)
    assert loaded.scans[0].qc is None


@pytest.mark.parametrize("value", [None, True, False])
def test_is_foil_tri_state_preserved(tmp_path, value):
    manifest = _minimal_manifest()
    manifest.scans[0].is_foil = value
    out = tmp_path / "session.json"
    save_manifest(manifest, out)
    loaded = load_manifest(out)
    assert loaded.scans[0].is_foil is value


def test_session_metadata_preserved(tmp_path):
    manifest = _minimal_manifest()
    manifest.session = SessionMeta(
        name="SOLARIS_K_kedge_20260513",
        operator="Cem Celikutku",
        beamline="SOLARIS / ASTRA",
        notes="K K-edge operando",
        started_iso="2026-05-13T09:00:00",
        ended_iso="2026-05-14T18:00:00",
    )
    out = tmp_path / "session.json"
    save_manifest(manifest, out)
    loaded = load_manifest(out)
    assert loaded.session.name == "SOLARIS_K_kedge_20260513"
    assert loaded.session.beamline == "SOLARIS / ASTRA"
    assert loaded.session.ended_iso == "2026-05-14T18:00:00"


# ---------------------------------------------------------------------------
# Loading errors
# ---------------------------------------------------------------------------

def test_load_nonexistent_file_raises(tmp_path):
    with pytest.raises(ManifestError):
        load_manifest(tmp_path / "missing.json")


def test_load_invalid_json_raises(tmp_path):
    p = tmp_path / "bad.json"
    p.write_text("{ not valid json", encoding="utf-8")
    with pytest.raises(ManifestError):
        load_manifest(p)


def test_load_non_object_top_level_raises(tmp_path):
    p = tmp_path / "list.json"
    p.write_text("[]", encoding="utf-8")
    with pytest.raises(ManifestValidationError):
        load_manifest(p)


def test_load_collects_multiple_errors(tmp_path):
    """Validation reports all problems in one error, not just the first."""
    p = tmp_path / "bad.json"
    p.write_text(
        json.dumps(
            {
                "schema_version": 99,
                "input": {"base_dir": ""},
                "config": {"source": "magic"},
                "scans": [{"filename": "a.xasd", "path": "", "status": "weird"}],
                "groups": [],
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ManifestValidationError) as exc:
        load_manifest(p)
    assert len(exc.value.errors) >= 4


# ---------------------------------------------------------------------------
# save_manifest plumbing
# ---------------------------------------------------------------------------

def test_save_is_atomic_no_tmp_remains(tmp_path):
    out = tmp_path / "session.json"
    save_manifest(_minimal_manifest(), out)
    assert list(tmp_path.glob("*.tmp")) == []


def test_save_without_path_uses_source_path(tmp_path):
    out = tmp_path / "session.json"
    save_manifest(_minimal_manifest(), out)
    loaded = load_manifest(out)
    loaded.session.notes = "edited"
    save_manifest(loaded)
    assert load_manifest(out).session.notes == "edited"


def test_save_without_path_or_source_path_raises():
    with pytest.raises(ManifestError):
        save_manifest(_minimal_manifest())


def test_save_creates_parent_directories(tmp_path):
    out = tmp_path / "deep" / "nested" / "session.json"
    save_manifest(_minimal_manifest(), out)
    assert out.exists()


# ---------------------------------------------------------------------------
# validate_manifest standalone
# ---------------------------------------------------------------------------

def test_validate_manifest_returns_empty_list_for_valid():
    assert validate_manifest(_minimal_manifest()) == []


def test_validate_manifest_returns_errors_without_raising():
    manifest = _minimal_manifest()
    manifest.schema_version = 99
    manifest.scans[0].status = "weird"
    errors = validate_manifest(manifest)
    assert len(errors) == 2
