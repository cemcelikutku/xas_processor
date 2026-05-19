"""AstraXAS session manifest (Phase 1).

A manifest is the curated record of which scans to process, how to group
them, and what comments apply. It is produced by Beamtime Mode (Phase 2),
consumed by the offline pipeline (Phase 1 wires this up), and editable by
users in between.

This module owns the schema, dataclasses, JSON load/save, and validation.
It deliberately does NOT consume the manifest in the processing pipeline;
that integration lives in ``processor.process_folder`` and is added in a
follow-up patch.

Schema overview (v1):

  - ``schema_version`` (int, required): currently ``1``.
  - ``astra_version`` (str): AstraXAS version that wrote the file.
  - ``created_iso``, ``updated_iso`` (str): set by ``save_manifest``.
  - ``session``: free-form metadata (operator, beamline, notes, dates).
  - ``input.base_dir``: directory containing the .xasd files. May be
    absolute or relative; relative paths resolve against the manifest
    file's parent directory.
  - ``config.source``: ``"path"`` or ``"inline"``.
        * ``"path"`` requires ``config.path``, a JSON file path resolved
          like ``input.base_dir``.
        * ``"inline"`` requires ``config.inline``, a JSON object containing
          a full AstraConfig. v1 validates the shape; Phase 2 implements
          actual consumption.
  - ``scans[]``: list of scan entries. Each has ``filename``, ``path``
    (relative to ``input.base_dir``), ``status``, optional ``is_foil``
    override, optional ``group`` reference, optional ``assigned_foil``
    (reserved for v2), ``comments[]``, and an optional ``qc`` snapshot.

    The order of ``scans[]`` is significant: within each group, scans
    are processed (and merged) in the order they appear in this array.
    This overrides the default folder-mode behaviour of sorting by the
    ``_NN`` replicate suffix.

    ``is_foil`` overrides ``AstraConfig.foil_keyword`` filename inference
    when non-null. ``true`` forces foil treatment; ``false`` forces
    sample treatment; ``null`` falls back to the filename heuristic.

  - ``groups[]``: list of group definitions with stable ``id`` and
    human-facing ``label``/``comment``. **Output filenames in manifest
    mode derive from ``id``, not ``label``.** Choose ``id`` to be a
    safe filename token (alphanumerics, dashes, underscores); use
    ``label`` for the human-facing display name.

The ``status`` field is *human curation*:
    - ``"unreviewed"`` — default, watcher-produced, awaiting review.
    - ``"good"``       — include in processing.
    - ``"bad"``        — data is wrong (failed, glitched, etc.).
    - ``"excluded"``   — fine data, not part of this analysis.

The ``qc`` field is the most recent *automated* check from Beamtime Mode
(parsing, jumps, finite-channel checks). ``qc`` and ``status`` are
independent; they can disagree. A scan can be ``qc.auto_status == "warn"``
yet ``status == "good"`` (user reviewed and confirmed) or vice versa.
"""

from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any


SCHEMA_VERSION = 1

VALID_STATUSES = ("unreviewed", "good", "bad", "excluded")
VALID_CONFIG_SOURCES = ("path", "inline")
# QC auto_status values are advisory; manifest validation does not enforce
# the vocabulary so that future watcher revisions can introduce new ones.
CANONICAL_AUTO_STATUSES = ("ok", "warn", "reject", "pending")


class ManifestError(Exception):
    """Base error for manifest load/save operations."""


class ManifestValidationError(ManifestError):
    """Raised when a manifest fails schema validation."""

    def __init__(self, errors: list[str]):
        self.errors = list(errors)
        super().__init__(
            "Manifest validation failed:\n" + "\n".join(f"- {e}" for e in self.errors)
        )


@dataclass
class Comment:
    iso: str = ""
    author: str = ""
    text: str = ""


@dataclass
class QCSummary:
    timestamp_iso: str = ""
    auto_status: str = ""
    n_warnings: int = 0
    n_jumps: int = 0
    notes: str = ""


@dataclass
class ScanEntry:
    filename: str = ""
    path: str = ""
    status: str = "unreviewed"
    is_foil: bool | None = None
    group: str | None = None
    # Reserved for v2 (per-scan foil override). v1 accepts and preserves
    # this field but the offline pipeline ignores it.
    assigned_foil: str | None = None
    comments: list[Comment] = field(default_factory=list)
    qc: QCSummary | None = None


@dataclass
class GroupEntry:
    id: str = ""
    label: str = ""
    comment: str = ""


@dataclass
class SessionMeta:
    name: str = ""
    operator: str = ""
    beamline: str = ""
    notes: str = ""
    started_iso: str = ""
    ended_iso: str = ""


@dataclass
class InputSpec:
    base_dir: str = ""


@dataclass
class ConfigRef:
    source: str = "path"
    path: str | None = None
    inline: dict[str, Any] | None = None


@dataclass
class Manifest:
    schema_version: int = SCHEMA_VERSION
    astra_version: str = ""
    created_iso: str = ""
    updated_iso: str = ""
    session: SessionMeta = field(default_factory=SessionMeta)
    input: InputSpec = field(default_factory=InputSpec)
    config: ConfigRef = field(default_factory=ConfigRef)
    scans: list[ScanEntry] = field(default_factory=list)
    groups: list[GroupEntry] = field(default_factory=list)
    # Not serialized. Set by load_manifest() and save_manifest() so that
    # relative paths in the manifest can be resolved against the file's
    # parent directory.
    source_path: Path | None = field(default=None, repr=False, compare=False)

    def manifest_dir(self) -> Path:
        if self.source_path is None:
            raise ManifestError(
                "Manifest has no source_path; cannot resolve relative paths. "
                "Set manifest.source_path or load via load_manifest()."
            )
        return Path(self.source_path).expanduser().resolve().parent

    def resolve_base_dir(self) -> Path:
        return _resolve_against(self.input.base_dir, self.manifest_dir())

    def resolve_config_path(self) -> Path | None:
        if self.config.source != "path" or not self.config.path:
            return None
        return _resolve_against(self.config.path, self.manifest_dir())

    def resolve_scan_path(self, scan: ScanEntry) -> Path:
        return _resolve_against(scan.path, self.resolve_base_dir())


def _resolve_against(value: str, anchor: Path) -> Path:
    p = Path(value).expanduser()
    if p.is_absolute():
        return p.resolve()
    return (anchor / p).resolve()


def _now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


# ---------------------------------------------------------------------------
# Parsing (dict -> dataclasses). Tolerant of missing optional fields and
# wrong types: best-effort extraction, with semantic checks deferred to
# validate_manifest().
# ---------------------------------------------------------------------------

def _str(value: Any, default: str = "") -> str:
    if value is None:
        return default
    return str(value)


def _str_or_none(value: Any) -> str | None:
    if value is None:
        return None
    return str(value)


def _bool_or_none(value: Any) -> bool | None:
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    return bool(value)


def _int(value: Any, default: int = 0) -> int:
    if value is None:
        return default
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _dict_or_empty(value: Any) -> dict:
    return value if isinstance(value, dict) else {}


def _list_or_empty(value: Any) -> list:
    return value if isinstance(value, list) else []


def _comment_from_dict(d: Any) -> Comment:
    d = _dict_or_empty(d)
    return Comment(
        iso=_str(d.get("iso")),
        author=_str(d.get("author")),
        text=_str(d.get("text")),
    )


def _qc_from_dict(d: Any) -> QCSummary | None:
    if d is None:
        return None
    d = _dict_or_empty(d)
    return QCSummary(
        timestamp_iso=_str(d.get("timestamp_iso")),
        auto_status=_str(d.get("auto_status")),
        n_warnings=_int(d.get("n_warnings")),
        n_jumps=_int(d.get("n_jumps")),
        notes=_str(d.get("notes")),
    )


def _scan_from_dict(d: Any) -> ScanEntry:
    d = _dict_or_empty(d)
    return ScanEntry(
        filename=_str(d.get("filename")),
        path=_str(d.get("path")),
        status=_str(d.get("status"), default="unreviewed"),
        is_foil=_bool_or_none(d.get("is_foil")),
        group=_str_or_none(d.get("group")),
        assigned_foil=_str_or_none(d.get("assigned_foil")),
        comments=[_comment_from_dict(c) for c in _list_or_empty(d.get("comments"))],
        qc=_qc_from_dict(d.get("qc")) if "qc" in d else None,
    )


def _group_from_dict(d: Any) -> GroupEntry:
    d = _dict_or_empty(d)
    return GroupEntry(
        id=_str(d.get("id")),
        label=_str(d.get("label")),
        comment=_str(d.get("comment")),
    )


def _session_from_dict(d: Any) -> SessionMeta:
    d = _dict_or_empty(d)
    return SessionMeta(
        name=_str(d.get("name")),
        operator=_str(d.get("operator")),
        beamline=_str(d.get("beamline")),
        notes=_str(d.get("notes")),
        started_iso=_str(d.get("started_iso")),
        ended_iso=_str(d.get("ended_iso")),
    )


def _input_from_dict(d: Any) -> InputSpec:
    d = _dict_or_empty(d)
    return InputSpec(base_dir=_str(d.get("base_dir")))


def _config_from_dict(d: Any) -> ConfigRef:
    d = _dict_or_empty(d)
    inline = d.get("inline")
    return ConfigRef(
        source=_str(d.get("source"), default="path"),
        path=_str_or_none(d.get("path")),
        inline=inline if isinstance(inline, dict) else None,
    )


def _from_dict(data: dict) -> Manifest:
    return Manifest(
        schema_version=_int(data.get("schema_version")),
        astra_version=_str(data.get("astra_version")),
        created_iso=_str(data.get("created_iso")),
        updated_iso=_str(data.get("updated_iso")),
        session=_session_from_dict(data.get("session")),
        input=_input_from_dict(data.get("input")),
        config=_config_from_dict(data.get("config")),
        scans=[_scan_from_dict(s) for s in _list_or_empty(data.get("scans"))],
        groups=[_group_from_dict(g) for g in _list_or_empty(data.get("groups"))],
    )


# ---------------------------------------------------------------------------
# Serialisation (dataclasses -> dict). Null values are written explicitly
# so the schema stays visible to readers of a hand-edited manifest.
# ---------------------------------------------------------------------------

def _config_to_dict(c: ConfigRef) -> dict[str, Any]:
    out: dict[str, Any] = {"source": c.source}
    # Always include the alternate field as null so the schema is visible.
    out["path"] = c.path
    out["inline"] = dict(c.inline) if c.inline is not None else None
    return out


def _scan_to_dict(s: ScanEntry) -> dict[str, Any]:
    return {
        "filename": s.filename,
        "path": s.path,
        "status": s.status,
        "is_foil": s.is_foil,
        "group": s.group,
        "assigned_foil": s.assigned_foil,
        "comments": [asdict(c) for c in s.comments],
        "qc": asdict(s.qc) if s.qc is not None else None,
    }


def _to_dict(manifest: Manifest) -> dict[str, Any]:
    return {
        "schema_version": manifest.schema_version,
        "astra_version": manifest.astra_version,
        "created_iso": manifest.created_iso,
        "updated_iso": manifest.updated_iso,
        "session": asdict(manifest.session),
        "input": asdict(manifest.input),
        "config": _config_to_dict(manifest.config),
        "scans": [_scan_to_dict(s) for s in manifest.scans],
        "groups": [asdict(g) for g in manifest.groups],
    }


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_manifest(manifest: Manifest) -> list[str]:
    """Return a list of validation error messages. Empty list means valid."""
    errors: list[str] = []

    if manifest.schema_version != SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {SCHEMA_VERSION}, got {manifest.schema_version!r}"
        )

    if not manifest.input.base_dir:
        errors.append("input.base_dir is required and must be non-empty")

    source = manifest.config.source
    if source not in VALID_CONFIG_SOURCES:
        errors.append(
            f"config.source must be one of {VALID_CONFIG_SOURCES}, got {source!r}"
        )
    elif source == "path":
        if not manifest.config.path:
            errors.append(
                "config.source='path' requires config.path to be a non-empty string"
            )
        if manifest.config.inline is not None:
            errors.append("config.inline must be omitted when config.source='path'")
    elif source == "inline":
        if manifest.config.inline is None:
            errors.append(
                "config.source='inline' requires config.inline to be a JSON object"
            )
        if manifest.config.path:
            errors.append("config.path must be omitted when config.source='inline'")

    group_ids: dict[str, int] = {}
    for i, group in enumerate(manifest.groups):
        if not group.id:
            errors.append(f"groups[{i}].id is required and must be non-empty")
            continue
        if group.id in group_ids:
            errors.append(
                f"groups[{i}].id={group.id!r} duplicates groups[{group_ids[group.id]}].id"
            )
        else:
            group_ids[group.id] = i

    seen_filenames: dict[str, int] = {}
    for i, scan in enumerate(manifest.scans):
        if not scan.filename:
            errors.append(f"scans[{i}].filename is required and must be non-empty")
        elif scan.filename in seen_filenames:
            errors.append(
                f"scans[{i}].filename={scan.filename!r} duplicates "
                f"scans[{seen_filenames[scan.filename]}].filename"
            )
        else:
            seen_filenames[scan.filename] = i

        if not scan.path:
            errors.append(f"scans[{i}].path is required and must be non-empty")

        if scan.status not in VALID_STATUSES:
            errors.append(
                f"scans[{i}].status={scan.status!r} is not one of {VALID_STATUSES}"
            )

        if scan.group is not None and scan.group not in group_ids:
            errors.append(
                f"scans[{i}].group={scan.group!r} does not match any groups[].id"
            )

        if scan.qc is not None:
            if scan.qc.n_warnings < 0:
                errors.append(f"scans[{i}].qc.n_warnings must be >= 0")
            if scan.qc.n_jumps < 0:
                errors.append(f"scans[{i}].qc.n_jumps must be >= 0")

    return errors


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def load_manifest(path: str | Path) -> Manifest:
    """Read a session manifest from disk and validate it.

    Raises:
        ManifestError if the file is missing or not valid JSON.
        ManifestValidationError if the contents fail schema validation.
    """
    path = Path(path).expanduser().resolve()
    if not path.exists():
        raise ManifestError(f"Manifest file not found: {path}")
    try:
        with path.open("r", encoding="utf-8") as f:
            data = json.load(f)
    except json.JSONDecodeError as exc:
        raise ManifestError(f"Manifest is not valid JSON ({path}): {exc}") from exc

    if not isinstance(data, dict):
        raise ManifestValidationError(["top-level JSON value must be an object"])

    manifest = _from_dict(data)
    manifest.source_path = path

    errors = validate_manifest(manifest)
    if errors:
        raise ManifestValidationError(errors)
    return manifest


def save_manifest(manifest: Manifest, path: str | Path | None = None) -> Path:
    """Validate the manifest and write it to disk atomically.

    Sets ``created_iso`` on first save and ``updated_iso`` on every save.
    If ``path`` is omitted, writes back to ``manifest.source_path``.

    Returns the resolved path that was written.
    """
    if path is None:
        if manifest.source_path is None:
            raise ManifestError(
                "save_manifest requires a path; manifest has no source_path."
            )
        path = manifest.source_path
    path = Path(path).expanduser().resolve()

    now = _now_iso()
    if not manifest.created_iso:
        manifest.created_iso = now
    manifest.updated_iso = now

    errors = validate_manifest(manifest)
    if errors:
        raise ManifestValidationError(errors)

    data = _to_dict(manifest)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    try:
        with tmp.open("w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, sort_keys=True)
            f.write("\n")
        os.replace(tmp, path)
    except Exception:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass
        raise

    manifest.source_path = path
    return path
