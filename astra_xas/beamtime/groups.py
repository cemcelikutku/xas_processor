from __future__ import annotations

import json
import os
import threading
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from astra_xas.config import AstraConfig
from astra_xas.export import save_two_col
from astra_xas.io import load_xasd, natural_key, sanitize_name, split_replicate_suffix
from astra_xas.processor import _run_pre_edge, interpolate_to_grid
from astra_xas.single_scan import _entry_from_scan


PLOT_DPI = 120


@dataclass
class GroupState:
    base_name: str
    accepted_filenames: list[str] = field(default_factory=list)
    accepted_paths: list[str] = field(default_factory=list)
    n_accepted: int = 0
    last_updated_iso: str = ""
    last_merge_status: str = "pending"
    last_merge_error: str = ""
    output_files: dict[str, str] = field(default_factory=dict)


def _timestamp() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _sort_pairs(filenames: list[str], paths: list[str]) -> tuple[list[str], list[str]]:
    pairs = list(zip(filenames, paths))

    def key(pair):
        filename = pair[0]
        _, replicate_id = split_replicate_suffix(filename)
        if replicate_id is not None:
            return (0, replicate_id, natural_key(filename))
        return (1, natural_key(filename))

    pairs.sort(key=key)
    return [name for name, _ in pairs], [path for _, path in pairs]


def _state_to_dict(state: GroupState) -> dict:
    return asdict(state)


def _state_from_dict(data: dict) -> GroupState:
    return GroupState(
        base_name=str(data.get("base_name", "")),
        accepted_filenames=[str(v) for v in data.get("accepted_filenames", [])],
        accepted_paths=[str(v) for v in data.get("accepted_paths", [])],
        n_accepted=int(data.get("n_accepted", 0)),
        last_updated_iso=str(data.get("last_updated_iso", "")),
        last_merge_status=str(data.get("last_merge_status", "pending")),
        last_merge_error=str(data.get("last_merge_error", "")),
        output_files={str(k): str(v) for k, v in dict(data.get("output_files", {})).items()},
    )


def _relative(path: Path, output_dir: Path) -> str:
    return path.relative_to(output_dir).as_posix()


def _summary_path(output_dir: Path, base_name: str) -> Path:
    return output_dir / "groups" / f"{sanitize_name(base_name)}_group_summary.json"


def _write_group_summary_atomic(path: Path, state_or_data) -> None:
    data = _state_to_dict(state_or_data) if isinstance(state_or_data, GroupState) else dict(state_or_data)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)
        f.write("\n")
    os.replace(tmp, path)


def _signal_from_entry(entry: dict, config: AstraConfig):
    mode = getattr(config, "analysis_mode", "fluo")
    key = {
        "fluo": "mu_fluo",
        "trans": "mu_trans",
        "ref": "mu_ref",
    }.get(mode, "mu_fluo")
    return np.asarray(entry[key], dtype=float)


def merge_and_normalize_group(
    group_snapshot: dict,
    config: AstraConfig,
    log=print,
) -> dict:
    paths = [Path(path) for path in group_snapshot.get("accepted_paths", [])]
    filenames = list(group_snapshot.get("accepted_filenames", []))
    if len(paths) < 2:
        raise ValueError("At least two accepted paths are required for group merging.")

    entries = []
    for path in paths:
        scan = load_xasd(path)
        entries.append(_entry_from_scan(scan, config, path=path))

    master_energy = np.asarray(entries[0]["energy"], dtype=float)
    replicate_signals = []
    for entry in entries:
        signal = _signal_from_entry(entry, config)
        interp = interpolate_to_grid(entry["energy"], signal, master_energy, kind=config.interp_kind)
        valid = np.isfinite(interp)
        if np.count_nonzero(valid) < 0.9 * len(master_energy):
            raise ValueError(f"Insufficient interpolation overlap for {entry.get('filename', 'unknown')}")
        if not np.all(valid):
            interp = np.interp(master_energy, master_energy[valid], interp[valid])
        replicate_signals.append(np.asarray(interp, dtype=float))

    replicate_array = np.asarray(replicate_signals, dtype=float)
    mu_avg = np.nanmean(replicate_array, axis=0)
    pe = _run_pre_edge(master_energy, mu_avg, config)
    return {
        "energy": master_energy,
        "mu_avg": mu_avg,
        "norm": np.asarray(pe["norm"], dtype=float),
        "flat": np.asarray(pe["flat"], dtype=float),
        "edge_step": float(pe["edge_step"]),
        "e0": float(pe["e0"]),
        "n_used": len(replicate_signals),
        "replicate_signals": replicate_signals,
        "replicate_filenames": filenames,
    }


def _header(base_name: str, filenames: list[str], generated: str, column: str) -> str:
    return (
        "AstraXAS Beamtime Mode (live)\n"
        f"group: {base_name}\n"
        f"n_replicates: {len(filenames)}\n"
        f"source files: {', '.join(filenames)}\n"
        f"generated: {generated}\n"
        f"energy_eV   {column}"
    )


def _write_two_col_atomic(path: Path, energy, values, header: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    save_two_col(tmp, energy, values, header=header, comments="# ")
    os.replace(tmp, path)


def _render_group_qc_plot(
    path: Path,
    base_name: str,
    config: AstraConfig,
    energy,
    replicate_signals: list[np.ndarray],
    replicate_filenames: list[str],
    mu_avg,
    log=print,
) -> None:
    fig = None
    tmp = path.with_name(path.name + ".tmp.png")
    try:
        energy = np.asarray(energy, dtype=float)
        mu_avg = np.asarray(mu_avg, dtype=float)
        fig, ax = plt.subplots(figsize=(9, 5), dpi=PLOT_DPI)
        for signal, filename in zip(replicate_signals, replicate_filenames):
            signal = np.asarray(signal, dtype=float)
            mask = np.isfinite(energy) & np.isfinite(signal)
            if mask.any():
                ax.plot(energy[mask], signal[mask], linewidth=0.8, alpha=0.75, label=filename)
        avg_mask = np.isfinite(energy) & np.isfinite(mu_avg)
        if avg_mask.any():
            ax.plot(energy[avg_mask], mu_avg[avg_mask], color="black", linewidth=2.0, label="merged average")
        lo = getattr(config, "plot_energy_min", None)
        hi = getattr(config, "plot_energy_max", None)
        try:
            lo = float(lo)
            hi = float(hi)
            if np.isfinite(lo) and np.isfinite(hi) and lo < hi:
                ax.set_xlim(lo, hi)
        except (TypeError, ValueError):
            pass
        ax.set_title(f"Group QC: {base_name} ({len(replicate_signals)} replicates)")
        ax.set_xlabel("Energy / eV")
        ax.set_ylabel("Analysis signal μ(E)")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=7, loc="best")
        fig.tight_layout()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(tmp, dpi=PLOT_DPI, bbox_inches="tight")
        os.replace(tmp, path)
    except Exception as exc:
        log(f"WARNING: could not render group QC plot for {base_name}: {exc}")
        try:
            if tmp.exists():
                tmp.unlink()
        except OSError:
            pass
    finally:
        if fig is not None:
            plt.close(fig)


def _write_group_outputs(output_dir: Path, snapshot: dict, result: dict, generated: str, log=print) -> dict[str, str]:
    base_name = snapshot["base_name"]
    sanitized = sanitize_name(base_name)
    filenames = list(snapshot["accepted_filenames"])
    groups_dir = output_dir / "groups"
    plot_dir = output_dir / "plots" / "group_qc"
    processed_path = groups_dir / f"{sanitized}_processed.dat"
    norm_path = groups_dir / f"{sanitized}_norm.dat"
    flat_path = groups_dir / f"{sanitized}_flat.dat"
    plot_path = plot_dir / f"{sanitized}_replicate_qc.png"

    _write_two_col_atomic(
        processed_path,
        result["energy"],
        result["mu_avg"],
        _header(base_name, filenames, generated, "mu"),
    )
    _write_two_col_atomic(
        norm_path,
        result["energy"],
        result["norm"],
        _header(base_name, filenames, generated, "norm"),
    )
    _write_two_col_atomic(
        flat_path,
        result["energy"],
        result["flat"],
        _header(base_name, filenames, generated, "flat"),
    )
    _render_group_qc_plot(
        plot_path,
        base_name,
        snapshot["config"],
        result["energy"],
        result["replicate_signals"],
        result["replicate_filenames"],
        result["mu_avg"],
        log=log,
    )
    return {
        "processed": _relative(processed_path, output_dir),
        "norm": _relative(norm_path, output_dir),
        "flat": _relative(flat_path, output_dir),
        "qc_plot": _relative(plot_path, output_dir),
    }


def update_group_with_entry(
    entry: dict,
    path: Path,
    status: str,
    output_dir: Path,
    config: AstraConfig,
    registry: dict,
    registry_lock: threading.Lock,
    log=print,
) -> None:
    try:
        if status not in {"ok", "warn"}:
            return
        output_dir = Path(output_dir)
        abs_path = str(Path(path).expanduser().resolve())
        filename = str(entry.get("filename", Path(path).name))
        base_name, _ = split_replicate_suffix(filename)
        persist_state = None
        merge_snapshot = None

        with registry_lock:
            state = registry.get(base_name)
            if state is None:
                state = GroupState(base_name=base_name)
                registry[base_name] = state
            if abs_path in state.accepted_paths or filename in state.accepted_filenames:
                return
            state.accepted_paths.append(abs_path)
            state.accepted_filenames.append(filename)
            state.accepted_filenames, state.accepted_paths = _sort_pairs(
                state.accepted_filenames,
                state.accepted_paths,
            )
            state.n_accepted = len(state.accepted_filenames)
            state.last_updated_iso = _timestamp()
            if state.n_accepted < 2:
                state.last_merge_status = "pending"
                state.last_merge_error = ""
                persist_state = _state_to_dict(state)
            else:
                merge_snapshot = {
                    "base_name": state.base_name,
                    "accepted_paths": list(state.accepted_paths),
                    "accepted_filenames": list(state.accepted_filenames),
                    "n_accepted": state.n_accepted,
                    "last_updated_iso": state.last_updated_iso,
                    "config": config,
                }

        if persist_state is not None:
            _write_group_summary_atomic(_summary_path(output_dir, base_name), persist_state)
            return

        if merge_snapshot is None:
            return

        try:
            result = merge_and_normalize_group(merge_snapshot, config, log=log)
            generated = _timestamp()
            output_files = _write_group_outputs(output_dir, merge_snapshot, result, generated, log=log)
            with registry_lock:
                state = registry[base_name]
                state.last_merge_status = "ready"
                state.last_merge_error = ""
                state.output_files = output_files
                state.last_updated_iso = generated
                persist_state = _state_to_dict(state)
        except Exception as exc:
            log(f"WARNING: live group merge failed for {base_name}: {type(exc).__name__}: {exc}")
            with registry_lock:
                state = registry[base_name]
                state.last_merge_status = "error"
                state.last_merge_error = f"{type(exc).__name__}: {exc}"
                persist_state = _state_to_dict(state)

        if persist_state is not None:
            _write_group_summary_atomic(_summary_path(output_dir, base_name), persist_state)
    except Exception as exc:
        log(f"WARNING: live group update failed for {entry.get('filename', 'unknown')}: {type(exc).__name__}: {exc}")


def restore_group_registry(
    output_dir: Path,
    registry: dict,
    registry_lock: threading.Lock,
    log=print,
) -> None:
    try:
        output_dir = Path(output_dir)
        with registry_lock:
            registry.clear()
            for summary_path in sorted((output_dir / "groups").glob("*_group_summary.json")):
                try:
                    with summary_path.open("r", encoding="utf-8") as f:
                        data = json.load(f)
                    state = _state_from_dict(data)
                except Exception as exc:
                    log(f"WARNING: could not restore group summary {summary_path.name}: {exc}")
                    continue
                filenames = []
                paths = []
                changed = False
                for filename, path in zip(state.accepted_filenames, state.accepted_paths):
                    if Path(path).exists():
                        filenames.append(filename)
                        paths.append(path)
                    else:
                        changed = True
                        log(f"WARNING: restored group {state.base_name} dropped missing file {filename}")
                if changed:
                    filenames, paths = _sort_pairs(filenames, paths)
                    state.accepted_filenames = filenames
                    state.accepted_paths = paths
                    state.n_accepted = len(filenames)
                    if state.n_accepted < 2 and state.last_merge_status == "ready":
                        state.last_merge_status = "pending"
                    _write_group_summary_atomic(summary_path, state)
                registry[state.base_name] = state
    except Exception as exc:
        log(f"WARNING: could not restore beamtime group registry: {exc}")
