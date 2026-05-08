from __future__ import annotations

from dataclasses import dataclass

from .config import AstraConfig


@dataclass(frozen=True)
class EdgePreset:
    key: str
    label: str
    element: str
    edge: str
    e0_ref: float
    family: str
    pre1_rel: float
    pre2_rel: float
    norm1_rel: float
    norm2_rel: float
    align_min_rel: float
    align_max_rel: float
    plot_min_rel: float
    plot_max_rel: float
    notes: str = ""


_TENDER_K = {
    "family": "tender_K",
    "pre1_rel": -80.0,
    "pre2_rel": -20.0,
    "norm1_rel": 80.0,
    "norm2_rel": 700.0,
    "align_min_rel": -20.0,
    "align_max_rel": 40.0,
    "plot_min_rel": -60.0,
    "plot_max_rel": 250.0,
}

_TRANSITION_METAL_K = {
    "family": "transition_metal_K",
    "pre1_rel": -200.0,
    "pre2_rel": -50.0,
    "norm1_rel": 150.0,
    "norm2_rel": 800.0,
    "align_min_rel": -20.0,
    "align_max_rel": 40.0,
    "plot_min_rel": -50.0,
    "plot_max_rel": 200.0,
}

_L3_EDGE = {
    "family": "l3_edge",
    "pre1_rel": -150.0,
    "pre2_rel": -40.0,
    "norm1_rel": 100.0,
    "norm2_rel": 600.0,
    "align_min_rel": -20.0,
    "align_max_rel": 40.0,
    "plot_min_rel": -50.0,
    "plot_max_rel": 180.0,
}

_COMMON_NOTE = (
    "Approximate reference/start E0 and editable fit-window template; "
    "verify E0, fit windows, alignment window, and plot range manually."
)

_PRESETS: tuple[EdgePreset, ...] = (
    EdgePreset(
        key="p_k",
        label="P K",
        element="P",
        edge="K",
        e0_ref=2145.5,
        notes=_COMMON_NOTE,
        **_TENDER_K,
    ),
    EdgePreset(
        key="s_k",
        label="S K",
        element="S",
        edge="K",
        e0_ref=2472.0,
        notes=_COMMON_NOTE,
        **_TENDER_K,
    ),
    EdgePreset(
        key="fe_k",
        label="Fe K",
        element="Fe",
        edge="K",
        e0_ref=7112.0,
        notes=_COMMON_NOTE,
        **_TRANSITION_METAL_K,
    ),
    EdgePreset(
        key="co_k",
        label="Co K",
        element="Co",
        edge="K",
        e0_ref=7709.0,
        notes=_COMMON_NOTE,
        **_TRANSITION_METAL_K,
    ),
    EdgePreset(
        key="ni_k",
        label="Ni K",
        element="Ni",
        edge="K",
        e0_ref=8333.0,
        notes=_COMMON_NOTE,
        **_TRANSITION_METAL_K,
    ),
    EdgePreset(
        key="ce_l3",
        label="Ce L3",
        element="Ce",
        edge="L3",
        e0_ref=5723.0,
        notes=_COMMON_NOTE,
        **_L3_EDGE,
    ),
)

_PRESET_MAP = {preset.key: preset for preset in _PRESETS}


def list_edge_presets() -> tuple[EdgePreset, ...]:
    return _PRESETS


def get_edge_preset(key: str) -> EdgePreset | None:
    if key == "custom":
        return None
    return _PRESET_MAP.get(key)


def apply_edge_preset_to_config(config: AstraConfig, key: str) -> AstraConfig:
    preset = get_edge_preset(key)
    if preset is None:
        config.edge_preset_key = "custom"
        config.edge_preset_label = "Custom"
        config.edge_preset_applied = False
        config.edge_preset_note = ""
        return config

    config.e0 = preset.e0_ref
    config.pre1 = preset.pre1_rel
    config.pre2 = preset.pre2_rel
    config.norm1 = preset.norm1_rel
    config.norm2 = preset.norm2_rel
    config.align_window_min = preset.e0_ref + preset.align_min_rel
    config.align_window_max = preset.e0_ref + preset.align_max_rel
    config.plot_energy_min = preset.e0_ref + preset.plot_min_rel
    config.plot_energy_max = preset.e0_ref + preset.plot_max_rel
    config.edge_preset_key = preset.key
    config.edge_preset_label = preset.label
    config.edge_preset_applied = True
    config.edge_preset_note = preset.notes
    return config

