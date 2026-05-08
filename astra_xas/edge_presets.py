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


_LIGHT_TENDER_K = {
    "family": "light_tender_K",
    "pre1_rel": -80.0,
    "pre2_rel": -20.0,
    "norm1_rel": 80.0,
    "norm2_rel": 700.0,
    "align_min_rel": -20.0,
    "align_max_rel": 40.0,
    "plot_min_rel": -60.0,
    "plot_max_rel": 250.0,
}

_MEDIUM_K = {
    "family": "medium_K",
    "pre1_rel": -120.0,
    "pre2_rel": -30.0,
    "norm1_rel": 100.0,
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

_HEAVIER_K = {
    "family": "heavier_K",
    "pre1_rel": -200.0,
    "pre2_rel": -50.0,
    "norm1_rel": 150.0,
    "norm2_rel": 900.0,
    "align_min_rel": -20.0,
    "align_max_rel": 40.0,
    "plot_min_rel": -50.0,
    "plot_max_rel": 250.0,
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
    "Approximate starting values. Verify E0 and fit windows manually. "
    "Also verify alignment and plot windows for the dataset."
)

_PRESETS: tuple[EdgePreset, ...] = (
    EdgePreset(
        key="mg_k",
        label="Mg K",
        element="Mg",
        edge="K",
        e0_ref=1303.0,
        notes=_COMMON_NOTE,
        **_LIGHT_TENDER_K,
    ),
    EdgePreset(
        key="al_k",
        label="Al K",
        element="Al",
        edge="K",
        e0_ref=1559.6,
        notes=_COMMON_NOTE,
        **_LIGHT_TENDER_K,
    ),
    EdgePreset(
        key="si_k",
        label="Si K",
        element="Si",
        edge="K",
        e0_ref=1839.0,
        notes=_COMMON_NOTE,
        **_LIGHT_TENDER_K,
    ),
    EdgePreset(
        key="p_k",
        label="P K",
        element="P",
        edge="K",
        e0_ref=2145.5,
        notes=_COMMON_NOTE,
        **_LIGHT_TENDER_K,
    ),
    EdgePreset(
        key="s_k",
        label="S K",
        element="S",
        edge="K",
        e0_ref=2472.0,
        notes=_COMMON_NOTE,
        **_LIGHT_TENDER_K,
    ),
    EdgePreset(
        key="cl_k",
        label="Cl K",
        element="Cl",
        edge="K",
        e0_ref=2822.4,
        notes=_COMMON_NOTE,
        **_LIGHT_TENDER_K,
    ),
    EdgePreset(
        key="k_k",
        label="K K",
        element="K",
        edge="K",
        e0_ref=3608.4,
        notes=_COMMON_NOTE,
        **_MEDIUM_K,
    ),
    EdgePreset(
        key="ca_k",
        label="Ca K",
        element="Ca",
        edge="K",
        e0_ref=4038.5,
        notes=_COMMON_NOTE,
        **_MEDIUM_K,
    ),
    EdgePreset(
        key="ti_k",
        label="Ti K",
        element="Ti",
        edge="K",
        e0_ref=4966.0,
        notes=_COMMON_NOTE,
        **_TRANSITION_METAL_K,
    ),
    EdgePreset(
        key="v_k",
        label="V K",
        element="V",
        edge="K",
        e0_ref=5465.0,
        notes=_COMMON_NOTE,
        **_TRANSITION_METAL_K,
    ),
    EdgePreset(
        key="cr_k",
        label="Cr K",
        element="Cr",
        edge="K",
        e0_ref=5989.0,
        notes=_COMMON_NOTE,
        **_TRANSITION_METAL_K,
    ),
    EdgePreset(
        key="mn_k",
        label="Mn K",
        element="Mn",
        edge="K",
        e0_ref=6539.0,
        notes=_COMMON_NOTE,
        **_TRANSITION_METAL_K,
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
        key="cu_k",
        label="Cu K",
        element="Cu",
        edge="K",
        e0_ref=8979.0,
        notes=_COMMON_NOTE,
        **_TRANSITION_METAL_K,
    ),
    EdgePreset(
        key="zn_k",
        label="Zn K",
        element="Zn",
        edge="K",
        e0_ref=9659.0,
        notes=_COMMON_NOTE,
        **_TRANSITION_METAL_K,
    ),
    EdgePreset(
        key="ga_k",
        label="Ga K",
        element="Ga",
        edge="K",
        e0_ref=10367.0,
        notes=_COMMON_NOTE,
        **_HEAVIER_K,
    ),
    EdgePreset(
        key="ge_k",
        label="Ge K",
        element="Ge",
        edge="K",
        e0_ref=11103.0,
        notes=_COMMON_NOTE,
        **_HEAVIER_K,
    ),
    EdgePreset(
        key="as_k",
        label="As K",
        element="As",
        edge="K",
        e0_ref=11867.0,
        notes=_COMMON_NOTE,
        **_HEAVIER_K,
    ),
    EdgePreset(
        key="se_k",
        label="Se K",
        element="Se",
        edge="K",
        e0_ref=12658.0,
        notes=_COMMON_NOTE,
        **_HEAVIER_K,
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
