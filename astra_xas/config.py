from dataclasses import dataclass
from pathlib import Path

@dataclass
class AstraConfig:
    version: str = "ASTRA XAS Processor 0.1.0"
    foil_keyword: str = "foil"
    analysis_mode: str = "fluo"      # fluo, trans, ref
    alignment_source: str = "separate_foil"  # separate_foil, inline_ref
    foil_alignment_mode: str = "trans"  # trans, ref, fluo
    fluo_multiplicative_constant: float = 1e-11

    pre1: float = -229.740
    pre2: float = -49.980
    norm1: float = 55.070
    norm2: float = 227.220
    nnorm: int = 1
    e0: float = 7121.030
    step: float | None = None
    nvict: int = 0
    make_flat: bool = True

    align_window_min: float = 7100.0
    align_window_max: float = 7140.0
    shift_bound_min: float = -5.0
    shift_bound_max: float = 5.0
    interp_kind: str = "linear"
    warn_shift_abs_eV: float = 2.0

    exclude_filenames: tuple[str, ...] | str = ()
    exclude_filename_contains: tuple[str, ...] | str = ()
    enable_shift_rejection: bool = False
    reject_shift_abs_eV: float = 3.0
    enable_auto_outlier_detection: bool = False
    outlier_rms_threshold: float = 0.08

    # Automatic plots. Avoid using "raw" for processed μ(E); true raw detector channels are separate.
    save_detector_raw_overview_plot: bool = False
    save_processed_overview_plot: bool = True
    save_bkgcorr_overview_plot: bool = False
    save_norm_overview_plot: bool = True
    save_replicate_qc_plots: bool = True
    save_raw_overview_plot: bool = False  # legacy alias only; do not use in new GUI
    save_drift_plot: bool = False
    save_foil_alignment_plots: bool = False

    plot_energy_min: float = 7100.0
    plot_energy_max: float = 7160.0

    @property
    def align_window(self) -> tuple[float, float]:
        return (self.align_window_min, self.align_window_max)

    @property
    def shift_bounds(self) -> tuple[float, float]:
        return (self.shift_bound_min, self.shift_bound_max)
