from __future__ import annotations
import numpy as np
from scipy.interpolate import interp1d, CubicSpline
from scipy.optimize import minimize_scalar


def zscore(y):
    y = np.asarray(y, dtype=float)
    s = np.std(y)
    if s == 0:
        return y * 0
    return (y - np.mean(y)) / s


def find_best_shift(
    E_ref,
    mu_ref,
    E_mov,
    mu_mov,
    window=(7100, 7140),
    bounds=(-5, 5),
    grid_points=50,
):
    """Find shift to apply to moving spectrum: shifted_energy = E_mov + shift."""
    grid_points = max(int(grid_points), 5)

    E_ref = np.asarray(E_ref, dtype=float)
    mu_ref = np.asarray(mu_ref, dtype=float)
    mask_ref = (E_ref >= window[0]) & (E_ref <= window[1])
    if np.count_nonzero(mask_ref) < 10:
        raise ValueError("Reference alignment window has too few points")
    Er_raw = E_ref[mask_ref]
    mu_ref_raw = mu_ref[mask_ref]
    finite_ref = np.isfinite(Er_raw) & np.isfinite(mu_ref_raw)
    Er = Er_raw[finite_ref]
    mu_ref_window = mu_ref_raw[finite_ref]
    if len(Er) < 10:
        raise ValueError("Reference alignment window has too few finite points")
    sort_ref = np.argsort(Er)
    Er = Er[sort_ref]
    mu_ref_window = mu_ref_window[sort_ref]

    E_mov = np.asarray(E_mov, dtype=float)
    mu_mov = np.asarray(mu_mov, dtype=float)
    finite_mov = np.isfinite(E_mov) & np.isfinite(mu_mov)
    E_mov_clean = E_mov[finite_mov]
    mu_mov_clean = mu_mov[finite_mov]
    sort_idx = np.argsort(E_mov_clean)
    E_mov_clean = E_mov_clean[sort_idx]
    mu_mov_clean = mu_mov_clean[sort_idx]
    _, unique_idx = np.unique(E_mov_clean, return_index=True)
    unique_idx = np.sort(unique_idx)
    E_mov_clean = E_mov_clean[unique_idx]
    mu_mov_clean = mu_mov_clean[unique_idx]
    if len(E_mov_clean) < 10:
        return 0.0, np.nan, 0.0

    raw_dref = np.gradient(mu_ref_window, Er)
    signal_scale = np.nanmax(mu_ref_window) - np.nanmin(mu_ref_window)
    amplitude_range = np.nanmax(raw_dref) - np.nanmin(raw_dref)
    floor = max(signal_scale * 1e-6, 1e-15)
    if (
        not np.isfinite(amplitude_range)
        or not np.isfinite(signal_scale)
        or signal_scale < 1e-15
        or amplitude_range < floor
    ):
        return 0.0, np.nan, 0.0
    dref = zscore(raw_dref)
    cs = CubicSpline(E_mov_clean, mu_mov_clean, extrapolate=False)

    def objective(shift):
        y = cs(Er - shift)
        if not np.all(np.isfinite(y)):
            return 1e12
        dmov = zscore(np.gradient(y, Er))
        return np.mean((dref - dmov) ** 2)

    grid = np.linspace(bounds[0], bounds[1], grid_points)
    grid_scores = [objective(s) for s in grid]
    best_grid_shift = float(grid[np.argmin(grid_scores)])

    half_width = (bounds[1] - bounds[0]) / max(grid_points - 1, 1)
    lo = max(bounds[0], best_grid_shift - half_width)
    hi = min(bounds[1], best_grid_shift + half_width)
    try:
        res = minimize_scalar(objective, bounds=(lo, hi), method="bounded")
        optimal_shift = float(res.x)
        fit_residual = float(res.fun)
    except Exception:
        optimal_shift = best_grid_shift
        fit_residual = float(objective(best_grid_shift))
        return optimal_shift, fit_residual, 0.0

    y_opt = cs(Er - optimal_shift)
    if not np.all(np.isfinite(y_opt)):
        return float(optimal_shift), float(fit_residual), 0.0
    dmov = zscore(np.gradient(y_opt, Er))
    raw_quality = float(np.corrcoef(dref, dmov)[0, 1])
    if not np.isfinite(raw_quality):
        quality = 0.0
    else:
        quality = max(0.0, min(1.0, raw_quality))
    return float(optimal_shift), float(fit_residual), quality
