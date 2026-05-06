from __future__ import annotations
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar


def zscore(y):
    y = np.asarray(y, dtype=float)
    s = np.std(y)
    if s == 0:
        return y * 0
    return (y - np.mean(y)) / s


def find_best_shift(E_ref, mu_ref, E_mov, mu_mov, window=(7100, 7140), bounds=(-5, 5)):
    """Find shift to apply to moving spectrum: shifted_energy = E_mov + shift."""
    mask_ref = (E_ref >= window[0]) & (E_ref <= window[1])
    if np.count_nonzero(mask_ref) < 10:
        raise ValueError("Reference alignment window has too few points")
    Er = E_ref[mask_ref]
    dref = zscore(np.gradient(mu_ref[mask_ref], Er))
    mov_interp = interp1d(E_mov, mu_mov, kind="linear", bounds_error=False, fill_value=np.nan)

    def objective(shift):
        y = mov_interp(Er - shift)
        if np.isnan(y).any():
            return 1e12
        dmov = zscore(np.gradient(y, Er))
        return np.mean((dref - dmov) ** 2)

    res = minimize_scalar(objective, bounds=bounds, method="bounded")
    return float(res.x), float(res.fun)
