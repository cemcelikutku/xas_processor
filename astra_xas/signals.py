from __future__ import annotations
import numpy as np
from .config import AstraConfig


def safe_log_ratio(a, b, eps: float = 1e-30):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return np.log(np.clip(a, eps, None) / np.clip(b, eps, None))


def compute_signals(scan: dict, config: AstraConfig) -> dict:
    I0, I1, I2, IF = scan["I0"], scan["I1"], scan["I2"], scan["IF"]
    return {
        "energy": scan["energy"],
        "mu_trans": safe_log_ratio(I0, I1),
        "mu_ref": safe_log_ratio(I1, I2),
        "mu_fluo": config.fluo_multiplicative_constant * IF / np.clip(I0, 1e-30, None),
    }


def get_signal(entry: dict, mode: str):
    if mode == "trans":
        return entry["mu_trans"]
    if mode == "ref":
        return entry["mu_ref"]
    if mode == "fluo":
        return entry["mu_fluo"]
    raise ValueError(f"Unknown signal mode: {mode}")
