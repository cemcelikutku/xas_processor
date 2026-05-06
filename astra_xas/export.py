from __future__ import annotations
from pathlib import Path
import numpy as np


def save_two_col(path: str | Path, x, y, header: str, comments: str = ""):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, np.column_stack([x, y]), header=header, fmt="%.8e", comments=comments)
