from __future__ import annotations
from collections import OrderedDict
from .io import natural_key


def group_samples(sample_entries: list[dict]) -> OrderedDict:
    groups = OrderedDict()
    for sample in sample_entries:
        group_key = (sample["base_name"], sample["assigned_foil"])
        groups.setdefault(group_key, []).append(sample)
    for _, scans in groups.items():
        scans.sort(key=lambda x: (
            x["replicate_id"] is None,
            x["replicate_id"] if x["replicate_id"] is not None else 999999,
            natural_key(x["filename"]),
        ))
    return groups
