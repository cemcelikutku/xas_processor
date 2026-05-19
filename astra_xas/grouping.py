from __future__ import annotations
from collections import OrderedDict
from .io import natural_key


def _manifest_sort_key(scan: dict):
    return scan["_order_in_manifest"]


def _folder_sort_key(scan: dict):
    rid = scan["replicate_id"]
    return (
        rid is None,
        rid if rid is not None else 999999,
        natural_key(scan["filename"]),
    )


def group_samples(sample_entries: list[dict]) -> OrderedDict:
    groups = OrderedDict()
    for sample in sample_entries:
        group_key = (sample["base_name"], sample["assigned_foil"])
        groups.setdefault(group_key, []).append(sample)
    for _, scans in groups.items():
        # When entries were produced from a manifest, their explicit
        # position drives the replicate order. Folder mode falls back to
        # the replicate-suffix + natural-key heuristic.
        if scans and "_order_in_manifest" in scans[0]:
            scans.sort(key=_manifest_sort_key)
        else:
            scans.sort(key=_folder_sort_key)
    return groups
