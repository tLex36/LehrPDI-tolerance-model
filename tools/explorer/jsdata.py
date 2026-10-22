"""Shared writer for the explorer's data files: JSON wrapped in a function call, e.g. TOL.gene({...});"""
import json
import numpy as np


def clean(o):
    """NaN / inf -> None, recursively, so the JSON is valid."""
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def write_js(path, call, obj):
    """Write obj as compact JSON inside call(...)."""
    text = json.dumps(clean(obj), separators=(",", ":"), allow_nan=False)
    with open(path, "w") as f:
        f.write(f"{call}({text});\n")


def enc(a, scale):
    """Per-residue values -> one compact string of integers: round(value * scale), comma-separated, empty when
    missing. The page decodes it with value = integer / scale (see SCALES in docs/assets/data.js)."""
    a = np.asarray(a, float)
    return ",".join("" if not np.isfinite(x) else str(int(round(x * scale))) for x in a)
