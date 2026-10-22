"""Publication-style HTML tables for pasting into Word / PowerPoint.

    from tables import styled_table
    styled_table(density, inferred, caption="gnomAD missense per residue")

Render in a notebook cell, then select the output and copy -- the clipboard carries HTML, so
Word and PowerPoint keep the fonts, rules and cell shading.
"""
import numpy as np, pandas as pd
import matplotlib as mpl
from IPython.display import HTML

# single-hue sequential ramps (light -> dark), used for cell shading
RAMPS = {
    "blue":   ["#eef4fc", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95"],
    "purple": ["#f3eef8", "#e0d3ee", "#c4addd", "#a486c9", "#8560b4", "#6a4396", "#512d76"],
    "green":  ["#eef7f0", "#cfe9d6", "#a7d7b4", "#7bc292", "#4fa972", "#348a58", "#226b42"],
    "warm":   ["#fdf1ec", "#fadbcd", "#f5bda5", "#ec9a78", "#dd744e", "#c05531", "#983f22"],
}

def _shade(v, vmin, vmax, ramp):
    if pd.isna(v): return "transparent", "#999"
    cmap = mpl.colors.LinearSegmentedColormap.from_list("r", ramp)
    t = 0 if vmax <= vmin else (v - vmin) / (vmax - vmin)
    r, g, b, _ = cmap(t)
    ink = "#ffffff" if (0.299 * r + 0.587 * g + 0.114 * b) < 0.55 else "#1a1a1a"
    return mpl.colors.to_hex((r, g, b)), ink

def styled_table(df, inferred=None, caption=None, note=None, ramp="blue", shade=True,
                 decimals=2, vmin=None, vmax=None, index_name=None, font="Helvetica Neue",
                 fontsize=13, na="—"):
    """df: numeric DataFrame (rows = genes, cols = regions).
    inferred: same-shaped boolean frame; True cells get an asterisk and italic type.
    shade: colour cells by value (one hue, light->dark). Set False for a plain journal table.
    Returns an IPython HTML object; use .data for the raw string."""
    vals = df.to_numpy(dtype=float)
    vmin = np.nanmin(vals) if vmin is None else vmin
    vmax = np.nanmax(vals) if vmax is None else vmax
    ramp = RAMPS.get(ramp, ramp if isinstance(ramp, list) else RAMPS["blue"])
    idx_name = index_name if index_name is not None else (df.index.name or "")

    css = f"""
    <style>
    .pubtbl {{ border-collapse: collapse; font-family: "{font}", Helvetica, Arial, sans-serif;
               font-size: {fontsize}px; color: #1a1a1a; margin: 6px 0 4px 0; }}
    .pubtbl caption {{ caption-side: top; text-align: left; font-weight: 600;
                       font-size: {fontsize + 1}px; padding: 0 0 8px 2px; color: #111; }}
    .pubtbl th {{ font-weight: 600; padding: 7px 14px; text-align: right; white-space: nowrap;
                  border-bottom: 1.5px solid #333; }}
    .pubtbl th.rowhdr, .pubtbl td.rowhdr {{ text-align: left; font-weight: 600;
                  border-right: 1px solid #e3e3e3; }}
    .pubtbl td {{ padding: 7px 14px; text-align: right; font-variant-numeric: tabular-nums;
                  border-bottom: 1px solid #eee; }}
    .pubtbl tr:last-child td {{ border-bottom: 1.5px solid #333; }}
    .pubtbl .inf {{ font-style: italic; }}
    .pubnote {{ font-family: "{font}", Helvetica, Arial, sans-serif; font-size: {fontsize - 2}px;
                color: #666; margin: 2px 0 14px 2px; max-width: 640px; }}
    </style>"""

    head = "".join(f"<th>{c}</th>" for c in df.columns)
    rows = []
    for r in df.index:
        cells = []
        for c in df.columns:
            v = df.loc[r, c]
            is_inf = bool(inferred.loc[r, c]) if inferred is not None else False
            if pd.isna(v):
                cells.append(f'<td style="color:#bbb">{na}</td>'); continue
            bg, ink = _shade(v, vmin, vmax, ramp) if shade else ("transparent", "#1a1a1a")
            star = "*" if is_inf else ""
            cls = ' class="inf"' if is_inf else ""
            cells.append(f'<td{cls} style="background:{bg};color:{ink}">{v:.{decimals}f}{star}</td>')
        rows.append(f'<tr><td class="rowhdr">{r}</td>{"".join(cells)}</tr>')

    cap = f"<caption>{caption}</caption>" if caption else ""
    html = (css + f'<table class="pubtbl">{cap}<thead><tr>'
            f'<th class="rowhdr">{idx_name}</th>{head}</tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table>')
    if note:
        html += f'<div class="pubnote">{note}</div>'
    return HTML(html)

def save_table_html(html_obj, path, title=None):
    """Write a standalone .html file (open in a browser, then copy into Word/PowerPoint)."""
    body = html_obj.data if hasattr(html_obj, "data") else str(html_obj)
    doc = (f"<!doctype html><meta charset='utf-8'>"
           f"<title>{title or 'table'}</title>"
           f"<div style='padding:28px'>{body}</div>")
    open(path, "w").write(doc)
    return path
