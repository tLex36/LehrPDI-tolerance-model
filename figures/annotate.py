"""General domain/region annotation for protein alignments.

    from annotate import receptor_domain_tracks, annotated_alignment_html
    tracks = receptor_domain_tracks({"GRIN2A": "Q12879", "GRIN2B": "Q13224"})
    annotated_alignment_html(aln, tracks)

Regions are supplied as {sequence: {label: (start, end)}} in each sequence's own 1-based
numbering, so any set of features can be added -- domains, motifs, custom annotations.
"""
from html import escape
import pandas as pd
from msa_print import residue_columns, conservation_char

# categorical palette, assigned in fixed order (never cycled); pale enough for black text
PALETTE = ["#cde2fb", "#ffe0b3", "#cfe9d6", "#e6d7f0", "#f9d3d3",
           "#d4ecec", "#ede3c8", "#dcdcdc", "#f6d9ea", "#d9e4c8"]

# sensible defaults for an ionotropic glutamate receptor
IGLUR_COLORS = {"signal": "#dcdcdc", "ATD": "#cde2fb", "S1": "#cfe9d6", "S2": "#cfe9d6",
                "M1": "#ffb38a", "M2": "#ff8a5c", "M3": "#ffb38a", "M4": "#ffb38a",
                "CTD": "#e6d7f0", "LBD": "#cfe9d6", "TM": "#ffb38a"}

def receptor_domain_tracks(accessions, merge_lbd=False, drop=("signal",)) -> dict:
    """{name: {domain: (start, end)}} from conservation.receptor_domains() for each accession.
    merge_lbd=True reports S1+S2 as a single 'LBD' span (discontinuous, shown as its extent)."""
    from conservation import receptor_domains
    out = {}
    for name, acc in accessions.items():
        d = receptor_domains(acc)
        d = {k: v for k, v in d.items() if k not in (drop or ())}
        if merge_lbd and "S1" in d and "S2" in d:
            d = {k: v for k, v in d.items() if k not in ("S1", "S2")}
            d["LBD"] = None                      # placeholder; S1/S2 kept separately below
        out[name] = d
    return out

def assign_colors(tracks, colors=None, palette=PALETTE) -> dict:
    """{label: colour} -- explicit `colors` win; anything unnamed takes the next palette slot."""
    labels = []
    for d in tracks.values():
        for k in d:
            if k not in labels: labels.append(k)
    out = dict(colors or {})
    nxt = 0
    for k in labels:
        if k in out: continue
        out[k] = palette[nxt % len(palette)]; nxt += 1
    return out

def region_columns(gapped_seq, regions):
    """{alignment column: label} for one sequence's regions (later entries win on overlap)."""
    idx = residue_columns(gapped_seq)
    out = {}
    for label, span in regions.items():
        if span is None: continue
        s, e = span
        for p in range(int(s), int(e) + 1):
            if 0 < p <= len(idx):
                out[idx[p - 1]] = label
    return out

def annotated_alignment_html(aln, tracks, colors=None, palette=PALETTE, width=60, region=None,
                             offsets=None, label_row=True, extra=None, font="Courier New",
                             fontsize=11, title=None, show_conservation=True, legend=True):
    """Alignment with named regions highlighted.

    tracks:  {sequence: {label: (start, end)}} in each sequence's own numbering
    colors:  {label: colour}; unspecified labels get palette colours in order
    extra:   {(sequence, column): style string} for any additional per-residue styling
    label_row: write the region name above the first residue of each span
    region:  (start, end) window in the FIRST sequence's numbering
    """
    from IPython.display import HTML
    colors = assign_colors(tracks, colors, palette)
    cols = {r.id: region_columns(str(r.seq), tracks.get(r.id, {})) for r in aln}
    offsets = offsets or {}
    L = aln.get_alignment_length()
    lo, hi = 0, L
    if region:
        idx = residue_columns(str(aln[0].seq))
        lo, hi = idx[region[0] - 1], idx[min(region[1], len(idx)) - 1] + 1
    names = [r.id for r in aln]
    pad = max(len(n) for n in names) + 2
    cons = "".join(conservation_char(aln[:, i]) for i in range(L))
    count = {r.id: offsets.get(r.id, 0) + sum(1 for c in str(r.seq)[:lo] if c != "-") for r in aln}

    def label_line(mapping, start, stop):
        row = [" "] * (stop - start)
        prev = None
        for c in range(start, stop):
            lab = mapping.get(c)
            if lab and lab != prev:                        # first column of a span
                i = c - start
                if all(row[j] == " " for j in range(i, min(i + len(lab) + 1, len(row)))):
                    for j, ch in enumerate(lab):
                        if i + j < len(row): row[i + j] = ch
            prev = lab
        return "".join(row)

    lines = []
    for start in range(lo, hi, width):
        stop = min(start + width, hi)
        for rec in aln:
            m = cols.get(rec.id, {})
            if label_row and m:
                lr = label_line(m, start, stop)
                if lr.strip():
                    lines.append(f'{"":<{pad}}<span style="color:#555;font-size:90%">{escape(lr)}</span>')
            chunk = str(rec.seq[start:stop])
            count[rec.id] += len(chunk.replace("-", ""))
            cells = []
            for i, ch in enumerate(chunk):
                col = start + i
                style = []
                lab = m.get(col)
                if lab: style.append(f"background:{colors[lab]}")
                ex = (extra or {}).get((rec.id, col))
                if ex: style.append(ex)
                cells.append(f'<span style="{";".join(style)}">{escape(ch)}</span>' if style else escape(ch))
            cells.append(" " * (width - len(chunk)))
            lines.append(f"{escape(rec.id + ':'):<{pad}}{''.join(cells)} {count[rec.id]}")
        if show_conservation:
            lines.append(f"{'':<{pad}}{escape(cons[start:stop])}")
        lines.append("")

    key = ""
    if legend:
        used = [k for k in colors if any(k in d for d in tracks.values())]
        key = ('<div style="font-family:sans-serif;font-size:10pt;margin:8px 0">'
               + " &nbsp; ".join(f'<span style="background:{colors[k]};padding:1px 8px">&nbsp;</span> {escape(k)}'
                                 for k in used) + "</div>")
    head = (f'<div style="font-family:sans-serif;font-size:12pt;font-weight:600">{escape(title)}</div>'
            if title else "")
    return HTML(head + key + f'<pre style="font-family:\'{font}\',monospace;font-size:{fontsize}pt;'
                             f'line-height:1.5;white-space:pre">{chr(10).join(lines)}</pre>')

def domain_table(tracks) -> pd.DataFrame:
    """Tidy view of the regions being drawn."""
    rows = [{"sequence": n, "region": k, "start": v[0], "end": v[1], "length": v[1] - v[0] + 1}
            for n, d in tracks.items() for k, v in d.items() if v]
    return pd.DataFrame(rows)
