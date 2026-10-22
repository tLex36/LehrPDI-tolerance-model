"""Scan neuroligin C-tails for the beta-catenin binding motif D-x-θ-θ-x-φ-x(2-7)-E (and its
short core D-x-θ-θ-x-φ) and draw the hits on a multiple alignment.

    from ctail_motifs import scan_ctail_motifs, ctail_motif_html
    long_  = scan_ctail_motifs(seqs, topologies)                 # one row per hit
    short_ = scan_ctail_motifs(seqs, topologies, last=None)      # D-x-θ-θ-x-φ
    ctail_motif_html(aln, pd.concat([long_, short_]), topologies, ref="Human_NLGN1")

Motif positions (1-based within the match):
    1 D   2 x   3 θ   4 θ   5 x   6 φ   7..6+k x (k = 2-7)   7+k E      (short motif: 1-6)
Only the five defined positions (D, θ, θ, φ, E) are scored; a near match is a window in which
up to `max_mismatch` of them disagree.
"""
from html import escape
import pandas as pd
from msa_print import residue_columns, conservation_char

HYDROPHOBIC = "AVLIMFWYC"
AROMATIC    = "FWY"

def motif_positions(k, first="D", theta=HYDROPHOBIC, phi=HYDROPHOBIC, last="E"):
    """[(offset in window, allowed residues, label)] for a spacer of length k.
    last=None gives the short core motif D-x-θ-θ-x-φ (no spacer, no terminal E)."""
    core = [(0, set(first), "D"), (2, set(theta), "θ"), (3, set(theta), "θ"), (5, set(phi), "φ")]
    return core if last is None else core + [(6 + k, set(last), "E")]

def scan_motif(seq, start, end, name=None, spacer=(2, 7), max_mismatch=1, **motif):
    """All windows in seq[start..end] (1-based, inclusive) matching the motif with at most
    `max_mismatch` defined positions wrong.  Overlapping windows are collapsed to the best one
    (fewest mismatches, then shortest), so each hit is one distinct site."""
    rows = []
    short = motif.get("last", "E") is None
    for i in range(start - 1, end):
        for k in ([0] if short else range(spacer[0], spacer[1] + 1)):
            n = 6 if short else 7 + k
            if i + n > end:
                break
            w = seq[i:i + n]
            pos = motif_positions(k, **motif)
            bad = [(o, lab) for o, allowed, lab in pos if w[o] not in allowed]
            if len(bad) <= max_mismatch:
                rows.append({"protein": name, "start": i + 1, "end": i + n, "peptide": w,
                             "spacer": None if short else k, "n_mismatch": len(bad),
                             "mismatch": ", ".join(f"{lab}{o + 1}={w[o]}" for o, lab in bad),
                             "key_positions": [i + 1 + o for o, _, _ in pos],
                             "mismatch_positions": [i + 1 + o for o, _ in bad]})
    rows.sort(key=lambda r: (r["n_mismatch"], r["start"], r["end"]))
    kept = []
    for r in rows:                                   # greedy: best, then earliest, wins its span
        if all(r["end"] < q["start"] or r["start"] > q["end"] for q in kept):
            kept.append(r)
    return sorted(kept, key=lambda r: r["start"])

LONG_MOTIF  = "D-x-θ-θ-x-φ-x(2-7)-E"
SHORT_MOTIF = "D-x-θ-θ-x-φ"

def scan_ctail_motifs(seqs, topologies, relax_first="E", max_mismatch=0, label=None,
                      **motif) -> pd.DataFrame:
    """Scan each sequence from the residue after its (last) TM helix to the C-terminus.

    tier 'exact'          : the motif as written (D at position 1)
    tier 'E at position 1': matches only once `relax_first` residues are accepted at position 1
                            (D->E is the conservative acidic swap); position 1 is flagged red
    max_mismatch > 0 additionally reports windows with that many defined positions wrong
    (tier 'near'), for exploration -- most such windows are noise."""
    first = motif.pop("first", "D")
    rows = []
    for name, seq in seqs.items():
        tm_end = max(e for _, e in topologies[name]["transmembrane"])
        a, b = tm_end + 1, len(seq)
        strict = scan_motif(seq, a, b, name=name, max_mismatch=max_mismatch, first=first, **motif)
        for r in strict:
            r["tier"] = "exact" if r["n_mismatch"] == 0 else "near"
        found = list(strict)
        if relax_first:
            for r in scan_motif(seq, a, b, name=name, max_mismatch=0,
                                first=first + relax_first, **motif):
                if seq[r["start"] - 1] in first:
                    continue                         # already an exact hit
                if any(not (r["end"] < q["start"] or r["start"] > q["end"]) and q["tier"] == "exact"
                       for q in found):
                    continue
                found = [q for q in found if q["tier"] == "exact"
                         or r["end"] < q["start"] or r["start"] > q["end"]]
                r.update(tier=f"{relax_first} at position 1", n_mismatch=1,
                         mismatch=f"D1={seq[r['start'] - 1]}", mismatch_positions=[r["start"]])
                found.append(r)
        rows += sorted(found, key=lambda r: r["start"])
    label = label or (SHORT_MOTIF if motif.get("last", "E") is None else LONG_MOTIF)
    for r in rows:
        r["motif"] = label
    cols = ["motif", "protein", "start", "end", "peptide", "spacer", "tier", "mismatch",
            "n_mismatch", "key_positions", "mismatch_positions"]
    return pd.DataFrame(rows, columns=cols)

# colours
C_EXACT   = "#b7e1b0"      # whole span, exact match
C_NEAR    = "#fde9a9"      # whole span, near match
C_KEY     = "#1b5e20"      # D/θ/θ/φ/E residues that fit
C_MISS    = "#c00000"      # defined position that does not fit
C_PDZ     = "#d9c8ee"      # C-terminal PDZ-binding motif
U_EXACT   = "#08519c"      # short-motif underline, exact match
U_NEAR    = "#6baed6"      # short-motif underline, near match

def _swatch(bg=None, underline=None):
    css = f"background:{bg};" if bg else f"box-shadow: inset 0 -3px 0 {underline};"
    return f'<span style="{css}padding:1px 7px">&nbsp;</span>'

def ctail_motif_html(aln, hits, topologies, ref=None, width=70, pdz_len=4, title=None,
                     font="Courier New", fontsize=11, near_label="E in place of D at position 1"):
    """Alignment of the C-tails (residue after the reference TM to the end of the alignment).

    Long-motif hits (with spacer and terminal E) are shaded: green = exact, yellow = near.
    Short-motif hits (D-x-θ-θ-x-φ, spacer is None) are underlined: dark blue = exact,
    light blue = near.  Both can sit on the same residues.  Defined motif positions are bold,
    green if they fit and red if they do not; the last `pdz_len` residues (PDZ ligand) are
    boxed in lavender."""
    from IPython.display import HTML
    names = [r.id for r in aln]
    ref = ref or names[0]
    cols = {r.id: residue_columns(str(r.seq)) for r in aln}
    is_short = hits["spacer"].isna() if len(hits) else pd.Series(dtype=bool)

    style = {}                                       # (name, column) -> {css}
    for h, short in zip(hits.itertuples(), is_short):
        idx = cols[h.protein]
        for p in range(h.start, h.end + 1):
            s = style.setdefault((h.protein, idx[p - 1]), {})
            if short:
                s["box-shadow"] = f"inset 0 -3px 0 {U_EXACT if h.n_mismatch == 0 else U_NEAR}"
            else:
                s["background"] = C_EXACT if h.n_mismatch == 0 else C_NEAR
        for p in h.key_positions:
            s = style[(h.protein, idx[p - 1])]
            s["font-weight"] = "700"
            s["color"] = C_MISS if p in h.mismatch_positions else C_KEY
    for n, idx in cols.items():                      # PDZ ligand at the very C-terminus
        for p in range(len(idx) - pdz_len + 1, len(idx) + 1):
            style.setdefault((n, idx[p - 1]), {}).setdefault("background", C_PDZ)

    ref_tm_end = max(e for _, e in topologies[ref]["transmembrane"])
    L = aln.get_alignment_length()
    lo, hi = cols[ref][ref_tm_end], L                # first C-tail column of the reference
    pad = max(len(n) for n in names) + 2
    cons = "".join(conservation_char(aln[:, i]) for i in range(L))
    count = {r.id: sum(1 for c in str(r.seq)[:lo] if c != "-") for r in aln}

    lines = []
    for start in range(lo, hi, width):
        stop = min(start + width, hi)
        for rec in aln:
            chunk = str(rec.seq[start:stop])
            count[rec.id] += len(chunk.replace("-", ""))
            cells = []
            for i, ch in enumerate(chunk):
                s = style.get((rec.id, start + i))
                css = ";".join(f"{k}:{v}" for k, v in s.items()) if s else ""
                cells.append(f'<span style="{css}">{escape(ch)}</span>' if css else escape(ch))
            cells.append(" " * (width - len(chunk)))
            lines.append(f"{escape(rec.id + ':'):<{pad}}{''.join(cells)} {count[rec.id]}")
        lines.append(f"{'':<{pad}}{escape(cons[start:stop])}")
        lines.append("")

    rows = []
    for short, lab, ex, nr in [(False, LONG_MOTIF, dict(bg=C_EXACT), dict(bg=C_NEAR)),
                               (True, SHORT_MOTIF, dict(underline=U_EXACT), dict(underline=U_NEAR))]:
        sub = hits[is_short == short] if len(hits) else hits
        if len(sub):
            lab = sub["motif"].iloc[0]
            rows.append(f'<b>{escape(lab)}</b>&nbsp; {_swatch(**ex)} exact &nbsp; '
                        f'{_swatch(**nr)} one position off ({escape(near_label)})')
    rows.append(f'<b style="color:{C_KEY}">X</b> defined position fits &nbsp; '
                f'<b style="color:{C_MISS}">X</b> defined position does not fit &nbsp; '
                f'{_swatch(bg=C_PDZ)} C-terminal PDZ ligand')
    key = ('<div style="font-family:sans-serif;font-size:10pt;margin:8px 0;line-height:1.9">'
           + "<br>".join(rows) + '</div>')
    head = (f'<div style="font-family:sans-serif;font-size:12pt;font-weight:600">{escape(title)}</div>'
            if title else "")
    body = "\n".join(lines)
    return HTML(head + key + f'<pre style="font-family:\'{font}\',monospace;font-size:{fontsize}pt;'
                             f'line-height:1.6;white-space:pre">{body}</pre>')
