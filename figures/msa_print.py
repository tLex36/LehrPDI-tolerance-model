import subprocess, tempfile, os
from Bio import AlignIO

# Clustal's residue similarity groups
STRONG = ["STA","NEQK","NHQK","NDEQ","QHRK","MILV","MILF","HY","FYW"]
WEAK   = ["CSA","ATV","SAG","STNK","STPA","SGND","SNDEQK","NDEQHK","NEQHRK","FVLIM","HFY"]

def conservation_char(col):
    res = set(col) - {"-"}
    if "-" in col or not res:      # any gap -> blank
        return " "
    if len(res) == 1:
        return "*"
    if any(res <= set(g) for g in STRONG):
        return ":"
    if any(res <= set(g) for g in WEAK):
        return "."
    return " "

def align(seqs: dict, order=None):
    """seqs: {name: sequence}. Returns Biopython MultipleSeqAlignment in `order`."""
    order = order or list(seqs)
    with tempfile.TemporaryDirectory() as d:
        fa, out = os.path.join(d, "in.fa"), os.path.join(d, "out.fa")
        # FASTA ids end at the first space, so write safe numeric ids and map back after
        names = list(seqs)
        with open(fa, "w") as f:
            for i, n in enumerate(names):
                f.write(f">s{i}\n{seqs[n]}\n")
        subprocess.run(["muscle", "-align", fa, "-output", out],
                       check=True, capture_output=True)
        aln = AlignIO.read(out, "fasta")
    recs = {}
    for r in aln:
        r.id = r.name = r.description = names[int(r.id[1:])]
        recs[r.id] = r
    aln = AlignIO.MultipleSeqAlignment([recs[n] for n in order])
    return aln

def format_clustal(aln, width=60, offsets=None):
    """Clustal-style text: wrapped blocks, residue numbers, conservation line.
    offsets: {name: first_residue_number - 1} so numbering matches the full protein."""
    offsets = offsets or {}
    names = [r.id for r in aln]
    pad = max(len(n) for n in names) + 2
    L = aln.get_alignment_length()
    cons = "".join(conservation_char(aln[:, i]) for i in range(L))
    count = {n: offsets.get(n, 0) for n in names}
    lines = []
    for start in range(0, L, width):
        for r in aln:
            chunk = str(r.seq[start:start+width])
            count[r.id] += len(chunk.replace("-", ""))
            lines.append(f"{r.id + ':':<{pad}}{chunk:<{width}} {count[r.id]}")
        lines.append(f"{'':<{pad}}{cons[start:start+width]}")
        lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------- motif highlighting
import re as _re
from html import escape as _escape

def residue_columns(gapped_seq):
    """[alignment column of residue 0, residue 1, ...] for one aligned sequence."""
    return [c for c, ch in enumerate(gapped_seq) if ch != "-"]

def find_motif_columns(gapped_seq, motif, regex=False, ignore_case=False):
    """Alignment columns covered by every occurrence of `motif` in this sequence.

    The search runs on the UNGAPPED sequence, so a motif split by alignment gaps is still
    found; the returned column set spans the gaps it straddles.
    Returns [set(columns), ...] -- one set per occurrence.
    """
    cols = residue_columns(gapped_seq)
    plain = "".join(ch for ch in gapped_seq if ch != "-")
    pat = motif if regex else _re.escape(motif)
    flags = _re.IGNORECASE if ignore_case else 0
    out = []
    for m in _re.finditer(pat, plain, flags):
        out.append({cols[i] for i in range(m.start(), m.end())})
    return out

def motif_column_map(aln, motifs, regex=False, ignore_case=False, sequences=None):
    """{(seq_id, column): colour} for every motif occurrence.

    motifs: {pattern: colour} -- searched in every sequence, or restrict with
            `sequences=[ids]`. Earlier keys win where occurrences overlap.
    """
    colored = {}
    for pattern, color in motifs.items():
        for rec in aln:
            if sequences and rec.id not in sequences:
                continue
            for occ in find_motif_columns(str(rec.seq), pattern, regex, ignore_case):
                for c in occ:
                    colored.setdefault((rec.id, c), color)
    return colored

def format_clustal_html(aln, motifs=None, width=60, offsets=None, regex=False, ignore_case=False,
                        sequences=None, font="Courier New", fontsize=11, show_conservation=None,
                        extra=None):
    """Clustal-style alignment as HTML, with motif occurrences highlighted.

    motifs:  {pattern: colour}, e.g. {"TTRV": "#ffe066", "PDYT": "#a8d8ff"}
    regex:   treat the patterns as regular expressions (for degenerate motifs like "S[TS].V")
    extra:   {(seq_id, residue_number): colour} for individual residues, in the same numbering
             as `offsets` (i.e. real protein positions, not alignment columns)
    Returns an IPython HTML object; render it, select the output and copy into Word/PowerPoint.
    """
    from IPython.display import HTML
    offsets = offsets or {}
    colored = motif_column_map(aln, motifs or {}, regex, ignore_case, sequences)
    if extra:                                          # residue numbers -> alignment columns
        for rec in aln:
            cols = residue_columns(str(rec.seq))
            off = offsets.get(rec.id, 0)
            for (sid, pos), color in extra.items():
                if sid != rec.id: continue
                i = pos - off - 1
                if 0 <= i < len(cols):
                    colored[(rec.id, cols[i])] = color

    names = [r.id for r in aln]
    pad = max(len(n) for n in names) + 2
    L = aln.get_alignment_length()
    if show_conservation is None:            # meaningless for a single sequence
        show_conservation = len(aln) > 1
    cons = "".join(conservation_char(aln[:, i]) for i in range(L))
    count = {n: offsets.get(n, 0) for n in names}

    lines = []
    for start in range(0, L, width):
        for rec in aln:
            chunk = str(rec.seq[start:start + width])
            count[rec.id] += len(chunk.replace("-", ""))
            cells = []
            for i, ch in enumerate(chunk):
                color = colored.get((rec.id, start + i))
                ch = _escape(ch)
                cells.append(f'<span style="background:{color}">{ch}</span>' if color else ch)
            pad_right = " " * (width - len(chunk))
            lines.append(f"{_escape(rec.id + ':'):<{pad}}{''.join(cells)}{pad_right} {count[rec.id]}")
        if show_conservation:
            lines.append(f"{'':<{pad}}{_escape(cons[start:start + width])}")
        lines.append("")
    body = "\n".join(lines)
    return HTML(f'<pre style="font-family:\'{font}\',monospace;font-size:{fontsize}pt;'
                f'line-height:1.35;white-space:pre">{body}</pre>')

def motif_report(aln, motifs, regex=False, ignore_case=False, offsets=None):
    """Where each motif was found: one row per occurrence, in real residue numbering."""
    import pandas as pd
    offsets = offsets or {}
    rows = []
    for pattern in motifs:
        for rec in aln:
            gapped = str(rec.seq)
            cols = residue_columns(gapped)
            plain = "".join(ch for ch in gapped if ch != "-")
            pat = pattern if regex else _re.escape(pattern)
            for m in _re.finditer(pat, plain, _re.IGNORECASE if ignore_case else 0):
                off = offsets.get(rec.id, 0)
                rows.append({"motif": pattern, "sequence": rec.id, "match": m.group(0),
                             "start": m.start() + 1 + off, "end": m.end() + off,
                             "aln_col_start": cols[m.start()] + 1, "aln_col_end": cols[m.end() - 1] + 1,
                             "split_by_gap": (cols[m.end() - 1] - cols[m.start() ]) != (m.end() - m.start() - 1)})
    return pd.DataFrame(rows)
