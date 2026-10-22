"""Exon structure of a transcript, mapped onto protein residue positions.

    ex = exon_table("ENST00000609686")      # GRIN2B MANE
    parse_variant("c.2362G>T"); parse_variant("Glu788Ter")

Junctions are reported in protein coordinates so they can be drawn on a protein alignment.
"""
from netutil import pkg_path
import json, re
import pandas as pd
from netutil import urlopen_retry

ENS = "https://rest.ensembl.org"

def transcript_structure(transcript_id) -> dict:
    """Ensembl transcript with its exons and translation bounds."""
    tid = transcript_id.split(".")[0]
    return json.loads(urlopen_retry(
        f"{ENS}/lookup/id/{tid}?expand=1;content-type=application/json").decode())

def exon_table(transcript_id, cache_dir=None) -> pd.DataFrame:
    """One row per exon: exon number, genomic span, coding span, cumulative CDS position and
    the protein residues it encodes.

    Columns: exon, exon_id, genomic_start/end, cds_start/end (1-based, coding only),
             aa_start/aa_end, junction_aa (residue the 3' junction falls in), phase
             (0 = junction between codons, 1/2 = junction splits a codon).
    """
    import os
    cache_dir = cache_dir or pkg_path("exon_cache")
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, f"{transcript_id.split('.')[0]}.csv")
    if os.path.exists(path):
        return pd.read_csv(path)
    t = transcript_structure(transcript_id)
    strand = t["strand"]
    tl = t["Translation"]
    cds_lo, cds_hi = tl["start"], tl["end"]          # genomic bounds of the CDS
    exons = sorted(t["Exon"], key=lambda e: e["start"], reverse=(strand == -1))

    rows, cds_pos = [], 0
    for i, e in enumerate(exons, start=1):
        lo, hi = max(e["start"], cds_lo), min(e["end"], cds_hi)   # coding part of this exon
        coding = max(0, hi - lo + 1)
        row = {"exon": i, "exon_id": e["id"],
               "genomic_start": e["start"], "genomic_end": e["end"],
               "coding_bases": coding,
               # genomic bounds of the CODING part only (differs from the exon at the ends)
               "coding_genomic_start": lo if coding else None,
               "coding_genomic_end":   hi if coding else None}
        if coding:
            row["cds_start"] = cds_pos + 1
            cds_pos += coding
            row["cds_end"] = cds_pos
            row["aa_start"] = (row["cds_start"] - 1) // 3 + 1
            row["aa_end"] = (row["cds_end"] - 1) // 3 + 1
            row["junction_aa"] = row["aa_end"]         # residue containing the 3' junction
            row["phase"] = cds_pos % 3                 # 0 = clean break between codons
        else:
            row.update({k: None for k in ("cds_start","cds_end","aa_start","aa_end",
                                          "junction_aa","phase")})
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)
    return df

# ---------------------------------------------------------------- variant notation
AA3 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H",
       "Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W",
       "Tyr":"Y","Val":"V","Ter":"*","Sec":"U","Xaa":"X"}

def parse_variant(text, cds_offset=0):
    """Parse HGVS-ish variant notation into a protein position.

    Accepts cDNA ('c.2362G>T') -> position inferred as (n-1)//3 + 1,
    and protein ('Glu788Ter', 'p.Glu788Ter', 'p.E788*', 'E788*', 'p.T652A').
    Returns {'kind', 'aa_pos', 'ref', 'alt', 'label', 'raw'}; aa_pos is None if unparseable.
    """
    s = str(text).strip()
    m = re.fullmatch(r"c\.(\d+)([ACGT])>([ACGT])", s, re.I)
    if m:
        n = int(m.group(1)) + cds_offset
        return {"kind": "cdna", "aa_pos": (n - 1) // 3 + 1, "cds_pos": n,
                "ref": None, "alt": None, "label": s, "raw": s}
    m = re.fullmatch(r"(?:p\.)?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", s)
    if m:
        return {"kind": "protein", "aa_pos": int(m.group(2)),
                "ref": AA3.get(m.group(1)), "alt": AA3.get(m.group(3)),
                "label": s, "raw": s}
    m = re.fullmatch(r"(?:p\.)?([A-Z])(\d+)([A-Z*=])", s)
    if m:
        return {"kind": "protein", "aa_pos": int(m.group(2)),
                "ref": m.group(1), "alt": m.group(3), "label": s, "raw": s}
    return {"kind": None, "aa_pos": None, "ref": None, "alt": None, "label": s, "raw": s}

def variant_positions(variants, cds_offset=0) -> pd.DataFrame:
    """Parse a list of variant strings into a table of protein positions."""
    return pd.DataFrame([parse_variant(v, cds_offset) for v in variants])


# ---------------------------------------------------------------- rendering on an alignment
def junction_columns(gapped_seq, exon_df):
    """{alignment column: exon number starting at that column} for each exon's 5' boundary."""
    from msa_print import residue_columns
    idx = residue_columns(gapped_seq)
    out = {}
    coding = exon_df.dropna(subset=["aa_start"])
    for _, r in coding.iterrows():
        p = int(r["aa_start"])
        if 0 < p <= len(idx):
            out[idx[p - 1]] = int(r["exon"])
    return out

# default palette for ionotropic glutamate receptor domains; override via `domain_colors`
IGLUR_DOMAIN_COLORS = {"signal": "#dcdcdc", "ATD": "#cde2fb", "S1": "#cfe9d6", "S2": "#cfe9d6",
                       "LBD": "#cfe9d6", "M1": "#ffc9a3", "M2": "#ff9e6b", "M3": "#ffc9a3",
                       "M4": "#ffc9a3", "TM": "#ffc9a3", "CTD": "#e6d7f0"}
DOMAIN_PALETTE = ["#cde2fb", "#ffe0b3", "#cfe9d6", "#e6d7f0", "#f9d3d3",
                  "#d4ecec", "#ede3c8", "#dcdcdc", "#f6d9ea", "#d9e4c8"]

def receptor_domain_tracks(accessions, drop=("signal",)) -> dict:
    """{name: {domain: (start, end)}} from conservation.receptor_domains(), ready to pass as
    `domains=`. Edit the returned dict to add, rename or remove regions."""
    from conservation import receptor_domains
    return {n: {k: v for k, v in receptor_domains(a).items() if k not in (drop or ())}
            for n, a in accessions.items()}

def _domain_columns(gapped_seq, regions):
    """{alignment column: label} for one sequence's regions (later entries win on overlap)."""
    from msa_print import residue_columns
    idx = residue_columns(gapped_seq)
    out = {}
    for label, span in (regions or {}).items():
        if not span: continue
        for p in range(int(span[0]), int(span[1]) + 1):
            if 0 < p <= len(idx): out[idx[p - 1]] = label
    return out

def exon_alignment_html(aln, exon_tables, variants=None, width=60, offsets=None,
                        junction_mark="⌄", font="Courier New", fontsize=11,
                        exon_color="#1a5fb4", variant_color="#c00000", title=None,
                        show_conservation=None, domains=None, domain_colors=None,
                        domain_label_row=True, region=None, legend=True):
    """Alignment with an exon-junction row above each sequence and a variant row below.

    exon_tables: {seq_id: exon_table DataFrame}
    variants:    {seq_id: [variant strings]} -- parsed with parse_variant()
    domains:     {seq_id: {label: (start, end)}} drawn as background colour, in each sequence's
                 own numbering. Any labels work -- receptor_domain_tracks() gives the standard
                 iGluR set; add your own entries for further annotations.
    domain_colors: {label: colour}; unnamed labels take DOMAIN_PALETTE in order
    region:      (start, end) window in the FIRST sequence's numbering
    domains:     {seq_id: {label: (start, end)}} highlighted as background colour, in each
                 sequence's own numbering. Any labels work -- use receptor_domain_tracks()
                 for the standard iGluR set, or write your own to add more annotations.
    domain_colors: {label: colour}; unnamed labels take DOMAIN_PALETTE in order
    region:      (start, end) window in the FIRST sequence's numbering
    Junctions are marked with `junction_mark` at the first residue of each exon, labelled with
    the exon number; exons absent from a transcript simply never appear in its row.
    """
    from IPython.display import HTML
    from html import escape
    from msa_print import residue_columns, conservation_char
    offsets = offsets or {}
    L = aln.get_alignment_length()
    if show_conservation is None:            # meaningless for a single sequence
        show_conservation = len(aln) > 1

    dom_cols, dcolors = {}, dict(domain_colors or {})
    if domains:
        labels = [k for d in domains.values() for k in d if k not in dcolors]
        for i, k in enumerate(dict.fromkeys(labels)):
            dcolors.setdefault(k, IGLUR_DOMAIN_COLORS.get(k, DOMAIN_PALETTE[i % len(DOMAIN_PALETTE)]))
        dom_cols = {r.id: _domain_columns(str(r.seq), domains.get(r.id)) for r in aln}
    names = [r.id for r in aln]
    pad = max(len(n) for n in names) + 2
    cons = "".join(conservation_char(aln[:, i]) for i in range(L))
    count = {n: offsets.get(n, 0) for n in names}

    junc = {r.id: junction_columns(str(r.seq), exon_tables[r.id])
            for r in aln if r.id in exon_tables}
    vmap = {}
    for sid, vs in (variants or {}).items():
        rec = next((r for r in aln if r.id == sid), None)
        if rec is None: continue
        idx = residue_columns(str(rec.seq))
        for v in vs:
            p = parse_variant(v)
            if p["aa_pos"] and 0 < p["aa_pos"] <= len(idx):
                vmap.setdefault(sid, {}).setdefault(idx[p["aa_pos"] - 1], []).append(p["label"])

    def marker_row(mapping, start, stop, fmt):
        """Place labels at their columns without overlapping."""
        row = [" "] * (stop - start)
        for col in sorted(k for k in mapping if start <= k < stop):
            lab = fmt(mapping[col])
            i = col - start
            if all(row[j] == " " for j in range(i, min(i + len(lab), len(row)))):
                for j, ch in enumerate(lab):
                    if i + j < len(row): row[i + j] = ch
        return "".join(row)

    lo, hi = 0, L
    if region:
        idx0 = residue_columns(str(aln[0].seq))
        lo, hi = idx0[region[0] - 1], idx0[min(region[1], len(idx0)) - 1] + 1
    for rec in aln:
        count[rec.id] += sum(1 for c in str(rec.seq)[:lo] if c != "-")

    lines = []
    for start in range(lo, hi, width):
        stop = min(start + width, hi)
        for rec in aln:
            if rec.id in junc:
                r = marker_row(junc[rec.id], start, stop,
                               lambda n: f"{junction_mark}{n}")
                if r.strip():
                    lines.append(f'{"":<{pad}}<span style="color:{exon_color}">{escape(r)}</span>')
            dm = dom_cols.get(rec.id, {})
            if domains and domain_label_row and dm:
                row = [" "] * (stop - start)
                # contiguous runs of one label within this block -> label centred on the run
                runs, cur, run_start = [], None, None
                for c in range(start, stop + 1):
                    lab = dm.get(c) if c < stop else None
                    if lab != cur:
                        if cur is not None: runs.append((cur, run_start, c - 1))
                        cur, run_start = lab, c
                for lab, a, b in runs:
                    i = (a + b) // 2 - start - len(lab) // 2          # centre of the run
                    i = max(0, min(i, len(row) - len(lab)))           # keep inside the block
                    if all(row[j] == " " for j in range(i, min(i + len(lab), len(row)))):
                        for j, ch in enumerate(lab):
                            if i + j < len(row): row[i + j] = ch
                if "".join(row).strip():
                    lines.append(f'{"":<{pad}}<span style="color:#555;font-size:90%">'
                                 f'{escape("".join(row))}</span>')
            chunk = str(rec.seq[start:stop])
            count[rec.id] += len(chunk.replace("-", ""))
            cells = []
            for i, ch in enumerate(chunk):
                lab = dm.get(start + i)
                cells.append(f'<span style="background:{dcolors[lab]}">{escape(ch)}</span>'
                             if lab else escape(ch))
            lines.append(f"{escape(rec.id + ':'):<{pad}}{''.join(cells)}"
                         f"{' ' * (width - len(chunk))} {count[rec.id]}")
            if rec.id in vmap:
                r = marker_row(vmap[rec.id], start, stop, lambda labs: "|" + ",".join(labs))
                if r.strip():
                    lines.append(f'{"":<{pad}}<span style="color:{variant_color}">{escape(r)}</span>')
        if show_conservation:
            lines.append(f"{'':<{pad}}{escape(cons[start:stop])}")
        lines.append("")

    key = (f'<div style="font-family:sans-serif;font-size:10pt;margin:6px 0">'
           f'<span style="color:{exon_color}">{junction_mark}N</span> = first residue of exon N'
           + (f' &nbsp;&nbsp; <span style="color:{variant_color}">|variant</span>' if variants else '')
           + (('<br>' + " &nbsp; ".join(
                 f'<span style="background:{dcolors[k]};padding:1px 8px">&nbsp;</span> {escape(k)}'
                 for k in dict.fromkeys(k for d in domains.values() for k in d)))
              if (domains and legend) else '')
           + '</div>')
    head = f'<div style="font-family:sans-serif;font-size:12pt;font-weight:600">{escape(title)}</div>' if title else ""
    return HTML(head + key + f'<pre style="font-family:\'{font}\',monospace;font-size:{fontsize}pt;'
                             f'line-height:1.45;white-space:pre">{chr(10).join(lines)}</pre>')


# ---------------------------------------------------------------- codon-aligned figure
def cds_sequence(transcript_id, cache_dir=None) -> str:
    """Coding sequence of a transcript (Ensembl), cached alongside the exon tables."""
    import os
    cache_dir = cache_dir or pkg_path("exon_cache")
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, f"{transcript_id.split('.')[0]}_cds.txt")
    if os.path.exists(path):
        return open(path).read().strip()
    fa = urlopen_retry(f"{ENS}/sequence/id/{transcript_id.split('.')[0]}"
                       f"?type=cds;content-type=text/x-fasta").decode()
    cds = "".join(fa.split("\n")[1:]).strip()
    open(path, "w").write(cds)
    return cds

def junction_cds_positions(exon_df):
    """{cds position of the first base of each exon: exon number} (coding exons only)."""
    d = exon_df.dropna(subset=["cds_start"])
    return {int(r.cds_start): int(r.exon) for _, r in d.iterrows()}

# ---- figure colours: change these, or pass the matching keyword argument ----------------
JUNCTION_COLOR = "#e0218a"      # exon-exon junction rule (pink)   -> junction_color=
VARIANT_COLOR  = "#c00000"      # variant text and labels          -> variant_color=
TRUNCATION_COLOR = "#f28b82"    # truncating-variant labels (light) -> truncation_color=
MISSENSE_COLOR   = "#e8a0a0"    # missense-variant labels (lightest) -> missense_color=
DOMAIN_LABEL_COLOR = "#555555"  # the domain-name row              -> not exposed


def codon_alignment_html(protein, cds, exons=None, domains=None, domain_colors=None,
                         variants=(), name="", region=None, width=40, offsets=0,
                         font="Courier New", fontsize=11, variant_color=VARIANT_COLOR,
                         junction_color=JUNCTION_COLOR, domain_label_row=True, legend=True,
                         title=None, show_exon_numbers=True, codon_spacing=True,
                         left_index=False, exon_label_offset=2.5, cds_variant_above=True,
                         truncation_variants=(), truncation_color=TRUNCATION_COLOR,
                         missense_variants=(), missense_color=MISSENSE_COLOR,
                         variant_angle=45, legend_only=False):
    """Protein sequence with its codons aligned underneath, exon junctions marked by a vertical
    rule, and variants highlighted.

    protein:  amino-acid sequence (no stop)
    cds:      coding sequence; 3 x len(protein), a retained stop codon is trimmed
    exons:    exon_table() DataFrame -- junctions drawn at each exon's first CDS base
    domains:  {label: (start, end)} in protein numbering; shading is applied to the amino-acid
              row only, so the nucleotide row stays legible
    variants: HGVS strings; 'c.NNNN X>Y' colours the base and is labelled below, protein
              notation ('p.E788*', 'Glu788Ter') colours the residue and is labelled above
    region:   (start, end) in protein numbering
    width:    codons per block
    codon_spacing: insert a space between codons (4 character cells each instead of 3)
    left_index:    also print the first residue / CDS position at the start of each line
    exon_label_offset: how far (in codons) the "exon N" labels sit either side of a junction
    cds_variant_above: put the cDNA label on the row above, merged with the protein change on
              one line ("c.2362G>T, p.E788*"), its first character over the affected residue.
              False keeps the cDNA label on its own row below the nucleotides.
    truncation_variants: a second list (e.g. ClinVar truncating variants) drawn in
              truncation_color -- a lighter red -- so they read as a separate set
    missense_variants:   a third list drawn in missense_color (lighter still)
    variant_angle:  degrees to rotate the variant labels counter-clockwise; the first
              character's bottom-left corner sits over the residue it affects. 0 = horizontal.
    legend_only:    return just the colour key (used to put one legend above several panels)
    """
    from IPython.display import HTML
    from html import escape
    n_aa = len(protein)
    if len(cds) >= 3 * n_aa + 3:
        cds = cds[:3 * n_aa]
    if len(cds) != 3 * n_aa:
        raise ValueError(f"cds length {len(cds)} != 3 x protein length {3*n_aa}")

    cw = 4 if codon_spacing else 3          # character cells per codon

    dcolors = dict(domain_colors or {})
    if domains:
        for i, k in enumerate(domains):
            dcolors.setdefault(k, IGLUR_DOMAIN_COLORS.get(k, DOMAIN_PALETTE[i % len(DOMAIN_PALETTE)]))
    dom_at = {}
    for lab, span in (domains or {}).items():
        if span:
            for p in range(int(span[0]), int(span[1]) + 1):
                dom_at[p] = lab

    junc = junction_cds_positions(exons) if exons is not None else {}

    def collect(vlist):
        """-> (by_residue, hit_residues, hit_bases).

        by_residue maps residue -> [(cdna_part, protein_part), ...]; a variant string may give
        either or both ('c.2362G>T, p.E788*'), so the two halves stay separate and the figure
        can place them above and/or below."""
        by_res, res, base = {}, set(), set()
        for v in vlist:
            parts = [x.strip() for x in str(v).split(",") if x.strip()]
            parsed = [parse_variant(x) for x in parts]
            cdna = next((q for q in parsed if q["kind"] == "cdna"), None)
            prot = next((q for q in parsed if q["kind"] == "protein" and q["aa_pos"]), None)
            pos = prot["aa_pos"] if prot else ((cdna["cds_pos"] - 1) // 3 + 1 if cdna else None)
            if pos is None: continue
            if cdna: base.add(cdna["cds_pos"])
            res.add(pos)
            by_res.setdefault(pos, []).append((cdna["label"] if cdna else None,
                                               prot["label"] if prot else None))
        return by_res, res, base

    v_by, hit_res, hit_base = collect(variants)
    t_by, t_res, t_base = collect(truncation_variants)
    m_by, m_res, m_base = collect(missense_variants)



    lo, hi = (region or (1, n_aa))
    lo, hi = max(1, lo), min(n_aa, hi)
    pad = max(len(name) + 2, 10)
    iw = 7 if left_index else 0             # width of the left index column

    def blank(block_lo, block_hi):
        return [" "] * ((block_hi - block_lo + 1) * cw)

    def put(row, i, lab):
        i = max(0, min(i, len(row) - len(lab)))
        if all(row[j] == " " for j in range(i, min(i + len(lab), len(row)))):
            for j, ch in enumerate(lab): row[i + j] = ch

    def prefix(idx):
        return f"{idx:<{iw}}" if left_index else ""

    def build_key():
        parts = [f'<span style="background:{dcolors[k]};padding:1px 8px">&nbsp;</span> {escape(k)}'
                 for k in (domains or {})]
        parts.append(f'<span style="border-left:3px solid {junction_color};padding-left:5px">'
                     f'&nbsp;</span> exon junction')
        parts.append(f'<span style="color:{variant_color};font-weight:700">variant</span>')
        if len(truncation_variants):
            parts.append(f'<span style="color:{truncation_color};font-weight:700">truncating</span>')
        if len(missense_variants):
            parts.append(f'<span style="color:{missense_color};font-weight:700">missense</span>')
        return ('<div style="font-family:sans-serif;font-size:10pt;margin:8px 0">'
                + " &nbsp; ".join(parts) + "</div>")

    if legend_only:
        return HTML(build_key())

    lines = []
    for blo in range(lo, hi + 1, width):
        bhi = min(blo + width - 1, hi)
        span = range(blo, bhi + 1)

        if domains and domain_label_row:
            row, runs, cur, rs = blank(blo, bhi), [], None, None
            for p in list(span) + [bhi + 1]:
                lab = dom_at.get(p)
                if lab != cur:
                    if cur is not None: runs.append((cur, rs, p - 1))
                    cur, rs = lab, p
            for lab, a, b in runs:
                put(row, ((a + b) // 2 - blo) * cw + 1 - len(lab) // 2, lab)
            if "".join(row).strip():
                lines.append(f'{"":<{pad+iw}}<span style="color:{DOMAIN_LABEL_COLOR};font-size:90%">'
                             f'{escape("".join(row))}</span>')

        # ---- variant labels -------------------------------------------------------------
        # merge every tier by residue first, so two variants at one position become a list
        # rather than two spans drawn on top of each other
        tiers = ((v_by, variant_color), (t_by, truncation_color), (m_by, missense_color))
        above_at, below_at = {}, {}
        for by_res, colour in tiers:
            for pos, comps in by_res.items():
                if cds_variant_above:                      # both halves on one line above
                    lab = ", ".join(", ".join(x for x in c if x) for c in comps)
                    if lab: above_at.setdefault(pos, []).append((lab, colour))
                else:                                      # protein above, cDNA below
                    up = ", ".join(c[1] for c in comps if c[1])
                    dn = ", ".join(c[0] for c in comps if c[0])
                    if up: above_at.setdefault(pos, []).append((up, colour))
                    if dn: below_at.setdefault(pos, []).append((dn, colour))

        def rotated_row(mapping, clockwise):
            """Absolutely positioned labels anchored at the residue; components keep their own
            colour. clockwise=False rises above the row, True hangs below it."""
            import math
            items = {p_: v for p_, v in mapping.items() if blo <= p_ <= bhi}
            if not items: return None
            longest = max(sum(len(l) for l, _ in v) + 2 * (len(v) - 1) for v in items.values())
            rad = math.radians(variant_angle)
            height = round(longest * 0.60 * math.sin(rad) * fontsize * 1.34) + 8
            spans = []
            for pos, comps in sorted(items.items()):
                col = pad + iw + (pos - blo) * cw + 1
                inner = "; ".join(f'<span style="color:{c}">{escape(l)}</span>' for l, c in comps)
                anchor = "top" if clockwise else "bottom"
                spans.append(
                    f'<span style="position:absolute;left:{col}ch;{anchor}:0;'
                    f'transform:rotate({"" if clockwise else "-"}{variant_angle}deg);'
                    f'transform-origin:left {anchor};white-space:pre">{inner}</span>')
            return f'<div style="position:relative;height:{height}px">' + "".join(spans) + "</div>"

        def flat_row(mapping, align_first):
            """Plain text row; a merged label takes the colour of its highest-priority tier."""
            order = [variant_color, truncation_color, missense_color]
            items = {p_: v for p_, v in mapping.items() if blo <= p_ <= bhi}
            if not items: return None
            row, colour_at = blank(blo, bhi), {}
            for pos, comps in sorted(items.items()):
                lab = "; ".join(l for l, _ in comps)
                colour_at[pos] = min((c for _, c in comps),
                                     key=lambda c: order.index(c) if c in order else 99)
                put(row, (pos - blo) * cw + 1 if align_first
                         else (pos - blo) * cw + 1 - len(lab) // 2, lab)
            txt = "".join(row)
            if not txt.strip(): return None
            colour = min(colour_at.values(), key=lambda c: order.index(c) if c in order else 99)
            return f'{"":<{pad+iw}}<span style="color:{colour}">{escape(txt)}</span>'

        row_html = (rotated_row(above_at, clockwise=False) if variant_angle
                    else flat_row(above_at, align_first=cds_variant_above))
        if row_html: lines.append(row_html)

        cells = []                                                # amino acids (domain shading here)
        for p in span:
            aa, lab = protein[p - 1], dom_at.get(p)
            style = []
            if lab: style.append(f"background:{dcolors[lab]}")
            if p in hit_res: style.append(f"color:{variant_color};font-weight:700")
            elif p in t_res: style.append(f"color:{truncation_color};font-weight:700")
            elif p in m_res: style.append(f"color:{missense_color};font-weight:700")
            body = f'<span style="{";".join(style)}">&nbsp;{escape(aa)}&nbsp;</span>' if style \
                   else f"&nbsp;{escape(aa)}&nbsp;"
            if codon_spacing:
                # shade the gap too when the run continues, so the block reads as one band
                same = lab and p < bhi and dom_at.get(p + 1) == lab
                body += (f'<span style="background:{dcolors[lab]}">&nbsp;</span>' if same
                         else "&nbsp;")
            cells.append(body)
        lines.append(f"{escape(name + ':'):<{pad}}{prefix(blo)}{''.join(cells)}   {bhi + offsets}")

        cells = []                                                # nucleotides (no shading)
        for p in span:
            for k in range(3):
                c = (p - 1) * 3 + k + 1
                style = []
                if c in hit_base: style.append(f"color:{variant_color};font-weight:700")
                elif c in t_base: style.append(f"color:{truncation_color};font-weight:700")
                elif c in m_base: style.append(f"color:{missense_color};font-weight:700")
                if c in junc:    style.append(f"border-left:2px solid {junction_color}")
                base = escape(cds[c - 1])
                cells.append(f'<span style="{";".join(style)}">{base}</span>' if style else base)
            if codon_spacing: cells.append(" ")
        nt_first = (blo - 1) * 3 + 1
        lines.append(f"{'':<{pad}}{prefix(nt_first)}{''.join(cells)}   {bhi * 3 + offsets * 3}")

        if show_exon_numbers and junc:
            on_line = sorted(c for c in junc if (blo - 1) * 3 < c <= bhi * 3)
            if on_line:                       # label only lines that contain a junction
                row = blank(blo, bhi)
                gap = int(round(exon_label_offset * cw))   # characters either side of the bar
                for c in on_line:
                    n = junc[c]
                    res, k = (c - 1) // 3 + 1, (c - 1) % 3
                    bar = (res - blo) * cw + k             # character offset of the junction rule
                    prev_n = n - 1
                    if prev_n >= 1:                        # exon ending at this junction, to the left
                        lab = f"exon {prev_n}"
                        put(row, bar - gap - len(lab) // 2, lab)
                    lab = f"exon {n}"                      # exon starting here, to the right
                    put(row, bar + gap - len(lab) // 2, lab)
                if "".join(row).strip():
                    lines.append(f'{"":<{pad+iw}}<span style="color:{junction_color};font-size:90%">'
                                 f'{escape("".join(row))}</span>')

        if below_at:                       # cDNA labels below, first character under the residue
            row_html = (rotated_row(below_at, clockwise=True) if variant_angle
                        else flat_row(below_at, align_first=True))
            if row_html: lines.append(row_html)
        lines.append("")

    key = build_key() if legend else ""

    head = (f'<div style="font-family:sans-serif;font-size:12pt;font-weight:600">{escape(title)}</div>'
            if title else "")
    body = "".join(l if l.startswith("<div") else f'<div style="white-space:pre">{l}</div>'
                   for l in lines)
    return HTML(head + key + f'<div style="font-family:\'{font}\',monospace;'
                             f'font-size:{fontsize}pt;line-height:1.5">{body}</div>')


# ---------------------------------------------------------------- genomic -> CDS, HGVS labels
COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]

def genomic_to_cds(genomic_pos, exon_df, strand):
    """CDS position (1-based) for a genomic coordinate, or None if not in the CDS."""
    p = int(genomic_pos)
    for _, r in exon_df.dropna(subset=["cds_start"]).iterrows():
        c_lo, c_hi = int(r.coding_genomic_start), int(r.coding_genomic_end)
        if not (c_lo <= p <= c_hi):
            continue
        # CDS runs 5'->3': with genomic position for + strand, against it for -
        return int(r.cds_start) + ((p - c_lo) if strand == 1 else (c_hi - p))
    return None

def hgvs_c_from_variant_id(variant_id, exon_df, strand, cds=None):
    """'12-13753873-G-A' -> 'c.454C>T'. Handles SNVs, and gives an approximate del/dup/ins
    notation for indels (no 3'-shifting, so a repeat-region indel may differ from ClinVar)."""
    parts = str(variant_id).split("-")
    if len(parts) < 4: return None
    pos, ref, alt = int(parts[1]), parts[2], parts[3]
    c = genomic_to_cds(pos, exon_df, strand)
    if c is None: return None
    if strand == -1:
        ref, alt = revcomp(ref), revcomp(alt)
    if len(ref) == 1 and len(alt) == 1:
        return f"c.{c}{ref}>{alt}"
    if len(ref) > len(alt):                                   # deletion
        n = len(ref) - len(alt)
        a, b = (c + 1, c + n) if strand == 1 else (c - n, c - 1)
        a, b = min(a, b), max(a, b)
        return f"c.{a}del" if n == 1 else f"c.{a}_{b}del"
    if len(alt) > len(ref):                                   # insertion / duplication
        ins = alt[len(ref):] if strand == 1 else revcomp(alt)[:len(alt) - len(ref)]
        if cds and len(ins) == 1 and cds[c - 1:c - 1 + len(ins)] == ins:
            return f"c.{c}dup"
        return f"c.{c}_{c + 1}ins{ins}"
    return f"c.{c}{ref}>{alt}"

def label_variants(items, exon_df, strand, cds=None, joiner=", "):
    """[{variant_id: effect}] -> ['c.454C>T, p.Q152*', ...] using the effect's
    short_description for the protein half. Falls back to the protein change alone."""
    out = []
    for item in items:
        vid = list(item.keys())[0]
        eff = item[vid]
        p = getattr(eff, "short_description", str(eff))
        c = hgvs_c_from_variant_id(vid, exon_df, strand, cds)
        out.append(f"{c}{joiner}{p}" if c else p)
    return out


def junction_panels_html(protein, cds, exons, flank=10, width=None, title=None,
                         panel_titles=True, name="", **kwargs):
    """One codon_alignment_html panel per exon-exon junction, stacked into a single document.

    Each panel shows `flank` residues either side of the junction. Any other keyword argument
    (domains, variants, truncation_variants, missense_variants, variant_angle, left_index,
    cds_variant_above, ...) is passed straight through to codon_alignment_html.

        html = junction_panels_html(prot, cds, ex, domains=dom, variants=v, flank=10)
        open("junctions.html", "w").write(html.data)
    """
    from IPython.display import HTML
    from html import escape
    coding = exons.dropna(subset=["cds_start"]).sort_values("cds_start")
    n_aa = len(protein)
    width = width or (2 * flank + 1)

    panels, rows = [], []
    for _, r in coding.iloc[1:].iterrows():            # every junction after the first exon
        exon_n = int(r.exon)
        junction_aa = (int(r.cds_start) - 1) // 3 + 1
        lo, hi = max(1, junction_aa - flank), min(n_aa, junction_aa + flank)
        h = codon_alignment_html(protein, cds, exons=exons, name=name,
                                 region=(lo, hi), width=width,
                                 title=None, legend=False, **kwargs)
        head = ""
        if panel_titles:
            head = (f'<div style="font-family:sans-serif;font-size:11pt;font-weight:600;'
                    f'margin:18px 0 4px 0;color:#333">'
                    f'exon {exon_n - 1} | exon {exon_n} &nbsp;&mdash;&nbsp; '
                    f'residues {lo}&ndash;{hi} (junction at {junction_aa}, CDS {int(r.cds_start)})'
                    f'</div>')
        panels.append(head + h.data)
        rows.append({"exon_5prime": exon_n - 1, "exon_3prime": exon_n,
                     "junction_cds": int(r.cds_start), "junction_aa": junction_aa,
                     "phase": int(r.phase) if pd.notna(r.phase) else None,
                     "region_start": lo, "region_end": hi})

    # one legend for the whole document
    legend = codon_alignment_html(protein, cds, exons=exons, region=(1, min(width, n_aa)),
                                  width=width, title=None, legend_only=True, **kwargs).data
    head = (f'<div style="font-family:sans-serif;font-size:13pt;font-weight:700;'
            f'margin-bottom:6px">{escape(title)}</div>' if title else "")
    doc = head + legend + "".join(panels)
    out = HTML(doc)
    out.junctions = pd.DataFrame(rows)                 # summary table alongside the figure
    return out
