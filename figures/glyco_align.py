"""Align the neuroligins and highlight predicted O-glycosylation sites, signal peptide and TM.

    from glyco_align import parse_netoglyc, glyco_alignment_html
    pred = parse_netoglyc("nlgn_o_glyc_pred.txt")
    glyco_alignment_html(pred, accessions, threshold=0.5)
"""
from netutil import pkg_path
import re
import pandas as pd
from conservation import fetch_uniprot, uniprot_topology
from msa_print import align, residue_columns, conservation_char
from html import escape

def parse_netoglyc(path) -> pd.DataFrame:
    """Parse one or more concatenated NetOGlyc 4 GFF outputs.
    Blocks are separated by a bare 'GENE:' line. Returns gene, position, score, positive."""
    rows, gene = [], None
    for line in open(path):
        line = line.rstrip("\n")
        m = re.fullmatch(r"([A-Za-z0-9_.-]+):\s*", line)
        if m:
            gene = m.group(1); continue
        if line.startswith("#") or "CARBOHYD" not in line:
            continue
        f = line.split("\t") if "\t" in line else line.split()
        pos, score = int(f[3]), float(f[5])
        rows.append({"gene": gene, "position": pos, "score": score,
                     "positive": "POSITIVE" in line})
    return pd.DataFrame(rows)

# feature colours
C_GLYCO_HI  = "#e8811a"      # predicted site, score >= threshold
C_GLYCO_LO  = "#fde2c0"      # sub-threshold site (shown faintly)
C_SIGNAL    = "#c9c9c9"      # signal peptide
C_TM        = "#9ecae1"      # transmembrane helix

def feature_columns(aln, accessions, pred, threshold=0.5, show_subthreshold=True,
                    topologies=None):
    """{(seq_id, alignment column): colour} for glyco sites, signal peptide and TM."""
    colored = {}
    topologies = topologies or {}
    for rec in aln:
        g = rec.id
        cols = residue_columns(str(rec.seq))
        topo = topologies.get(g) or uniprot_topology(accessions[g])
        def paint(start, end, color, overwrite=True):
            for p in range(start, end + 1):
                if 0 < p <= len(cols):
                    key = (g, cols[p - 1])
                    if overwrite or key not in colored:
                        colored[key] = color
        if topo.get("signal"): paint(*topo["signal"], C_SIGNAL)
        for s, e in topo.get("transmembrane", []): paint(s, e, C_TM)
        sub = pred[pred.gene == g]
        if show_subthreshold:
            for p in sub[sub.score < threshold].position: paint(p, p, C_GLYCO_LO)
        for p in sub[sub.score >= threshold].position: paint(p, p, C_GLYCO_HI)
    return colored

def glyco_alignment_html(pred, accessions, threshold=0.5, width=60, show_subthreshold=True,
                         region=None, font="Courier New", fontsize=10.5, seqs=None, aln=None,
                         topologies=None):
    """Clustal-style HTML alignment with features highlighted.
    region: (start, end) in the FIRST sequence's numbering, to show only that window."""
    from IPython.display import HTML
    seqs = seqs or fetch_uniprot(accessions)
    aln = aln if aln is not None else align(seqs)
    topologies = topologies or {g: topology_in(seqs[g], a) for g, a in accessions.items()}
    colored = feature_columns(aln, accessions, pred, threshold, show_subthreshold, topologies)

    L = aln.get_alignment_length()
    lo, hi = 0, L
    if region:
        ref = aln[0]
        cols = residue_columns(str(ref.seq))
        lo, hi = cols[region[0] - 1], cols[min(region[1], len(cols)) - 1] + 1

    names = [r.id for r in aln]
    pad = max(len(n) for n in names) + 2
    show_cons = len(aln) > 1                 # meaningless for a single sequence
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
                col = colored.get((rec.id, start + i))
                cells.append(f'<span style="background:{col}">{escape(ch)}</span>' if col else escape(ch))
            cells.append(" " * (width - len(chunk)))
            lines.append(f"{escape(rec.id + ':'):<{pad}}{''.join(cells)} {count[rec.id]}")
        if show_cons:
            lines.append(f"{'':<{pad}}{escape(cons[start:stop])}")
        lines.append("")

    key = (f'<div style="font-family:sans-serif;font-size:10pt;margin:8px 0">'
           f'<b>Key</b>&nbsp; '
           f'<span style="background:{C_GLYCO_HI};padding:1px 6px">&nbsp;</span> predicted O-glycosylation (score &ge; {threshold}) &nbsp; '
           + (f'<span style="background:{C_GLYCO_LO};padding:1px 6px">&nbsp;</span> sub-threshold site &nbsp; ' if show_subthreshold else '')
           + f'<span style="background:{C_SIGNAL};padding:1px 6px">&nbsp;</span> signal peptide &nbsp; '
           f'<span style="background:{C_TM};padding:1px 6px">&nbsp;</span> transmembrane helix</div>')
    body = "\n".join(lines)
    return HTML(key + f'<pre style="font-family:\'{font}\',monospace;font-size:{fontsize}pt;'
                      f'line-height:1.35;white-space:pre">{body}</pre>')

def glyco_summary(pred, accessions, threshold=0.5, topologies=None, seqs=None) -> pd.DataFrame:
    """Predicted sites per gene, split by domain (stalk = between the globular domain and TM)."""
    topologies = topologies or {g: topology_in(seqs[g], a) if seqs else uniprot_topology(a)
                                for g, a in accessions.items()}
    rows = []
    for g, topo in topologies.items():
        tm = min(topo["transmembrane"])[0] if topo["transmembrane"] else None
        sub = pred[pred.gene == g]
        sig_end = topo["signal"][1] if topo["signal"] else 0
        def count(lo, hi):
            s = sub[(sub.position > lo) & (sub.position <= hi)]
            return len(s), int((s.score >= threshold).sum()), round(s.score.max(), 2) if len(s) else None
        stalk_lo = tm - 80 if tm else 0
        for label, (lo, hi) in {"signal peptide": (0, sig_end),
                                "ectodomain": (sig_end, stalk_lo),
                                "stalk (80 aa before TM)": (stalk_lo, tm - 1),
                                "TM + cytoplasmic": (tm - 1, 10**6)}.items():
            n, npos, mx = count(lo, hi)
            rows.append({"gene": g, "region": label, "n_sites": n, "n_positive": npos, "max_score": mx})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- MANE-based inputs
NLGN_MANE = {"NLGN1": "ENSP00000511841", "NLGN2": "ENSP00000305288", "NLGN3": "ENSP00000351591",
             "NLGN4X": "ENSP00000370485", "NLGN4Y": "ENSP00000510011"}

def mane_seqs(ids=NLGN_MANE, cache=None):
    """{gene: MANE Select protein sequence}, cached. NetOGlyc is usually run on these, and for
    NLGN1 / NLGN4Y they differ in length from the UniProt canonical entry."""
    import json, os
    from netutil import urlopen_retry
    cache = cache or pkg_path("mane_proteins_nlgn.json")
    seqs = json.load(open(cache)) if os.path.exists(cache) else {}
    for g, p in ids.items():
        if g in seqs: continue
        fa = urlopen_retry(f"https://rest.ensembl.org/sequence/id/{p}?type=protein;content-type=text/x-fasta").decode()
        seqs[g] = "".join(fa.split("\n")[1:]).strip().rstrip("*")
    json.dump(seqs, open(cache, "w"))
    return {g: seqs[g] for g in ids}

def topology_in(seq, accession):
    """UniProt topology transferred onto `seq`'s numbering by alignment (no-op if identical)."""
    from structure_pipeline import map_all_positions
    from structure import fetch_uniprot_seq
    topo = uniprot_topology(accession)
    ref = fetch_uniprot_seq(accession)
    if ref == seq:
        return topo
    m = map_all_positions(ref, seq)
    def tr(rng):
        a, b = m.get(rng[0]), m.get(rng[1])
        return (a, b) if a and b else None
    out = {"signal": tr(topo["signal"]) if topo["signal"] else None,
           "transmembrane": [r for r in (tr(x) for x in topo["transmembrane"]) if r],
           "extracellular": [r for r in (tr(x) for x in topo["extracellular"]) if r],
           "cytoplasmic": [r for r in (tr(x) for x in topo["cytoplasmic"]) if r]}
    out["ctd"] = max(out["cytoplasmic"], key=lambda r: r[0]) if out["cytoplasmic"] else None
    out["length"] = len(seq)
    return out

def set_tm_start(topologies, seqs, motif="LSVT"):
    """Start each sequence's first TM helix at `motif` (UniProt starts the neuroligin TM two
    residues later, at VT).  The TM end is kept.  Raises if the motif is missing or repeated,
    so a sequence never silently keeps the old boundary."""
    import re, copy
    out = copy.deepcopy(topologies)
    for name, topo in out.items():
        hits = [m.start() + 1 for m in re.finditer(motif, seqs[name])]
        if len(hits) != 1:
            raise ValueError(f"{name}: expected one {motif}, found {len(hits)} at {hits}")
        tms = sorted(topo["transmembrane"])
        tms[0] = (hits[0], tms[0][1])
        topo["transmembrane"] = tms
    return out

def drop_tm_sites(pred, topologies):
    """Remove O-glycosylation predictions inside a TM helix -- the luminal glycosyltransferases
    cannot reach them, and drawn on top they would hide the TM shading."""
    keep = pd.Series(True, index=pred.index)
    for name, topo in topologies.items():
        for s, e in topo.get("transmembrane") or []:
            keep &= ~((pred.gene == name) & pred.position.between(s, e))
    return pred[keep]

def check_predictions(pred, seqs) -> pd.DataFrame:
    """Sanity check: NetOGlyc only scores Ser/Thr, so every predicted position must be S or T in
    the sequence the prediction was run on. A low match rate means the wrong isoform."""
    rows = []
    for g, s in seqs.items():
        sub = pred[pred.gene == g]
        if not len(sub): continue
        ok = sum(1 for p in sub.position if p <= len(s) and s[p - 1] in "ST")
        rows.append({"gene": g, "seq_len": len(s), "max_predicted_pos": int(sub.position.max()),
                     "n_sites": len(sub), "matching_ST": ok,
                     "pct_match": round(100 * ok / len(sub), 1),
                     "ok": ok == len(sub)})
    return pd.DataFrame(rows)


def glyco_alignment_rtf(pred, accessions, out="alignment.rtf", threshold=0.5, width=60,
                        show_subthreshold=True, seqs=None, aln=None, topologies=None,
                        font="Courier New", fontsize=10, title=None):
    """Same figure as glyco_alignment_html, written as RTF -- editable text that keeps its
    highlighting when pasted into Word or PowerPoint."""
    from export_figure import alignment_rtf
    seqs = seqs or fetch_uniprot(accessions)
    aln = aln if aln is not None else align(seqs)
    topologies = topologies or {g: topology_in(seqs[g], a) for g, a in accessions.items()}
    colored = feature_columns(aln, accessions, pred, threshold, show_subthreshold, topologies)
    key = [(f"predicted O-glycosylation (score >= {threshold})", C_GLYCO_HI)]
    if show_subthreshold: key.append(("sub-threshold site", C_GLYCO_LO))
    key += [("signal peptide", C_SIGNAL), ("transmembrane helix", C_TM)]
    return alignment_rtf(aln, colored, out=out, width=width, font=font, fontsize=fontsize,
                         key=key, title=title)


# ---------------------------------------------------------------- cleavage-site overlay
PROTEASE_COLORS = {"MMP3": "#d62728", "MMP9": "#1f77b4", "ADAM10": "#2ca02c",
                   "ADAM17": "#9467bd", "MMP2": "#ff7f0e", "MMP7": "#8c564b"}

def cleavage_overlay(aln, sites, protease_colors=PROTEASE_COLORS):
    """{(seq_id, column): {'underline': colour, 'cut': colour}} from a scan_juxtamembrane table.
    The matched motif is underlined; the scissile bond gets a vertical rule before P1'."""
    from msa_print import residue_columns
    cols = {r.id: residue_columns(str(r.seq)) for r in aln}
    ov = {}
    for _, s in sites.iterrows():
        name = s["protein"]
        if name not in cols: continue
        c = protease_colors.get(s["protease"], "#333333")
        idx = cols[name]
        for p in range(int(s["start"]), int(s["end"]) + 1):
            if 0 < p <= len(idx):
                ov.setdefault((name, idx[p - 1]), {}).setdefault("underline", c)
        p1prime = int(str(s["cut_between"]).split("|")[1])
        if 0 < p1prime <= len(idx):
            ov.setdefault((name, idx[p1prime - 1]), {})["cut"] = c
    return ov

def glyco_cleavage_html(pred, accessions, sites, threshold=0.5, width=60, region=None,
                        show_subthreshold=True, seqs=None, aln=None, topologies=None,
                        protease_colors=PROTEASE_COLORS, font="Courier New", fontsize=11,
                        title=None, missense=None, clinvar_color="#c00000"):
    """Alignment with O-glycosylation / signal peptide / TM as background colours, plus
    protease cleavage sites as coloured underlines with a rule at the scissile bond.

    sites: the DataFrame from mmp_sites.scan_juxtamembrane (or any table with columns
           protein, protease, start, end, cut_between).
    missense: optional table (variants.load_missense) -- residues carrying a ClinVar variant are
           drawn in clinvar_color, residues with a gnomAD variant in bold, both if both.
           Positions must use the same numbering as `seqs` (MANE)."""
    from IPython.display import HTML
    from msa_print import residue_columns, conservation_char
    seqs = seqs or fetch_uniprot(accessions)
    aln = aln if aln is not None else align(seqs)
    topologies = topologies or {g: topology_in(seqs[g], a) for g, a in accessions.items()}
    bg = feature_columns(aln, accessions, pred, threshold, show_subthreshold, topologies)
    ov = cleavage_overlay(aln, sites, protease_colors)
    var = variant_overlay(aln, missense) if missense is not None else {}

    L = aln.get_alignment_length()
    lo, hi = 0, L
    if region:
        idx = residue_columns(str(aln[0].seq))
        lo, hi = idx[region[0] - 1], idx[min(region[1], len(idx)) - 1] + 1
    names = [r.id for r in aln]
    pad = max(len(n) for n in names) + 2
    show_cons = len(aln) > 1                 # meaningless for a single sequence
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
                col = start + i
                style = []
                b = bg.get((rec.id, col))
                if b: style.append(f"background:{b}")
                o = ov.get((rec.id, col), {})
                if "underline" in o: style.append(f"box-shadow: inset 0 -3px 0 {o['underline']}")
                if "cut" in o:       style.append(f"border-left:2px solid {o['cut']}")
                v = var.get((rec.id, col), {})
                if v.get("clinvar"): style.append(f"color:{clinvar_color}")
                if v.get("gnomad"):  style.append("font-weight:700")
                cells.append(f'<span style="{";".join(style)}">{escape(ch)}</span>' if style else escape(ch))
            cells.append(" " * (width - len(chunk)))
            lines.append(f"{escape(rec.id + ':'):<{pad}}{''.join(cells)} {count[rec.id]}")
        if show_cons:
            lines.append(f"{'':<{pad}}{escape(cons[start:stop])}")
        lines.append("")

    used = [p for p in protease_colors if p in set(sites.protease)]
    key = (f'<div style="font-family:sans-serif;font-size:10pt;margin:8px 0;line-height:1.9">'
           f'<b>Features</b>&nbsp; '
           f'<span style="background:{C_GLYCO_HI};padding:1px 7px">&nbsp;</span> predicted O-glycosylation (&ge;{threshold}) &nbsp; '
           + (f'<span style="background:{C_GLYCO_LO};padding:1px 7px">&nbsp;</span> sub-threshold &nbsp; ' if show_subthreshold else '')
           + f'<span style="background:{C_SIGNAL};padding:1px 7px">&nbsp;</span> signal peptide &nbsp; '
             f'<span style="background:{C_TM};padding:1px 7px">&nbsp;</span> transmembrane<br>'
             f'<b>Cleavage</b>&nbsp; '
           + " &nbsp; ".join(
               f'<span style="box-shadow: inset 0 -3px 0 {protease_colors[p]};padding:1px 8px">&nbsp;</span> {p}'
               for p in used)
           + ' &nbsp;&nbsp; <span style="border-left:2px solid #333;padding-left:5px">|</span> scissile bond'
           + (f'<br><b>Variants</b>&nbsp; <span style="color:{clinvar_color}">ClinVar</span> &nbsp; '
              f'<b>gnomAD</b> &nbsp; <b style="color:{clinvar_color}">both</b>' if missense is not None else '')
           + '</div>')
    body = "\n".join(lines)
    head = f'<div style="font-family:sans-serif;font-size:12pt;font-weight:600">{escape(title)}</div>' if title else ""
    return HTML(head + key + f'<pre style="font-family:\'{font}\',monospace;font-size:{fontsize}pt;'
                             f'line-height:1.6;white-space:pre">{body}</pre>')


# ---------------------------------------------------------------- variant overlay
def variant_overlay(aln, missense, genes=None):
    """{(seq_id, column): {'clinvar': bool, 'gnomad': bool}} from the missense table.
    Positions are assumed to be in the same numbering as the aligned sequences (MANE)."""
    from msa_print import residue_columns
    cols = {r.id: residue_columns(str(r.seq)) for r in aln}
    ov = {}
    for g, idx in cols.items():
        if genes and g not in genes: continue
        sub = missense[missense.gene == g]
        for src, flag in (("ClinVar", "clinvar"), ("gnomAD", "gnomad")):
            for p in sub.position[sub.source == src].unique():
                if 0 < p <= len(idx):
                    ov.setdefault((g, idx[int(p) - 1]), {})[flag] = True
    return ov

def variants_in_region(missense, seqs, region_by_gene) -> pd.DataFrame:
    """Variants falling inside a per-gene residue window.
    region_by_gene: {gene: (start, end)} inclusive, in the sequence's own numbering."""
    rows = []
    for g, (lo, hi) in region_by_gene.items():
        sub = missense[(missense.gene == g) & (missense.position >= lo) & (missense.position <= hi)]
        if not len(sub): continue
        grp = (sub.groupby(["gene", "position", "reference", "alternate", "protein_consequence"])
                  .agg(sources=("source", lambda s: "+".join(sorted(set(s)))),
                       n_rows=("source", "size"),
                       gnomad_alleles=("total_count", "sum"),
                       blosum62=("blosum62", "first"))
                  .reset_index())
        rows.append(grp)
    if not rows: return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    out["in_both"] = out.sources.str.contains(r"ClinVar\+gnomAD")
    return out.sort_values(["gene", "position"]).reset_index(drop=True)
