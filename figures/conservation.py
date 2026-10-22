"""Conservation pipeline: fetch -> select -> compare.

    seqs   = fetch_uniprot({"Human": "Q8N2Q7", "Mouse": "Q99K10"})   # {label: full sequence}
    subset = {k: v[695:] for k, v in seqs.items()}                    # slice to a domain if wanted
    result = compare(subset)                                          # scores + printed MSA
"""
import urllib.request
from netutil import urlopen_retry as _urlopen_retry, pkg_path
from itertools import combinations
import pandas as pd
from msa_print import align, format_clustal, conservation_char

# ---------------------------------------------------------------- fetch
def fetch_uniprot(accessions: dict) -> dict:
    """{label: accession} -> {label: sequence string}."""
    seqs = {}
    for label, acc in accessions.items():
        fa = _urlopen_retry(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta").decode()
        seqs[label] = "".join(fa.split("\n")[1:]).strip()
    return seqs

# ---------------------------------------------------------------- score
def pairwise_identity(aln) -> pd.DataFrame:
    """% identical residues per pair, over columns where neither has a gap."""
    names = [r.id for r in aln]
    mat = pd.DataFrame(100.0, index=names, columns=names)
    for a, b in combinations(aln, 2):
        pairs = [(x, y) for x, y in zip(a.seq, b.seq) if x != "-" and y != "-"]
        pid = 100 * sum(x == y for x, y in pairs) / len(pairs) if pairs else float("nan")
        mat.loc[a.id, b.id] = mat.loc[b.id, a.id] = pid
    return mat

def column_conservation(aln) -> dict:
    """Fraction of alignment columns that are identical / similar across ALL sequences."""
    L = aln.get_alignment_length()
    symbols = [conservation_char(aln[:, i]) for i in range(L)]
    n_ident = symbols.count("*")
    n_sim = n_ident + symbols.count(":") + symbols.count(".")
    n_ungapped = sum("-" not in aln[:, i] for i in range(L))
    return {
        "columns": L,
        "ungapped_columns": n_ungapped,
        "pct_identical": 100 * n_ident / L,
        "pct_similar": 100 * n_sim / L,                      # identical + strong + weak
        "pct_identical_of_ungapped": 100 * n_ident / n_ungapped if n_ungapped else float("nan"),
    }

# ---------------------------------------------------------------- run
def compare(seqs: dict, width=60, offsets=None, show=True) -> dict:
    """Align the given sequences (or substrings) and report conservation."""
    aln = align(seqs)
    text = format_clustal(aln, width=width, offsets=offsets)
    out = {
        "alignment": aln,
        "text": text,
        "identity_matrix": pairwise_identity(aln).round(1),
        "overall": column_conservation(aln),
    }
    if show:
        o = out["overall"]
        print(f"{o['columns']} alignment columns | identical: {o['pct_identical']:.1f}% | "
              f"similar: {o['pct_similar']:.1f}%\n")
        print(out["identity_matrix"], "\n")
        print(text)
    return out

# ---------------------------------------------------------------- true pairwise scores
from Bio.Align import PairwiseAligner, substitution_matrices
_BLOSUM62 = substitution_matrices.load("BLOSUM62")

def _needle_aligner():
    """Global alignment with EMBOSS Needle's defaults: BLOSUM62, gap open 10, extend 0.5."""
    a = PairwiseAligner()
    a.substitution_matrix = _BLOSUM62
    a.open_gap_score, a.extend_gap_score = -10, -0.5
    a.mode = "global"
    return a

def pairwise_scores(seqs: dict, count_gaps=True) -> dict:
    """Align every pair on its own and score it (independent of any other sequence).

    Returns {"identity": DataFrame, "similarity": DataFrame, "score": DataFrame}.
      identity   = 100 * identical positions / length
      similarity = 100 * positions whose BLOSUM62 score > 0 / length   (EMBOSS convention:
                   identical + conservative substitutions)
      length     = full alignment length incl. gap columns if count_gaps else only
                   columns where both have a residue
    """
    aligner = _needle_aligner()
    names = list(seqs)
    ident = pd.DataFrame(100.0, index=names, columns=names)
    simil = pd.DataFrame(100.0, index=names, columns=names)
    score = pd.DataFrame(float("nan"), index=names, columns=names)
    for a, b in combinations(names, 2):
        aln = aligner.align(seqs[a], seqs[b])[0]
        s1, s2 = aln[0], aln[1]                      # gapped strings
        n_id = n_sim = n_res = 0
        for x, y in zip(s1, s2):
            if x == "-" or y == "-":
                continue
            n_res += 1
            if x == y:
                n_id += 1; n_sim += 1
            elif _BLOSUM62[x][y] > 0:
                n_sim += 1
        length = len(s1) if count_gaps else n_res
        ident.loc[a, b] = ident.loc[b, a] = 100 * n_id / length
        simil.loc[a, b] = simil.loc[b, a] = 100 * n_sim / length
        score.loc[a, b] = score.loc[b, a] = aln.score
    return {"identity": ident.round(1), "similarity": simil.round(1), "score": score}

# ---------------------------------------------------------------- species x gene vs human
def species_vs_reference(gene_seqs: dict, reference="Human", metric="identity",
                         count_gaps=True, collapse=True) -> pd.DataFrame:
    """gene_seqs: {gene: {species: seq}} (already sliced to the region of interest).
    Returns DataFrame rows=species, cols=genes, cell = % `metric` of species vs `reference`
    for that gene (direct pairwise alignment). Missing species/gene -> NaN.
    collapse: fold paralog labels like "Zebrafish (nlgn2a)" into "Zebrafish", keeping the
    higher-scoring paralog."""
    aligner = _needle_aligner()
    rows = {}
    for gene, seqs in gene_seqs.items():
        ref = seqs.get(reference)
        if ref is None:
            continue
        for sp, seq in seqs.items():
            aln = aligner.align(ref, seq)[0]
            s1, s2 = aln[0], aln[1]
            n_hit = n_res = 0
            for x, y in zip(s1, s2):
                if x == "-" or y == "-": continue
                n_res += 1
                if x == y or (metric == "similarity" and _BLOSUM62[x][y] > 0):
                    n_hit += 1
            length = len(s1) if count_gaps else n_res
            label = sp.split(" (")[0] if collapse else sp
            val = 100 * n_hit / length
            rows.setdefault(label, {})[gene] = max(val, rows.get(label, {}).get(gene, -1))
    return pd.DataFrame(rows).T.round(1)

# ---------------------------------------------------------------- signal peptide / mature ECD
import json as _json

def uniprot_signal_end(accession) -> int | None:
    """Last residue (1-based) of the UniProt 'Signal' feature, or None if not annotated."""
    j = _json.loads(_urlopen_retry(f"https://rest.uniprot.org/uniprotkb/{accession}.json"))
    for f in j.get("features", []):
        if f["type"] == "Signal":
            return f["location"]["end"]["value"]
    return None

def map_position(ref_seq, seq, ref_pos) -> int | None:
    """Transfer a 0-based index in ref_seq to the equivalent index in seq via global alignment.
    Returns None if ref_pos falls opposite a gap (no equivalent residue)."""
    aln = _needle_aligner().align(ref_seq, seq)[0]
    r, s = aln[0], aln[1]
    i_ref = i_seq = 0
    for x, y in zip(r, s):
        if x != "-" and i_ref == ref_pos:
            return i_seq if y != "-" else None
        if x != "-": i_ref += 1
        if y != "-": i_seq += 1
    return None

def mature_start_indices(seqs: dict, ref="Human", ref_signal_end=None, ref_accession=None) -> dict:
    """{species: 0-based index of the first mature residue}, transferred from the reference.
    Give ref_signal_end (1-based last signal residue) or ref_accession to look it up."""
    if ref_signal_end is None:
        ref_signal_end = uniprot_signal_end(ref_accession)
    ref_seq = seqs[ref]
    out = {}
    for sp, s in seqs.items():
        idx = ref_signal_end if sp == ref else map_position(ref_seq, s, ref_signal_end)
        out[sp] = idx
    return out

# ---------------------------------------------------------------- per-residue conservation profile
def conservation_profile(seqs: dict, ref="Human", window=50, sliding=False) -> pd.DataFrame:
    """Per-residue conservation along the reference sequence, from the MSA of `seqs`.
    Returns a DataFrame indexed by 0-based position in the ref string with columns:
      residue     ref amino acid
      identical   1 if the column is identical across all sequences
      similar     1 if identical or all in one Clustal strong/weak group
      win_identical / win_similar   windowed fraction (blocked, or centred-sliding if sliding=True)
    """
    aln = align(seqs)
    r = next(x for x in aln if x.id == ref)
    L = aln.get_alignment_length()
    sym = [conservation_char(aln[:, c]) for c in range(L)]
    rows = []
    for c, ch in enumerate(str(r.seq)):
        if ch == "-": continue
        rows.append({"residue": ch, "identical": int(sym[c] == "*"), "similar": int(sym[c] in "*:.")})
    df = pd.DataFrame(rows)
    for k in ("identical", "similar"):
        if sliding:
            df[f"win_{k}"] = df[k].rolling(window, center=True, min_periods=1).mean()
        else:
            df[f"win_{k}"] = df.groupby(df.index // window)[k].transform("mean")
    return df

# ---------------------------------------------------------------- topology-based domains
def uniprot_topology(accession) -> dict:
    """Domain boundaries from UniProt's curated topology features.

    Returns {"signal": (s,e) or None, "extracellular": [(s,e),...], "transmembrane": [...],
             "cytoplasmic": [...], "ctd": (s,e) or None}  -- 1-based inclusive.
    `ctd` is the final cytoplasmic stretch (the C-terminal domain of a multi-pass receptor);
    for a single-pass protein it is simply the intracellular tail.
    """
    j = _json.loads(_urlopen_retry(f"https://rest.uniprot.org/uniprotkb/{accession}.json"))
    out = {"signal": None, "extracellular": [], "transmembrane": [], "cytoplasmic": [], "ctd": None}
    for f in j.get("features", []):
        s, e = f["location"]["start"]["value"], f["location"]["end"]["value"]
        t, d = f["type"], f.get("description", "")
        if t == "Signal":                     out["signal"] = (s, e)
        elif t in ("Transmembrane", "Intramembrane"): out["transmembrane"].append((s, e))
        elif t == "Topological domain":
            if d.startswith("Extra"):  out["extracellular"].append((s, e))
            elif d.startswith("Cyto"): out["cytoplasmic"].append((s, e))
    if out["cytoplasmic"]:
        out["ctd"] = max(out["cytoplasmic"], key=lambda r: r[0])
    out["length"] = j["sequence"]["length"]
    return out

def domain_slices(seq, topo, mature=True):
    """{'ECD': str, 'TM': str, 'CTD': str} for one sequence given its own topology dict.
    ECD = all extracellular stretches concatenated; TM = all membrane segments; CTD = final
    cytoplasmic stretch. mature=True drops the signal peptide from the ECD."""
    def cat(ranges):
        return "".join(seq[s - 1:e] for s, e in sorted(ranges))
    ecd = list(topo["extracellular"])
    if mature and topo["signal"]:
        ecd = [(max(s, topo["signal"][1] + 1), e) for s, e in ecd if e > topo["signal"][1]]
    return {"ECD": cat(ecd), "TM": cat(topo["transmembrane"]),
            "CTD": seq[topo["ctd"][0] - 1:topo["ctd"][1]] if topo["ctd"] else ""}

def transfer_domain_slices(seqs: dict, ref, ref_topo, mature=True) -> dict:
    """{domain: {label: subsequence}} for every sequence in `seqs`, using the reference's
    topology transferred by global alignment. Labels whose boundaries cannot be mapped are
    dropped from that domain."""
    ref_seq = seqs[ref]
    out = {"ECD": {}, "TM": {}, "CTD": {}}
    for lab, s in seqs.items():
        if lab == ref:
            d = domain_slices(s, ref_topo, mature)
        else:
            m = map_position                     # reuse the single-position mapper
            idx = {}
            topo = {"extracellular": [], "transmembrane": [], "cytoplasmic": [], "signal": None, "ctd": None}
            def tr(r):
                a = map_position(ref_seq, s, r[0] - 1); b = map_position(ref_seq, s, r[1] - 1)
                return (a + 1, b + 1) if a is not None and b is not None else None
            topo["extracellular"] = [x for x in (tr(r) for r in ref_topo["extracellular"]) if x]
            topo["transmembrane"] = [x for x in (tr(r) for r in ref_topo["transmembrane"]) if x]
            if ref_topo["signal"]: topo["signal"] = tr(ref_topo["signal"])
            if ref_topo["ctd"]:
                c = tr((ref_topo["ctd"][0], ref_topo["ctd"][0]))      # map the start only
                topo["ctd"] = (c[0], len(s)) if c else None           # ...tail runs to the end
            d = domain_slices(s, topo, mature)
        for k, v in d.items():
            if v: out[k][lab] = v
    return out

# ---------------------------------------------------------------- iGluR domain architecture
INTERPRO_LIG_CHAN = "IPR001320"        # "Ionotropic glutamate receptor, C-terminal" = S1..S2

def interpro_domains(accession) -> dict:
    """{InterPro accession: [(start, end), ...]} for one UniProt entry."""
    d = _json.loads(_urlopen_retry(
        f"https://www.ebi.ac.uk/interpro/api/entry/InterPro/protein/uniprot/{accession}/?page_size=100"))
    out = {}
    for r in d["results"]:
        acc = r["metadata"]["accession"]
        for e in r.get("proteins", []):
            for loc in e.get("entry_protein_locations", []):
                out.setdefault(acc, []).extend((f["start"], f["end"]) for f in loc["fragments"])
    return out

def receptor_domains(accession, topo=None) -> dict:
    """Domain boundaries for an ionotropic glutamate receptor subunit, 1-based inclusive.

    Combines UniProt topology (authoritative for the membrane segments and the CTD) with the
    InterPro Lig_chan domain start, which marks where the ATD ends and S1 begins.

    Returns an ordered dict: signal, ATD, S1, M1, M2 (pore loop), M3, S2, M4, CTD.
    Membrane helices are numbered in sequence order: the Intramembrane feature is M2.
    """
    topo = topo or uniprot_topology(accession)
    L = topo["length"]
    tms = sorted(topo["transmembrane"])                     # M1, M2(pore), M3, M4 in order
    ipr = interpro_domains(accession)
    s1_start = min(s for s, _ in ipr.get(INTERPRO_LIG_CHAN, [(431, 0)]))
    mature = (topo["signal"][1] + 1) if topo["signal"] else 1

    d = {}
    if topo["signal"]: d["signal"] = topo["signal"]
    d["ATD"] = (mature, s1_start - 1)
    d["S1"]  = (s1_start, tms[0][0] - 1)
    names = ["M1", "M2", "M3", "M4"]
    for n, (s, e) in zip(names, tms):
        d[n] = (s, e)
    d["S2"]  = (tms[2][1] + 1, tms[3][0] - 1)               # between M3 and M4
    d["CTD"] = topo["ctd"] if topo["ctd"] else (tms[3][1] + 1, L)
    return d

def domain_table(accessions: dict, labels=None) -> pd.DataFrame:
    """Tidy table of receptor domain boundaries for several entries.
    accessions: {name: uniprot accession}."""
    rows = []
    for name, acc in accessions.items():
        d = receptor_domains(acc)
        for dom, (s, e) in d.items():
            rows.append({"protein": name, "accession": acc, "domain": dom,
                         "start": s, "end": e, "length": e - s + 1})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- general feature retrieval
# grouped so a table can be filtered to the level of detail you want
FEATURE_GROUPS = {
    "processing":   ["Signal", "Transit peptide", "Propeptide", "Chain", "Peptide"],
    "topology":     ["Topological domain", "Transmembrane", "Intramembrane"],
    "domain":       ["Domain", "Repeat", "Zinc finger", "DNA binding", "Coiled coil"],
    "region":       ["Region", "Motif", "Compositional bias"],
    "site":         ["Active site", "Binding site", "Site", "Metal binding"],
    "ptm":          ["Modified residue", "Glycosylation", "Disulfide bond", "Lipidation",
                     "Cross-link"],
    "variant":      ["Natural variant", "Mutagenesis", "Sequence conflict", "Alternative sequence"],
    "secondary":    ["Helix", "Beta strand", "Turn"],
}
DEFAULT_GROUPS = ("processing", "topology", "domain", "region", "site")

def uniprot_features(accession, groups=DEFAULT_GROUPS, types=None,
                     include_secondary=False) -> pd.DataFrame:
    """Every annotated feature of a UniProt entry as a tidy table.

    groups: which FEATURE_GROUPS to keep ('all' for everything)
    types:  explicit list of feature types, overriding `groups`
    Returns: accession, entry, group, type, start, end, length, description, evidence, note
    """
    j = _json.loads(_urlopen_retry(f"https://rest.uniprot.org/uniprotkb/{accession}.json"))
    entry = j.get("uniProtkbId", accession)
    wanted = None
    if types is not None:
        wanted = set(types)
    elif groups != "all":
        wanted = {t for g in groups for t in FEATURE_GROUPS.get(g, [])}
        if include_secondary:
            wanted |= set(FEATURE_GROUPS["secondary"])
    type_group = {t: g for g, ts in FEATURE_GROUPS.items() for t in ts}
    rows = []
    for f in j.get("features", []):
        t = f["type"]
        if wanted is not None and t not in wanted:
            continue
        loc = f["location"]
        s, e = loc["start"].get("value"), loc["end"].get("value")
        rows.append({"accession": accession, "entry": entry,
                     "group": type_group.get(t, "other"), "type": t,
                     "start": s, "end": e,
                     "length": (e - s + 1) if (s and e) else None,
                     "description": f.get("description", ""),
                     "evidence": ";".join(sorted({ev.get("evidenceCode", "")
                                                  for ev in f.get("evidences", [])})),
                     "note": (f.get("featureCrossReferences") or [{}])[0].get("id", "")})
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values(["start", "end"], na_position="last").reset_index(drop=True)
    return df

def feature_table(accessions, groups=DEFAULT_GROUPS, types=None, sequence=False,
                  include_secondary=False) -> pd.DataFrame:
    """uniprot_features() for several entries at once.
    accessions: {name: accession} or a list of accessions.
    sequence=True adds the residues each feature spans (truncated for long spans)."""
    if not isinstance(accessions, dict):
        accessions = {a: a for a in accessions}
    out = []
    for name, acc in accessions.items():
        df = uniprot_features(acc, groups=groups, types=types,
                              include_secondary=include_secondary)
        df.insert(0, "name", name)
        if sequence and len(df):
            seq = fetch_uniprot({"s": acc})["s"]
            df["residues"] = [
                (seq[int(s) - 1:int(e)] if e - s < 40 else
                 f"{seq[int(s)-1:int(s)+14]}...{seq[int(e)-15:int(e)]}")
                if pd.notna(s) and pd.notna(e) else ""
                for s, e in zip(df.start, df.end)]
        out.append(df)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()

def domain_tracks_from_features(accession, groups=("topology", "domain"), keep_types=None,
                                label_from="description") -> dict:
    """{label: (start, end)} ready to pass to exons.exon_alignment_html(domains=...).
    Labels come from the feature description, falling back to its type."""
    df = uniprot_features(accession, groups=groups, types=keep_types)
    out, seen = {}, {}
    for _, r in df.iterrows():
        if pd.isna(r.start) or pd.isna(r.end): continue
        lab = (r[label_from] or r["type"]).strip() or r["type"]
        seen[lab] = seen.get(lab, 0) + 1
        if seen[lab] > 1: lab = f"{lab} {seen[lab]}"
        out[lab] = (int(r.start), int(r.end))
    return out


# ---------------------------------------------------------------- MANE lookup by gene symbol
MANE_CSV = "../mane_transcripts.csv"        # relative to the package; override if it moves

def mane_lookup(symbols, csv_path=None) -> pd.DataFrame:
    """MANE Select transcript / protein identifiers for one or more gene symbols."""
    csv_path = csv_path or pkg_path(MANE_CSV)
    df = pd.read_csv(csv_path, dtype=str)
    df.columns = [c.lstrip("#") for c in df.columns]
    if isinstance(symbols, str): symbols = [symbols]
    hit = df[df.symbol.isin(symbols) & df.MANE_status.str.contains("Select", na=False)]
    cols = ["symbol", "Ensembl_nuc", "Ensembl_prot", "RefSeq_nuc", "RefSeq_prot",
            "Ensembl_Gene", "GRCh38_chr", "chr_start", "chr_end", "chr_strand"]
    return hit[[c for c in cols if c in hit.columns]].reset_index(drop=True)

def mane_protein(symbols, csv_path=None) -> dict:
    """{gene symbol: MANE Select protein sequence}, fetched from Ensembl."""
    pass
    out = {}
    for _, r in mane_lookup(symbols, csv_path).iterrows():
        pid = r.Ensembl_prot.split(".")[0]
        fa = _urlopen_retry(f"https://rest.ensembl.org/sequence/id/{pid}"
                            f"?type=protein;content-type=text/x-fasta").decode()
        out[r.symbol] = "".join(fa.split("\n")[1:]).strip().rstrip("*")
    return out

def mane_cds(symbols, csv_path=None) -> dict:
    """{gene symbol: MANE Select coding sequence}."""
    pass
    out = {}
    for _, r in mane_lookup(symbols, csv_path).iterrows():
        tid = r.Ensembl_nuc.split(".")[0]
        fa = _urlopen_retry(f"https://rest.ensembl.org/sequence/id/{tid}"
                            f"?type=cds;content-type=text/x-fasta").decode()
        out[r.symbol] = "".join(fa.split("\n")[1:]).strip()
    return out
