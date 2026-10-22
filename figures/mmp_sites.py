"""Scan a protein sequence for MMP cleavage-site motifs.

    from mmp_sites import scan_motifs, MMP_MOTIFS
    scan_motifs(seq, MMP_MOTIFS)

Motifs are written as space-separated positions; '/' lists the residues allowed at that
position; 'X' (or '.') means any residue. Positions are P5 P4 P3 P2 P1 | P1' P2' P3', i.e.
the scissile bond falls after position 5 (cut_after=5).
"""
from netutil import pkg_path
import re
import pandas as pd

MMP_MOTIFS = {
    "MMP9": "V V P/V L/Y S L/M/I/Y/F R/T/Y/V/I S/A/G",
    "MMP3": "N/I K/V/I/R P/V/I F/Y/L/M/A S/E M/I/L/Y/F M/K/I/R M/A",
}

def motif_regex(motif: str) -> str:
    """'P/V L/Y S' -> '[PV][LY][S]'"""
    out = []
    for pos in motif.split():
        res = "".join(pos.split("/"))
        out.append("." if res.upper() in ("X", ".") else f"[{res}]")
    return "".join(out)

def scan_motifs(seq, motifs=MMP_MOTIFS, cut_after=5, offset=0, max_mismatch=0) -> pd.DataFrame:
    """Find motif matches in `seq`.
    cut_after: number of motif positions before the scissile bond (P1|P1').
    offset:    added to reported positions, if `seq` is a slice of a larger protein.
    max_mismatch: allow this many positions to disagree (0 = exact motif match).
    Returns start/end (1-based, inclusive), the matched peptide, and the cut position.
    """
    rows = []
    for name, motif in motifs.items():
        pats = [re.compile(p) for p in [motif_regex(motif)]]
        n = len(motif.split())
        if max_mismatch == 0:
            spans = [(m.start(), m.end()) for p in pats
                     for m in re.finditer(f"(?=({p.pattern}))", seq)]
            spans = [(s, s + n) for s, _ in spans]
        else:                                        # sliding window, count position mismatches
            sets = [set(p.split("/")) for p in motif.split()]
            spans = []
            for i in range(len(seq) - n + 1):
                bad = sum(1 for j, s in enumerate(sets)
                          if "X" not in s and "." not in s and seq[i + j] not in s)
                if bad <= max_mismatch:
                    spans.append((i, i + n))
        for s, e in spans:
            pep = seq[s:e]
            rows.append({"enzyme": name, "start": s + 1 + offset, "end": e + offset,
                         "peptide": pep,
                         "cut_between": f"{s + cut_after + offset}|{s + cut_after + 1 + offset}",
                         "P1": pep[cut_after - 1], "P1_prime": pep[cut_after],
                         "mismatches": sum(1 for j, p in enumerate(motif.split())
                                           if p.upper() not in ("X", ".") and pep[j] not in p.split("/"))})
    cols = ["enzyme", "start", "end", "peptide", "cut_between", "P1", "P1_prime", "mismatches"]
    if not rows:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(rows).sort_values(["enzyme", "start"]).reset_index(drop=True)[cols]


# ---------------------------------------------------------------- MEROPS observed cleavages
import os, urllib.request
import numpy as np

MEROPS_IDS = {"MMP3": "M10.005", "MMP9": "M10.004", "ADAM10": "M12.221",
              "MMP2": "M10.003", "MMP7": "M10.008", "ADAM17": "M12.217"}
POSITIONS = ["P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'"]
AAS = "ACDEFGHIKLMNPQRSTVWY"
THREE_TO_ONE = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G",
                "His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S",
                "Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}
# UniProt/Swiss-Prot average amino-acid composition, used as the null model
BACKGROUND = {"A":.0825,"C":.0137,"D":.0545,"E":.0675,"F":.0386,"G":.0707,"H":.0227,"I":.0596,
              "K":.0584,"L":.0966,"M":.0242,"N":.0406,"P":.0470,"Q":.0393,"R":.0553,"S":.0656,
              "T":.0534,"V":.0687,"W":.0108,"Y":.0292}

def fetch_merops_substrates(protease, cache_dir=None, merops_id=None) -> pd.DataFrame:
    """Observed cleavages for one protease from MEROPS, as a table of P4..P4' residues.
    Cached locally so the site is queried once per protease."""
    mid = merops_id or MEROPS_IDS[protease]
    cache_dir = cache_dir or pkg_path("merops_cache")
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, f"{protease}_{mid}.csv")
    if os.path.exists(path):
        return pd.read_csv(path, keep_default_na=False, na_values=[""])
    url = f"https://www.ebi.ac.uk/merops/cgi-bin/substrates?id={mid}"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (research use)"})
    html = urllib.request.urlopen(req, timeout=60).read().decode("utf8", "ignore")
    from io import StringIO
    tables = pd.read_html(StringIO(html))
    df = max(tables, key=len)
    keep = [c for c in df.columns if str(c) in ("Substrate", "Uniprot", "Cleavage type",
                                                "Evidence", *POSITIONS)]
    df = df[keep]
    df.to_csv(path, index=False)
    return df

def build_pssm(df, positions=POSITIONS, pseudocount=0.5, background=BACKGROUND,
               physiological_only=True) -> pd.DataFrame:
    """Log-odds position-specific scoring matrix from observed cleavages.
    physiological_only: keep only rows MEROPS marks as physiological cleavages."""
    if physiological_only and "Cleavage type" in df.columns:
        # MEROPS codes this column P (physiological), S (synthetic), N (non-physiological)
        t = df["Cleavage type"].astype(str).str.strip().str.upper()
        df = df[t.str.startswith("P")]
        if df.empty:
            raise ValueError("no physiological cleavages for this protease; "
                             "use physiological_only=False")
    mat = {}
    for p in positions:
        if p not in df.columns: continue
        col = df[p].astype(str).str.strip()
        counts = {a: pseudocount for a in AAS}
        for v in col:                                  # MEROPS uses three-letter codes
            a = THREE_TO_ONE.get(v.capitalize(), v.upper() if len(v) == 1 else None)
            if a in counts: counts[a] += 1
        total = sum(counts.values())
        mat[p] = {a: np.log2((counts[a] / total) / background[a]) for a in AAS}
    return pd.DataFrame(mat)[list(mat)]

def score_sequence(seq, pssm, cut_after=4) -> pd.DataFrame:
    """Slide the PSSM along `seq`. Returns one row per window with its log-odds score.
    cut_after: motif positions before the scissile bond (4 for a P4..P4' matrix)."""
    n = len(pssm.columns)
    rows = []
    for i in range(len(seq) - n + 1):
        win = seq[i:i + n]
        if any(a not in AAS for a in win): continue
        s = sum(pssm.loc[a, p] for a, p in zip(win, pssm.columns))
        rows.append({"start": i + 1, "end": i + n, "peptide": win, "score": round(s, 2),
                     "cut_between": f"{i + cut_after}|{i + cut_after + 1}",
                     "P1": win[cut_after - 1], "P1_prime": win[cut_after]})
    return pd.DataFrame(rows)

def scan_pssm(seq, proteases=("MMP3", "MMP9", "ADAM10"), offset=0, top=10,
              physiological_only=False, min_score=None) -> pd.DataFrame:
    """Score `seq` against the MEROPS-derived matrix for each protease.
    Returns the top-scoring windows per protease (or everything above min_score)."""
    out = []
    for p in proteases:
        df = fetch_merops_substrates(p)
        pssm = build_pssm(df, physiological_only=physiological_only)
        s = score_sequence(seq, pssm)
        s["protease"] = p
        s["n_observed_cleavages"] = len(df)
        s["start"] += offset; s["end"] += offset
        s["cut_between"] = s.cut_between.map(
            lambda c: "|".join(str(int(x) + offset) for x in c.split("|")))
        s = s.sort_values("score", ascending=False)
        out.append(s[s.score >= min_score] if min_score is not None else s.head(top))
    cols = ["protease", "start", "end", "peptide", "cut_between", "P1", "P1_prime",
            "score", "n_observed_cleavages"]
    return pd.concat(out, ignore_index=True)[cols]


def scan_juxtamembrane(seqs, accessions, window=40, proteases=("MMP3", "MMP9", "ADAM10"),
                       top=None, min_score=None, physiological_only=None,
                       topologies=None) -> pd.DataFrame:
    """Score the `window` residues immediately before the TM helix in each sequence.

    seqs:        {name: sequence}
    accessions:  {name: uniprot accession} -- used for the TM boundary
    window:      residues before the TM start to scan (the scored window is extended by the
                 motif length so sites spanning the boundary are still found)
    top:         keep this many best-scoring windows per (sequence, protease); None = all
    min_score:   alternatively, keep everything at or above this log-odds score
    Returns one row per candidate site with positions in full-protein numbering.
    """
    from conservation import uniprot_topology
    rows = []
    for name, seq in seqs.items():
        topo = (topologies or {}).get(name) or uniprot_topology(accessions[name])
        tms = sorted(topo.get("transmembrane") or [])
        if not tms:
            print(f"  {name}: no transmembrane annotation - skipped"); continue
        tm_start = tms[0][0]
        lo = max(0, tm_start - 1 - window)                  # 0-based slice start
        sub = seq[lo:tm_start - 1]
        if len(sub) < 8: continue
        phys = {p: (physiological_only if physiological_only is not None
                    else p not in ("ADAM10", "ADAM17")) for p in proteases}
        for p in proteases:
            df = fetch_merops_substrates(p)
            pssm = build_pssm(df, physiological_only=phys[p])
            s = score_sequence(sub, pssm)
            if s.empty: continue
            s["start"] += lo; s["end"] += lo
            s["cut_between"] = s.cut_between.map(
                lambda c: "|".join(str(int(x) + lo) for x in c.split("|")))
            s["protein"] = name; s["protease"] = p
            s["tm_start"] = tm_start
            s["dist_to_TM"] = tm_start - s.end          # residues between site end and the TM
            s["physiological_only"] = phys[p]
            s = s.sort_values("score", ascending=False)
            if min_score is not None: s = s[s.score >= min_score]
            elif top is not None:     s = s.head(top)
            rows.append(s)
    cols = ["protein", "protease", "start", "end", "peptide", "cut_between", "P1", "P1_prime",
            "score", "dist_to_TM", "tm_start", "physiological_only"]
    if not rows: return pd.DataFrame(columns=cols)
    out = pd.concat(rows, ignore_index=True)[cols]
    return out.sort_values(["protein", "protease", "score"], ascending=[True, True, False]).reset_index(drop=True)
