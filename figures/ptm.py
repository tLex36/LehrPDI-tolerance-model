"""Post-translational modification sites: UniProt annotations + what is actually modeled in
the deposited structures.

    tab = ptm_table("Q9NZ94", [("7CEE","A"),("8GS3","A")], results=results, gene="NLGN3")
    view_ptms(results, "NLGN3", pdb="7CEE", accession="Q9NZ94").show()
"""
import json
import pandas as pd
import py3Dmol
from netutil import urlopen_retry
from structure import fetch_pdb, sifts_chains, fetch_uniprot_seq
from structure_pipeline import map_all_positions

PTM_TYPES = ("Modified residue", "Glycosylation", "Disulfide bond", "Lipidation", "Cross-link")
PTM_COLORS = {"Glycosylation": "#2ca02c", "Disulfide bond": "#e8c547",
              "Modified residue": "#ff7f0e", "Lipidation": "#9467bd", "Cross-link": "#8c564b"}

def uniprot_ptms(accession) -> pd.DataFrame:
    """Annotated PTM features for one UniProt entry (1-based positions)."""
    j = json.loads(urlopen_retry(f"https://rest.uniprot.org/uniprotkb/{accession}.json"))
    rows = []
    for f in j.get("features", []):
        if f["type"] not in PTM_TYPES: continue
        s, e = f["location"]["start"]["value"], f["location"]["end"]["value"]
        rows.append({"type": f["type"], "start": s, "end": e,
                     "description": f.get("description", ""),
                     "evidence": f.get("evidences", [{}])[0].get("evidenceCode", "")})
    return pd.DataFrame(rows)

def observed_glycans(pdb_id, chain=None, text=None) -> pd.DataFrame:
    """Glycans modeled on a structure, read from LINK records (sugar -> Asn/Ser/Thr)."""
    txt = text or fetch_pdb(pdb_id)
    SUGARS = {"NAG", "BMA", "MAN", "FUC", "GAL", "NDG", "GLC", "XYL", "SIA"}
    rows = []
    for l in txt.split("\n"):
        if not l.startswith("LINK"): continue
        a1, r1, c1, n1 = l[12:16].strip(), l[17:20].strip(), l[21], l[22:26].strip()
        a2, r2, c2, n2 = l[42:46].strip(), l[47:50].strip(), l[51], l[52:56].strip()
        for (ra, aa, ca, na), (rb, ab, cb, nb) in ((( r1,a1,c1,n1),(r2,a2,c2,n2)),
                                                   ((r2,a2,c2,n2),(r1,a1,c1,n1))):
            if rb in SUGARS and ra not in SUGARS:
                if chain and ca != chain: continue
                rows.append({"pdb": pdb_id, "chain": ca, "residue": ra, "pdb_resnum": int(na),
                             "atom": aa, "sugar": rb, "sugar_resnum": int(nb)})
    return pd.DataFrame(rows).drop_duplicates()

def observed_disulfides(pdb_id, chain=None, text=None) -> pd.DataFrame:
    """Disulfide bonds modeled in a structure, from SSBOND records."""
    txt = text or fetch_pdb(pdb_id)
    rows = []
    for l in txt.split("\n"):
        if not l.startswith("SSBOND"): continue
        c1, n1, c2, n2 = l[15], l[17:21].strip(), l[29], l[31:35].strip()
        if chain and chain not in (c1, c2): continue
        rows.append({"pdb": pdb_id, "chain_1": c1, "resnum_1": int(n1),
                     "chain_2": c2, "resnum_2": int(n2), "intra_chain": c1 == c2})
    return pd.DataFrame(rows)

def ptm_table(accession, structures=(), results=None, gene=None) -> pd.DataFrame:
    """Annotated PTM sites, plus whether each one is resolved and modified in each structure.

    structures: [(pdb, chain), ...] -- the chain's own UniProt entry may be another species;
      positions are mapped back to `accession` numbering by alignment.
    results/gene: optional consolidated structure annotation, to add rsa and class.
    """
    ref = fetch_uniprot_seq(accession)
    ann = uniprot_ptms(accession)
    sites = []
    for _, r in ann.iterrows():
        positions = [r.start] if r.start == r.end else [r.start, r.end]
        for p in positions:
            sites.append({"type": r.type, "position": p, "residue": ref[p - 1],
                          "description": r.description or ("S-S partner" if r.type == "Disulfide bond" else ""),
                          "partner": (r.end if p == r.start else r.start) if r.type == "Disulfide bond" else None,
                          "evidence": r.evidence})
    df = pd.DataFrame(sites)

    for pdb, chain in structures:
        txt = fetch_pdb(pdb)
        info = sifts_chains(pdb).loc[chain]
        chain_seq = fetch_uniprot_seq(info.accession)
        to_chain = map_all_positions(ref, chain_seq) if info.accession != accession else {i: i for i in range(1, len(ref) + 1)}
        gly = observed_glycans(pdb, chain, txt)
        ss = observed_disulfides(pdb, chain, txt)
        glyco_sites = set(gly.pdb_resnum) if len(gly) else set()
        ss_sites = set(ss.resnum_1) | set(ss.resnum_2) if len(ss) else set()
        resolved = {int(l[22:26]) for l in txt.split("\n")
                    if l.startswith("ATOM") and l[21] == chain}
        col_r, col_m = f"{pdb}:{chain} resolved", f"{pdb}:{chain} modified"
        df[col_r], df[col_m] = False, ""
        for i, row in df.iterrows():
            p = to_chain.get(row.position)
            if p is None: continue
            df.at[i, col_r] = p in resolved
            marks = []
            if p in glyco_sites: marks.append("glycan")
            if p in ss_sites:    marks.append("S-S")
            df.at[i, col_m] = "+".join(marks)
    if results is not None and gene in results:
        t = results[gene]["table"]
        df["rsa"] = [t.rsa.get(p) for p in df.position]
        df["class"] = [t["class"].get(p) for p in df.position]
    return df.sort_values("position").reset_index(drop=True)

def view_ptms(results, gene, pdb, chain="A", accession=None, ptm_colors=PTM_COLORS,
              cartoon_color="#d9d9d9", show_glycans=True, glycan_color="#2ca02c",
              label=True, stick_radius=0.2, highlight=None, highlight_color="red",
              highlight_sidechain=True, highlight_label=True, width=850, height=650):
    """Cartoon of the chain with annotated PTM sites drawn as coloured side-chain sticks,
    and any glycans modeled in the structure shown as sticks too.
    Sites outside the resolved region (e.g. cytoplasmic-tail phosphosites) are skipped --
    use ptm_table() to see those.
    highlight: extra residues to mark, in `accession` numbering -- a list (all take
      highlight_color) or a {position: colour} dict. Side chains are drawn as sticks."""
    txt = fetch_pdb(pdb)
    info = sifts_chains(pdb).loc[chain]
    ref = fetch_uniprot_seq(accession)
    chain_seq = fetch_uniprot_seq(info.accession)
    to_chain = map_all_positions(ref, chain_seq) if info.accession != accession else {i: i for i in range(1, len(ref) + 1)}

    resolved = {int(l[22:26]) for l in txt.split("\n") if l.startswith("ATOM") and l[21] == chain}

    v = py3Dmol.view(width=width, height=height)
    v.addModel(txt, "pdb")
    v.setStyle({}, {})
    v.setStyle({"chain": chain}, {"cartoon": {"color": cartoon_color}})

    for _, r in uniprot_ptms(accession).iterrows():
        color = ptm_colors.get(r.type, "red")
        for p in ({r.start, r.end} if r.type == "Disulfide bond" else {r.start}):
            q = to_chain.get(p)
            if q is None or q not in resolved:     # site not present in this structure
                continue
            sel = {"chain": chain, "resi": [int(q)]}
            v.addStyle(sel, {"stick": {"color": color, "radius": stick_radius}})
            v.addStyle(sel, {"cartoon": {"color": color}})
            if label:
                v.addLabel(f"{ref[p-1]}{p}", {"backgroundColor": color, "fontColor": "white",
                                              "fontSize": 11, "backgroundOpacity": 0.85}, sel)
    if highlight:
        hl = {int(p): highlight_color for p in highlight} if not isinstance(highlight, dict) else \
             {int(k): v for k, v in highlight.items()}
        for pos, color in hl.items():
            q = to_chain.get(pos)
            if q is None or q not in resolved:
                print(f"  residue {pos} is not resolved in {pdb}:{chain} - skipped")
                continue
            sel = {"chain": chain, "resi": [int(q)]}
            v.addStyle(sel, {"cartoon": {"color": color}})
            if highlight_sidechain:                    # backbone N/C/O dropped so the R group reads
                v.addStyle({**sel, "not": {"atom": ["N", "C", "O"]}},
                           {"stick": {"color": color, "radius": stick_radius * 1.3}})
            if highlight_label:
                v.addLabel(f"{ref[pos-1]}{pos}", {"backgroundColor": color, "fontColor": "white",
                                                  "fontSize": 12, "backgroundOpacity": 0.9}, sel)

    if show_glycans:
        gly = observed_glycans(pdb, chain, txt)
        if len(gly):
            v.addStyle({"chain": chain, "resn": sorted(set(gly.sugar))},
                       {"stick": {"color": glycan_color, "radius": stick_radius}})
    v.zoomTo({"chain": chain})
    return v
