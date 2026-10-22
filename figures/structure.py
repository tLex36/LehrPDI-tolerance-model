"""Core / surface classification from a structure via relative solvent accessibility (RSA).

    df = residue_accessibility("Q8N2Q7")            # AlphaFold model, UniProt numbering
    df = residue_accessibility(pdb_file="3B3Q.pdb", chain="A")
    df["class"]  ->  "core" (RSA<0.2) / "boundary" / "surface" (RSA>0.4)
"""
import io, urllib.request
import pandas as pd
from netutil import urlopen_retry as _urlopen_retry
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import three_to_index, index_to_one

# Max SASA per residue (Å²), Tien et al. 2013 (theoretical)
MAX_ASA = {"A":129,"R":274,"N":195,"D":193,"C":167,"Q":225,"E":223,"G":104,"H":224,"I":197,
           "L":201,"K":236,"M":224,"F":240,"P":159,"S":155,"T":172,"W":285,"Y":263,"V":174}

def fetch_alphafold(uniprot_acc, out=None):
    import json
    meta = json.load(urllib.request.urlopen(f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_acc}"))
    url = next(m["pdbUrl"] for m in meta if m["entryId"] == f"AF-{uniprot_acc}-F1")
    data = _urlopen_retry(url).decode()
    if out: open(out, "w").write(data)
    return data

def fetch_pdb(pdb_id, out=None):
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    data = _urlopen_retry(url).decode()
    if out: open(out, "w").write(data)
    return data

def residue_accessibility(uniprot_acc=None, pdb_file=None, pdb_text=None, chain=None,
                          core=0.20, surface=0.40) -> pd.DataFrame:
    """RSA per residue. Give one of uniprot_acc (AlphaFold), pdb_file, or pdb_text.
    Returns DataFrame: resnum (as numbered in the file), aa, sasa, rsa, class, bfactor
    (bfactor = pLDDT for AlphaFold models)."""
    if pdb_text is None:
        pdb_text = fetch_alphafold(uniprot_acc) if uniprot_acc else open(pdb_file).read()
    struct = PDBParser(QUIET=True).get_structure("s", io.StringIO(pdb_text))
    model = struct[0]
    if chain:                                       # drop other chains so SASA is for this chain alone
        for ch in list(model):
            if ch.id != chain: model.detach_child(ch.id)
    ShrakeRupley().compute(model, level="R")
    rows = []
    for ch in model:
        for res in ch:
            if res.id[0] != " ": continue           # skip waters/ligands
            try: aa = index_to_one(three_to_index(res.get_resname()))
            except KeyError: continue
            rsa = res.sasa / MAX_ASA[aa]
            rows.append({"chain": ch.id, "resnum": res.id[1], "aa": aa, "sasa": round(res.sasa, 1),
                         "rsa": round(rsa, 3),
                         "class": "core" if rsa < core else "surface" if rsa > surface else "boundary",
                         "bfactor": round(next(iter(res)).get_bfactor(), 1)})
    return pd.DataFrame(rows)

def interface_residues(complex_pdb_text, chain, partner_chains, min_drop=1.0) -> pd.DataFrame:
    """Residues of `chain` whose SASA drops by >= min_drop Å² when partner_chains are present."""
    alone = residue_accessibility(pdb_text=complex_pdb_text, chain=chain).set_index("resnum")
    struct = PDBParser(QUIET=True).get_structure("s", io.StringIO(complex_pdb_text))
    model = struct[0]
    keep = {chain, *partner_chains}
    for ch in list(model):
        if ch.id not in keep: model.detach_child(ch.id)
    ShrakeRupley().compute(model, level="R")
    bound = {res.id[1]: res.sasa for res in model[chain] if res.id[0] == " "}
    alone["sasa_bound"] = pd.Series(bound)
    alone["buried"] = (alone["sasa"] - alone["sasa_bound"]).round(1)
    return alone[alone["buried"] >= min_drop].reset_index()

# ---------------------------------------------------------------- crystal-structure annotation
import copy, json as _json
from Bio.Align import PairwiseAligner, substitution_matrices


def sifts_chains(pdb_id) -> pd.DataFrame:
    """Chain -> UniProt mapping from PDBe SIFTS: accession, entry name, species, UniProt range."""
    d = _json.loads(_urlopen_retry(f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"))
    rows = []
    for acc, v in d[pdb_id.lower()]["UniProt"].items():
        for m in v["mappings"]:
            rows.append({"chain": m["chain_id"], "accession": acc, "entry": v["identifier"],
                         "species": v["identifier"].split("_")[-1], "unp_start": m["unp_start"], "unp_end": m["unp_end"]})
    return pd.DataFrame(rows).drop_duplicates("chain").set_index("chain")

def fetch_uniprot_seq(acc):
    fa = _urlopen_retry(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta").decode()
    return "".join(fa.split("\n")[1:]).strip()

def _chain_residues(chain):
    """Ordered list of (pdb_resnum, one-letter aa, residue) for standard residues."""
    out = []
    for res in chain:
        if res.id[0] != " ": continue
        try: aa = index_to_one(three_to_index(res.get_resname()))
        except KeyError: continue
        out.append((res.id[1], aa, res))
    return out

def map_to_uniprot(observed, unp_seq):
    """observed: list of (pdb_resnum, aa). Returns {pdb_resnum: unp_resnum (1-based)} by aligning
    the observed chain sequence to the full UniProt sequence (handles offsets, missing loops, tags)."""
    a = PairwiseAligner(); a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score, a.extend_gap_score = -10, -0.5; a.mode = "local"
    obs_seq = "".join(aa for _, aa in observed)
    aln = a.align(obs_seq, unp_seq)[0]
    mapping = {}
    for (o0, o1), (u0, u1) in zip(*aln.aligned):          # aligned blocks (0-based, half-open)
        for k in range(o1 - o0):
            mapping[observed[o0 + k][0]] = u0 + k + 1
    return mapping

def _sasa_of(model, keep_chains, strip_hetero=True):
    m = copy.deepcopy(model)
    for ch in list(m):
        if ch.id not in keep_chains: m.detach_child(ch.id); continue
        if strip_hetero:                                   # waters, glycans, ligands occlude the surface
            for res in list(ch):
                if res.id[0] != " ": ch.detach_child(res.id)
    ShrakeRupley().compute(m, level="R")
    return {ch.id: {res.id[1]: res.sasa for res in ch if res.id[0] == " "} for ch in m}

def to_ranges(nums):
    """[1,2,3,7,8,10] -> [(1,3),(7,8),(10,10)]"""
    nums = sorted(set(nums)); out = []
    for n in nums:
        if out and n == out[-1][1] + 1: out[-1] = (out[-1][0], n)
        else: out.append((n, n))
    return out

def annotate_pdb(pdb_id, chain, partners=None, core=0.20, surface=0.40, min_buried=1.0):
    """Classify residues of `chain` in crystal structure `pdb_id`.

    partners: {label: [chain ids]} e.g. {"dimer": ["B"], "NRXN1B": ["E"]} -> for each label,
      residues whose SASA drops >= min_buried Å² when those chains are present are flagged.

    Returns dict with:
      accession, entry, species   which UniProt entry this chain is (use to fetch the sequence
                                  for the MSA step)
      sequence                    full UniProt sequence of that entry
      residues                    DataFrame: unp_resnum, pdb_resnum, aa, rsa, class,
                                  buried_<label>, iface_<label>
      ranges                      {"core": [(s,e),..], "boundary": .., "surface": ..,
                                   "iface_<label>": ..}   (UniProt numbering, 1-based inclusive)
    """
    partners = partners or {}
    info = sifts_chains(pdb_id).loc[chain]
    unp_seq = fetch_uniprot_seq(info.accession)
    model = PDBParser(QUIET=True).get_structure(pdb_id, io.StringIO(fetch_pdb(pdb_id)))[0]
    observed = [(n, aa) for n, aa, _ in _chain_residues(model[chain])]
    p2u = map_to_uniprot(observed, unp_seq)

    alone = _sasa_of(model, {chain})[chain]
    bound = {lab: _sasa_of(model, {chain, *chs})[chain] for lab, chs in partners.items()}

    rows = []
    for n, aa in observed:
        if n not in p2u: continue
        rsa = alone[n] / MAX_ASA[aa]
        r = {"unp_resnum": p2u[n], "pdb_resnum": n, "aa": aa, "sasa": round(alone[n], 1), "rsa": round(rsa, 3),
             "class": "core" if rsa < core else "surface" if rsa > surface else "boundary"}
        for lab in partners:
            d = alone[n] - bound[lab][n]
            r[f"buried_{lab}"] = round(d, 1); r[f"iface_{lab}"] = d >= min_buried
        rows.append(r)
    df = pd.DataFrame(rows)
    # sanity: mapped residues must match the UniProt sequence
    mism = df[df.apply(lambda r: unp_seq[r.unp_resnum - 1] != r.aa, axis=1)]
    if len(mism): print(f"warning: {len(mism)} residues differ from UniProt (mutations/engineered?): "
                        f"{mism.unp_resnum.tolist()[:10]}")
    ranges = {c: to_ranges(df.unp_resnum[df["class"] == c]) for c in ("core", "boundary", "surface")}
    for lab in partners:
        ranges[f"iface_{lab}"] = to_ranges(df.unp_resnum[df[f"iface_{lab}"]])
    return {"pdb": pdb_id, "chain": chain, "accession": info.accession, "entry": info.entry,
            "species": info.species, "sequence": unp_seq, "residues": df, "ranges": ranges}

def find_contacts(pdb_id, chain, min_buried=1.0, min_residues=5) -> pd.DataFrame:
    """Which other chains bury surface of `chain`? Returns entry/species per contacting chain
    and how many residues of `chain` they bury -- use it to choose `partners` for annotate_pdb."""
    info = sifts_chains(pdb_id)
    model = PDBParser(QUIET=True).get_structure(pdb_id, io.StringIO(fetch_pdb(pdb_id)))[0]
    alone = _sasa_of(model, {chain})[chain]
    rows = []
    for other in model:
        if other.id == chain: continue
        bound = _sasa_of(model, {chain, other.id})[chain]
        n = sum((alone[r] - bound[r]) >= min_buried for r in alone)
        if n >= min_residues:
            e = info.loc[other.id] if other.id in info.index else None
            rows.append({"partner_chain": other.id, "entry": e.entry if e is not None else "?",
                         "buried_residues": n})
    return pd.DataFrame(rows).sort_values("buried_residues", ascending=False)
