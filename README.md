# Missense pathogenicity in postsynaptic proteins from human variation and protein structure

Per-residue **tolerance scores** (smoothed gnomAD variant density minus smoothed ClinVar density) for 1,080
postsynaptic (SynGO) proteins, combined with ESM-1v scores, protein structure (RSA, interfaces, disorder) and
annotations (InterPro, UniProt, MobiDB, ELM, PDBe-KB) in a logistic regression that classifies ClinVar missense
variants as benign or pathogenic.

On 1,357 ClinVar variants in 85 genes held out from all training and model selection (pre-registered test):
AUROC **0.89** for the model vs **0.62** for ESM-1v alone.

## Contents
| Path | What it is |
|---|---|
| `figures/notebooks/05_tolerance_score_modelling.ipynb` | full pipeline: data, tolerance scores, structure and annotation features, models, final test (§16) |
| `figures/notebooks/06_poster_figures.ipynb` | recreates every poster figure (SVG) with explanations |
| `figures/*.py` | helper modules (`convolve.py`, `tolerance.py`, `structure.py`, `evolution_conservation.py`, `struct_render.py`, ...) |
| `figures/results/tolerance_model/final_test/` | pre-registered plan, test predictions and results |
| `figures/results/tolerance_model/poster/` | poster SVGs, PNG previews, captions (`poster_text.md`), figure scripts |
| `figures/results/tolerance_model/report/` | plain-language guide to how the model was trained (PDF / DOCX / HTML) |
| `tools/poster_layout/` | scripts used to lay out the PowerPoint poster and its references |

## Data not included
Large intermediate tables, structure/annotation caches and ESM representations (~3.4 GB) are not in the repository;
notebook 05 rebuilds them from the public sources (gnomAD v4.1, ClinVar, AlphaFold DB, PDBe, InterPro, UniProt,
MobiDB, ELM, Ensembl). Paths in the notebooks are absolute to the original machine and need adjusting.

## Environment
`pip install -r requirements.txt` (Python 3.10).
