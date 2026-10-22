# Missense pathogenicity in postsynaptic proteins from human variation and protein structure

Per-residue **tolerance scores** (smoothed gnomAD variant density minus smoothed ClinVar density) for 1,082
postsynaptic (SynGO) proteins, combined with ESM-1v scores, protein structure (RSA, interfaces, disorder) and
annotations (InterPro, UniProt, MobiDB, ELM, PDBe-KB) in a logistic regression that classifies ClinVar missense
variants as benign or pathogenic.

On 1,357 ClinVar variants in 85 genes held out from all training and model selection (pre-registered test):
AUROC **0.89** for the model vs **0.62** for ESM-1v alone.

## Interactive explorer
**[Open the interactive explorer](https://tlex36.github.io/LehrPDI-tolerance-model/)** to browse the data and the experiments:
the tolerance score of any of the 1,082 proteins (rebuilt in your browser, with a recipe you can change), the held-out
variants, the structural classes, how the model was chosen, the test on untouched genes, and what the model relies on.

The page is in `docs/` and is served by GitHub Pages. To view it without Pages, open `docs/index.html` in a browser
(the charting library loads from the internet). Its data files are rebuilt from notebook 05's saved results with
`tools/explorer/build_data.py` and `tools/explorer/build_model_data.py` (see the top of each script).

## Contents
| Path | What it is |
|---|---|
| `figures/notebooks/05_tolerance_score_modelling.ipynb` | full pipeline: data, tolerance scores, structure and annotation features, models, final test (§16) |
| `figures/notebooks/06_poster_figures.ipynb` | recreates every poster figure (SVG) with explanations |
| `figures/*.py` | helper modules (`convolve.py`, `tolerance.py`, `structure.py`, `evolution_conservation.py`, `struct_render.py`, ...) |
| `figures/results/tolerance_model/final_test/` | pre-registered plan, test predictions and results |
| `figures/results/tolerance_model/poster/` | poster SVGs, PNG previews, captions (`poster_text.md`), figure scripts |
| `figures/results/tolerance_model/report/` | plain-language guide to how the model was trained, plus Part II on the poster: each panel, a presentation script, likely questions and a pre-print checklist (PDF / DOCX / HTML) |
| `tools/poster_layout/` | scripts used to lay out the PowerPoint poster and its references |
| `docs/` | the interactive explorer (GitHub Pages site): `index.html`, `assets/` (page code), `data/` (one file per gene, plus results) |
| `tools/explorer/` | scripts that build the explorer's data files from notebook 05's results |

## Data not included
Large intermediate tables, structure/annotation caches and ESM representations (~3.4 GB) are not in the repository;
notebook 05 rebuilds them from the public sources (gnomAD v4.1, ClinVar, AlphaFold DB, PDBe, InterPro, UniProt,
MobiDB, ELM, Ensembl). Paths in the notebooks are absolute to the original machine and need adjusting.

## Environment
`pip install -r requirements.txt` (Python 3.10).
