# Explorer data files

Every file is JSON wrapped in a function call, for example `TOL.gene({...});`. The page loads them as `<script>` tags,
so it also works when `index.html` is opened straight from disk. To read one in Python:

```python
import json
text = open("docs/data/genes/GLRA1.js").read()
gene = json.loads(text[text.index("(") + 1 : text.rindex(")")])
```

Rebuild them with `tools/explorer/build_data.py` and `tools/explorer/build_model_data.py`.

## `genes/<GENE>.js` (one per protein, from `build_data.py`)
| Key | Meaning |
|---|---|
| `seq` | MANE Select protein sequence |
| `split` | `train`, `val`, `test` or `excluded` (gene-family split of notebook 05 §10) |
| `dis`, `bind`, `rsa`, `plddt`, `esm`, `cons` | one value per residue, as comma-separated integers: AIUPred disorder, AIUPred binding and RSA (× 100); AlphaFold pLDDT (× 1); ESM-1v mean LLR and conservation (× 10). Empty = missing |
| `ann` | annotation ranges `[[start, end], ...]` (1-based, inclusive): domains, transmembrane segments, ligand sites, interfaces, ... |
| `gn` | gnomAD variants, e.g. `E154Du`: reference, position, alternate, then `u` = used in the score (allele count ≥ 2 and the X rule), `s` = seen once, `x` = left out by the X rule |
| `cv` | ClinVar variants, e.g. `A300Pp1`: reference, position, alternate, class (`p` pathogenic, `b` benign, `v` VUS, `o` other), review stars |
| `check` | sum and sum of absolute values of the official `tolerance_v4` track, used to confirm the browser rebuild |
| `test` | test genes only: `[variant, benign (1/0), model P(benign), ESM-1v LLR, tolerance (own residue hidden)]` |

## Other files
| File | Contents |
|---|---|
| `genes.js` | the gene index: length, split, variant counts |
| `test.js` | the official predictions for the 1,357 test variants, and `final_test/results.json` |
| `summary.js` | result tables shown on the page: split sizes, cross-validation, feature groups, structural classes, Design A, conservation |
| `model.js` (`build_model_data.py`) | held-out histograms (§6g), the final model's weights and intercept |
| `model/<GENE>.js` | for each training or validation gene with labelled variants: each variant's P(benign), z and its 12 largest feature contributions (the final model, fitted on these variants) |
