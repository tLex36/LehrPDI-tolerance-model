# Poster: figures and copy-ready text

All figures are in `svg/` (editable: text stays text, font Arial). `png_preview/` is for viewing only, not for printing.
They were drawn for a column about 12–15 in wide; scale them proportionally (never stretch one axis).
The colours are the same in every figure:
- blue = gnomAD, benign, or the final model
- red = ClinVar pathogenic
- yellow = VUS
- grey = ESM-1v alone

Regenerate everything by running `figures/notebooks/06_poster_figures.ipynb` top to bottom in the LehrPDI environment (~4 min; the Ensembl conservation data are cached in `conservation_cache/`). `poster_figs.py` still recreates the original 11 figures (I1, M1–M4, R1–R5) without a notebook (~1.5 min).

---------------------------------------------------------------------------------------------------------

## Where each figure goes

| Poster section | File | Shows |
|---|---|---|
| Introduction | `I0_conservation_phylogeny.svg` | species tree + heatmap: how conserved postsynaptic, olfactory and ribosomal proteins, regulatory DNA and the whole genome are, from chimpanzee to *C. elegans* |
| Introduction | `I1_esm_by_variant_group.svg` | ESM-1v separates gnomAD from pathogenic variants, but rates VUS like gnomAD variants |
| Methods 1: tolerance | `M0_tolerance_diagram.svg` | schematic: severity-weighted lollipops, wide vs narrow rolling windows, the three aligned heatmaps (simulated protein) |
| Methods 1: tolerance | `M1_tolerance_score_build.svg` | how the score is built, on one gene (GLRA1) |
| Methods 1: tolerance | `M2_tolerance_by_variant_group.svg` | mean tolerance of held-out gnomAD, benign, VUS and pathogenic variants |
| Methods 2: structure | `M3a_glra1_structure_overview.svg` | GLRA1 AlphaFold model: structured vs unstructured, with AIUPred and pLDDT along the sequence |
| Methods 2: structure | `M3b_core_vs_surface.svg` | cross-section of the extracellular domain where glycine binds: core vs surface residues, and % pathogenic by class |
| Methods 2: structure | `M3c_interfaces.svg` | the same cross-section with both neighbours: the glycine pockets between subunits, interface residues, and % pathogenic |
| Methods 2: structure | `M3_structural_feature_tracks.svg` | disorder, RSA (surface vs core), interfaces, domains / TM / ligand sites along GLRA1 |
| Methods 2: structure | `M3d_annotation_sources.svg` | table: each database, what it looks for, what it means for the protein, the model features it feeds |
| Methods 3: assembling features | `M4_residue_to_variant_row.svg` | Figure 8 of the guide: per-residue data → one row per variant (ESM, conservation, ...) |
| Results 1 | `R1_model_auroc_table.svg` | every model's AUROC (table) |
| Results 2 | `R2_logistic_regression_curve.svg` + `R2b_roc_untouched_genes.svg` | how the model works (Figure 21) + that it beats ESM-1v on untouched genes |
| Results 3 | `R3_untouched_genes_vs_esm.svg` | the pre-registered test on untouched genes, with 95% CIs |
| Results 4 (optional) | `R4_feature_groups.svg`, `R5_model_weights.svg` | what each data type adds (Figure 23); what the model relies on (Figure 26) |

Figure 21 shows *how* the model scores variants and that test variants spread out along its curve. It does not compare with ESM.
To show the model is *better*, put R2b (ROC: 0.89 vs 0.62) or R3 next to it.

---------------------------------------------------------------------------------------------------------

## Captions and text (copy-ready)

### I0: Conservation across evolution
**Conservation of human sequence across evolution.**
- **Left:** % similarity of human proteins to their orthologues, i.e. identical or conservative amino acids (positive BLOSUM62). The bold number is the median over genes with an orthologue; the small text is the share of genes that have one. Three gene sets: postsynaptic proteins (1,080 genes), olfactory receptors (399) and ribosomal proteins (86).
- **Right:** % of human bases identical in Ensembl pairwise genome alignments, for 80 non-coding regulatory elements and 80 protein-coding exons.
- **Hatched:** too few orthologues, or no genome alignment available.

Postsynaptic proteins are highly conserved across vertebrates and ribosomal proteins across all animals. Olfactory receptors and regulatory DNA are conserved only among mammals, and coding DNA is far more conserved than regulatory DNA. Data: Ensembl 116.

Key numbers (similarity / identity):

| | Mouse | Chicken | Zebrafish | Fly | Worm |
|---|---|---|---|---|---|
| Postsynaptic protein | 97% | 90% | 81% | 59% (55% of genes have an orthologue) | 51% (53%) |
| Ribosomal protein | 100% | 98% | 96% | 84% | 80% |
| Coding DNA | 82% | 51% | 46% | no alignment | no alignment |
| Regulatory DNA | 33% | 3.5% | 0.2% | no alignment | no alignment |

The full explanation, with caveats, is in notebook 06 §3.

### I1: ESM-1v by variant group
**ESM-1v scores separate healthy-population from pathogenic variants.**
- Mean log-likelihood ratio (LLR) of the substitution, for variants held out from all score building:
  - gnomAD −5.77 (n = 32,517)
  - ClinVar benign −4.44 (n = 1,912)
  - ClinVar VUS −5.53 (n = 27,158)
  - ClinVar pathogenic −10.53 (n = 763)
- AUROC, gnomAD vs pathogenic: 0.80.
- ESM-1v cannot distinguish VUS from gnomAD variants (AUROC 0.48).

### M0: How the tolerance score is built (diagram)
**Schematic (simulated protein, real recipe).**
- ClinVar (above) and gnomAD (below) missense variants are weighted by substitution severity (lollipop height: BLOSUM62 in structured, Grantham in unstructured regions).
- They are averaged with a binomial rolling window: wide in structured regions (w = 31), narrow in unstructured regions (w = 5).
- The two densities are scaled 0–1 and subtracted: tolerance = gnomAD − ClinVar (middle lane; blue = tolerated, red = intolerant).

Wording note: a binomial window concentrates its weight near the centre. About 97% of a w = 31 window's weight lies within ±8 residues, and 98% of a w = 5 window's within ±3. If you write "±31 residues", that is the window's reach, not its effective width.

### M1: Building the tolerance score
**Tolerance = smoothed density of gnomAD variants − smoothed density of ClinVar variants.**
1. **Weighting:** each variant is weighted by the severity of the substitution (BLOSUM62 in ordered regions, Grantham distance in disordered regions).
2. **Smoothing:** the weight is spread along the protein with a binomial kernel (±31 residues in ordered regions, ±5 in disordered regions, using AIUPred disorder ≥ 0.5).
3. **Scaling and subtraction:** both densities are scaled 0–1 per gene and subtracted.

Above 0: more healthy-population variation than disease reports (tolerated). Below 0: intolerant.
gnomAD variants count only if seen at least twice (allele count ≥ 2; on chrX, only if carried by a male or a homozygous female).

### M2: Tolerance by variant group
**Tolerance, built without the scored variants, separates the variant groups.**
10% of variants (62,942) were held out, and tolerance was rebuilt from the remaining 90%. Mean tolerance at the held-out variants:
- gnomAD +0.20
- ClinVar benign +0.06
- VUS +0.04
- pathogenic −0.34

AUROC, gnomAD vs pathogenic: 0.90 (ESM-1v: 0.80).

**Wording caution:** VUS score lower than gnomAD variants (AUROC 0.67), but their mean is almost identical to ClinVar-benign and far from pathogenic.
The ClinVar side of the score counts every ClinVar submission, VUS included. So part of the VUS–gnomAD gap reflects VUS sitting near other ClinVar submissions.
Say "VUS fall between gnomAD and pathogenic variants, close to ClinVar-benign", not "VUS look pathogenic".

### M3a: GLRA1 structure, structured vs unstructured
**GLRA1 (glycine receptor α1) has structured and unstructured regions, and the unstructured ones still carry annotated sites.**
- **Left:** AlphaFold model coloured by AIUPred disorder (yellow, ≥ 0.5 = unstructured; signal peptide not shown).
- **Right:** domains, AIUPred disorder, AlphaFold confidence (pLDDT) and annotations along the sequence.
- **119 of 449 residues are unstructured:** the N-terminus, the intracellular M3–M4 loop and the end of M4. AlphaFold's confidence is low in the same regions (mean pLDDT 63 vs 93).
- **Annotation examples in the unstructured C-terminal half:**
  - a predicted binding stretch (AIUPred, MobiDB), E341–K353;
  - a curated binding stretch (MobiDB), the proline-rich T391–S398;
  - lipid contacts (PDBe-KB) along the end of M4, which UniProt annotates as a membrane helix.
- **Wording note:** the RFRRK core of E341–K353 matches a reported Gβγ-binding motif (Yevenes et al., *PNAS* 2008). That link is from the literature, not the model; check it before quoting.

### M3b: Core and surface
**Core residues hold the fold together; surface residues tolerate change.**
- **Panel:** the extracellular domain of one α1 subunit (PDB 8DN5, glycine bound), sliced at the level where glycine binds and coloured by relative solvent accessibility (RSA: core < 0.2, surface > 0.4).
- **What it shows:** core residues fill the space between the two β-sheets of the sandwich.
- **Held-out ClinVar variants that are pathogenic:** core 60% (n = 547), boundary 27% (n = 235), surface 15% (n = 344).

### M3c: Interfaces where glycine binds
**Glycine binds at interfaces between subunits, and interface residues are as constrained as the core.**
- **The structure:** the human α1β glycine receptor with glycine bound (PDB 8DN5).
- **The slice:** the same cross-section as M3b, widened to both neighbours.
  - Violet: exposed residues of subunit A that touch a neighbour (within 4.5 Å).
  - Pale violet: the neighbour's side.
  - Orange: glycine.
- **Each glycine pocket is formed by two subunits:**
  - F187, Y230 and F235 of one α1 subunit with F91, R93 and S157 of the next;
  - in mature numbering, F159, Y202, F207 / F63, R65, S129.
- **Held-out ClinVar variants that are pathogenic:** ligand-binding residues 86% (n = 99), interface 68% (n = 174), core 60%, other surface 15%.
- **Wording note:** the model's interface feature comes from PDBe-KB, which flags 101 of the 105 contact residues seen in 8DN5. Its ligand-site feature also counts sugar and lipid-like contacts.

### M3d: Where the features come from (table)
**Each database is queried per protein, mapped onto the MANE sequence and expanded to one value per residue.**
- AIUPred and MobiDB: structured or not.
- PDB (via PDBe) and AlphaFold DB: burial.
- PDBe-KB, MobiDB LIPs and ELM: interaction sites.
- InterPro and UniProt: domains and known sites.

UniProt natural variants and mutagenesis are excluded so clinical labels cannot leak into the features.

### M3: Structural features
**Per-residue structural context, mapped onto the MANE sequence:**
- disorder and binding (AIUPred);
- relative solvent accessibility (RSA: core < 0.2, surface > 0.4), from the best PDB chain, otherwise from AlphaFold (residues with pLDDT < 70 excluded);
- protein–protein interfaces (PDBe-KB);
- domains, transmembrane segments, ligand and PTM sites (InterPro, UniProt, MobiDB, ELM).

Dotted line: GLRA1 p.A300P (pathogenic).

### M4: Assembling the features (Figure 8)
**Every data source is aligned per residue, then looked up for each variant.**
- Each variant row combines its residue's features with features of the substitution itself (ESM-1v LLR, BLOSUM62 / Grantham, hydrophobicity change).
- It also includes human-variant features recomputed with the variant's own residue hidden, so no variant contributes to its own score.
- In total, 42 features per variant.

### R1: Models tested (table)
**Models were compared by 5-fold cross-validation grouped by gene family (15,339 ClinVar variants, 716 genes).**
- Logistic regression on 42 features matched more complex models (gradient-boosted trees, MLP, CNN and BiLSTM on raw variant tracks).
- It was pre-registered as the primary model for one test on untouched genes.
- AUROC: 0.5 = chance, 1.0 = perfect.

### R2 + R2b: The model, and its performance on untouched genes
**Logistic regression combines 42 standardised features into one probability that a variant is benign.**
On 1,357 ClinVar variants in 85 genes never used for training or model selection, it reached:

| | AUROC, all variants | AUROC where ESM-1v is unsure |
|---|---|---|
| Final model | 0.89 | 0.88 |
| Tolerance alone | 0.83 | 0.84 |
| ESM-1v alone | 0.62 | 0.54 |

"Where ESM-1v is unsure" = the 668 variants whose LLR is above the median of benign training variants.

### R3: Pre-registered test on untouched genes
**The model and comparisons were fixed before the test genes were scored.**
Differences in AUROC vs ESM-1v alone, with 95% CIs from 1,000 paired resamples of whole genes:

| Comparison | Δ AUROC | 95% CI | |
|---|---|---|---|
| All variants | +0.27 | +0.05 to +0.45 | significant |
| Where ESM-1v is unsure | +0.34 | +0.18 to +0.47 | significant |
| Within genes (only 8 genes testable) | +0.16 | −0.001 to +0.32 | not significant |
| ESM + tolerance vs ESM alone | +0.21 | +0.02 to +0.30 | significant |

### R4: What each data type adds (Figure 23)
**Human-variant data add information that ESM-1v lacks** (gradient-boosted trees, cross-validated). Added to ESM scores:
- human variants nearby: +0.029 AUROC (95% CI +0.018 to +0.040)
- structure: +0.036 (+0.014 to +0.066)
- disorder + conservation: +0.040 (+0.015 to +0.076)

The groups overlap: all together reach 0.958, vs 0.904 for ESM alone.

### R5: What the model relies on (Figure 26)
**Largest weights of the final model, per standard deviation of each feature.** Blue pushes towards benign, red towards pathogenic.
- The ESM-1v LLR dominates (+1.53).
- Next come local human variation:
  - gnomAD variants within ±15 residues: +1.07
  - ClinVar variants within ±15 residues: −0.99
- Then conservation (−0.90), interfaces (−0.49), solvent exposure (+0.44) and substitution severity (−0.40).

### Take-home message (suggested)
**Human population variation adds information that evolution alone (ESM-1v) misses.**
- A simple model combining gnomAD / ClinVar variant density with ESM-1v, structure and annotations classified ClinVar variants in untouched genes with AUROC 0.89, vs 0.62 for ESM-1v.
- The gain was largest where ESM-1v was uncertain: 0.88 vs 0.54.

### Small print worth including (honest limits)
- **The test set is dominated by three genes.** RYR1, ABL1 and VHL hold 52% of test variants. ESM-1v struggles in RYR1 and VHL, which inflates the all-variant gap. In the other 82 test genes: 0.96 vs 0.91.
- **The probabilities are not yet calibrated.** At P(benign) ≥ 0.5, 97% of benign but only 52% of pathogenic test variants were classified correctly.
- **Performance on VUS is not yet measured.** The test used confidently classified variants; VUS are harder.

---------------------------------------------------------------------------------------------------------

## Reading Figure 23 (R4) and Figure 26 (R5)

**Figure 23: feature groups.** It uses the gradient-boosted tree model, cross-validated on the 15,339 training/validation variants. Each feature group was refitted alone, added to the ESM scores, or removed from the full model.
- **Panel A, each group alone.** How much each data type knows by itself:
  - ESM 0.90; human variants nearby 0.89; disorder + conservation 0.89; structure 0.88; annotations 0.80; gene-wide rates 0.72; substitution chemistry 0.71.
  - Light bars are the within-gene AUROC. Gene-wide rates drop to 0.505 there, because they are the same for every variant in a gene.
- **Panel B, added to ESM.** The direct test of the hypothesis: does adding a group to ESM's scores improve them?
  - Human variants +0.029, structure +0.036, disorder + conservation +0.040. Their confidence intervals exclude 0.
  - Gains within genes are smaller (+0.009 to +0.017).
- **PDF-only panels.** The PDF version also shows the cost of removing each group from the full model (at most 0.005, because the groups overlap) and ROC curves.
- **Takeaway:** human-variant data carry real information beyond evolution. Structure and disorder carry a similar, overlapping amount, because all three mark constrained regions.

**Figure 26: the final model's weights.**
- **How to read a bar:** each bar is how much one feature moves the model's log-odds of "benign" when that feature rises by one standard deviation, with the others held fixed.
- **Intuitive directions:**
  - higher LLR (a less surprising swap) → benign
  - more gnomAD variants nearby → benign
  - more ClinVar variants nearby → pathogenic
  - higher conservation → pathogenic
  - interface residue → pathogenic
  - exposed residue → benign
  - drastic swap → pathogenic
  - disordered → benign
- **One counter-intuitive pair:**
  - ClinVar variants per residue across the whole gene get a *positive* weight (+0.85), next to the negative local count (−0.99).
  - Together they mean the model compares local ClinVar density with the gene's average: an excess *near* the variant is what signals pathogenicity.
  - The same holds for gnomAD (local +1.07, gene-wide −0.26).
- **The tolerance score itself gets a modest weight (+0.29).** Its information is shared with the raw counts, so the credit is split. Correlated features share weight, so don't read any single bar as an independent effect.
