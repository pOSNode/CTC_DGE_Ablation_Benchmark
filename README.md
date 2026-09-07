# Ablation benchmark for low-input CTC RNA-seq differential expression

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22646903.svg)](https://doi.org/10.5281/zenodo.22646903)

Ground-truth simulation and one-at-a-time ablation of every documented design
decision in [`Rare_CTC_RNA_DGE_Pipeline`](https://github.com/pOSNode/Rare_CTC_RNA_DGE_Pipeline).

The question: a low-input pipeline departs from standard parameters in nine
places at once, each justified by argument. Which of those departures actually
carry the benefit?

## Headline results

| Finding | Evidence |
|---|---|
| The pipeline beats edgeR defaults at every sample size | AP +9.0% / +17.8% / +32.5% at n = 3 / 5 / 8 |
| It controls FDR despite four permissive choices | Empirical FDR 0.026-0.068 vs nominal 0.10 |
| Additive batch modelling is the load-bearing decision | Pre-subtraction: FDR 0.21 balanced, 0.73 imbalanced |
| Pre-subtracting batch invents signal from noise | 17.1 false positives/experiment under a complete null at n = 3 |
| Exclusive retention works by an unintended mechanism | Rescues <=1.2 genes/experiment while admitting ~3,650 |
| Four decisions are inert | AP change < 2x10^-5 |
| Ignoring tumour purity breaks FDR control | 0.75 empirical FDR at a 0.40 purity gap; 0.10 covariate |
| Purity as a covariate fixes it, filtering does not | Covariate 0.04-0.10 throughout; filtering 0.68 |

## Real-data validation

`GSE67980` (Miyamoto et al. 2015), 169 prostate samples: 77 lineage-confirmed
CTCs, 45 candidates, 30 cell-line single cells, 12 primary tumours, 5
leukocytes. Real CTCs are 79% zeros at a 4.4M median library.

- **Purity score validation.** The signature score orders sample types exactly
  as biology predicts (cell line 0.845 > confirmed CTC 0.706 > primary tumour
  0.550 > candidate CTC 0.427 > leukocyte 0.288), with AUC 0.990 separating
  confirmed CTCs from leukocytes. It reproduces the original authors' lineage
  confirmation without being shown the labels.
- **Permutation null.** Lineage-confirmed CTCs from the six donors with >= 6
  cells are randomly relabelled into two fake conditions, balanced within
  donor. There is no true differential expression by construction, so every
  gene called is a false positive - an empirical FDR measured on real noise,
  with no simulator assumptions. An imbalanced variant reproduces partial
  donor/condition confounding.

Running the pipeline against this data is also what exposed the
upper-quartile normalisation crash: UQ is undefined when a sample's upper
quartile is zero, which never occurred in simulation but happens routinely at
real CTC sparsity.

## Layout

```
sim/simulate_ctc.R        NB + dropout simulator with retained ground truth,
                          leukocyte background and per-sample purity
realdata/                 GSE67980 loader and real-data experiments
bench/arms.R              the 12 arms (pipeline, defaults, 7 ablations, 2 batch variants)
bench/score.R             sensitivity / empirical FDR / average precision
bench/run_benchmark.R     driver:  Rscript bench/run_benchmark.R <experiment> [reps] [cores]
figures/make_figures.R    all seven publication figures
results/*.csv             per-replicate results
bench/verify_numbers.R    recomputes every value quoted in the paper from the
                          result tables and fails on any mismatch
pipeline/                 the pipeline under test (separate repository)
```

## Reproducing

```bash
Rscript -e 'BiocManager::install(c("edgeR","limma","statmod"))'
for e in main batch exclusive depth null purity; do
  Rscript bench/run_benchmark.R $e 50 8
done
for e in rd1 rd2 rd3; do
  Rscript realdata/run_realdata.R $e 50 8   # downloads GSE67980 on first use
done
Rscript figures/make_figures.R
Rscript bench/verify_numbers.R
```

Every replicate records its seed, so any individual simulation can be regenerated
in isolation.

## Method notes

- **Sensitivity is scored over all true positives**, including genes a filter
  discarded before testing. Scoring only tested genes would make aggressive
  filtering appear free.
- **Average precision** is computed over the full truth set with untested genes
  ranked last, giving a threshold-free comparison between arms operating at
  different nominal FDR cutoffs.
- **The benchmark and the pipeline share one implementation**
  (`pipeline/scripts/R/ctc_dge_core.R`), so they cannot drift apart. An ablation
  is a one-field change to a parameter list, not a reimplementation.
- **The simulator was validated under a complete null** before use: p-values
  uniform (KS 0.022), type-I error 0.045 at nominal 0.05, 0.1 false positives per
  null experiment at FDR < 0.05. It does not manufacture the effects reported.

## Scope

Simulation is count-level only. Eight of the nine design decisions act on the count matrix;
the upstream ones (trimming, bootstrap count, quasi-mapping) would need
read-level simulation and are out of scope. No public CTC RNA-seq dataset
carries known differential expression, so real-data validation uses a
permutation null and a signature score with known sample types rather than a
labelled DE truth set.
