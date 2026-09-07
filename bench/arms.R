# Benchmark arms: the shipped pipeline, a standard-practice comparator,
# and one-at-a-time ablations of each documented design decision.

source(file.path(BENCH_ROOT, "pipeline/scripts/R/ctc_purity.R"))
source(file.path(BENCH_ROOT, "pipeline/scripts/R/ctc_dge_core.R"))

ctc_v1 <- ctc_default_params()          # pipeline as shipped at v1.0.0

edger_default <- modifyList(ctc_v1, list(
  filter_mode = "default", keep_exclusive = FALSE, norm_method = "TMM",
  low_input_mode = FALSE, ql_shrink_floor = FALSE, fdr_cutoff = 0.05,
  prior_count = 2, priority_override = FALSE))

benchmark_arms <- list(
  # --- the two headline arms ---------------------------------------------
  ctc_pipeline       = ctc_v1,
  edger_default      = edger_default,

  # --- one-at-a-time ablations from ctc_pipeline --------------------------
  ablate_filter      = modifyList(ctc_v1, list(filter_mode = "default")),
  ablate_exclusive   = modifyList(ctc_v1, list(keep_exclusive = FALSE)),
  ablate_tmmwsp      = modifyList(ctc_v1, list(norm_method = "TMM")),
  ablate_priorcount  = modifyList(ctc_v1, list(prior_count = 2)),
  ablate_fdr         = modifyList(ctc_v1, list(fdr_cutoff = 0.05)),
  ablate_shrinkfloor = modifyList(ctc_v1, list(ql_shrink_floor = FALSE)),
  ablate_lowinput    = modifyList(ctc_v1, list(low_input_mode = FALSE)),

  # --- batch handling: the pipeline's central methodological claim --------
  batch_removed      = modifyList(ctc_v1, list(batch_mode = "removed")),
  batch_ignored      = modifyList(ctc_v1, list(batch_mode = "ignored")),

  # --- the corrected filter call ------------------------------------------
  fix_filter         = modifyList(ctc_v1, list(filter_mode = "fixed"))
)

# Purity arms are kept separate: they need compartment marker sets, which the
# simulator supplies per replicate, so they are attached at run time.
purity_arms <- list(
  purity_ignored   = ctc_v1,
  purity_covariate = modifyList(ctc_v1, list(purity_mode = "covariate")),
  purity_filtered  = modifyList(ctc_v1, list(purity_mode = "filter",
                                             purity_min = 0.30))
)
