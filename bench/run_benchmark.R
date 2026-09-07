#!/usr/bin/env Rscript
# ============================================================
# Benchmark driver.
#   Rscript bench/run_benchmark.R <experiment> [n_reps] [n_cores]
# Experiments: main, batch, exclusive, depth, null
# Writes results/<experiment>.csv
# ============================================================

BENCH_ROOT <- Sys.getenv("BENCH_ROOT",
                         "/Users/owaissiddiqi/Documents/methods_paper")
suppressPackageStartupMessages({library(edgeR); library(parallel)})
source(file.path(BENCH_ROOT, "sim/simulate_ctc.R"))
source(file.path(BENCH_ROOT, "bench/arms.R"))
source(file.path(BENCH_ROOT, "bench/external_arms.R"))
source(file.path(BENCH_ROOT, "bench/score.R"))

args    <- commandArgs(trailingOnly = TRUE)
exp_nm  <- if (length(args) >= 1) args[1] else "main"
n_reps  <- if (length(args) >= 2) as.integer(args[2]) else 50L
n_cores <- if (length(args) >= 3) as.integer(args[3]) else 8L

# ---- scenario grids ------------------------------------------------------
grid_main <- expand.grid(
  n_per_group = c(3L, 5L, 8L), stringsAsFactors = FALSE)
grid_main$arms <- I(rep(list(names(benchmark_arms)), nrow(grid_main)))

grid_batch <- expand.grid(
  batch_sd = c(0, 0.35, 0.70), batch_balanced = c(TRUE, FALSE),
  stringsAsFactors = FALSE)
grid_batch$n_per_group <- 5L
grid_batch$arms <- I(rep(list(c("ctc_pipeline", "batch_removed",
                                "batch_ignored")), nrow(grid_batch)))

grid_excl <- data.frame(exclusive_cpm = c(0.1, 0.5, 1, 2, 5, 10, 20, 50))
grid_excl$n_per_group <- 5L
grid_excl$arms <- I(rep(list(c("ctc_pipeline", "ablate_exclusive",
                               "edger_default")), nrow(grid_excl)))

# depth sweep: does the filter's effective CPM cutoff drift with library size?
grid_depth <- data.frame(lib_lo = c(1e5, 3e5, 1e6, 5e6, 2e7))
grid_depth$lib_hi        <- grid_depth$lib_lo * 3
grid_depth$n_per_group   <- 5L
grid_depth$arms <- I(rep(list(c("ctc_pipeline", "fix_filter",
                                "ablate_filter")), nrow(grid_depth)))

grid_null <- data.frame(n_per_group = c(3L, 5L, 8L))
grid_null$prop_de <- 0; grid_null$n_exclusive <- 0L
grid_null$arms <- I(rep(list(names(benchmark_arms)), nrow(grid_null)))

# purity: leukocyte contamination varying between conditions
grid_purity <- data.frame(purity_cond_diff = c(0, 0.10, 0.20, 0.30, 0.40))
grid_purity$purity_mean <- 0.60
grid_purity$purity_sd   <- 0.15
grid_purity$n_per_group <- 5L
grid_purity$arms <- I(rep(list(names(purity_arms)), nrow(grid_purity)))

# comparators: the CTC pipeline against the three most-used bulk DE tools,
# each at its own documented defaults
grid_comparators <- data.frame(n_per_group = c(3L, 5L, 8L))
grid_comparators$arms <- I(rep(list(c("ctc_pipeline", "edger_default",
                                      "limma_voom", "deseq2")),
                               nrow(grid_comparators)))

grids <- list(comparators = grid_comparators, main = grid_main, batch = grid_batch, exclusive = grid_excl,
              depth = grid_depth, null = grid_null, purity = grid_purity)
stopifnot(exp_nm %in% names(grids))
grid <- grids[[exp_nm]]

# ---- one replicate of one scenario ---------------------------------------
run_scenario_rep <- function(row, rep_i) {
  sim_args <- list(seed = 10000L * row + rep_i, n_per_group = grid$n_per_group[row])
  for (nm in c("batch_sd", "batch_balanced", "exclusive_cpm", "prop_de",
               "n_exclusive", "purity_mean", "purity_sd", "purity_cond_diff"))
    if (!is.null(grid[[nm]])) sim_args[[nm]] <- grid[[nm]][row]
  if (!is.null(grid$lib_lo))
    sim_args$lib_size_range <- c(grid$lib_lo[row], grid$lib_hi[row])

  sim <- do.call(simulate_ctc, sim_args)
  arm_set <- if (exp_nm == "purity") purity_arms else benchmark_arms
  out <- lapply(grid$arms[[row]], function(a) {
    # arms implemented by another package run through their own entry point;
    # every arm returns the same structure, so scoring is identical
    if (a %in% names(external_arms)) {
      r <- tryCatch(suppressMessages(
             external_arms[[a]](sim$counts, sim$samples,
                                list(fdr_cutoff = 0.05))),
           error = function(e) NULL)
      if (is.null(r)) return(NULL)
      s <- score_arm(r, sim$truth)
      s$arm <- a; s$rep <- rep_i; s$scenario_row <- row
      for (nm in setdiff(names(grid), "arms")) s[[nm]] <- grid[[nm]][row]
      s$median_lib <- median(colSums(sim$counts))
      return(s)
    }
    pars <- arm_set[[a]]
    # purity is estimated from the counts; the simulator supplies only the
    # compartment marker sets, never the true purity values
    if (exp_nm == "purity") pars <- modifyList(pars,
                              list(purity_signatures = sim$signatures))
    r <- tryCatch(suppressMessages(
           ctc_dge_core(sim$counts, sim$samples, pars)),
         error = function(e) NULL)
    if (is.null(r)) return(NULL)
    s <- score_arm(r, sim$truth)
    s$arm <- a; s$rep <- rep_i; s$scenario_row <- row
    for (nm in setdiff(names(grid), "arms")) s[[nm]] <- grid[[nm]][row]
    s$median_lib <- median(colSums(sim$counts))
    s
  })
  do.call(rbind, out)
}

cat(sprintf("[%s] experiment '%s': %d scenarios x %d reps on %d cores\n",
            format(Sys.time(), "%H:%M:%S"), exp_nm, nrow(grid), n_reps, n_cores))

jobs <- expand.grid(row = seq_len(nrow(grid)), rep_i = seq_len(n_reps))
res  <- mclapply(seq_len(nrow(jobs)), function(k)
          run_scenario_rep(jobs$row[k], jobs$rep_i[k]),
          mc.cores = n_cores, mc.preschedule = TRUE)

failed <- sum(vapply(res, is.null, logical(1)))
df <- do.call(rbind, Filter(Negate(is.null), res))
dir.create(file.path(BENCH_ROOT, "results"), showWarnings = FALSE)
out_f <- file.path(BENCH_ROOT, "results", paste0(exp_nm, ".csv"))
write.csv(df, out_f, row.names = FALSE)
cat(sprintf("[%s] wrote %s : %d rows (%d failed jobs)\n",
            format(Sys.time(), "%H:%M:%S"), out_f, nrow(df), failed))
