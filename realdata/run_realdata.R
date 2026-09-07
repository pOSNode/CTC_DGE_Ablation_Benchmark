#!/usr/bin/env Rscript
# ============================================================
# Real-data validation on GSE67980 (Miyamoto et al. 2015),
# prostate CTC RNA-seq. No simulator assumptions are used here.
#
#   RD1  permutation null, donor-balanced   -> false positives on real noise
#   RD2  permutation null, donor-imbalanced -> batch handling on real structure
#   RD3  CTC vs leukocyte                   -> known-marker recovery
#
# Under a permutation null there is by construction no true differential
# expression, so every gene called significant is a false positive.
# ============================================================

ROOT <- Sys.getenv("BENCH_ROOT", "/Users/owaissiddiqi/Documents/methods_paper")
suppressPackageStartupMessages({library(edgeR); library(parallel)})
source(file.path(ROOT, "realdata/load_gse67980.R"))
BENCH_ROOT <- ROOT
source(file.path(ROOT, "bench/arms.R"))

args    <- commandArgs(trailingOnly = TRUE)
which_e <- if (length(args) >= 1) args[1] else "rd1"
n_perm  <- if (length(args) >= 2) as.integer(args[2]) else 50L
n_cores <- if (length(args) >= 3) as.integer(args[3]) else 8L

d <- load_gse67980()
lc <- d$samples$source == "lineage-confirmed PCa CTC"
donor_n <- table(d$samples$donor[lc])
keep_donors <- names(donor_n)[donor_n >= 6]
idx <- which(lc & d$samples$donor %in% keep_donors)

counts_ctc <- d$counts[, idx, drop = FALSE]
donors     <- d$samples$donor[idx]
# drop genes with no reads anywhere - they carry no information in any arm
counts_ctc <- counts_ctc[rowSums(counts_ctc) > 0, , drop = FALSE]

cat(sprintf("GSE67980: %d lineage-confirmed CTCs from %d donors, %d genes\n",
            ncol(counts_ctc), length(keep_donors), nrow(counts_ctc)))
cat(sprintf("  donors: %s\n", paste(sprintf("%s(%d)", names(table(donors)),
                                            table(donors)), collapse = ", ")))
cat(sprintf("  median library %s, zero fraction %.0f%%\n\n",
            format(round(median(colSums(counts_ctc))), big.mark = ","),
            100 * mean(counts_ctc == 0)))

# ---- assign a fake condition ------------------------------------------
# balanced : within each donor, split cells evenly between the two labels
# imbalanced: each donor contributes ~75% of its cells to one label,
#             reproducing partial donor/condition confounding
assign_labels <- function(donors, mode, seed) {
  set.seed(seed)
  lab <- character(length(donors))
  ud  <- unique(donors)
  for (k in seq_along(ud)) {
    i <- which(donors == ud[k])
    i <- sample(i)
    if (mode == "balanced") {
      lab[i] <- rep(c("A", "B"), length.out = length(i))
    } else {
      major <- if (k <= ceiling(length(ud) / 2)) "A" else "B"
      minor <- setdiff(c("A", "B"), major)
      nmaj  <- max(1, round(0.75 * length(i)))
      lab[i] <- c(rep(major, nmaj), rep(minor, length(i) - nmaj))
    }
  }
  lab
}

run_perm <- function(p, mode, arms) {
  lab <- assign_labels(donors, mode, seed = 20000L + p)
  if (min(table(lab)) < 3) return(NULL)
  samples <- data.frame(
    sample_id    = colnames(counts_ctc),
    condition    = factor(lab, levels = c("A", "B")),
    batch        = factor(donors),
    cell_count   = 1L, control_type = "", stringsAsFactors = FALSE)
  out <- lapply(arms, function(a) {
    r <- tryCatch(suppressMessages(
           ctc_dge_core(counts_ctc, samples, benchmark_arms[[a]])),
         error = function(e) NULL)
    if (is.null(r)) return(data.frame(arm = a, perm = p, mode = mode,
                                      n_tested = NA, n_significant = NA,
                                      batch_in_model = NA, failed = TRUE))
    data.frame(arm = a, perm = p, mode = mode,
               n_tested = r$n_tested, n_significant = r$n_significant,
               batch_in_model = r$batch_in_model, failed = FALSE)
  })
  do.call(rbind, out)
}

dir.create(file.path(ROOT, "results"), showWarnings = FALSE)

if (which_e %in% c("rd1", "rd2")) {
  mode <- if (which_e == "rd1") "balanced" else "imbalanced"
  arms <- names(benchmark_arms)
  cat(sprintf("[%s] %s permutation null: %d permutations x %d arms\n",
              format(Sys.time(), "%H:%M:%S"), mode, n_perm, length(arms)))
  res <- mclapply(seq_len(n_perm), function(p) run_perm(p, mode, arms),
                  mc.cores = n_cores)
  df <- do.call(rbind, Filter(Negate(is.null), res))
  f <- file.path(ROOT, "results", paste0("real_", which_e, ".csv"))
  write.csv(df, f, row.names = FALSE)
  cat(sprintf("[%s] wrote %s (%d rows)\n", format(Sys.time(), "%H:%M:%S"),
              f, nrow(df)))
  a <- aggregate(n_significant ~ arm, df, function(x) mean(x, na.rm = TRUE))
  a <- a[order(a$n_significant), ]
  cat("\nFalse positives per experiment (no true DE by construction):\n")
  print(a, row.names = FALSE)
}

if (which_e == "rd3") {
  # real biological contrast: CTC vs leukocyte
  sel <- d$samples$source %in% c("lineage-confirmed PCa CTC", "white blood cell")
  cm  <- d$counts[, sel, drop = FALSE]
  cm  <- cm[rowSums(cm) > 0, , drop = FALSE]
  samples <- data.frame(
    sample_id  = colnames(cm),
    condition  = factor(ifelse(d$samples$source[sel] == "white blood cell",
                               "WBC", "CTC"), levels = c("CTC", "WBC")),
    batch      = factor("b1"), cell_count = 1L, control_type = "",
    stringsAsFactors = FALSE)
  cat(sprintf("CTC vs WBC: %d vs %d samples, %d genes\n",
              sum(samples$condition == "CTC"), sum(samples$condition == "WBC"),
              nrow(cm)))
  epi <- c("KLK3","KLK2","AR","EPCAM","KRT8","KRT18","KRT19","TMPRSS2","FOLH1","STEAP1")
  wbc <- c("PTPRC","CD52","LYZ","CSF1R","ITGAM","CXCR4","HLA-DRA","LCP1")
  out <- lapply(names(benchmark_arms), function(a) {
    r <- tryCatch(suppressMessages(ctc_dge_core(cm, samples, benchmark_arms[[a]])),
                  error = function(e) NULL)
    if (is.null(r)) return(NULL)
    t <- r$table
    up_ctc <- t$gene_id[t$significant & t$logFC > 0]
    up_wbc <- t$gene_id[t$significant & t$logFC < 0]
    data.frame(arm = a, n_tested = r$n_tested, n_significant = r$n_significant,
               epi_tested = sum(epi %in% t$gene_id),
               epi_up_in_CTC = sum(epi %in% up_ctc),
               wbc_tested = sum(wbc %in% t$gene_id),
               wbc_up_in_WBC = sum(wbc %in% up_wbc))
  })
  df <- do.call(rbind, out)
  write.csv(df, file.path(ROOT, "results", "real_rd3.csv"), row.names = FALSE)
  print(df, row.names = FALSE)
}
