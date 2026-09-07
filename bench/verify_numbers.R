#!/usr/bin/env Rscript
# Recompute every figure quoted in the manuscript from the result tables and
# compare against the value as written. Run before submission.
ROOT <- Sys.getenv("BENCH_ROOT", "/Users/owaissiddiqi/Documents/methods_paper")
setwd(ROOT)
fails <- 0L; checks <- 0L
chk <- function(label, got, want, tol = 5e-4) {
  checks <<- checks + 1L
  ok <- is.finite(got) && abs(got - want) <= tol
  if (!ok) fails <<- fails + 1L
  cat(sprintf("%-58s manuscript %-10s computed %-10s %s\n", label,
              formatC(want, format="g", digits=6),
              formatC(got,  format="g", digits=6),
              if (ok) "ok" else "**MISMATCH**"))
}
m  <- function(f) read.csv(file.path("results", f))
ag <- function(d, v, by) aggregate(d[[v]], d[by], function(x) mean(x, na.rm=TRUE))
pick <- function(a, ...) { k <- list(...); r <- a
  for (n in names(k)) r <- r[r[[n]] == k[[n]], , drop=FALSE]; r$x }

cat("=== Table 1: pipeline vs edgeR defaults ===\n")
d <- m("main.csv")
a <- ag(d, "avg_precision", c("arm","n_per_group"))
f <- ag(d, "empirical_fdr", c("arm","n_per_group"))
t_ <- ag(d, "n_tested",      c("arm","n_per_group"))
for (n in c(3,5,8)) {
  chk(sprintf("  AP  pipeline n=%d", n), pick(a, arm="ctc_pipeline",  n_per_group=n),
      c(`3`=0.2015,`5`=0.2536,`8`=0.3032)[[as.character(n)]])
  chk(sprintf("  AP  defaults n=%d", n), pick(a, arm="edger_default", n_per_group=n),
      c(`3`=0.1849,`5`=0.2153,`8`=0.2289)[[as.character(n)]])
  chk(sprintf("  FDR pipeline n=%d", n), pick(f, arm="ctc_pipeline",  n_per_group=n),
      c(`3`=0.026,`5`=0.041,`8`=0.068)[[as.character(n)]], tol=6e-4)
  chk(sprintf("  tested pipeline n=%d", n), pick(t_, arm="ctc_pipeline", n_per_group=n),
      c(`3`=11962,`5`=12540,`8`=11269)[[as.character(n)]], tol=1)
  chk(sprintf("  tested defaults n=%d", n), pick(t_, arm="edger_default", n_per_group=n),
      c(`3`=4449,`5`=4515,`8`=4100)[[as.character(n)]], tol=1)
}
rel <- function(n) 100*(pick(a,arm="ctc_pipeline",n_per_group=n)/pick(a,arm="edger_default",n_per_group=n)-1)
for (n in c(3,5,8)) chk(sprintf("  relative AP gain n=%d (%%)", n), rel(n),
                        c(`3`=9.0,`5`=17.8,`8`=32.5)[[as.character(n)]], tol=0.06)

cat("\n=== Table 2: ablation deltas at n=5 ===\n")
b <- a[a$n_per_group==5,]; base <- pick(a, arm="ctc_pipeline", n_per_group=5)
want <- c(edger_default=-0.0384, batch_removed=-0.0239, ablate_filter=-0.0141,
          ablate_exclusive=-0.0136, batch_ignored=-0.0101, ablate_tmmwsp=-0.0004)
for (k in names(want)) chk(paste0("  delta AP ", k), pick(a, arm=k, n_per_group=5)-base,
                           want[[k]], tol=6e-5)

cat("\n=== Batch experiment ===\n")
d <- m("batch.csv"); a2 <- ag(d, "empirical_fdr", c("arm","batch_sd","batch_balanced"))
chk("  pipeline FDR range min (all cells)",
    min(a2$x[a2$arm=="ctc_pipeline"]), 0.037, tol=6e-4)
chk("  pipeline FDR range max (all cells)",
    max(a2$x[a2$arm=="ctc_pipeline"]), 0.044, tol=6e-4)
chk("  removed, balanced, sd=0",
    pick(a2, arm="batch_removed", batch_sd=0, batch_balanced=TRUE), 0.192, tol=6e-4)
chk("  removed, imbalanced, sd=0.35",
    pick(a2, arm="batch_removed", batch_sd=0.35, batch_balanced=FALSE), 0.722, tol=6e-4)
chk("  ignored, imbalanced, sd=0.70",
    pick(a2, arm="batch_ignored", batch_sd=0.70, batch_balanced=FALSE), 0.516, tol=6e-4)

cat("\n=== Null calibration ===\n")
d <- m("null.csv"); a3 <- ag(d, "n_significant", c("arm","n_per_group"))
for (n in c(3,5,8)) chk(sprintf("  batch_removed FP n=%d", n),
    pick(a3, arm="batch_removed", n_per_group=n),
    c(`3`=17.1,`5`=2.7,`8`=1.1)[[as.character(n)]], tol=0.06)

cat("\n=== Exclusive retention ===\n")
d <- m("exclusive.csv"); e <- d[d$arm=="ctc_pipeline",]
chk("  max true exclusives rescued (of 120)",
    max(ag(e,"sig_excl_rescued",c("exclusive_cpm"))$x), 1.20, tol=0.01)
chk("  genes admitted, mean across sweep",
    mean(e$n_exclusive_added), 3650, tol=90)

cat("\n=== Depth sweep ===\n")
d <- m("depth.csv"); a4 <- ag(d, "median_lib", c("arm","lib_lo"))
lo <- min(d$lib_lo); hi <- max(d$lib_lo)
ml_lo <- pick(a4, arm="ctc_pipeline", lib_lo=lo); ml_hi <- pick(a4, arm="ctc_pipeline", lib_lo=hi)
chk("  edgeR default cutoff at lowest depth (CPM)", 10/ml_lo*1e6, 54, tol=1.5)
chk("  legacy cutoff at lowest depth (CPM)",  0.1/ml_lo*1e6, 0.54, tol=0.02)
chk("  legacy cutoff at highest depth (CPM)", 0.1/ml_hi*1e6, 0.0026, tol=0.0002)
chk("  drift factor across depth range", (0.1/ml_lo*1e6)/(0.1/ml_hi*1e6), 208, tol=6)
t5 <- ag(d, "n_tested", c("arm","lib_lo"))
chk("  genes tested, defaults, lowest depth", pick(t5, arm="ablate_filter", lib_lo=lo), 5857, tol=1)
for (i in seq_along(sort(unique(d$lib_lo))[1:3])) {
  L <- sort(unique(d$lib_lo))[i]
  chk(sprintf("  legacy == corrected tested set, depth %d", i),
      pick(t5, arm="ctc_pipeline", lib_lo=L) - pick(t5, arm="fix_filter", lib_lo=L), 0, tol=0.5)
}

cat("\n=== Table 3: purity ===\n")
d <- m("purity.csv"); a6 <- ag(d, "empirical_fdr", c("arm","purity_cond_diff"))
w <- list("0"=c(0.037,0.041,0.030), "0.1"=c(0.104,0.039,0.062),
          "0.2"=c(0.232,0.043,0.167), "0.3"=c(0.564,0.064,0.449),
          "0.4"=c(0.746,0.099,0.680))
arms <- c("purity_ignored","purity_covariate","purity_filtered")
for (g in names(w)) for (j in 1:3)
  chk(sprintf("  gap %-3s %s", g, arms[j]),
      pick(a6, arm=arms[j], purity_cond_diff=as.numeric(g)), w[[g]][j], tol=6e-4)

cat("\n=== Real data: purity score ===\n")
d <- m("real_purity.csv"); a7 <- aggregate(purity_score ~ type, d, mean)
mp <- function(t) a7$purity_score[a7$type==t]
chk("  cell line",        mp("single cell from PCa cell line"), 0.845, tol=6e-4)
chk("  confirmed CTC",    mp("lineage-confirmed PCa CTC"),      0.706, tol=6e-4)
chk("  primary tumour",   mp("primary PCa tumor"),              0.550, tol=6e-4)
chk("  candidate CTC",    mp("candidate PCa CTC"),              0.427, tol=6e-4)
chk("  leukocyte",        mp("white blood cell"),               0.288, tol=6e-4)

cat("\n=== Real data: permutation null ===\n")
d <- m("real_rd1.csv"); d <- d[!d$failed,]
a8 <- aggregate(n_significant ~ arm, d, mean)
fp <- function(x) a8$n_significant[a8$arm==x]
chk("  ctc_pipeline FP/experiment",  fp("ctc_pipeline"),  1135, tol=1)
chk("  edger_default FP/experiment", fp("edger_default"), 243,  tol=1)
chk("  batch_removed FP/experiment", fp("batch_removed"), 2546, tol=1)
chk("  batch_ignored FP/experiment", fp("batch_ignored"), 2882, tol=1)

cat("\n=== Table 1: comparators ===\n")
d <- m("comparators.csv")
a <- ag(d, "avg_precision", c("arm","n_per_group"))
f <- ag(d, "empirical_fdr", c("arm","n_per_group"))
want_ap <- list("ctc_pipeline"=c(0.2047,0.2514,0.3025),
                "deseq2"=c(0.1845,0.2260,0.2760),
                "edger_default"=c(0.1880,0.2150,0.2286),
                "limma_voom"=c(0.1783,0.2006,0.2147))
for (arm in names(want_ap)) for (i in seq_along(c(3,5,8)))
  chk(sprintf("  AP %-14s n=%d", arm, c(3,5,8)[i]),
      pick(a, arm=arm, n_per_group=c(3,5,8)[i]), want_ap[[arm]][i], tol=6e-4)
for (i in seq_along(c(3,5,8)))
  chk(sprintf("  DESeq2 emp FDR n=%d", c(3,5,8)[i]),
      pick(f, arm="deseq2", n_per_group=c(3,5,8)[i]),
      c(0.225,0.113,0.057)[i], tol=6e-4)

cat("\n=== paired advantage over each comparator ===\n")
w <- reshape(d[,c("rep","n_per_group","arm","avg_precision")],
             idvar=c("rep","n_per_group"), timevar="arm", direction="wide")
for (cmp in c("edger_default","limma_voom","deseq2")) {
  diff <- w$avg_precision.ctc_pipeline - w[[paste0("avg_precision.", cmp)]]
  chk(sprintf("  mean AP advantage vs %s", cmp), mean(diff, na.rm=TRUE),
      c(edger_default=0.042, limma_voom=0.055, deseq2=0.024)[[cmp]], tol=6e-4)
  chk(sprintf("  replicates won vs %s", cmp), sum(diff > 0, na.rm=TRUE), 150, tol=0)
}

cat(sprintf("\n%d checks, %d mismatches\n", checks, fails))
quit(status = if (fails > 0) 1 else 0)
