#!/usr/bin/env Rscript
# Publication figures for the CTC DGE benchmark.
# Categorical slots 1-3 of the validated reference palette; never more than
# three series in a panel; one axis per panel; recessive grid.

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(scales)
})

ROOT <- Sys.getenv("BENCH_ROOT", "/Users/owaissiddiqi/Documents/methods_paper")
FIG  <- file.path(ROOT, "figures")

S1 <- "#2a78d6"; S2 <- "#eb6834"; S3 <- "#1baf7a"
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#9a9992"; GRID <- "#e6e5e1"

theme_ctc <- function(base = 9) {
  theme_minimal(base_size = base) +
    theme(
      text             = element_text(colour = INK),
      axis.text        = element_text(colour = INK2, size = base - 1),
      axis.title       = element_text(colour = INK2, size = base),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
      panel.background = element_rect(fill = NA, colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      strip.text       = element_text(colour = INK, face = "bold", size = base),
      plot.title       = element_text(face = "bold", size = base + 2),
      plot.subtitle    = element_text(colour = INK2, size = base),
      plot.caption     = element_text(colour = MUTED, size = base - 1, hjust = 0),
      legend.position  = "top",
      legend.title     = element_blank(),
      legend.key.size  = unit(9, "pt")
    )
}

arm_label <- c(
  ctc_pipeline = "CTC pipeline (all decisions)",
  edger_default = "edgeR defaults",
  ablate_filter = "- relaxed CPM filter",
  ablate_exclusive = "- exclusive retention",
  ablate_tmmwsp = "- TMMwsp (use TMM)",
  ablate_priorcount = "- prior count 0.5",
  ablate_fdr = "- FDR 0.1 (use 0.05)",
  ablate_shrinkfloor = "- QL shrinkage floor",
  ablate_lowinput = "- low-input mode",
  batch_removed = "batch pre-subtracted",
  batch_ignored = "batch not modelled",
  fix_filter = "corrected filter call")

arm_class <- function(a)
  ifelse(a == "ctc_pipeline", "CTC pipeline",
  ifelse(a == "edger_default", "edgeR defaults",
  ifelse(a %in% c("batch_removed", "batch_ignored"), "batch variant",
         "ablation")))

se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# =====================================================================
# Figure 1 - average precision and empirical FDR across all arms
# =====================================================================
fig1 <- function() {
  d <- read.csv(file.path(ROOT, "results/main.csv"))
  a <- aggregate(cbind(avg_precision, empirical_fdr, sensitivity) ~ arm + n_per_group,
                 d, mean, na.rm = TRUE)
  s <- aggregate(cbind(avg_precision, empirical_fdr) ~ arm + n_per_group, d, se)
  a$ap_se <- s$avg_precision; a$fdr_se <- s$empirical_fdr
  a$class <- arm_class(a$arm)
  a$label <- arm_label[a$arm]
  ord <- a$label[a$n_per_group == 5][order(a$avg_precision[a$n_per_group == 5])]
  a$label <- factor(a$label, levels = ord)
  a$facet <- factor(sprintf("n = %d per group", a$n_per_group))

  pal <- c("CTC pipeline" = S1, "edgeR defaults" = S2,
           "batch variant" = S3, "ablation" = MUTED)

  p1 <- ggplot(a, aes(avg_precision, label, colour = class)) +
    geom_errorbar(aes(xmin = avg_precision - ap_se, xmax = avg_precision + ap_se),
                  orientation = "y", width = 0, linewidth = 0.5, alpha = 0.7) +
    geom_point(size = 2.2) +
    facet_wrap(~facet, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = pal) +
    scale_x_continuous(expand = expansion(mult = 0.12),
                       labels = label_number(accuracy = 0.01)) +
    labs(x = "Average precision (higher is better)", y = NULL,
         title = "A  Ranking quality of every design decision",
         subtitle = "Mean over 50 simulated CTC experiments; bars are standard errors") +
    theme_ctc()

  p2 <- ggplot(a, aes(empirical_fdr, label, colour = class)) +
    geom_vline(xintercept = 0.10, linetype = "22", colour = INK2, linewidth = 0.4) +
    geom_errorbar(aes(xmin = empirical_fdr - fdr_se, xmax = empirical_fdr + fdr_se),
                  orientation = "y", width = 0, linewidth = 0.5, alpha = 0.7) +
    geom_point(size = 2.2) +
    facet_wrap(~facet, nrow = 1) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(labels = label_number(accuracy = 0.01),
                       expand = expansion(mult = 0.08)) +
    labs(x = "Empirical false discovery rate", y = NULL,
         title = "B  Realised FDR against the nominal 0.10 threshold (dashed)",
         caption = paste("Arms testing at nominal FDR 0.05 (edgeR defaults,",
                         "- FDR 0.1) are expected to sit lower.")) +
    theme_ctc()

  combined <- (p1 / p2) + plot_layout(guides = "collect") &
    theme(legend.position = "top")
  ggsave(file.path(FIG, "fig1_main_arms.png"), combined, width = 9.5,
         height = 8.5, dpi = 300, bg = "white")
  cat("fig1 written\n")
}

# =====================================================================
# Figure 2 - batch handling
# =====================================================================
fig2 <- function() {
  f <- file.path(ROOT, "results/batch.csv"); if (!file.exists(f)) return(invisible())
  d <- read.csv(f)
  a <- aggregate(cbind(empirical_fdr, sensitivity, avg_precision) ~
                   arm + batch_sd + batch_balanced, d, mean, na.rm = TRUE)
  s <- aggregate(empirical_fdr ~ arm + batch_sd + batch_balanced, d, se)
  a$fdr_se <- s$empirical_fdr
  a$arm <- factor(arm_label[a$arm],
                  levels = c("CTC pipeline (all decisions)",
                             "batch pre-subtracted", "batch not modelled"))
  a$facet <- factor(ifelse(a$batch_balanced, "Batch balanced across conditions",
                           "Batch imbalanced (75/25)"))

  p <- ggplot(a, aes(batch_sd, empirical_fdr, colour = arm, group = arm)) +
    geom_hline(yintercept = 0.10, linetype = "22", colour = INK2, linewidth = 0.4) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.4) +
    geom_errorbar(aes(ymin = empirical_fdr - fdr_se, ymax = empirical_fdr + fdr_se),
                  width = 0, linewidth = 0.5, alpha = 0.7) +
    facet_wrap(~facet) +
    scale_colour_manual(values = c(S1, S2, S3)) +
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    labs(x = "Batch effect size (SD of per-gene log batch factor)",
         y = "Empirical false discovery rate",
         title = "Pre-subtracting batch inflates the false discovery rate",
         subtitle = paste("Nominal FDR 0.10 (dashed). n = 5 per group,",
                          "50 simulations per point."),
         caption = paste("Additive batch modelling holds at or below nominal;",
                         "removeBatchEffect does not.")) +
    theme_ctc()
  ggsave(file.path(FIG, "fig2_batch.png"), p, width = 8, height = 4.2,
         dpi = 300, bg = "white")
  cat("fig2 written\n")
}

# =====================================================================
# Figure 3 - condition-exclusive retention
# =====================================================================
fig3 <- function() {
  f <- file.path(ROOT, "results/exclusive.csv"); if (!file.exists(f)) return(invisible())
  d <- read.csv(f)
  a <- aggregate(cbind(sens_exclusive, n_excl_kept, n_exclusive_added,
                       sig_excl_rescued, avg_precision) ~ arm + exclusive_cpm,
                 d, mean, na.rm = TRUE)
  a$armf <- factor(arm_label[a$arm],
                   levels = c("CTC pipeline (all decisions)",
                              "- exclusive retention", "edgeR defaults"))

  pA <- ggplot(a, aes(exclusive_cpm, sens_exclusive, colour = armf, group = armf)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
    scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 5, 10, 20, 50)) +
    scale_colour_manual(values = c(S1, S2, S3)) +
    labs(x = "Expression of exclusive genes (CPM, present condition)",
         y = "Fraction detected at FDR < 0.10",
         title = "A  Detection of condition-exclusive genes") +
    theme_ctc()

  b <- a[a$arm == "ctc_pipeline", ]
  bb <- rbind(
    data.frame(exclusive_cpm = b$exclusive_cpm, n = b$n_excl_kept,
               what = "True exclusive genes retained (max 120)"),
    data.frame(exclusive_cpm = b$exclusive_cpm, n = b$n_exclusive_added,
               what = "Genes admitted by the rule (total)"))
  pB <- ggplot(bb, aes(exclusive_cpm, n, colour = what, group = what)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
    scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 5, 10, 20, 50)) +
    scale_y_log10(labels = label_comma()) +
    scale_colour_manual(values = c(S2, S1)) +
    labs(x = "Expression of exclusive genes (CPM, present condition)",
         y = "Genes (log scale)",
         title = "B  What the retention rule actually admits",
         caption = paste("The rule's CPM thresholds (0.05 / 0.10) fall below one",
                         "read at CTC library sizes, so it admits\nthousands of",
                         "genes that are zero in one condition by chance.")) +
    theme_ctc()

  ggsave(file.path(FIG, "fig3_exclusive.png"), pA / pB, width = 8, height = 7,
         dpi = 300, bg = "white")
  cat("fig3 written\n")
}

# =====================================================================
# Figure 4 - depth dependence of the filter
# =====================================================================
fig4 <- function() {
  f <- file.path(ROOT, "results/depth.csv"); if (!file.exists(f)) return(invisible())
  d <- read.csv(f)
  a <- aggregate(cbind(n_tested, sensitivity, avg_precision, median_lib) ~
                   arm + lib_lo, d, mean, na.rm = TRUE)
  a$armf <- factor(c(ctc_pipeline = "v1.0.0 filter call (legacy)",
                     fix_filter = "corrected: true 0.1 CPM cutoff",
                     ablate_filter = "edgeR defaults")[a$arm],
                   levels = c("v1.0.0 filter call (legacy)",
                              "corrected: true 0.1 CPM cutoff", "edgeR defaults"))
  a$eff_cpm <- ifelse(a$arm == "ctc_pipeline", 0.1 / a$median_lib * 1e6,
               ifelse(a$arm == "fix_filter", 0.1, 10 / a$median_lib * 1e6))

  pA <- ggplot(a, aes(median_lib, eff_cpm, colour = armf, group = armf)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
    scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
    scale_y_log10(labels = label_number(accuracy = 0.001)) +
    scale_colour_manual(values = c(S1, S2, S3)) +
    labs(x = "Median library size (reads)", y = "Effective CPM cutoff",
         title = "A  The shipped filter's threshold moves with sequencing depth") +
    theme_ctc()

  pB <- ggplot(a, aes(median_lib, n_tested, colour = armf, group = armf)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
    scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
    scale_y_continuous(labels = label_comma()) +
    scale_colour_manual(values = c(S1, S2, S3)) +
    labs(x = "Median library size (reads)", y = "Genes tested",
         title = "B  Consequence for the number of genes carried into testing",
         caption = paste("min.count is a COUNT in filterByExpr; passing --min_cpm",
                         "to it yields a cutoff of min_cpm / medianLib x 1e6.")) +
    theme_ctc()

  ggsave(file.path(FIG, "fig4_depth.png"), pA / pB, width = 8, height = 7,
         dpi = 300, bg = "white")
  cat("fig4 written\n")
}

# =====================================================================
# Figure 5 - null calibration
# =====================================================================
fig5 <- function() {
  f <- file.path(ROOT, "results/null.csv"); if (!file.exists(f)) return(invisible())
  d <- read.csv(f)
  a <- aggregate(n_significant ~ arm + n_per_group, d, mean, na.rm = TRUE)
  a$label <- factor(arm_label[a$arm],
                    levels = arm_label[names(sort(tapply(a$n_significant, a$arm, mean)))])
  a$class <- arm_class(a$arm)
  a$facet <- factor(sprintf("n = %d per group", a$n_per_group))
  p <- ggplot(a, aes(n_significant, label, colour = class)) +
    geom_point(size = 2.2) + facet_wrap(~facet, nrow = 1) +
    scale_colour_manual(values = c("CTC pipeline" = S1, "edgeR defaults" = S2,
                                   "batch variant" = S3, "ablation" = MUTED)) +
    labs(x = "False positives per experiment (no true DE present)", y = NULL,
         title = "Calibration under a complete null",
         subtitle = "50 simulations with zero differentially expressed genes") +
    theme_ctc()
  ggsave(file.path(FIG, "fig5_null.png"), p, width = 9.5, height = 4.5,
         dpi = 300, bg = "white")
  cat("fig5 written\n")
}

# =====================================================================
# Figure 6 - purity confounding
# =====================================================================
fig6 <- function() {
  f <- file.path(ROOT, "results/purity.csv"); if (!file.exists(f)) return(invisible())
  d <- read.csv(f)
  a <- aggregate(cbind(empirical_fdr, sensitivity, avg_precision) ~
                   arm + purity_cond_diff, d, mean, na.rm = TRUE)
  s <- aggregate(empirical_fdr ~ arm + purity_cond_diff, d, se)
  a$fdr_se <- s$empirical_fdr
  a$armf <- factor(c(purity_ignored   = "Purity ignored (v1.1.0)",
                     purity_covariate = "Purity as covariate",
                     purity_filtered  = "Low-purity samples dropped")[a$arm],
                   levels = c("Purity ignored (v1.1.0)", "Purity as covariate",
                              "Low-purity samples dropped"))
  p <- ggplot(a, aes(purity_cond_diff, empirical_fdr, colour = armf, group = armf)) +
    geom_hline(yintercept = 0.10, linetype = "22", colour = INK2, linewidth = 0.4) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.4) +
    geom_errorbar(aes(ymin = empirical_fdr - fdr_se, ymax = empirical_fdr + fdr_se),
                  width = 0, linewidth = 0.5, alpha = 0.7) +
    scale_colour_manual(values = c(S2, S1, S3)) +
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    labs(x = "Difference in mean tumour purity between conditions",
         y = "Empirical false discovery rate",
         title = "Modelling purity restores FDR control under contamination confounding",
         subtitle = paste("Nominal FDR 0.10 (dashed). Purity is estimated from",
                          "the counts, never supplied. n = 5 per group, 50 simulations."),
         caption = paste("Dropping low-purity samples removes data without removing",
                         "the confound; entering purity in the\ndesign matrix holds",
                         "at nominal and costs almost nothing when no confounding exists.")) +
    theme_ctc()
  ggsave(file.path(FIG, "fig6_purity.png"), p, width = 8, height = 4.6,
         dpi = 300, bg = "white")
  cat("fig6 written\n")
}

# =====================================================================
# Figure 7 - purity score validated on real data (GSE67980)
# =====================================================================
fig7 <- function() {
  f <- file.path(ROOT, "results/real_purity.csv"); if (!file.exists(f)) return(invisible())
  d <- read.csv(f)
  lab <- c("white blood cell" = "Leukocyte",
           "candidate PCa CTC" = "Candidate CTC",
           "primary PCa tumor" = "Primary tumour",
           "lineage-confirmed PCa CTC" = "Confirmed CTC",
           "single cell from PCa cell line" = "Cell line")
  d$grp <- lab[d$type]
  d <- d[!is.na(d$grp) & !is.na(d$purity_score), ]
  ord <- names(sort(tapply(d$purity_score, d$grp, median)))
  d$grp <- factor(d$grp, levels = ord)
  med <- aggregate(purity_score ~ grp, d, median)

  p <- ggplot(d, aes(purity_score, grp)) +
    geom_jitter(height = 0.16, size = 1.5, alpha = 0.55, colour = S1) +
    geom_point(data = med, size = 3.2, shape = 18, colour = INK) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(x = "Signature-based tumour purity score", y = NULL,
         title = "Purity score against known sample type (GSE67980)",
         subtitle = paste("169 prostate samples. Diamonds are group medians.",
                          "AUC 0.990 separating confirmed CTCs from leukocytes."),
         caption = paste("The score was not shown the labels. It nevertheless",
                         "separates the authors' lineage-confirmed CTCs from\ntheir",
                         "unconfirmed candidates, and places primary tumours below",
                         "CTCs as stromal infiltrate predicts.")) +
    theme_ctc()
  ggsave(file.path(FIG, "fig7_real_purity.png"), p, width = 8, height = 4.4,
         dpi = 300, bg = "white")
  cat("fig7 written\n")
}

for (f in list(fig1, fig2, fig3, fig4, fig5, fig6, fig7))
  tryCatch(f(), error = function(e) cat("figure failed:", conditionMessage(e), "\n"))
