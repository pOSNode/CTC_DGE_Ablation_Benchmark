# ============================================================
# Comparator arms from other established DE tools, each run at its own
# documented defaults, returning the same structure as ctc_dge_core() so that
# score_arm() can score them identically.
#
# The point is to test whether the CTC pipeline's advantage holds against the
# tools an analyst would actually reach for, not only against edgeR's defaults.
# ============================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})

# shared shell so every arm reports the same fields
.as_core_result <- function(tab, kept, n_all, fdr_cutoff, extra = list()) {
  tab$significant <- tab$FDR < fdr_cutoff
  c(list(table = tab,
         n_tested = nrow(tab),
         n_significant = sum(tab$significant, na.rm = TRUE),
         genes_kept = kept,
         keep_base_ids = kept,
         excl_only_ids = character(0),
         n_exclusive_added = 0L,
         norm_method_used = NA_character_,
         filter_mode_used = "default",
         filter_cpm_cutoff = NA_real_,
         batch_in_model = TRUE,
         batch_confounded = FALSE,
         purity_in_model = FALSE,
         shrink_applied = FALSE,
         common_dispersion = NA_real_,
         params = list(fdr_cutoff = fdr_cutoff)),
    extra)
}

.design_of <- function(samples) {
  use_batch <- nlevels(droplevels(samples$batch)) > 1 &&
               !any(rowSums(table(samples$condition, samples$batch) > 0) == 1)
  d <- if (use_batch) model.matrix(~ 0 + condition + batch, data = samples)
       else           model.matrix(~ 0 + condition, data = samples)
  colnames(d) <- make.names(colnames(d))
  list(design = d, use_batch = use_batch)
}

.drop_controls <- function(counts, samples) {
  if (!is.null(samples$control_type)) {
    k <- !samples$control_type %in% c("blank", "positive")
    counts <- counts[, k, drop = FALSE]; samples <- samples[k, , drop = FALSE]
  }
  samples$condition <- droplevels(as.factor(samples$condition))
  samples$batch     <- droplevels(as.factor(samples$batch))
  list(counts = counts, samples = samples)
}

# ---- limma-voom, standard workflow ---------------------------------------
# filterByExpr defaults -> TMM -> voom -> lmFit -> eBayes
voom_arm <- function(counts, samples, params = list(), ...) {
  fdr <- if (is.null(params$fdr_cutoff)) 0.05 else params$fdr_cutoff
  d <- .drop_controls(counts, samples); counts <- d$counts; samples <- d$samples
  dd <- .design_of(samples); design <- dd$design

  y <- DGEList(counts = counts)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")

  v   <- voom(y, design)
  fit <- eBayes(lmFit(v, design))
  cond <- grep("^condition", colnames(design), value = TRUE)
  ct <- rep(0, ncol(design)); names(ct) <- colnames(design)
  ct[cond[1]] <- 1; ct[cond[2]] <- -1
  tt <- topTable(eBayes(contrasts.fit(fit, ct)), number = Inf, sort.by = "P")
  tt <- data.frame(gene_id = rownames(tt), logFC = tt$logFC,
                   logCPM = tt$AveExpr, PValue = tt$P.Value,
                   FDR = tt$adj.P.Val, stringsAsFactors = FALSE)
  .as_core_result(tt, rownames(y), nrow(counts), fdr,
                  list(norm_method_used = "TMM (voom)"))
}

# ---- DESeq2, standard workflow -------------------------------------------
# vignette pre-filter -> median-of-ratios -> Wald test with independent filtering
deseq2_arm <- function(counts, samples, params = list(), ...) {
  if (!requireNamespace("DESeq2", quietly = TRUE))
    stop("DESeq2 not installed")
  fdr <- if (is.null(params$fdr_cutoff)) 0.05 else params$fdr_cutoff
  d <- .drop_controls(counts, samples); counts <- d$counts; samples <- d$samples
  dd <- .design_of(samples)

  fml <- if (dd$use_batch) ~ batch + condition else ~ condition
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                        colData = samples, design = fml)
  # the light pre-filter the DESeq2 vignette recommends
  keep <- rowSums(DESeq2::counts(dds)) >= 10
  dds  <- dds[keep, ]
  dds  <- DESeq2::DESeq(dds, quiet = TRUE, fitType = "local")
  lv   <- levels(samples$condition)
  res  <- DESeq2::results(dds, contrast = c("condition", lv[1], lv[2]),
                          alpha = fdr)
  res  <- as.data.frame(res)
  tt <- data.frame(gene_id = rownames(res), logFC = res$log2FoldChange,
                   logCPM = log2(res$baseMean + 1), PValue = res$pvalue,
                   FDR = res$padj, stringsAsFactors = FALSE)
  # genes dropped by independent filtering have NA padj: not discoveries
  tt$FDR[is.na(tt$FDR)] <- 1
  tt$PValue[is.na(tt$PValue)] <- 1
  tt <- tt[order(tt$PValue), ]
  .as_core_result(tt, rownames(res), nrow(counts), fdr,
                  list(norm_method_used = "median-of-ratios (DESeq2)"))
}

external_arms <- list(limma_voom = voom_arm, deseq2 = deseq2_arm)
