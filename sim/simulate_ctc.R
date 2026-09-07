# ============================================================
# CTC low-input RNA-seq count simulator with known ground truth
#
# Generative model (per replicate):
#   1. Gene-level relative abundance      p_g ~ LogNormal, sum-normalised
#   2. Condition effect on a known DE set  mu_gB = mu_gA * 2^lfc_g
#   3. Condition-exclusive genes           mu_gA > 0, mu_gB = 0  (true DE)
#   4. Per-gene batch factor               exp(N(0, batch_sd)) per batch
#   5. Library size per sample             drawn from a CTC tier range
#   6. Counts                              NB(mu = lib * p, disp = BCV trend)
#   7. Low-input dropout                   Bernoulli zero-inflation, biased
#                                          toward low-abundance genes
#
# Every quantity a benchmark needs to score against is returned in `truth`.
# ============================================================

suppressPackageStartupMessages(library(stats))

# --- mean-dispersion trend -------------------------------------------------
# BCV^2 = disp = asymptotic_disp + trend_coef / mu
# Low-input CTC data has far higher BCV than bulk tissue (~0.6 vs ~0.15).
ctc_dispersion <- function(mu, bcv_asymptotic = 0.55, trend_coef = 3.0) {
  bcv_asymptotic^2 + trend_coef / pmax(mu, 0.05)
}

simulate_ctc <- function(n_genes          = 20000,
                         n_per_group      = 4,
                         conditions       = c("tumour", "control"),
                         prop_de          = 0.08,
                         n_exclusive      = 120,
                         exclusive_cpm    = 1.0,
                         lfc_min          = 0.5,
                         lfc_rate         = 1.0,
                         lib_size_range   = c(3e5, 1.5e6),
                         n_batches        = 2,
                         batch_sd         = 0.35,
                         batch_balanced   = TRUE,
                         bcv_asymptotic   = 0.45,
                         trend_coef       = 3.0,
                         dropout_max      = 0.60,
                         dropout_scale    = 15,
                         n_blanks         = 0,
                         purity_mean      = 1.0,
                         purity_sd        = 0.0,
                         purity_cond_diff = 0.0,
                         seed             = NULL) {

  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(conditions) == 2, n_per_group >= 2)

  gene_ids <- sprintf("ENSG%011d", seq_len(n_genes))

  # ---- 1. baseline relative abundance ------------------------------------
  # Log-normal on the log10 scale reproduces the usual RNA-seq abundance
  # spread: a few thousand well-expressed genes, a long low-abundance tail.
  log10_abund <- rnorm(n_genes, mean = -0.6, sd = 1.15)
  p_base      <- 10^log10_abund
  p_base      <- p_base / sum(p_base)

  # ---- 2. true DE set -----------------------------------------------------
  n_de <- round(prop_de * n_genes)
  # Exclusive genes are assigned an explicit abundance (expected CPM in the
  # condition where they are present) so that detectability can be swept
  # independently of the background abundance distribution. Drawing them from
  # the low-abundance tail instead would make most of them unobservable at CTC
  # depth, which tests sequencing depth rather than the retention rule.
  excl_idx     <- sample(seq_len(n_genes), n_exclusive)
  de_candidate <- setdiff(seq_len(n_genes), excl_idx)
  de_idx       <- sample(de_candidate, n_de)

  true_lfc          <- numeric(n_genes)
  true_lfc[de_idx]  <- (lfc_min + rexp(n_de, rate = lfc_rate)) *
                         sample(c(-1, 1), n_de, replace = TRUE)

  is_de              <- logical(n_genes)
  is_de[de_idx]      <- TRUE
  is_de[excl_idx]    <- TRUE          # exclusive genes are genuinely DE
  is_exclusive       <- logical(n_genes)
  is_exclusive[excl_idx] <- TRUE

  # ---- 3. per-condition abundance ----------------------------------------
  p_cond <- matrix(p_base, nrow = n_genes, ncol = 2,
                   dimnames = list(gene_ids, conditions))
  p_cond[, 2] <- p_base * 2^true_lfc
  # exclusive: present in condition 1 at a specified CPM, structurally absent
  # from condition 2 (e.g. an epithelial marker against a WBC background)
  p_cond[excl_idx, 1] <- exclusive_cpm / 1e6
  p_cond[excl_idx, 2] <- 0
  # true LFC for exclusive genes is undefined (-Inf); record as NA
  true_lfc[excl_idx] <- NA_real_

  # ---- 4. sample sheet ----------------------------------------------------
  n_real  <- 2 * n_per_group
  cond_v  <- rep(conditions, each = n_per_group)

  if (n_batches > 1) {
    if (batch_balanced) {
      # each condition split as evenly as possible across batches
      batch_v <- unlist(lapply(seq_along(conditions), function(i)
        rep_len(seq_len(n_batches), n_per_group)))
    } else {
      # condition 1 concentrated in batch 1 (75/25), condition 2 the reverse
      k <- max(1, round(0.75 * n_per_group))
      batch_v <- c(c(rep(1, k), rep_len(2:n_batches, n_per_group - k)),
                   c(rep(2, k), rep_len(c(1, seq_len(n_batches)[-(1:2)]),
                                        n_per_group - k)))
      batch_v <- batch_v[seq_len(n_real)]
    }
  } else {
    batch_v <- rep(1L, n_real)
  }
  batch_v <- paste0("batch", batch_v)

  samples <- data.frame(
    sample_id    = sprintf("S%02d", seq_len(n_real)),
    condition    = factor(cond_v, levels = conditions),
    batch        = factor(batch_v),
    cell_count   = 10L,
    control_type = "",
    stringsAsFactors = FALSE
  )
  # placeholder; filled in after purity is drawn below
  samples$purity_true <- 1

  # ---- 5. per-gene batch factors -----------------------------------------
  ub <- levels(samples$batch)
  batch_fac <- matrix(1, nrow = n_genes, ncol = length(ub),
                      dimnames = list(gene_ids, ub))
  if (length(ub) > 1 && batch_sd > 0) {
    for (b in seq_along(ub))
      batch_fac[, b] <- exp(rnorm(n_genes, 0, batch_sd))
  }

  # ---- 5b. leukocyte background and per-sample purity ---------------------
  # A CTC library is a mixture of tumour and co-enriched leukocyte RNA. The
  # leukocyte compartment has its own abundance profile and carries NO
  # condition effect: differences driven by varying purity are confounding,
  # not biology, and are scored as false positives.
  p_wbc <- 10^rnorm(n_genes, mean = -0.6, sd = 1.15)
  p_wbc <- p_wbc / sum(p_wbc)

  purity <- rep(1, n_real)
  if (purity_mean < 1 || purity_sd > 0 || purity_cond_diff != 0) {
    shift  <- ifelse(cond_v == conditions[1], purity_cond_diff / 2,
                     -purity_cond_diff / 2)
    purity <- pmin(0.99, pmax(0.01, rnorm(n_real, purity_mean + shift, purity_sd)))
  }

  # ---- 6. counts ----------------------------------------------------------
  lib_sizes <- runif(n_real, lib_size_range[1], lib_size_range[2])
  counts    <- matrix(0L, nrow = n_genes, ncol = n_real,
                      dimnames = list(gene_ids, samples$sample_id))

  for (j in seq_len(n_real)) {
    ci <- match(as.character(samples$condition[j]), conditions)
    bi <- match(as.character(samples$batch[j]), ub)

    p_j <- p_cond[, ci] * batch_fac[, bi]
    p_j <- p_j / sum(p_j)
    if (purity[j] < 1) p_j <- purity[j] * p_j + (1 - purity[j]) * p_wbc
    mu  <- lib_sizes[j] * p_j

    disp <- ctc_dispersion(mu, bcv_asymptotic, trend_coef)
    cnt  <- rnbinom(n_genes, mu = mu, size = 1 / disp)

    # low-input dropout: capture failure, concentrated on low-abundance genes
    if (dropout_max > 0) {
      p_drop <- dropout_max * exp(-mu / dropout_scale)
      cnt[runif(n_genes) < p_drop] <- 0L
    }
    counts[, j] <- as.integer(cnt)
  }

  # ---- 7. optional blank capture wells ------------------------------------
  if (n_blanks > 0) {
    blank_ids <- sprintf("BLANK%02d", seq_len(n_blanks))
    bcounts   <- matrix(0L, nrow = n_genes, ncol = n_blanks,
                        dimnames = list(gene_ids, blank_ids))
    for (j in seq_len(n_blanks)) {
      mu <- runif(1, 1e4, 4e4) * p_base
      bcounts[, j] <- as.integer(rnbinom(n_genes, mu = mu,
                                         size = 1 / ctc_dispersion(mu)))
    }
    counts  <- cbind(counts, bcounts)
    samples <- rbind(samples, data.frame(
      sample_id    = blank_ids,
      condition    = factor(conditions[1], levels = conditions),
      batch        = factor(levels(samples$batch)[1],
                            levels = levels(samples$batch)),
      cell_count   = 0L,
      control_type = "blank",
      stringsAsFactors = FALSE))
  }

  samples$purity_true <- c(purity, rep(NA_real_, nrow(samples) - n_real))

  truth <- data.frame(
    gene_id      = gene_ids,
    is_de        = is_de,
    is_exclusive = is_exclusive,
    true_lfc     = true_lfc,
    base_mean    = p_base * mean(lib_sizes),
    stringsAsFactors = FALSE
  )

  # Marker sets a purity estimator can use, defined exactly as real signatures
  # are: genes most specific to one compartment. Returned so the estimator can
  # be run on simulated data without ever seeing the true purity values.
  lr  <- log2((p_wbc + 1e-9) / (p_base + 1e-9))
  sig <- list(leukocyte  = gene_ids[order(lr, decreasing = TRUE)][1:60],
              epithelial = gene_ids[order(lr)][1:60])

  list(counts  = counts,
       samples = samples,
       truth   = truth,
       signatures = sig,
       params  = list(n_genes = n_genes, n_per_group = n_per_group,
                      prop_de = prop_de, n_exclusive = n_exclusive,
                      lib_size_range = lib_size_range, n_batches = n_batches,
                      batch_sd = batch_sd, batch_balanced = batch_balanced,
                      exclusive_cpm = exclusive_cpm,
                      purity_mean = purity_mean, purity_sd = purity_sd,
                      purity_cond_diff = purity_cond_diff,
                      bcv_asymptotic = bcv_asymptotic,
                      dropout_max = dropout_max, seed = seed))
}
