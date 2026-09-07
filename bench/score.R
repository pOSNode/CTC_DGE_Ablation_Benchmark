# Scoring against simulated ground truth.
#
# Sensitivity is computed over ALL true-DE genes in the simulation, including
# those a filter discarded. A gene removed before testing is a missed
# discovery, not an excluded one - otherwise aggressive filtering appears free.

score_arm <- function(res, truth, low_abundance_q = 0.25) {

  tt  <- res$table
  all_genes <- truth$gene_id
  sig <- setNames(rep(FALSE, length(all_genes)), all_genes)
  sig[tt$gene_id] <- tt$significant

  # p-values; genes never tested rank last
  pv <- setNames(rep(1, length(all_genes)), all_genes)
  pv[tt$gene_id] <- tt$PValue

  is_de   <- truth$is_de
  n_de    <- sum(is_de)
  TP      <- sum(sig &  is_de)
  FP      <- sum(sig & !is_de)
  n_sig   <- TP + FP

  # threshold-free average precision over the full truth set
  ord   <- order(pv, runif(length(pv)))
  hit   <- is_de[ord]
  prec  <- cumsum(hit) / seq_along(hit)
  ap    <- if (n_de > 0) sum(prec * hit) / n_de else NA_real_

  # stratified sensitivity
  lowq  <- truth$base_mean <= quantile(truth$base_mean, low_abundance_q)
  excl  <- truth$is_exclusive
  de_lo <- is_de & lowq & !excl

  data.frame(
    n_tested          = res$n_tested,
    n_significant     = n_sig,
    true_positives    = TP,
    false_positives   = FP,
    n_true_de         = n_de,
    sensitivity       = TP / n_de,
    precision         = if (n_sig > 0) TP / n_sig else NA_real_,
    empirical_fdr     = if (n_sig > 0) FP / n_sig else NA_real_,
    nominal_fdr       = res$params$fdr_cutoff,
    avg_precision     = ap,
    sens_exclusive    = if (sum(excl)  > 0) sum(sig & excl)  / sum(excl)  else NA_real_,
    sens_low_abund    = if (sum(de_lo) > 0) sum(sig & de_lo) / sum(de_lo) else NA_real_,
    # true exclusive genes that reached significance AND were admitted only by
    # the retention rule - i.e. discoveries the base filter would have lost
    sig_excl_rescued  = sum(sig & excl & (all_genes %in% res$excl_only_ids)),
    n_excl_kept       = sum(excl & (all_genes %in% res$genes_kept)),
    n_exclusive_added = res$n_exclusive_added,
    excl_admit_fpr    = if (res$n_exclusive_added > 0)
                          1 - sum(excl & (all_genes %in% res$excl_only_ids)) /
                              res$n_exclusive_added else NA_real_,
    norm_used         = res$norm_method_used,
    batch_in_model    = res$batch_in_model,
    shrink_applied    = res$shrink_applied,
    common_disp       = res$common_dispersion,
    stringsAsFactors  = FALSE
  )
}
