# adapt_crosssectional.R
#
# Ablated variant of ADAPT-ME with no subject random intercept.
# Every observation is treated as independent (cross-sectional Tobit).
# Used as a comparison baseline to isolate the contribution of the
# mixed-effects model in ADAPT-ME.
#
# ── ISOLATION GUARANTEE ──────────────────────────────────────────────────────
# This file introduces ONLY new symbols. It does NOT modify, override, or
# call into adapt(), count_ratio(), fit_tobit_mixed(), or preprocess().
# BUM_fit() and BUM_llk() are shared pure-utility functions that are
# called but not modified.
# ─────────────────────────────────────────────────────────────────────────────

# ── Fixed-effects Tobit log-likelihood (no random intercept) ─────────────────

cs_tobit_loglik <- function(par, y, censor_point, censored, X) {
  beta  <- par[seq_len(ncol(X))]
  sigma <- exp(par[ncol(X) + 1L])
  eta   <- drop(X %*% beta)
  ll    <- 0
  exact <- !censored
  if (any(exact))
    ll <- ll + sum(dnorm(y[exact], mean = eta[exact], sd = sigma, log = TRUE))
  if (any(censored))
    ll <- ll + sum(pnorm(censor_point[censored], mean = eta[censored],
                         sd = sigma, log.p = TRUE))
  if (!is.finite(ll)) return(.Machine$double.xmax^0.5)
  -ll
}

fit_tobit_cs <- function(y, censor_point, censored, X) {
  X <- as.matrix(X)
  if (ncol(X) == 0 || nrow(X) != length(y))
    stop("fit_tobit_cs: invalid input dimensions")
  if (qr(X)$rank < ncol(X))
    stop("fit_tobit_cs: model matrix is rank deficient")

  lm_fit     <- lm.fit(x = X, y = y)
  beta_start <- lm_fit$coefficients
  beta_start[is.na(beta_start)] <- 0
  beta_start <- beta_start[seq_len(ncol(X))]
  sigma_start <- sd(lm_fit$residuals)
  if (!is.finite(sigma_start) || sigma_start <= 0) sigma_start <- sd(y)
  if (!is.finite(sigma_start) || sigma_start <= 0) sigma_start <- 1

  start <- c(beta_start, log(sigma_start))
  fit   <- optim(start, cs_tobit_loglik,
                 y = y, censor_point = censor_point,
                 censored = censored, X = X,
                 method = "BFGS", control = list(maxit = 300))
  if (fit$convergence != 0 || !is.finite(fit$value))
    stop("fit_tobit_cs: optimisation did not converge")
  list(beta = fit$par[seq_len(ncol(X))], logLik = -fit$value)
}

# ── Count-ratio loop using fixed Tobit ───────────────────────────────────────

count_ratio_cs <- function(count_table, design_matrix,
                            reftaxa = NULL, censor = 1, test_all = FALSE) {
  alltaxa  <- colnames(count_table)
  if (is.null(reftaxa)) reftaxa <- alltaxa
  TBDtaxa  <- if (test_all) alltaxa else reftaxa

  refcounts   <- rowSums(count_table[, reftaxa, drop = FALSE])
  null_filter <- refcounts != 0
  refcounts   <- refcounts[null_filter]
  subset_dm   <- as.matrix(design_matrix[null_filter, , drop = FALSE])
  TBD_counts  <- count_table[null_filter, TBDtaxa, drop = FALSE]
  existence   <- 1L * (TBD_counts > 0)
  TBD_counts[TBD_counts == 0] <- censor
  prevalences <- colMeans(existence)

  CR_result <- data.frame(
    Taxa             = TBDtaxa,
    prevalence       = prevalences,
    log10foldchange  = NA_real_,
    teststat         = NA_real_,
    pval             = NA_real_,
    within_subject_sd = NA_real_,   # NA for CS; kept for DAresult compatibility
    stringsAsFactors = FALSE
  )

  condition_col <- match("condition", colnames(subset_dm))
  if (is.na(condition_col))
    stop("count_ratio_cs: design matrix must contain a 'condition' column")
  null_dm <- subset_dm[, -condition_col, drop = FALSE]

  for (i in seq_len(ncol(TBD_counts))) {
    y_vec        <- as.numeric(log(TBD_counts[, i]) - log(refcounts))
    censor_point <- log(censor) - log(refcounts)
    censored     <- existence[, i] == 0

    fit_result <- tryCatch({
      full_fit <- fit_tobit_cs(y_vec, censor_point, censored, subset_dm)
      null_fit <- fit_tobit_cs(y_vec, censor_point, censored, null_dm)
      lrt      <- max(0, 2 * (full_fit$logLik - null_fit$logLik))
      list(beta = full_fit$beta[condition_col], teststat = lrt)
    }, error = function(e) NULL)

    if (!is.null(fit_result)) {
      CR_result$log10foldchange[i] <- fit_result$beta / log(10)
      CR_result$teststat[i]        <- fit_result$teststat
      CR_result$pval[i]            <- 1 - pchisq(fit_result$teststat, 1)
    }
  }

  failed <- sum(is.na(CR_result$pval))
  if (failed > 0)
    warning(sprintf("count_ratio_cs: %d taxa failed to fit (NA p-values)", failed))
  CR_result
}

# ── Minimal preprocessing (no subject.var) ───────────────────────────────────

preprocess_cs <- function(input_data, cond.var, base.cond, adj.var,
                           prev.filter, depth.filter) {
  stopifnot("Input must be a phyloseq object" =
              is(input_data, "phyloseq"))
  stopifnot("prev.filter must be in [0, 1)" =
              prev.filter >= 0 && prev.filter < 1)
  stopifnot("depth.filter must be nonnegative" =
              depth.filter >= 0)

  subset_data <- filter_taxa(input_data,
                              function(x) mean(x > 0) > prev.filter, TRUE)
  subset_data <- prune_samples(sample_sums(subset_data) > depth.filter,
                                subset_data)

  metadata <- data.frame(sample_data(subset_data))
  allcols  <- colnames(metadata)

  stopifnot("cond.var must be a single character string" =
              is.character(cond.var) && length(cond.var) == 1)
  stopifnot("adj.var must be NULL or a character vector" =
              is.null(adj.var) || is.character(adj.var))

  selected_cols <- unique(c(cond.var, adj.var))
  missing_cols  <- selected_cols[!selected_cols %in% allcols]
  if (length(missing_cols) > 0)
    stop(sprintf("preprocess_cs: columns not found: %s",
                 paste(missing_cols, collapse = ", ")))

  subset_metadata <- metadata[, selected_cols, drop = FALSE]
  if (any(is.na(subset_metadata)))
    stop("preprocess_cs: no missing values allowed in metadata")

  main_variable <- subset_metadata[, cond.var]
  if (length(unique(main_variable)) == 1)
    stop("preprocess_cs: all samples share the same condition")

  if (!is.numeric(main_variable)) {
    if (is.null(base.cond) || !base.cond %in% main_variable)
      base.cond <- unique(main_variable)[1]
    cat(sprintf("Choose '%s' as the baseline condition\n", base.cond))
    others <- if (length(unique(main_variable)) == 2)
      setdiff(unique(main_variable), base.cond) else "others"
    main_variable <- as.integer(main_variable != base.cond)
    cond.var <- sprintf("%s (%s VS %s)", cond.var, others, base.cond)
  }

  adjustments <- NULL
  if (!is.null(adj.var)) {
    adj_df      <- subset_metadata[, adj.var, drop = FALSE]
    adjustments <- model.matrix(~., data = adj_df)[, -1, drop = FALSE]
  }

  design_matrix <- cbind(1, main_variable, adjustments)
  colnames(design_matrix)[1:2] <- c("(Intercept)", "condition")

  count_table <- otu_table(subset_data)
  if (taxa_are_rows(subset_data)) count_table <- t(count_table)
  count_table <- as(count_table, "matrix")
  if (any(is.na(count_table)))
    stop("preprocess_cs: no missing values allowed in count table")

  cat(sprintf("%d taxa and %d samples (cross-sectional)...\n",
              ncol(count_table), nrow(count_table)))
  list(count_table   = count_table,
       design_matrix = design_matrix,
       DAAname       = cond.var)
}

# ── Main function ─────────────────────────────────────────────────────────────

#' ADAPT-ME Cross-Sectional Ablation
#'
#' @description
#' Runs the full ADAPT-ME pipeline (BUM-guided reference selection +
#' Tobit count-ratio tests) but with **no subject random intercept**.
#' Every observation is treated as an independent draw, identical to the
#' behaviour of the original cross-sectional ADAPT framework.
#'
#' This function exists solely as an ablation baseline. It is intentionally
#' isolated: it shares no internal code paths with `adapt()` and does not
#' require a `subject.var`.
#'
#' The reference score uses only median-distance (original ADAPT criterion),
#' not the stability term added by ADAPT-ME, so the comparison isolates the
#' single methodological difference: **random effects vs. no random effects**.
#'
#' @param input_data a phyloseq object
#' @param cond.var condition variable name (character)
#' @param base.cond baseline level for categorical conditions
#' @param adj.var adjustment variable names (character vector or NULL)
#' @param censor censoring value for zero counts (default 1)
#' @param prev.filter prevalence filter threshold (default 0.05)
#' @param depth.filter minimum library size (default 1000)
#' @param alpha BH-adjusted p-value cutoff (default 0.05)
#'
#' @returns a `DAresult` object (same class as `adapt()`); use `summary()`
#'   and `plot()` normally
#' @export
adapt_crosssectional <- function(input_data, cond.var, base.cond = NULL,
                                  adj.var = NULL, censor = 1,
                                  prev.filter = 0.05, depth.filter = 1000,
                                  alpha = 0.05) {
  stopifnot("censor must be positive"           = censor > 0)
  stopifnot("alpha must be between 0 and 0.5"   = alpha > 0 && alpha < 0.5)

  pp          <- preprocess_cs(input_data, cond.var, base.cond, adj.var,
                                prev.filter, depth.filter)
  count_table   <- pp$count_table
  design_matrix <- pp$design_matrix
  taxa_names    <- colnames(count_table)

  # Reference selection — median fold-change only (original ADAPT criterion)
  reftaxa <- taxa_names
  cat("Selecting Reference Set (CS)... ")
  while (TRUE) {
    relabd <- count_ratio_cs(count_table, design_matrix,
                              reftaxa = reftaxa, censor = censor,
                              test_all = FALSE)
    estimated_effect <- relabd$log10foldchange
    pvals            <- relabd$pval
    names(pvals) <- names(estimated_effect) <- relabd$Taxa

    bumfit <- BUM_fit(pvals)
    loglik <- BUM_llk(bumfit$estim_params, pvals[!is.na(pvals)])
    if (2 * loglik > qchisq(0.95, 1)) {
      distance2med <- abs(estimated_effect -
                            median(estimated_effect, na.rm = TRUE))
      dist_scale   <- median(distance2med, na.rm = TRUE)
      if (!is.finite(dist_scale) || dist_scale <= 0) dist_scale <- 1
      ref_score <- distance2med / dist_scale
      ref_score[!is.finite(ref_score)] <- Inf
      sorted_score <- sort(ref_score)
      reftaxa <- names(sorted_score)[seq_len(length(sorted_score) / 2)]
    } else {
      break
    }
  }
  cat(sprintf("%d reference taxa selected\n", length(reftaxa)))

  # Final test pass across all taxa
  cat("\nRunning fixed-effects Tobit models across taxa...\n")
  all_results  <- count_ratio_cs(count_table, design_matrix,
                                  reftaxa = reftaxa, censor = censor,
                                  test_all = TRUE)
  all_pvals    <- all_results$pval
  names(all_pvals) <- all_results$Taxa
  all_adj      <- p.adjust(all_pvals, method = "BH")
  all_results$adjusted_pval <- all_adj

  sig      <- all_adj[!is.na(all_adj) & all_adj < alpha]
  DiffTaxa <- if (length(sig) > 0) names(sig) else ""
  n_da     <- if (identical(DiffTaxa, "")) 0L else length(DiffTaxa)
  cat(sprintf("%d differentially abundant taxa detected (cross-sectional)\n",
              n_da))

  new("DAresult",
      DAAname   = paste0(pp$DAAname, " [CS]"),
      reference = reftaxa,
      signal    = DiffTaxa,
      details   = all_results,
      input     = input_data)
}
