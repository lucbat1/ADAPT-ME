#' @importFrom stats aggregate as.formula ave plogis rbinom rnorm rpois runif setNames
NULL

#' ADAPT-ME Validation Scenarios
#'
#' @description
#' Standard scenario definitions for simulation-based validation of ADAPT-ME.
#'
#' @returns A named list of simulation parameter lists.
#' @export
adaptme_validation_scenarios <- function() {
  list(
    null = list(n_da = 0, effect_log10 = 0, zero_inflation = 0.05,
                unbalanced = FALSE, confounded = FALSE, condition_type = "binary"),
    strong_signal = list(n_da = 3, effect_log10 = 0.7, zero_inflation = 0.05,
                         unbalanced = FALSE, confounded = FALSE, condition_type = "binary"),
    sparse = list(n_da = 3, effect_log10 = 0.7, zero_inflation = 0.35,
                  unbalanced = FALSE, confounded = FALSE, condition_type = "binary"),
    unbalanced = list(n_da = 3, effect_log10 = 0.7, zero_inflation = 0.05,
                      unbalanced = TRUE, confounded = FALSE, condition_type = "binary"),
    confounded = list(n_da = 3, effect_log10 = 0.7, zero_inflation = 0.05,
                      unbalanced = FALSE, confounded = TRUE, condition_type = "binary"),
    continuous_time = list(n_da = 3, effect_log10 = 0.25, zero_inflation = 0.05,
                           unbalanced = FALSE, confounded = FALSE, condition_type = "continuous")
  )
}

#' Simulate Longitudinal Microbiome Data for ADAPT-ME
#'
#' @description
#' Creates a phyloseq object with repeated measurements, known subject-level
#' baselines, known condition effects, optional sparse zeros, unbalanced visits,
#' and optional confounding. The first `n_da` taxa are truly differential.
#'
#' @param seed random seed
#' @param n_subjects number of subjects
#' @param n_timepoints number of repeated measurements per subject before optional imbalance
#' @param n_taxa number of taxa
#' @param n_da number of truly differential taxa
#' @param effect_log10 condition effect for DA taxa on the log10 count-ratio scale
#' @param subject_sd standard deviation of taxon-specific subject baselines
#' @param confounder_log10 confounder effect on the log10 scale
#' @param zero_inflation additional probability of observed zero counts
#' @param unbalanced whether to randomly drop non-baseline visits
#' @param drop_fraction fraction of non-baseline visits to drop when `unbalanced = TRUE`
#' @param confounded whether the time-varying confounder is correlated with condition
#' @param condition_type either `"binary"` or `"continuous"`
#'
#' @returns A list with `phyobj`, `truth`, `metadata`, and `counts`.
#' @export
simulate_adaptme_phyloseq <- function(seed = 1, n_subjects = 20, n_timepoints = 2,
                                      n_taxa = 12, n_da = 3, effect_log10 = 0.7,
                                      subject_sd = 0.5, confounder_log10 = 0.25,
                                      zero_inflation = 0.05, unbalanced = FALSE,
                                      drop_fraction = 0.25, confounded = FALSE,
                                      condition_type = c("binary", "continuous")) {
  condition_type <- match.arg(condition_type)
  if (n_subjects < 2) stop("n_subjects must be at least 2")
  if (n_timepoints < 2) stop("n_timepoints must be at least 2")
  if (n_taxa < 4) stop("n_taxa must be at least 4")
  if (n_da < 0 || n_da >= n_taxa) stop("n_da must be nonnegative and smaller than n_taxa")
  if (zero_inflation < 0 || zero_inflation >= 1) stop("zero_inflation must be in [0, 1)")

  set.seed(seed)
  metadata <- expand.grid(
    time_index = seq_len(n_timepoints) - 1,
    subject = seq_len(n_subjects)
  )
  metadata <- metadata[order(metadata$subject, metadata$time_index), , drop = FALSE]

  if (unbalanced) {
    keep <- metadata$time_index == 0 | runif(nrow(metadata)) > drop_fraction
    metadata <- metadata[keep, , drop = FALSE]
    repeated <- ave(metadata$subject, metadata$subject, FUN = length) > 1
    metadata <- metadata[repeated, , drop = FALSE]
  }

  metadata$condition <- if (condition_type == "binary") {
    as.integer(metadata$time_index > 0)
  } else {
    metadata$time_index / max(metadata$time_index)
  }

  if (confounded) {
    # Confounded: binary covariate correlated with condition (0.75 / 0.25 probability)
    confounder_probability <- ifelse(metadata$condition > 0, 0.75, 0.25)
    metadata$confounder <- rbinom(nrow(metadata), size = 1, prob = confounder_probability)
  } else {
    # Non-confounded: observation-level continuous covariate, independent of both
    # condition and subject. Using a subject-level variable would be near-collinear
    # with the random intercept, creating multicollinearity that unfairly hurts
    # mixed-effects methods (ADAPT-ME, ANCOM-BC2, MaAsLin2) in the optimizer.
    metadata$confounder <- rnorm(nrow(metadata))
  }

  rownames(metadata) <- paste0("Sample", seq_len(nrow(metadata)))
  taxon_names <- paste0("Taxon", seq_len(n_taxa))
  da_taxa <- taxon_names[seq_len(n_da)]

  taxon_intercept <- runif(n_taxa, min = 2.7, max = 4.2)
  sparse_start <- max(n_da + 1, floor(n_taxa * 0.65))
  if (sparse_start <= n_taxa) {
    sparse_taxa <- seq.int(sparse_start, n_taxa)
    taxon_intercept[sparse_taxa] <- taxon_intercept[sparse_taxa] - 1.8
  }
  condition_effect <- rep(0, n_taxa)
  if (n_da > 0) {
    condition_effect[seq_len(n_da)] <- effect_log10 * log(10)
  }
  # Taxon-specific confounder effects so the adjustment doesn't cancel in log-ratio
  # space. A uniform effect on all taxa would vanish in any ratio-based method.
  confounder_effect <- rnorm(n_taxa, mean = 0, sd = confounder_log10 * log(10))
  subject_effect <- matrix(
    rnorm(n_subjects * n_taxa, sd = subject_sd),
    nrow = n_subjects,
    ncol = n_taxa
  )
  library_factor <- exp(rnorm(nrow(metadata), sd = 0.25))

  counts <- matrix(0, nrow = nrow(metadata), ncol = n_taxa,
                   dimnames = list(rownames(metadata), taxon_names))
  for (j in seq_len(n_taxa)) {
    eta <- taxon_intercept[j] +
      subject_effect[metadata$subject, j] +
      condition_effect[j] * metadata$condition +
      confounder_effect[j] * metadata$confounder +
      log(library_factor)
    lambda <- pmax(exp(eta), .Machine$double.eps)
    counts[, j] <- rpois(nrow(metadata), lambda = lambda)
  }

  zero_probability <- zero_inflation + (1 - zero_inflation) * plogis(1.6 - log1p(counts))
  zero_mask <- matrix(runif(length(counts)), nrow = nrow(counts)) < zero_probability
  counts[zero_mask] <- 0

  phyobj <- phyloseq::phyloseq(
    phyloseq::otu_table(counts, taxa_are_rows = FALSE),
    phyloseq::sample_data(metadata)
  )
  truth <- data.frame(
    Taxa = taxon_names,
    is_da = taxon_names %in% da_taxa,
    log10_effect = condition_effect / log(10),
    stringsAsFactors = FALSE
  )
  rownames(truth) <- truth$Taxa

  list(phyobj = phyobj, truth = truth, metadata = metadata, counts = counts)
}

#' Summarize ADAPT-ME Simulation Performance
#'
#' @param result a `DAresult` object from `adapt`
#' @param truth truth table returned by `simulate_adaptme_phyloseq`
#' @param alpha adjusted p-value cutoff
#'
#' @returns A one-row data frame with FPR, FDR, power, bias, and fit counts.
#' @export
adaptme_simulation_metrics <- function(result, truth, alpha = 0.05) {
  details <- result@details
  details <- details[details$Taxa %in% truth$Taxa, , drop = FALSE]
  truth <- truth[details$Taxa, , drop = FALSE]
  called <- !is.na(details$adjusted_pval) & details$adjusted_pval < alpha
  is_da <- truth$is_da
  null_taxa <- !is_da
  tp <- sum(called & is_da)
  fp <- sum(called & null_taxa)
  discoveries <- sum(called)
  fitted <- !is.na(details$pval)
  da_bias <- if (any(is_da & fitted)) {
    mean(details$log10foldchange[is_da & fitted] - truth$log10_effect[is_da & fitted])
  } else {
    NA_real_
  }

  data.frame(
    n_taxa = nrow(details),
    n_da = sum(is_da),
    n_fit = sum(fitted),
    discoveries = discoveries,
    true_positives = tp,
    false_positives = fp,
    fpr = if (any(null_taxa)) fp / sum(null_taxa) else NA_real_,
    fdr = if (discoveries > 0) fp / discoveries else 0,
    power = if (any(is_da)) tp / sum(is_da) else NA_real_,
    mean_da_bias = da_bias,
    stringsAsFactors = FALSE
  )
}

#' Run ADAPT-ME Simulation Validation
#'
#' @description
#' Repeats one or more validation scenarios and returns replicate-level
#' performance metrics. The defaults are intentionally small; increase
#' `n_replicates`, `n_subjects`, and `n_taxa` for manuscript-scale validation.
#'
#' @param scenarios named list of scenario parameter lists
#' @param n_replicates number of replicates per scenario
#' @param seed base random seed
#' @param n_subjects number of subjects passed to each simulation
#' @param n_timepoints number of timepoints passed to each simulation
#' @param n_taxa number of taxa passed to each simulation
#' @param alpha adjusted p-value cutoff
#' @param prev.filter prevalence filter passed to `adapt`
#' @param depth.filter depth filter passed to `adapt`
#'
#' @returns A data frame with one row per scenario replicate.
#' @export
run_adaptme_simulation_validation <- function(scenarios = adaptme_validation_scenarios(),
                                              n_replicates = 3, seed = 1,
                                              n_subjects = 20, n_timepoints = 2,
                                              n_taxa = 12, alpha = 0.05,
                                              prev.filter = 0,
                                              depth.filter = 0) {
  rows <- list()
  row_id <- 1
  for (scenario_name in names(scenarios)) {
    scenario <- scenarios[[scenario_name]]
    for (replicate_id in seq_len(n_replicates)) {
      sim_seed <- seed + row_id - 1
      sim <- do.call(simulate_adaptme_phyloseq, c(
        list(seed = sim_seed, n_subjects = n_subjects,
             n_timepoints = n_timepoints, n_taxa = n_taxa),
        scenario
      ))

      elapsed <- system.time({
        result <- tryCatch(
          suppressMessages(suppressWarnings(adapt(
            input_data = sim$phyobj,
            cond.var = "condition",
            adj.var = "confounder",
            subject.var = "subject",
            prev.filter = prev.filter,
            depth.filter = depth.filter,
            alpha = alpha
          ))),
          error = function(e) e
        )
      })[["elapsed"]]

      if (inherits(result, "error")) {
        metrics <- data.frame(
          n_taxa = n_taxa, n_da = scenario$n_da, n_fit = 0,
          discoveries = NA_integer_, true_positives = NA_integer_,
          false_positives = NA_integer_, fpr = NA_real_, fdr = NA_real_,
          power = NA_real_, mean_da_bias = NA_real_,
          stringsAsFactors = FALSE
        )
        error_message <- conditionMessage(result)
      } else {
        metrics <- adaptme_simulation_metrics(result, sim$truth, alpha = alpha)
        error_message <- NA_character_
      }

      rows[[row_id]] <- cbind(
        data.frame(
          scenario = scenario_name,
          replicate = replicate_id,
          seed = sim_seed,
          elapsed_sec = as.numeric(elapsed),
          error = error_message,
          stringsAsFactors = FALSE
        ),
        metrics
      )
      row_id <- row_id + 1
    }
  }
  do.call(rbind, rows)
}

#' Aggregate ADAPT-ME Validation Metrics
#'
#' @param validation_results output from `run_adaptme_simulation_validation`
#'
#' @returns A scenario-level summary data frame.
#' @export
summarize_adaptme_validation <- function(validation_results) {
  metric_cols <- c("elapsed_sec", "n_fit", "discoveries", "fpr", "fdr", "power", "mean_da_bias")
  aggregate(validation_results[, metric_cols, drop = FALSE],
            by = list(scenario = validation_results$scenario),
            FUN = function(x) mean(x, na.rm = TRUE))
}

#' ADAPT-ME Stress-Test Scenarios
#'
#' @description
#' Scenarios designed to stress-test ADAPT-ME under challenging real-world
#' conditions: extreme zero-inflation, strong confounding correlated with
#' condition, and unbalanced designs with sparse taxa.
#'
#' @returns A named list of simulation parameter lists compatible with
#'   `run_adaptme_simulation_validation`.
#' @export
adaptme_stress_scenarios <- function() {
  list(
    extreme_zeros = list(
      n_da = 2, effect_log10 = 0.7, zero_inflation = 0.75,
      unbalanced = FALSE, confounded = FALSE, condition_type = "binary"
    ),
    strong_confound = list(
      n_da = 2, effect_log10 = 0.7, zero_inflation = 0.1,
      unbalanced = FALSE, confounded = TRUE, condition_type = "binary"
    ),
    unbalanced_sparse = list(
      n_da = 2, effect_log10 = 0.7, zero_inflation = 0.4,
      unbalanced = TRUE, drop_fraction = 0.4,
      confounded = FALSE, condition_type = "binary"
    )
  )
}

#' Naive LM Baseline Metrics for Simulation Comparison
#'
#' @description
#' Fits a simple pseudocount log-proportion linear model—ignoring
#' repeated measures entirely—to simulation output from
#' `simulate_adaptme_phyloseq`. Intended for benchmarking: this naive
#' approach does not account for subject-level random effects and is
#' expected to show inflated FPR under longitudinal designs with strong
#' subject baselines, where ADAPT-ME should improve.
#'
#' @param sim output list from `simulate_adaptme_phyloseq`
#' @param alpha BH-adjusted p-value cutoff for calling discoveries
#'
#' @returns A one-row data frame with columns `method`, `n_da`,
#'   `discoveries`, `true_positives`, `false_positives`, `fpr`, `fdr`,
#'   and `power`.
#' @export
naive_lm_metrics <- function(sim, alpha = 0.05) {
  counts <- sim$counts
  meta   <- sim$metadata
  truth  <- sim$truth

  counts_pc <- counts + 0.5
  lib_size  <- rowSums(counts_pc)
  log_prop  <- log(counts_pc) - log(lib_size)

  has_confounder <- "confounder" %in% colnames(meta)

  pvals <- vapply(seq_len(ncol(log_prop)), function(j) {
    df <- data.frame(
      y         = log_prop[, j],
      condition = meta$condition,
      stringsAsFactors = FALSE
    )
    if (has_confounder) df$confounder <- meta$confounder
    formula_obj <- if (has_confounder) {
      y ~ condition + confounder
    } else {
      y ~ condition
    }
    fit <- tryCatch(stats::lm(formula_obj, data = df), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    coefs <- summary(fit)$coefficients
    if ("condition" %in% rownames(coefs)) coefs["condition", "Pr(>|t|)"] else NA_real_
  }, numeric(1))
  names(pvals) <- colnames(counts)

  adj_pvals  <- p.adjust(pvals, method = "BH")
  called     <- !is.na(adj_pvals) & adj_pvals < alpha
  is_da      <- truth$is_da[match(colnames(counts), truth$Taxa)]
  null_taxa  <- !is_da
  tp         <- sum(called & is_da,    na.rm = TRUE)
  fp         <- sum(called & null_taxa, na.rm = TRUE)
  discoveries <- sum(called, na.rm = TRUE)

  data.frame(
    method          = "naive_lm",
    n_da            = sum(is_da),
    discoveries     = discoveries,
    true_positives  = tp,
    false_positives = fp,
    fpr  = if (any(null_taxa)) fp / sum(null_taxa) else NA_real_,
    fdr  = if (discoveries > 0) fp / discoveries else 0,
    power = if (any(is_da)) tp / sum(is_da) else NA_real_,
    stringsAsFactors = FALSE
  )
}

#' Compare ADAPT-ME Against Naive LM Baseline
#'
#' @description
#' Runs both ADAPT-ME (mixed-effects Tobit with subject random intercepts)
#' and a naive pseudocount log-proportion linear model on replicated
#' simulated datasets. This demonstrates that ADAPT-ME provides better
#' false-positive-rate control under longitudinal designs where the naive
#' model ignores within-subject correlation.
#'
#' @param n_sims number of simulation replicates
#' @param seed base random seed
#' @param n_subjects number of subjects
#' @param n_timepoints number of timepoints per subject
#' @param n_taxa total number of taxa
#' @param n_da number of truly differential taxa
#' @param effect_log10 condition effect for DA taxa on the log10 scale
#' @param zero_inflation additional zero-inflation probability
#' @param alpha BH-adjusted p-value cutoff
#'
#' @returns A data frame with two rows per replicate (one per method) and
#'   columns `replicate`, `method`, `n_da`, `discoveries`,
#'   `true_positives`, `false_positives`, `fpr`, `fdr`, and `power`.
#' @export
compare_adaptme_vs_naive <- function(n_sims = 5, seed = 42,
                                     n_subjects = 20, n_timepoints = 3,
                                     n_taxa = 15, n_da = 3,
                                     effect_log10 = 0.5, zero_inflation = 0.1,
                                     alpha = 0.05) {
  rows <- list()
  for (i in seq_len(n_sims)) {
    sim <- simulate_adaptme_phyloseq(
      seed        = seed + i - 1,
      n_subjects  = n_subjects,
      n_timepoints = n_timepoints,
      n_taxa      = n_taxa,
      n_da        = n_da,
      effect_log10 = effect_log10,
      zero_inflation = zero_inflation
    )

    me_result <- tryCatch(
      suppressMessages(suppressWarnings(adapt(
        input_data   = sim$phyobj,
        cond.var     = "condition",
        adj.var      = "confounder",
        subject.var  = "subject",
        prev.filter  = 0,
        depth.filter = 0,
        alpha        = alpha
      ))),
      error = function(e) e
    )
    me_row <- if (inherits(me_result, "error")) {
      data.frame(
        replicate = i, method = "adaptme", n_da = n_da,
        discoveries = NA_integer_, true_positives = NA_integer_,
        false_positives = NA_integer_, fpr = NA_real_, fdr = NA_real_,
        power = NA_real_, stringsAsFactors = FALSE
      )
    } else {
      m <- adaptme_simulation_metrics(me_result, sim$truth, alpha)
      data.frame(
        replicate = i, method = "adaptme", n_da = m$n_da,
        discoveries = m$discoveries, true_positives = m$true_positives,
        false_positives = m$false_positives, fpr = m$fpr, fdr = m$fdr,
        power = m$power, stringsAsFactors = FALSE
      )
    }

    naive_row <- tryCatch(
      {
        nm <- naive_lm_metrics(sim, alpha)
        data.frame(
          replicate = i, method = nm$method, n_da = nm$n_da,
          discoveries = nm$discoveries, true_positives = nm$true_positives,
          false_positives = nm$false_positives, fpr = nm$fpr, fdr = nm$fdr,
          power = nm$power, stringsAsFactors = FALSE
        )
      },
      error = function(e) {
        data.frame(
          replicate = i, method = "naive_lm", n_da = n_da,
          discoveries = NA_integer_, true_positives = NA_integer_,
          false_positives = NA_integer_, fpr = NA_real_, fdr = NA_real_,
          power = NA_real_, stringsAsFactors = FALSE
        )
      }
    )

    rows[[2 * i - 1]] <- me_row
    rows[[2 * i]]     <- naive_row
  }
  do.call(rbind, rows)
}

# ── Internal helper: extract standard metrics from a called/truth pair ────────

.benchmark_metrics <- function(method, called, taxa_names, truth, n_da_scenario) {
  is_da     <- truth$is_da[match(taxa_names, truth$Taxa)]
  null_taxa <- !is_da
  tp <- sum(called & is_da,     na.rm = TRUE)
  fp <- sum(called & null_taxa, na.rm = TRUE)
  discoveries <- sum(called, na.rm = TRUE)
  data.frame(
    method          = method,
    n_da            = sum(is_da, na.rm = TRUE),
    discoveries     = discoveries,
    true_positives  = tp,
    false_positives = fp,
    fpr   = if (any(null_taxa)) fp / sum(null_taxa) else NA_real_,
    fdr   = if (discoveries > 0) fp / discoveries else 0,
    power = if (any(is_da)) tp / sum(is_da) else NA_real_,
    stringsAsFactors = FALSE
  )
}

#' MaAsLin2 Metrics for Simulation Benchmarking
#'
#' @description
#' Runs MaAsLin2 (linear mixed model on TSS+LOG transformed counts) on
#' a simulated dataset and returns the same metric columns as
#' `naive_lm_metrics`. Requires the Maaslin2 package.
#'
#' @param sim output list from `simulate_adaptme_phyloseq`
#' @param alpha BH adjusted p-value cutoff
#' @returns A one-row data frame with method, FPR, FDR, power, etc.
#' @export
maaslin2_metrics <- function(sim, alpha = 0.05) {
  if (!requireNamespace("Maaslin2", quietly = TRUE))
    stop("Maaslin2 is not installed. Install with: BiocManager::install('Maaslin2')")

  counts <- sim$counts
  meta   <- sim$metadata
  truth  <- sim$truth
  taxa_names <- colnames(counts)

  meta_df <- as.data.frame(meta)
  meta_df$subject <- as.character(meta_df$subject)

  adj_cols <- intersect(c("condition", "confounder"), colnames(meta_df))
  rand_cols <- "subject"

  output_dir <- tempfile(pattern = "maaslin2_")
  dir.create(output_dir, showWarnings = FALSE)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

  fit <- tryCatch(
    suppressMessages(suppressWarnings(
      Maaslin2::Maaslin2(
        input_data      = as.data.frame(counts),
        input_metadata  = meta_df[, c(adj_cols, rand_cols), drop = FALSE],
        output          = output_dir,
        fixed_effects   = adj_cols,
        random_effects  = rand_cols,
        normalization   = "TSS",
        transform       = "LOG",
        analysis_method = "LM",
        min_abundance   = 0,
        min_prevalence  = 0,
        max_significance = 1,
        plot_heatmap    = FALSE,
        plot_scatter    = FALSE
      )
    )),
    error = function(e) NULL
  )

  if (is.null(fit))
    return(.benchmark_metrics("maaslin2", rep(FALSE, length(taxa_names)),
                               taxa_names, truth, sum(truth$is_da)))

  res <- fit$results
  res <- res[!is.na(res$metadata) & res$metadata == "condition", , drop = FALSE]

  called <- setNames(rep(FALSE, length(taxa_names)), taxa_names)
  if (nrow(res) > 0) {
    sig_feats <- res$feature[!is.na(res$qval) & res$qval < alpha]
    # MaAsLin2 may sanitise feature names — match loosely
    clean <- function(x) gsub("[^A-Za-z0-9]", ".", x)
    for (feat in sig_feats) {
      idx <- which(clean(taxa_names) == clean(feat))
      if (length(idx) > 0) called[idx] <- TRUE
    }
  }

  .benchmark_metrics("maaslin2", called, taxa_names, truth, sum(truth$is_da))
}

#' ANCOM-BC2 Metrics for Simulation Benchmarking
#'
#' @description
#' Runs ANCOM-BC2 with a subject random intercept on a simulated dataset
#' and returns the same metric columns as `naive_lm_metrics`.
#' Requires the ANCOMBC package (>= 2.0.0).
#'
#' @param sim output list from `simulate_adaptme_phyloseq`
#' @param alpha significance threshold passed to `ancombc2`
#' @returns A one-row data frame with method, FPR, FDR, power, etc.
#' @export
ancombc2_metrics <- function(sim, alpha = 0.05) {
  if (!requireNamespace("ANCOMBC", quietly = TRUE))
    stop("ANCOMBC is not installed. Install with: BiocManager::install('ANCOMBC')")

  counts <- sim$counts
  meta   <- sim$metadata
  truth  <- sim$truth
  taxa_names <- colnames(counts)

  meta_df <- as.data.frame(meta)
  meta_df$subject <- as.character(meta_df$subject)
  rownames(meta_df) <- rownames(counts)

  phyobj <- phyloseq::phyloseq(
    phyloseq::otu_table(counts, taxa_are_rows = FALSE),
    phyloseq::sample_data(meta_df)
  )

  fix_vars <- intersect(c("condition", "confounder"), colnames(meta_df))
  fix_formula <- paste(fix_vars, collapse = " + ")

  fit <- tryCatch(
    suppressMessages(suppressWarnings(
      ANCOMBC::ancombc2(
        data         = phyobj,
        fix_formula  = fix_formula,
        rand_formula = "(1|subject)",
        p_adj_method = "BH",
        prv_cut      = 0,
        lib_cut      = 0,
        alpha        = alpha,
        verbose      = FALSE
      )
    )),
    error = function(e) NULL
  )

  if (is.null(fit))
    return(.benchmark_metrics("ancombc2", rep(FALSE, length(taxa_names)),
                               taxa_names, truth, sum(truth$is_da)))

  res <- fit$res
  called <- setNames(rep(FALSE, length(taxa_names)), taxa_names)

  # Prefer diff_condition column; fall back to q_condition < alpha
  if ("diff_condition" %in% colnames(res) && "taxon" %in% colnames(res)) {
    sig_taxa <- res$taxon[!is.na(res$diff_condition) & res$diff_condition == TRUE]
    called[intersect(sig_taxa, taxa_names)] <- TRUE
  } else if ("q_condition" %in% colnames(res) && "taxon" %in% colnames(res)) {
    sig_taxa <- res$taxon[!is.na(res$q_condition) & res$q_condition < alpha]
    called[intersect(sig_taxa, taxa_names)] <- TRUE
  }

  .benchmark_metrics("ancombc2", called, taxa_names, truth, sum(truth$is_da))
}

#' ZINQ-L Metrics for Simulation Benchmarking
#'
#' @description
#' Runs ZINQ-L (zero-inflated quantile approach for longitudinal data) on a
#' simulated dataset one taxon at a time and returns the same metric columns as
#' `naive_lm_metrics`. Raw p-values are BH-corrected across all taxa.
#' Requires the ZINQL package (install with
#' `devtools::install_github("AlbertSL98/ZINQ-L")`).
#'
#' @param sim output list from `simulate_adaptme_phyloseq`
#' @param alpha BH adjusted p-value cutoff
#' @param taus quantile levels to test (default: c(0.1, 0.25, 0.5, 0.75, 0.9))
#' @returns A one-row data frame with method, FPR, FDR, power, etc.
#' @export
zinql_metrics <- function(sim, alpha = 0.05,
                          taus = c(0.1, 0.25, 0.5, 0.75, 0.9)) {
  if (!requireNamespace("ZINQL", quietly = TRUE))
    stop("ZINQL is not installed. Install with: devtools::install_github('AlbertSL98/ZINQ-L')")

  counts <- sim$counts
  meta   <- sim$metadata
  truth  <- sim$truth
  taxa_names <- colnames(counts)

  meta_df <- as.data.frame(meta)
  meta_df$subject <- as.character(meta_df$subject)

  # Build formula: condition (+ confounder if present) + subject random intercept
  fix_vars <- intersect(c("condition", "confounder"), colnames(meta_df))
  formula_str <- paste("y ~", paste(fix_vars, collapse = " + "), "+ (1|subject)")
  fmla <- as.formula(formula_str)

  # One ZINQL_fit call per taxon; collect raw p-values then BH-correct
  raw_pvals <- vapply(seq_len(ncol(counts)), function(j) {
    y_j <- counts[, j]
    res <- tryCatch(
      suppressMessages(suppressWarnings(
        ZINQL::ZINQL_fit(
          y       = y_j,
          meta    = meta_df,
          formula = fmla,
          C       = "condition",
          taus    = taus,
          method  = "MinP",
          seed    = 1L
        )
      )),
      error = function(e) NULL
    )
    if (is.null(res) || is.null(res$Final_P_value)) return(NA_real_)
    pv <- res$Final_P_value
    # method='MinP' returns a named scalar: ZINQL_MinP
    pv_val <- if (!is.null(names(pv)) && "ZINQL_MinP" %in% names(pv)) {
      pv[["ZINQL_MinP"]]
    } else {
      pv[[1]]
    }
    if (is.numeric(pv_val) && length(pv_val) == 1) pv_val else NA_real_
  }, numeric(1))

  # BH correction across all taxa; treat NA as non-significant
  qvals  <- p.adjust(raw_pvals, method = "BH")
  called <- !is.na(qvals) & qvals < alpha

  .benchmark_metrics("zinql", called, taxa_names, truth, sum(truth$is_da))
}

#' Run Multi-Method Benchmark Comparison
#'
#' @description
#' Runs ADAPT-ME, MaAsLin2, ANCOM-BC2, and naive LM on identical simulated
#' datasets across multiple scenarios and replicates. Returns per-replicate
#' FPR, FDR, and power for each method.
#'
#' @param scenarios named list of scenario parameter lists (defaults to null +
#'   strong_signal + sparse)
#' @param n_replicates replicates per scenario
#' @param seed base random seed
#' @param n_subjects subjects per simulation
#' @param n_timepoints timepoints per subject
#' @param n_taxa taxa per simulation
#' @param alpha BH threshold
#' @param methods character vector of methods to run; subset of
#'   `c("adaptme","maaslin2","ancombc2","naive_lm")`
#'
#' @returns A data frame with one row per scenario × replicate × method.
#' @export
run_benchmark_comparison <- function(
    scenarios    = NULL,
    n_replicates = 10,
    seed         = 42,
    n_subjects   = 20,
    n_timepoints = 2,
    n_taxa       = 40,
    alpha        = 0.05,
    methods      = c("adaptme", "maaslin2", "ancombc2", "naive_lm", "zinql")) {

  if (is.null(scenarios)) {
    scenarios <- list(
      null = list(n_da = 0, effect_log10 = 0, zero_inflation = 0.05,
                  unbalanced = FALSE, confounded = FALSE,
                  condition_type = "binary"),
      strong_signal = list(n_da = 3, effect_log10 = 0.7, zero_inflation = 0.05,
                           unbalanced = FALSE, confounded = FALSE,
                           condition_type = "binary"),
      sparse = list(n_da = 3, effect_log10 = 0.7, zero_inflation = 0.35,
                    unbalanced = FALSE, confounded = FALSE,
                    condition_type = "binary")
    )
  }

  rows <- list()
  row_id <- 1

  for (scenario_name in names(scenarios)) {
    scenario <- scenarios[[scenario_name]]
    cat(sprintf("\nScenario: %s\n", scenario_name))

    for (rep_id in seq_len(n_replicates)) {
      sim_seed <- seed + (row_id - 1)
      sim <- tryCatch(
        do.call(simulate_adaptme_phyloseq, c(
          list(seed = sim_seed, n_subjects = n_subjects,
               n_timepoints = n_timepoints, n_taxa = n_taxa),
          scenario
        )),
        error = function(e) NULL
      )
      if (is.null(sim)) { row_id <- row_id + 1; next }

      for (method in methods) {
        elapsed <- system.time({
          result <- tryCatch({
            if (method == "adaptme") {
              res <- suppressMessages(suppressWarnings(adapt(
                input_data   = sim$phyobj,
                cond.var     = "condition",
                adj.var      = "confounder",
                subject.var  = "subject",
                prev.filter  = 0,
                depth.filter = 0,
                alpha        = alpha
              )))
              m <- adaptme_simulation_metrics(res, sim$truth, alpha)
              data.frame(method = "adaptme", n_da = m$n_da,
                         discoveries = m$discoveries,
                         true_positives = m$true_positives,
                         false_positives = m$false_positives,
                         fpr = m$fpr, fdr = m$fdr, power = m$power,
                         error = NA_character_, stringsAsFactors = FALSE)
            } else if (method == "maaslin2") {
              m <- maaslin2_metrics(sim, alpha)
              cbind(m, error = NA_character_, stringsAsFactors = FALSE)
            } else if (method == "ancombc2") {
              m <- ancombc2_metrics(sim, alpha)
              cbind(m, error = NA_character_, stringsAsFactors = FALSE)
            } else if (method == "naive_lm") {
              m <- naive_lm_metrics(sim, alpha)
              cbind(m, error = NA_character_, stringsAsFactors = FALSE)
            } else if (method == "zinql") {
              m <- zinql_metrics(sim, alpha)
              cbind(m, error = NA_character_, stringsAsFactors = FALSE)
            }
          }, error = function(e) {
            data.frame(method = method, n_da = scenario$n_da,
                       discoveries = NA_integer_, true_positives = NA_integer_,
                       false_positives = NA_integer_, fpr = NA_real_,
                       fdr = NA_real_, power = NA_real_,
                       error = conditionMessage(e), stringsAsFactors = FALSE)
          })
          result
        })[["elapsed"]]

        rows[[row_id * length(methods) + match(method, methods)]] <- cbind(
          data.frame(scenario = scenario_name, replicate = rep_id,
                     seed = sim_seed, elapsed_sec = as.numeric(elapsed),
                     stringsAsFactors = FALSE),
          result
        )
      }
      cat(sprintf("  rep %d/%d done\n", rep_id, n_replicates))
      row_id <- row_id + 1
    }
  }

  do.call(rbind, Filter(Negate(is.null), rows))
}
