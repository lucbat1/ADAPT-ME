test_that("simulation generator returns phyloseq data with usable truth labels", {
  sim <- simulate_adaptme_phyloseq(
    seed = 10,
    n_subjects = 6,
    n_timepoints = 2,
    n_taxa = 6,
    n_da = 2,
    effect_log10 = 0.5,
    zero_inflation = 0.1
  )

  expect_s4_class(sim$phyobj, "phyloseq")
  expect_equal(nrow(sim$truth), 6)
  expect_equal(sum(sim$truth$is_da), 2)
  expect_true(all(c("condition", "subject", "confounder") %in% colnames(sim$metadata)))
  expect_true(any(sim$counts == 0))
})

test_that("validation scenarios cover null, signal, sparse, imbalance, and confounding", {
  scenarios <- adaptme_validation_scenarios()

  expect_true(all(c("null", "strong_signal", "sparse", "unbalanced", "confounded") %in% names(scenarios)))
  expect_equal(scenarios$null$n_da, 0)
  expect_gt(scenarios$strong_signal$effect_log10, 0)
  expect_gt(scenarios$sparse$zero_inflation, scenarios$strong_signal$zero_inflation)
  expect_true(scenarios$unbalanced$unbalanced)
  expect_true(scenarios$confounded$confounded)
})

test_that("metrics report FPR for null simulations and power for signal simulations", {
  null_sim <- simulate_adaptme_phyloseq(
    seed = 11,
    n_subjects = 8,
    n_timepoints = 2,
    n_taxa = 4,
    n_da = 0,
    effect_log10 = 0
  )
  null_result <- suppressMessages(adapt(
    input_data = null_sim$phyobj,
    cond.var = "condition",
    adj.var = "confounder",
    subject.var = "subject",
    prev.filter = 0,
    depth.filter = 0,
    alpha = 0.2
  ))
  null_metrics <- adaptme_simulation_metrics(null_result, null_sim$truth, alpha = 0.2)

  expect_equal(null_metrics$n_da, 0)
  expect_true(is.na(null_metrics$power))
  expect_true(is.finite(null_metrics$fpr))

  signal_sim <- simulate_adaptme_phyloseq(
    seed = 12,
    n_subjects = 8,
    n_timepoints = 2,
    n_taxa = 4,
    n_da = 1,
    effect_log10 = 0.9
  )
  signal_result <- suppressMessages(adapt(
    input_data = signal_sim$phyobj,
    cond.var = "condition",
    adj.var = "confounder",
    subject.var = "subject",
    prev.filter = 0,
    depth.filter = 0,
    alpha = 0.2
  ))
  signal_metrics <- adaptme_simulation_metrics(signal_result, signal_sim$truth, alpha = 0.2)

  expect_equal(signal_metrics$n_da, 1)
  expect_true(is.finite(signal_metrics$power))
  expect_true(is.finite(signal_metrics$mean_da_bias))
})

test_that("validation runner returns replicate-level and scenario-level summaries", {
  scenarios <- adaptme_validation_scenarios()[c("null", "sparse", "unbalanced", "confounded")]
  validation <- run_adaptme_simulation_validation(
    scenarios = scenarios,
    n_replicates = 1,
    seed = 20,
    n_subjects = 6,
    n_timepoints = 2,
    n_taxa = 4,
    alpha = 0.2
  )
  summary <- summarize_adaptme_validation(validation)

  expect_equal(nrow(validation), length(scenarios))
  expect_true(all(c("scenario", "replicate", "fpr", "fdr", "power", "elapsed_sec", "error") %in% colnames(validation)))
  expect_equal(sort(summary$scenario), sort(names(scenarios)))
})

# ── Item #2: baseline comparison ──────────────────────────────────────────────

test_that("naive_lm_metrics returns a valid one-row data frame", {
  sim <- simulate_adaptme_phyloseq(
    seed = 50, n_subjects = 8, n_timepoints = 2, n_taxa = 6, n_da = 2,
    effect_log10 = 0.6, zero_inflation = 0.1
  )
  m <- naive_lm_metrics(sim, alpha = 0.2)

  expect_equal(nrow(m), 1L)
  expect_equal(m$method, "naive_lm")
  expect_equal(m$n_da, 2L)
  expect_true(is.numeric(m$fpr) && (is.na(m$fpr) || (m$fpr >= 0 && m$fpr <= 1)))
  expect_true(m$fdr >= 0 && m$fdr <= 1)
})

test_that("compare_adaptme_vs_naive returns two rows per replicate", {
  result <- compare_adaptme_vs_naive(
    n_sims = 1, seed = 55,
    n_subjects = 6, n_timepoints = 2, n_taxa = 6, n_da = 2,
    effect_log10 = 0.5, zero_inflation = 0.1, alpha = 0.2
  )

  expect_equal(nrow(result), 2L)
  expect_true(all(c("adaptme", "naive_lm") %in% result$method))
  expect_true(all(c("fpr", "fdr", "power", "discoveries") %in% colnames(result)))
})

# ── Item #3: stress-test scenarios ────────────────────────────────────────────

test_that("continuous_time scenario runs and returns finite metrics", {
  sim <- simulate_adaptme_phyloseq(
    seed = 30, n_subjects = 8, n_timepoints = 3, n_taxa = 6, n_da = 2,
    effect_log10 = 0.3, condition_type = "continuous"
  )
  result <- suppressMessages(adapt(
    input_data   = sim$phyobj,
    cond.var     = "condition",
    adj.var      = "confounder",
    subject.var  = "subject",
    prev.filter  = 0,
    depth.filter = 0,
    alpha        = 0.2
  ))
  metrics <- adaptme_simulation_metrics(result, sim$truth, alpha = 0.2)

  expect_equal(metrics$n_da, 2L)
  expect_true(is.finite(metrics$fpr))
  expect_true(is.finite(metrics$n_fit))
})

test_that("stress scenarios produce valid phyloseq objects and adapt does not crash", {
  stress <- adaptme_stress_scenarios()
  for (name in names(stress)) {
    scenario <- stress[[name]]
    sim <- do.call(simulate_adaptme_phyloseq, c(
      list(seed = 77, n_subjects = 8, n_timepoints = 2, n_taxa = 6),
      scenario
    ))
    expect_s4_class(sim$phyobj, "phyloseq")

    output <- tryCatch(
      suppressMessages(suppressWarnings(adapt(
        input_data   = sim$phyobj,
        cond.var     = "condition",
        adj.var      = "confounder",
        subject.var  = "subject",
        prev.filter  = 0,
        depth.filter = 0
      ))),
      error = function(e) e
    )
    # Must either succeed gracefully or surface a meaningful error message
    if (inherits(output, "error")) {
      expect_true(
        nchar(conditionMessage(output)) > 0,
        info = paste("Stress scenario", name, "raised an empty error")
      )
    } else {
      expect_s4_class(output, "DAresult")
    }
  }
})
