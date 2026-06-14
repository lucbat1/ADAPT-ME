#!/usr/bin/env Rscript
# run-benchmark-comparison.R
#
# Benchmarks ADAPT-ME against MaAsLin2, ANCOM-BC2, and a naive LM baseline
# across three simulation scenarios (null, strong signal, sparse data).
#
# Requires:
#   BiocManager::install("Maaslin2")
#   BiocManager::install("ANCOMBC")   # version >= 2.0.0 for rand_formula
#
# Run from the package root:
#   Rscript inst/scripts/run-benchmark-comparison.R
#
# Scale up for manuscript (default is 10 reps x 3 scenarios, ~45 min):
#   ADAPTME_BENCH_REPS=25 ADAPTME_BENCH_TAXA=40 Rscript inst/scripts/run-benchmark-comparison.R
#
# Output goes to: adaptme-benchmark-results/

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})

# ── Config ────────────────────────────────────────────────────────────────────
read_env_int <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (identical(v, "")) default else as.integer(v)
}

n_reps       <- read_env_int("ADAPTME_BENCH_REPS",     10)
n_subjects   <- read_env_int("ADAPTME_BENCH_SUBJECTS",  20)
n_timepoints <- read_env_int("ADAPTME_BENCH_TIMEPOINTS", 2)
n_taxa       <- read_env_int("ADAPTME_BENCH_TAXA",       40)
output_dir   <- Sys.getenv("ADAPTME_BENCH_DIR", unset = "adaptme-benchmark-results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Sample sizes to benchmark. The primary run uses n_subjects (default 20).
# An additional run at n_large (default 50) is performed for the three key
# scenarios to show how power scales with sample size.
n_large <- read_env_int("ADAPTME_BENCH_SUBJECTS_LARGE", 50)

cat("ADAPT-ME Benchmark Comparison\n")
cat(sprintf("  replicates: %d  |  subjects: %d (+ %d large)  |  timepoints: %d  |  taxa: %d\n",
            n_reps, n_subjects, n_large, n_timepoints, n_taxa))

# ── Check required packages ───────────────────────────────────────────────────
missing_bioc <- c()
if (!requireNamespace("Maaslin2", quietly = TRUE))
  missing_bioc <- c(missing_bioc, "Maaslin2")
if (!requireNamespace("ANCOMBC", quietly = TRUE))
  missing_bioc <- c(missing_bioc, "ANCOMBC")

if (length(missing_bioc) > 0) {
  cat("\nInstalling missing Bioconductor packages:", paste(missing_bioc, collapse = ", "), "\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

# ZINQ-L (GitHub package — install separately if wanted)
if (!requireNamespace("ZINQL", quietly = TRUE)) {
  cat("\nZINQL not found. To include ZINQ-L in the benchmark run:\n")
  cat("  devtools::install_github('AlbertSL98/ZINQ-L')\n")
  cat("Continuing without ZINQ-L.\n\n")
}

methods_to_run <- c("adaptme", "naive_lm")
if (requireNamespace("Maaslin2", quietly = TRUE)) methods_to_run <- c(methods_to_run, "maaslin2")
if (requireNamespace("ANCOMBC",  quietly = TRUE)) methods_to_run <- c(methods_to_run, "ancombc2")
if (requireNamespace("ZINQL",    quietly = TRUE)) methods_to_run <- c(methods_to_run, "zinql")
cat(sprintf("  methods: %s\n\n", paste(methods_to_run, collapse = ", ")))

# ── Scenarios ─────────────────────────────────────────────────────────────────
# Effect sizes are chosen to reflect realistic microbiome studies:
#   moderate = 0.3 log10 (~2-fold) — typical in observational cohorts
#   strong   = 0.7 log10 (~5-fold) — large, used to verify power ceiling
# Running both makes the comparison more discriminating: methods that hit
# power=1.000 on strong_signal may diverge substantially on moderate_signal.
scenarios <- list(
  null = list(
    n_da = 0, effect_log10 = 0, zero_inflation = 0.05,
    unbalanced = FALSE, confounded = FALSE, condition_type = "binary"
  ),
  moderate_signal = list(
    n_da = 3, effect_log10 = 0.3, zero_inflation = 0.05,
    unbalanced = FALSE, confounded = FALSE, condition_type = "binary"
  ),
  strong_signal = list(
    n_da = 3, effect_log10 = 0.7, zero_inflation = 0.05,
    unbalanced = FALSE, confounded = FALSE, condition_type = "binary"
  ),
  moderate_sparse = list(
    n_da = 3, effect_log10 = 0.3, zero_inflation = 0.35,
    unbalanced = FALSE, confounded = FALSE, condition_type = "binary"
  ),
  sparse = list(
    n_da = 3, effect_log10 = 0.7, zero_inflation = 0.35,
    unbalanced = FALSE, confounded = FALSE, condition_type = "binary"
  ),
  unbalanced = list(
    n_da = 3, effect_log10 = 0.3, zero_inflation = 0.05,
    unbalanced = TRUE, drop_fraction = 0.25,
    confounded = FALSE, condition_type = "binary"
  ),
  confounded = list(
    n_da = 3, effect_log10 = 0.3, zero_inflation = 0.05,
    unbalanced = FALSE, confounded = TRUE, condition_type = "binary"
  )
)

# Rough runtime estimate
# ZINQ-L is taxon-by-taxon (quantile mixed models), ~2–5 s/taxon
n_adapt_taxa_passes  <- n_taxa * 2        # reference selection + final
est_sec_adaptme      <- n_adapt_taxa_passes * 0.08
est_sec_zinql        <- if ("zinql" %in% methods_to_run) n_taxa * 3 else 0
est_total_min <- round(
  n_reps * length(scenarios) * (est_sec_adaptme + est_sec_zinql + 15) / 60
)
cat(sprintf("Estimated runtime: ~%d minutes\n\n", est_total_min))

# ── Run benchmark — primary (n = n_subjects) ──────────────────────────────────
t_start <- proc.time()[["elapsed"]]
cat(sprintf("Pass 1: all scenarios, n_subjects = %d\n", n_subjects))

results_primary <- run_benchmark_comparison(
  scenarios    = scenarios,
  n_replicates = n_reps,
  seed         = 42,
  n_subjects   = n_subjects,
  n_timepoints = n_timepoints,
  n_taxa       = n_taxa,
  alpha        = 0.05,
  methods      = methods_to_run
)
results_primary$n_subjects <- n_subjects

# ── Run benchmark — sample-size comparison (n = n_large) ──────────────────────
# Only the three most informative scenarios are re-run at the larger sample
# size: moderate_signal (power bottleneck), strong_signal (anchor), and
# sparse (zero-inflation stress test).
key_scenarios <- scenarios[intersect(
  c("moderate_signal", "strong_signal", "sparse"),
  names(scenarios)
)]

cat(sprintf("\nPass 2: key scenarios, n_subjects = %d\n", n_large))

results_large <- run_benchmark_comparison(
  scenarios    = key_scenarios,
  n_replicates = n_reps,
  seed         = 9999,          # different seed to avoid overlap with pass 1
  n_subjects   = n_large,
  n_timepoints = n_timepoints,
  n_taxa       = n_taxa,
  alpha        = 0.05,
  methods      = methods_to_run
)
results_large$n_subjects <- n_large

results <- rbind(results_primary, results_large)

elapsed_total <- proc.time()[["elapsed"]] - t_start
cat(sprintf("\nBoth passes completed in %.1f minutes\n", elapsed_total / 60))

# ── Save replicate-level results ──────────────────────────────────────────────
write.csv(results,
          file.path(output_dir, "benchmark-replicate-results.csv"),
          row.names = FALSE)

# ── Aggregate summary ─────────────────────────────────────────────────────────
summary_tbl <- aggregate(
  cbind(fpr, fdr, power) ~ method + scenario + n_subjects,
  data = results[is.na(results$error), ],
  FUN  = function(x) round(mean(x, na.rm = TRUE), 4)
)
summary_tbl <- summary_tbl[order(summary_tbl$n_subjects,
                                  summary_tbl$scenario,
                                  summary_tbl$method), ]

write.csv(summary_tbl,
          file.path(output_dir, "benchmark-summary.csv"),
          row.names = FALSE)

# ── Print results ─────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════\n")
cat("BENCHMARK SUMMARY (mean across replicates)\n")
cat("══════════════════════════════════════════════════════════════════\n\n")

for (ns in sort(unique(summary_tbl$n_subjects))) {
  cat(sprintf("── n_subjects = %d ──────────────────────────────────────\n", ns))
  for (sc in unique(summary_tbl$scenario[summary_tbl$n_subjects == ns])) {
    cat(sprintf("Scenario: %s\n", sc))
    sub <- summary_tbl[summary_tbl$scenario == sc &
                       summary_tbl$n_subjects == ns,
                       c("method","fpr","fdr","power")]
    print(sub, row.names = FALSE)
    cat("\n")
  }
}

# ── Key checks ────────────────────────────────────────────────────────────────
cat("── Key checks ──────────────────────────────────────────────────────\n")

for (sc in c("null", "moderate_signal", "strong_signal")) {
  if (!sc %in% summary_tbl$scenario) next
  # Use n_subjects (primary run) for key checks
  sub <- summary_tbl[summary_tbl$scenario == sc &
                     summary_tbl$n_subjects == n_subjects, ]
  adaptme_row  <- sub[sub$method == "adaptme",  ]
  maaslin_row  <- sub[sub$method == "maaslin2", ]
  ancombc_row  <- sub[sub$method == "ancombc2", ]
  zinql_row    <- sub[sub$method == "zinql",    ]

  if (sc == "null" && nrow(adaptme_row) > 0) {
    cat(sprintf("[null FPR]  ADAPT-ME  = %.4f  (target <= 0.05): %s\n",
                adaptme_row$fpr,
                if (adaptme_row$fpr <= 0.05) "PASS" else "REVIEW"))
    if (nrow(maaslin_row) > 0)
      cat(sprintf("            MaAsLin2  = %.4f\n", maaslin_row$fpr))
    if (nrow(ancombc_row) > 0)
      cat(sprintf("            ANCOM-BC2 = %.4f\n", ancombc_row$fpr))
    if (nrow(zinql_row) > 0)
      cat(sprintf("            ZINQ-L    = %.4f\n", zinql_row$fpr))
  }

  if (sc %in% c("moderate_signal", "strong_signal") && nrow(adaptme_row) > 0) {
    cat(sprintf("[power]     ADAPT-ME  = %.4f  (target > 0.5): %s\n",
                adaptme_row$power,
                if (adaptme_row$power > 0.5) "PASS" else "REVIEW"))
    if (nrow(maaslin_row) > 0)
      cat(sprintf("            MaAsLin2  = %.4f\n", maaslin_row$power))
    if (nrow(ancombc_row) > 0)
      cat(sprintf("            ANCOM-BC2 = %.4f\n", ancombc_row$power))
    if (nrow(zinql_row) > 0)
      cat(sprintf("            ZINQ-L    = %.4f\n", zinql_row$power))
    cat(sprintf("[FDR]       ADAPT-ME  = %.4f  (target <= 0.05): %s\n",
                adaptme_row$fdr,
                if (adaptme_row$fdr <= 0.05) "PASS" else "REVIEW"))
    if (nrow(maaslin_row) > 0)
      cat(sprintf("            MaAsLin2  = %.4f\n", maaslin_row$fdr))
    if (nrow(ancombc_row) > 0)
      cat(sprintf("            ANCOM-BC2 = %.4f\n", ancombc_row$fdr))
    if (nrow(zinql_row) > 0)
      cat(sprintf("            ZINQ-L    = %.4f\n", zinql_row$fdr))
  }
}

# ── Error summary ─────────────────────────────────────────────────────────────
errors <- results[!is.na(results$error), ]
if (nrow(errors) > 0) {
  cat(sprintf("\n%d replicate(s) had errors:\n", nrow(errors)))
  print(errors[, c("scenario","replicate","method","error")])
} else {
  cat("\nAll replicates completed without error.\n")
}

cat(sprintf("\nResults written to: %s/\n", output_dir))
cat("Files:\n")
for (f in list.files(output_dir, full.names = FALSE))
  cat(sprintf("  %s\n", f))
