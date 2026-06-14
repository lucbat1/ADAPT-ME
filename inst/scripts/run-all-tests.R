#!/usr/bin/env Rscript
# run-all-tests.R
#
# Full ADAPT-ME test pipeline. Run from the package root:
#
#   Rscript inst/scripts/run-all-tests.R
#
# Or with a real dataset (steps 1-4 are simulation-based; step 5 is optional):
#
#   ADAPTME_REAL_DATA=path/to/mydata.rds Rscript inst/scripts/run-all-tests.R
#
# ── Environment variables ──────────────────────────────────────────────────────
#   ADAPTME_VALIDATION_DIR   output folder               [default: adaptme-validation-results]
#   ADAPTME_REPLICATES       replicates per scenario     [default: 25]
#   ADAPTME_SUBJECTS         subjects per simulation     [default: 40]
#   ADAPTME_TIMEPOINTS       timepoints per subject      [default: 3]
#   ADAPTME_TAXA             taxa per simulation         [default: 40]
#   ADAPTME_ALPHA            FDR threshold               [default: 0.05]
#   ADAPTME_REAL_DATA        path to a .rds phyloseq     [optional]
#   ADAPTME_COND_VAR         condition variable name     [required with real data]
#   ADAPTME_SUBJECT_VAR      subject variable name       [required with real data]
#   ADAPTME_ADJ_VAR          comma-separated adj vars    [optional, e.g. "Age,BMI"]
#   ADAPTME_BASE_COND        baseline category           [optional]

suppressPackageStartupMessages({
  if (!requireNamespace("ADAPTME", quietly = TRUE)) {
    message("ADAPTME not installed — loading from source via devtools::load_all()")
    devtools::load_all(".")
  } else {
    library(ADAPTME)
  }
})

# ── helpers ────────────────────────────────────────────────────────────────────
read_env_int <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (identical(v, "")) default else as.integer(v)
}
read_env_dbl <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (identical(v, "")) default else as.numeric(v)
}
read_env_str <- function(name, default = "") {
  v <- Sys.getenv(name, unset = "")
  if (identical(v, "")) default else v
}
section <- function(title) {
  bar <- strrep("─", 60)
  cat(sprintf("\n%s\n  %s\n%s\n", bar, title, bar))
}

# ── config ─────────────────────────────────────────────────────────────────────
output_dir    <- read_env_str("ADAPTME_VALIDATION_DIR", "adaptme-validation-results")
n_replicates  <- read_env_int("ADAPTME_REPLICATES",   5)   # use 25+ for manuscript
n_subjects    <- read_env_int("ADAPTME_SUBJECTS",    20)   # use 40+ for manuscript
n_timepoints  <- read_env_int("ADAPTME_TIMEPOINTS",   2)   # use 3+  for manuscript
n_taxa        <- read_env_int("ADAPTME_TAXA",        15)   # use 40+ for manuscript
alpha         <- read_env_dbl("ADAPTME_ALPHA",       0.05)
real_data_path <- read_env_str("ADAPTME_REAL_DATA")
cond_var      <- read_env_str("ADAPTME_COND_VAR")
subject_var   <- read_env_str("ADAPTME_SUBJECT_VAR")
adj_var_raw   <- read_env_str("ADAPTME_ADJ_VAR")
base_cond     <- read_env_str("ADAPTME_BASE_COND")
adj_var       <- if (nchar(adj_var_raw) > 0) trimws(strsplit(adj_var_raw, ",")[[1]]) else NULL

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("ADAPT-ME full test pipeline\n")
cat(sprintf("  output dir:  %s\n", output_dir))
cat(sprintf("  replicates:  %d  |  subjects: %d  |  timepoints: %d  |  taxa: %d\n",
            n_replicates, n_subjects, n_timepoints, n_taxa))
cat(sprintf("  alpha:       %.3f\n\n", alpha))


# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 – Unit tests
# ══════════════════════════════════════════════════════════════════════════════
section("STEP 1: Unit tests")

unit_ok <- tryCatch({
  result <- testthat::test_dir(
    system.file("../tests/testthat", package = "ADAPTME"),
    reporter = "progress",
    stop_on_failure = FALSE
  )
  failed <- sum(as.data.frame(result)$failed)
  cat(sprintf("\nUnit tests: %d failed\n", failed))
  failed == 0
}, error = function(e) {
  # fallback when running from source tree
  result <- testthat::test_dir("tests/testthat", reporter = "progress",
                                stop_on_failure = FALSE)
  failed <- sum(as.data.frame(result)$failed)
  cat(sprintf("\nUnit tests: %d failed\n", failed))
  failed == 0
})

if (!unit_ok) {
  cat("WARNING: unit tests have failures — review before relying on downstream results.\n")
}


# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 – Simulation validation (all 6 scenarios)
# ══════════════════════════════════════════════════════════════════════════════
section("STEP 2: Simulation validation (6 scenarios)")

validation <- run_adaptme_simulation_validation(
  scenarios    = adaptme_validation_scenarios(),
  n_replicates = n_replicates,
  seed         = 1,
  n_subjects   = n_subjects,
  n_timepoints = n_timepoints,
  n_taxa       = n_taxa,
  alpha        = alpha
)
val_summary <- summarize_adaptme_validation(validation)

write.csv(validation,  file.path(output_dir, "step2-replicate-metrics.csv"),  row.names = FALSE)
write.csv(val_summary, file.path(output_dir, "step2-scenario-summary.csv"),   row.names = FALSE)

cat("\nScenario-level summary:\n")
print(val_summary)

null_fpr <- val_summary$fpr[val_summary$scenario == "null"]
if (length(null_fpr) > 0 && !is.na(null_fpr)) {
  cat(sprintf("\n[check] Null FPR = %.3f  (target ≤ %.3f): %s\n",
              null_fpr, alpha,
              if (null_fpr <= alpha + 0.02) "PASS" else "REVIEW"))
}
sig_power <- val_summary$power[val_summary$scenario == "strong_signal"]
if (length(sig_power) > 0 && !is.na(sig_power)) {
  cat(sprintf("[check] Strong-signal power = %.3f  (target > 0.5): %s\n",
              sig_power,
              if (sig_power > 0.5) "PASS" else "REVIEW (may need more subjects)"))
}


# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 – Baseline comparison (ADAPT-ME vs naive LM)
# ══════════════════════════════════════════════════════════════════════════════
section("STEP 3: Baseline comparison — ADAPT-ME vs naive LM")

cat("Running null condition (n_da=0) to compare FPR ...\n")
cmp_null <- compare_adaptme_vs_naive(
  n_sims       = n_replicates,
  seed         = 100,
  n_subjects   = n_subjects,
  n_timepoints = n_timepoints,
  n_taxa       = n_taxa,
  n_da         = 0,            # pure null — only FPR matters
  effect_log10 = 0,
  alpha        = alpha
)

cat("Running signal condition (n_da=3) to compare power ...\n")
cmp_signal <- compare_adaptme_vs_naive(
  n_sims       = n_replicates,
  seed         = 200,
  n_subjects   = n_subjects,
  n_timepoints = n_timepoints,
  n_taxa       = n_taxa,
  n_da         = 3,
  effect_log10 = 0.5,
  alpha        = alpha
)

cmp_all <- rbind(
  cbind(condition = "null",   cmp_null),
  cbind(condition = "signal", cmp_signal)
)
write.csv(cmp_all, file.path(output_dir, "step3-baseline-comparison.csv"), row.names = FALSE)

cmp_summary <- aggregate(
  cbind(fpr, fdr, power) ~ method + condition,
  data = cmp_all,
  FUN  = function(x) round(mean(x, na.rm = TRUE), 3)
)
cat("\nComparison summary (mean across replicates):\n")
print(cmp_summary)

null_adaptme <- cmp_summary$fpr[cmp_summary$method == "adaptme"  & cmp_summary$condition == "null"]
null_naive   <- cmp_summary$fpr[cmp_summary$method == "naive_lm" & cmp_summary$condition == "null"]
if (length(null_adaptme) > 0 && length(null_naive) > 0 &&
    !is.na(null_adaptme) && !is.na(null_naive)) {
  cat(sprintf("\n[check] ADAPT-ME null FPR = %.3f, naive_lm null FPR = %.3f: %s\n",
              null_adaptme, null_naive,
              if (null_adaptme <= null_naive) "ADAPT-ME ≤ naive (expected)" else "REVIEW"))
}


# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 – Stress tests
# ══════════════════════════════════════════════════════════════════════════════
section("STEP 4: Stress tests")

stress_validation <- run_adaptme_simulation_validation(
  scenarios    = adaptme_stress_scenarios(),
  n_replicates = max(5L, as.integer(n_replicates / 5)),
  seed         = 300,
  n_subjects   = n_subjects,
  n_timepoints = n_timepoints,
  n_taxa       = n_taxa,
  alpha        = alpha
)
stress_summary <- summarize_adaptme_validation(stress_validation)

write.csv(stress_validation, file.path(output_dir, "step4-stress-metrics.csv"),  row.names = FALSE)
write.csv(stress_summary,    file.path(output_dir, "step4-stress-summary.csv"),  row.names = FALSE)

cat("\nStress-test summary:\n")
print(stress_summary)

error_rows <- stress_validation[!is.na(stress_validation$error), ]
if (nrow(error_rows) > 0) {
  cat(sprintf("\n[check] %d/%d stress replicates errored (non-zero is acceptable for extreme scenarios):\n",
              nrow(error_rows), nrow(stress_validation)))
  print(error_rows[, c("scenario", "replicate", "error")])
} else {
  cat("\n[check] All stress replicates completed without error: PASS\n")
}


# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 – Real data (optional)
# ══════════════════════════════════════════════════════════════════════════════
section("STEP 5: Real-data analysis (optional)")

# ── Required phyloseq format ───────────────────────────────────────────────────
#
#  phyloseq(
#    otu_table(counts_matrix, taxa_are_rows = FALSE),   # raw integer counts
#    sample_data(metadata_dataframe)                    # one row per sample
#  )
#
#  metadata must contain:
#    • subject variable  – one stable ID per person (e.g. "SubjectID")
#                          the SAME ID must appear in EVERY row for that person
#    • condition variable – numeric (continuous/binary) OR categorical (2+ levels)
#                           e.g. "Timepoint" = c("baseline","week4","week8")
#                               "BMI_zscore" = c(-0.3, 1.2, ...)
#    • adjustment vars    – optional; any numeric or categorical columns
#                           (no missing values; do NOT include the subject var)
#
#  count table rules:
#    • raw counts (integers), NOT relative abundances or CLR-transformed values
#    • NAs are not allowed
#    • taxa_are_rows can be TRUE or FALSE (handled automatically)
#
#  typical call:
#    adapt(
#      input_data  = my_phyloseq,
#      cond.var    = "Timepoint",   base.cond = "baseline",
#      adj.var     = c("Age", "Sex"),
#      subject.var = "SubjectID",
#      prev.filter = 0.10,          # drop taxa present in < 10 % of samples
#      depth.filter = 1000,         # drop samples with < 1000 reads
#      alpha       = 0.05
#    )
# ──────────────────────────────────────────────────────────────────────────────

if (nchar(real_data_path) == 0) {
  cat("No real dataset provided. Set ADAPTME_REAL_DATA to a .rds phyloseq path to run step 5.\n")
} else if (nchar(cond_var) == 0 || nchar(subject_var) == 0) {
  cat("ADAPTME_COND_VAR and ADAPTME_SUBJECT_VAR must be set when providing real data.\n")
} else {
  cat(sprintf("Loading real data from: %s\n", real_data_path))
  phyobj <- readRDS(real_data_path)

  cat(sprintf("  %d samples, %d taxa\n",
              phyloseq::nsamples(phyobj), phyloseq::ntaxa(phyobj)))
  cat(sprintf("  condition: %s  |  subject: %s\n", cond_var, subject_var))
  if (!is.null(adj_var)) cat(sprintf("  adjustment: %s\n", paste(adj_var, collapse = ", ")))

  real_result <- tryCatch(
    adapt(
      input_data   = phyobj,
      cond.var     = cond_var,
      base.cond    = if (nchar(base_cond) > 0) base_cond else NULL,
      adj.var      = adj_var,
      subject.var  = subject_var,
      prev.filter  = 0.10,
      depth.filter = 1000,
      alpha        = alpha
    ),
    error = function(e) {
      cat("ERROR during real-data analysis:\n")
      cat(" ", conditionMessage(e), "\n")
      NULL
    }
  )

  if (!is.null(real_result)) {
    real_summary <- summary(real_result, select = "da")
    out_path <- file.path(output_dir, "step5-real-data-results.csv")
    if (!is.null(real_summary)) {
      write.csv(real_summary, out_path, row.names = TRUE)
      cat(sprintf("DA results written to: %s\n", out_path))
    }
    plot_path <- file.path(output_dir, "step5-real-data-volcano.pdf")
    grDevices::pdf(plot_path, width = 7, height = 5)
    print(plot(real_result, n.label = 10))
    grDevices::dev.off()
    cat(sprintf("Volcano plot written to: %s\n", plot_path))
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# Done
# ══════════════════════════════════════════════════════════════════════════════
section("All steps complete")
cat(sprintf("Results written to: %s/\n", output_dir))
cat("Files:\n")
for (f in list.files(output_dir, full.names = FALSE)) {
  cat(sprintf("  %s\n", f))
}
