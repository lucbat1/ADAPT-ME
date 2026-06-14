#!/usr/bin/env Rscript
# run-hmp2-analysis.R
#
# Real-data validation of ADAPT-ME using the iHMP 16S stool data (EH1037).
#
# Dataset:  Human Microbiome Project Phase 2 (iHMP), 16S V1-3 hypervariable
#           region, stool body site. Available on ExperimentHub without
#           needing the HMP2Data package (which has a dplyr incompatibility).
#
# Design:   42 healthy adults sampled at Visit 1 (baseline) and Visit 2
#           (follow-up). Condition = visit number (1→2). Since this is a
#           healthy cohort with no intervention between visits, we expect
#           mostly stable microbiome — a real-data FPR validation.
#           Sex is included as an adjustment covariate.
#
# Run from the package root:
#   Rscript inst/scripts/run-hmp2-analysis.R
#
# Or override the prevalence filter (trade speed for resolution):
#   ADAPTME_HMP2_PREV=0.30 Rscript inst/scripts/run-hmp2-analysis.R  # ~15 min
#   ADAPTME_HMP2_PREV=0.20 Rscript inst/scripts/run-hmp2-analysis.R  # ~30 min (default)
#   ADAPTME_HMP2_PREV=0.10 Rscript inst/scripts/run-hmp2-analysis.R  # ~60 min
#
# Output goes to: adaptme-hmp2-results/

suppressPackageStartupMessages({
  if (!requireNamespace("ADAPTME", quietly = TRUE)) {
    devtools::load_all(".")
  } else {
    library(ADAPTME)
  }
  library(ExperimentHub)
  library(SummarizedExperiment)
  library(phyloseq)
})

output_dir  <- "adaptme-hmp2-results"
prev_filter <- as.numeric(Sys.getenv("ADAPTME_HMP2_PREV", unset = "0.20"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load raw data from ExperimentHub ───────────────────────────────────────
cat("Loading iHMP 16S V1-3 data from ExperimentHub (EH1037)...\n")
hub <- ExperimentHub()
se  <- hub[["EH1037"]]
cat(sprintf("  Raw dimensions: %d OTUs x %d samples\n", nrow(se), ncol(se)))

# ── 2. Filter to stool, two-visit subjects ────────────────────────────────────
cat("\nSubsetting to stool body site...\n")
meta_all <- as.data.frame(colData(se))

stool_idx <- which(meta_all$hmp_body_subsite == "Stool")
meta_stool <- meta_all[stool_idx, , drop = FALSE]

# Identify subjects with exactly visits 1 AND 2 (drop the one subject with 3)
visits_per_subject <- tapply(meta_stool$visit_number,
                              meta_stool$run_sample_id,
                              function(v) sort(unique(v)))
two_visit_subjects <- names(Filter(function(v) 1 %in% v && 2 %in% v,
                                    visits_per_subject))

cat(sprintf("  Subjects with both Visit 1 and Visit 2: %d\n",
            length(two_visit_subjects)))

# Keep only those subjects, visits 1 and 2 only
keep_mask <- meta_stool$run_sample_id %in% two_visit_subjects &
  meta_stool$visit_number %in% c(1L, 2L)
meta_final <- meta_stool[keep_mask, , drop = FALSE]

# Rename columns to clean names for ADAPT-ME
meta_clean <- data.frame(
  subject      = as.character(meta_final$run_sample_id),
  visit        = as.integer(meta_final$visit_number),
  sex          = as.character(meta_final$sex),
  row.names    = rownames(meta_final),
  stringsAsFactors = FALSE
)

cat(sprintf("  Final sample set: %d samples (%d subjects x ~2 visits)\n",
            nrow(meta_clean), length(unique(meta_clean$subject))))
cat(sprintf("  Sex breakdown: %s\n",
            paste(names(table(meta_clean$sex)), table(meta_clean$sex),
                  sep = "=", collapse = ", ")))

# ── 3. Build count matrix ─────────────────────────────────────────────────────
cat("\nExtracting count matrix...\n")
counts_raw <- assay(se, "16SrRNA")[, rownames(meta_clean), drop = FALSE]
# samples as rows, OTUs as columns  (phyloseq taxa_are_rows = FALSE)
counts_mat <- t(counts_raw)
cat(sprintf("  Count matrix: %d samples x %d OTUs\n",
            nrow(counts_mat), ncol(counts_mat)))
cat(sprintf("  Zero fraction: %.1f%%\n", 100 * mean(counts_mat == 0)))

# ── 4. Parse taxonomy ─────────────────────────────────────────────────────────
cat("\nParsing GreenGenes taxonomy strings...\n")
lineages     <- as.character(rowData(se)$consensus_lineage)
names(lineages) <- rownames(se)
lineages     <- lineages[colnames(counts_mat)]   # keep only our OTUs

tax_levels   <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")

tax_mat <- t(vapply(seq_along(lineages), function(i) {
  parts  <- trimws(strsplit(lineages[i], ";")[[1]])
  result <- setNames(rep(NA_character_, length(tax_levels)), tax_levels)
  for (j in seq_along(tax_prefixes)) {
    hit <- grep(paste0("^", tax_prefixes[j]), parts, value = TRUE)
    if (length(hit) > 0) {
      val <- sub(tax_prefixes[j], "", hit[1])
      if (nchar(trimws(val)) > 0) result[tax_levels[j]] <- trimws(val)
    }
  }
  result
}, character(length(tax_levels))))
rownames(tax_mat) <- names(lineages)

cat(sprintf("  Phyla represented: %s\n",
            paste(sort(unique(na.omit(tax_mat[, "Phylum"]))), collapse = ", ")))

# ── 5. Assemble phyloseq object ───────────────────────────────────────────────
cat("\nAssembling phyloseq object...\n")
phyobj <- phyloseq(
  otu_table(counts_mat, taxa_are_rows = FALSE),
  sample_data(meta_clean),
  tax_table(tax_mat)
)
cat(sprintf("  phyloseq: %d taxa, %d samples\n",
            ntaxa(phyobj), nsamples(phyobj)))

# ── 6. Sanity checks before running ADAPT-ME ──────────────────────────────────
cat("\nPre-flight sanity checks...\n")

# Check subject × visit structure
visit_tbl <- table(meta_clean$subject, meta_clean$visit)
cat(sprintf("  All subjects have exactly 2 samples: %s\n",
            all(rowSums(visit_tbl > 0) == 2)))
cat(sprintf("  Min library size: %d reads\n", min(sample_sums(phyobj))))
cat(sprintf("  Samples with < 1000 reads: %d\n",
            sum(sample_sums(phyobj) < 1000)))

# Estimate usable taxa after prevalence filter
prev_10pct <- mean(colMeans(counts_mat > 0) > 0.10)
cat(sprintf("  OTUs passing 10%% prevalence filter: ~%d of %d\n",
            round(prev_10pct * ncol(counts_mat)), ncol(counts_mat)))

# ── 7. Run ADAPT-ME ───────────────────────────────────────────────────────────
cat("\n── Running ADAPT-ME ────────────────────────────────────────────────────\n")
cat(sprintf("  cond.var    = 'visit'       (1 = baseline, 2 = follow-up)\n"))
cat(sprintf("  subject.var = 'subject'     (subject random intercept)\n"))
cat(sprintf("  adj.var     = 'sex'         (adjustment covariate)\n"))
cat(sprintf("  prev.filter = %.2f          (present in >= %.0f%% of samples)\n",
            prev_filter, prev_filter * 100))
cat(sprintf("  depth.filter = 1000         (remove very shallow samples)\n"))

# Rough runtime estimate (1.6 sec/taxon/pass, 2 passes for reference selection)
n_taxa_est <- sum(colMeans(counts_mat > 0) > prev_filter)
est_min    <- round(n_taxa_est * 1.6 * 2 / 60)
cat(sprintf("  ~%d taxa after filtering; estimated runtime: ~%d minutes\n\n",
            n_taxa_est, est_min))

t_start <- proc.time()[["elapsed"]]
result <- tryCatch(
  adapt(
    input_data   = phyobj,
    cond.var     = "visit",
    adj.var      = "sex",
    subject.var  = "subject",
    prev.filter  = prev_filter,
    depth.filter = 1000,
    alpha        = 0.05
  ),
  error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    NULL
  }
)
elapsed <- proc.time()[["elapsed"]] - t_start

if (!is.null(result)) {
  cat(sprintf("\nCompleted in %.1f seconds\n", elapsed))

  # ── 8. Save and report results ───────────────────────────────────────────────
  details    <- summary(result, select = "all")
  da_summary <- summary(result, select = "da")

  write.csv(details,
            file.path(output_dir, "hmp2-all-taxa-results.csv"),
            row.names = TRUE)

  cat("\n── Results ─────────────────────────────────────────────────────────────\n")
  cat(sprintf("  Reference taxa:  %d\n", length(result@reference)))
  n_da <- length(result@signal)
  if (n_da == 1 && result@signal[1] == "") n_da <- 0
  cat(sprintf("  DA taxa (FDR<5%%): %d\n", n_da))

  if (!is.null(da_summary) && nrow(da_summary) > 0) {
    cat("\n  Differentially abundant taxa:\n")
    show_cols <- intersect(c("Taxa", "log10foldchange", "pval", "adjusted_pval",
                              "prevalence", "Phylum", "Genus"),
                           colnames(da_summary))
    print(da_summary[, show_cols, drop = FALSE])
    write.csv(da_summary,
              file.path(output_dir, "hmp2-da-taxa.csv"),
              row.names = TRUE)
  } else {
    cat("  No taxa reached FDR < 5% — consistent with a stable healthy microbiome.\n")
    cat("  (This is the expected result and validates FPR control on real data.)\n")
  }

  # ── 9. FPR estimate on real null data ────────────────────────────────────────
  n_tested <- nrow(result@details)
  n_fit    <- sum(!is.na(result@details$pval))
  pvals    <- result@details$pval[!is.na(result@details$pval)]

  cat(sprintf("\n  Taxa tested / fit: %d / %d\n", n_tested, n_fit))

  # Rough calibration check: p-value histogram
  p_bins <- c("< 0.05" = mean(pvals < 0.05),
              "0.05–0.2" = mean(pvals >= 0.05 & pvals < 0.2),
              ">= 0.2" = mean(pvals >= 0.2))
  cat("  P-value distribution (should be roughly uniform if null):\n")
  for (nm in names(p_bins)) {
    cat(sprintf("    %-12s  %.1f%%\n", nm, 100 * p_bins[nm]))
  }

  # Save volcano plot
  plot_path <- file.path(output_dir, "hmp2-volcano.pdf")
  grDevices::pdf(plot_path, width = 7, height = 5)
  print(plot(result, n.label = 5))
  grDevices::dev.off()
  cat(sprintf("\n  Volcano plot saved: %s\n", plot_path))

  cat(sprintf("\nAll results written to: %s/\n", output_dir))
}
