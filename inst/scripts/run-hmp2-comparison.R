#!/usr/bin/env Rscript
# run-hmp2-comparison.R
#
# Side-by-side comparison of ADAPT-ME vs ADAPT-CS (cross-sectional ablation)
# on the iHMP HMP2 16S stool dataset (EH1037).
#
# Builds the same phyloseq object as run-hmp2-analysis.R, then runs:
#   1. adapt()               — ADAPT-ME (mixed-effects Tobit, with subject RE)
#   2. adapt_crosssectional() — ADAPT-CS ablation (fixed-effects Tobit, no subject RE)
#
# The phyloseq object is cached to adaptme-hmp2-results/hmp2-phyloseq.rds
# so that subsequent runs skip the ExperimentHub download.
#
# Run from the package root:
#   Rscript inst/scripts/run-hmp2-comparison.R
#
# Override prevalence filter:
#   ADAPTME_HMP2_PREV=0.30 Rscript inst/scripts/run-hmp2-comparison.R

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library(ExperimentHub)
  library(SummarizedExperiment)
  library(phyloseq)
})

output_dir  <- "adaptme-hmp2-results"
prev_filter <- as.numeric(Sys.getenv("ADAPTME_HMP2_PREV", unset = "0.20"))
phyobj_path <- file.path(output_dir, "hmp2-phyloseq.rds")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load or build phyloseq ─────────────────────────────────────────────────

if (file.exists(phyobj_path)) {
  cat(sprintf("Loading cached phyloseq from %s ...\n", phyobj_path))
  phyobj <- readRDS(phyobj_path)
  cat(sprintf("  phyloseq: %d taxa, %d samples\n", ntaxa(phyobj), nsamples(phyobj)))
} else {
  cat("Loading iHMP 16S V1-3 data from ExperimentHub (EH1037)...\n")
  hub <- ExperimentHub()
  se  <- hub[["EH1037"]]
  cat(sprintf("  Raw dimensions: %d OTUs x %d samples\n", nrow(se), ncol(se)))

  meta_all   <- as.data.frame(colData(se))
  stool_idx  <- which(meta_all$hmp_body_subsite == "Stool")
  meta_stool <- meta_all[stool_idx, , drop = FALSE]

  visits_per_subject <- tapply(meta_stool$visit_number,
                                meta_stool$run_sample_id,
                                function(v) sort(unique(v)))
  two_visit_subjects <- names(Filter(function(v) 1 %in% v && 2 %in% v,
                                      visits_per_subject))
  keep_mask  <- meta_stool$run_sample_id %in% two_visit_subjects &
                meta_stool$visit_number  %in% c(1L, 2L)
  meta_final <- meta_stool[keep_mask, , drop = FALSE]

  meta_clean <- data.frame(
    subject = as.character(meta_final$run_sample_id),
    visit   = as.integer(meta_final$visit_number),
    sex     = as.character(meta_final$sex),
    row.names = rownames(meta_final),
    stringsAsFactors = FALSE
  )

  counts_raw <- assay(se, "16SrRNA")[, rownames(meta_clean), drop = FALSE]
  counts_mat <- t(counts_raw)

  lineages     <- as.character(rowData(se)$consensus_lineage)
  names(lineages) <- rownames(se)
  lineages     <- lineages[colnames(counts_mat)]
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

  phyobj <- phyloseq(
    otu_table(counts_mat, taxa_are_rows = FALSE),
    sample_data(meta_clean),
    tax_table(tax_mat)
  )
  cat(sprintf("  phyloseq: %d taxa, %d samples — caching to %s\n",
              ntaxa(phyobj), nsamples(phyobj), phyobj_path))
  saveRDS(phyobj, phyobj_path)
  cat(sprintf("  Subjects: %d  |  Sex: %s\n",
              length(unique(meta_clean$subject)),
              paste(names(table(meta_clean$sex)), table(meta_clean$sex),
                    sep = "=", collapse = ", ")))
}

cat(sprintf("\nPrevalence filter: %.2f  |  Depth filter: 1000\n", prev_filter))

# ── 2. Run ADAPT-ME ───────────────────────────────────────────────────────────

cat("\n── Running ADAPT-ME (mixed-effects Tobit) ──────────────────────────────\n")
t0 <- proc.time()[["elapsed"]]
result_me <- tryCatch(
  adapt(
    input_data   = phyobj,
    cond.var     = "visit",
    adj.var      = "sex",
    subject.var  = "subject",
    prev.filter  = prev_filter,
    depth.filter = 1000,
    alpha        = 0.05
  ),
  error = function(e) { cat("ADAPT-ME ERROR:", conditionMessage(e), "\n"); NULL }
)
t_me <- proc.time()[["elapsed"]] - t0
cat(sprintf("\nADAPT-ME completed in %.1f seconds\n", t_me))

# ── 3. Run ADAPT-CS ───────────────────────────────────────────────────────────

cat("\n── Running ADAPT-CS (cross-sectional ablation) ─────────────────────────\n")
t0 <- proc.time()[["elapsed"]]
result_cs <- tryCatch(
  adapt_crosssectional(
    input_data   = phyobj,
    cond.var     = "visit",
    adj.var      = "sex",
    prev.filter  = prev_filter,
    depth.filter = 1000,
    alpha        = 0.05
  ),
  error = function(e) { cat("ADAPT-CS ERROR:", conditionMessage(e), "\n"); NULL }
)
t_cs <- proc.time()[["elapsed"]] - t0
cat(sprintf("\nADAPT-CS completed in %.1f seconds\n", t_cs))

# ── 4. Build comparison table ─────────────────────────────────────────────────

safe_signal <- function(res) {
  if (is.null(res)) return(character(0))
  sig <- res@signal
  if (length(sig) == 1 && sig == "") character(0) else sig
}

me_sig <- safe_signal(result_me)
cs_sig <- safe_signal(result_cs)

n_me <- length(me_sig)
n_cs <- length(cs_sig)
both <- intersect(me_sig, cs_sig)
me_only <- setdiff(me_sig, cs_sig)
cs_only <- setdiff(cs_sig, me_sig)

cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("COMPARISON SUMMARY\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat(sprintf("  ADAPT-ME  DA taxa:   %d  (runtime: %.0f s)\n", n_me, t_me))
cat(sprintf("  ADAPT-CS  DA taxa:   %d  (runtime: %.0f s)\n", n_cs, t_cs))
cat(sprintf("  Shared:              %d\n", length(both)))
cat(sprintf("  ME only:             %d  %s\n", length(me_only),
            if (length(me_only) > 0) paste0("(", paste(me_only, collapse=", "), ")") else ""))
cat(sprintf("  CS only:             %d  %s\n", length(cs_only),
            if (length(cs_only) > 0) paste0("(", paste(cs_only, collapse=", "), ")") else ""))

# ── 5. Full per-taxon comparison table ───────────────────────────────────────

build_detail <- function(res, method_tag) {
  if (is.null(res)) return(NULL)
  d <- res@details
  d$method  <- method_tag
  d$da_flag <- d$Taxa %in% safe_signal(res)
  d
}

detail_me <- build_detail(result_me, "ADAPT-ME")
detail_cs <- build_detail(result_cs, "ADAPT-CS")

if (!is.null(detail_me) && !is.null(detail_cs)) {
  # Join on Taxa
  shared_taxa <- intersect(detail_me$Taxa, detail_cs$Taxa)
  me_sub <- detail_me[detail_me$Taxa %in% shared_taxa,
                       c("Taxa", "log10foldchange", "pval", "adjusted_pval",
                         "within_subject_sd", "da_flag")]
  cs_sub <- detail_cs[detail_cs$Taxa %in% shared_taxa,
                       c("Taxa", "log10foldchange", "pval", "adjusted_pval", "da_flag")]
  colnames(me_sub)[-1] <- paste0(colnames(me_sub)[-1], "_ME")
  colnames(cs_sub)[-1] <- paste0(colnames(cs_sub)[-1], "_CS")
  cmp_tbl <- merge(me_sub, cs_sub, by = "Taxa")
  cmp_tbl$da_either <- cmp_tbl$da_flag_ME | cmp_tbl$da_flag_CS

  write.csv(cmp_tbl,
            file.path(output_dir, "hmp2-comparison-table.csv"),
            row.names = FALSE)
  cat(sprintf("\nFull comparison table (%d shared taxa) written to:\n  %s\n",
              nrow(cmp_tbl), file.path(output_dir, "hmp2-comparison-table.csv")))

  # Print DA taxa from each method with taxonomy
  print_da <- function(res, label) {
    if (is.null(res)) return(invisible(NULL))
    sig <- safe_signal(res)
    if (length(sig) == 0) {
      cat(sprintf("\n%s: no DA taxa\n", label)); return(invisible(NULL))
    }
    d <- res@details[res@details$Taxa %in% sig, ]
    if (!is.null(res@input) && !is.null(tax_table(res@input, errorIfNULL = FALSE))) {
      tt <- as.data.frame(tax_table(res@input))
      tt$Taxa <- rownames(tt)
      d <- merge(d, tt[, c("Taxa", "Phylum", "Genus")], by = "Taxa", all.x = TRUE)
    }
    d <- d[order(d$adjusted_pval, na.last = TRUE), ]
    cat(sprintf("\n%s DA taxa (n=%d):\n", label, nrow(d)))
    show_cols <- intersect(c("Taxa", "Phylum", "Genus", "log10foldchange",
                              "pval", "adjusted_pval"), colnames(d))
    print(d[, show_cols, drop = FALSE], row.names = FALSE)
  }

  print_da(result_me, "ADAPT-ME")
  print_da(result_cs, "ADAPT-CS")
}

# ── 6. Side-by-side volcano plots (saved as PDF) ─────────────────────────────

if (!is.null(result_me) && !is.null(result_cs)) {
  plot_path <- file.path(output_dir, "hmp2-comparison-volcano.pdf")
  grDevices::pdf(plot_path, width = 14, height = 5)
  oldpar <- par(mfrow = c(1, 2))
  print(plot(result_me, n.label = 8))
  print(plot(result_cs, n.label = 8))
  par(oldpar)
  grDevices::dev.off()
  cat(sprintf("\nComparison volcano plots saved: %s\n", plot_path))
}

# ── 7. Save CS full results ───────────────────────────────────────────────────

if (!is.null(result_cs)) {
  cs_details <- result_cs@details
  write.csv(cs_details,
            file.path(output_dir, "hmp2-cs-all-taxa-results.csv"),
            row.names = TRUE)

  cs_da <- cs_details[cs_details$Taxa %in% cs_sig, ]
  if (nrow(cs_da) > 0) {
    write.csv(cs_da,
              file.path(output_dir, "hmp2-cs-da-taxa.csv"),
              row.names = TRUE)
  }
}

cat(sprintf("\nAll comparison outputs written to: %s/\n", output_dir))
