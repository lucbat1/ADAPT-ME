make_longitudinal_phyloseq <- function(effect = 5, categorical = FALSE) {
  set.seed(1)
  n_subjects <- 12
  subject <- rep(seq_len(n_subjects), each = 2)
  condition <- rep(c(0, 1), n_subjects)
  subject_baseline <- rep(rpois(n_subjects, 30) + 20, each = 2)
  counts <- cbind(
    TaxonA = round(subject_baseline * ifelse(condition == 1, effect, 1)),
    TaxonB = round(subject_baseline * 3),
    TaxonC = round(subject_baseline * 2),
    TaxonD = round(subject_baseline * 4)
  )
  rownames(counts) <- paste0("Sample", seq_len(nrow(counts)))
  metadata <- data.frame(
    condition = if (categorical) ifelse(condition == 1, "post", "pre") else condition,
    subject = subject,
    confounder = rep(rep(c(0, 1), length.out = n_subjects), each = 2),
    row.names = rownames(counts)
  )
  phyloseq::phyloseq(
    phyloseq::otu_table(counts, taxa_are_rows = FALSE),
    phyloseq::sample_data(metadata)
  )
}

test_that("longitudinal ADAPT detects a strong within-subject condition effect", {
  output <- adapt(
    input_data = make_longitudinal_phyloseq(),
    cond.var = "condition",
    adj.var = "confounder",
    subject.var = "subject",
    prev.filter = 0,
    depth.filter = 0,
    alpha = 0.1
  )
  results <- summary(output)

  expect_true("TaxonA" %in% output@signal)
  expect_true(all(c("within_subject_sd", "adjusted_pval") %in% colnames(results)))
  expect_false(any(is.na(results$pval)))
})

test_that("categorical conditions honor the requested baseline", {
  output <- adapt(
    input_data = make_longitudinal_phyloseq(categorical = TRUE),
    cond.var = "condition",
    base.cond = "pre",
    subject.var = "subject",
    prev.filter = 0,
    depth.filter = 0,
    alpha = 0.1
  )

  expect_match(output@DAAname, "post VS pre", fixed = TRUE)
  expect_true("TaxonA" %in% output@signal)
})

test_that("subject variable is required and must identify repeated measures", {
  phyobj <- make_longitudinal_phyloseq()
  phyloseq::sample_data(phyobj)$unique_subject <- seq_len(phyloseq::nsamples(phyobj))

  expect_error(
    adapt(phyobj, cond.var = "condition", prev.filter = 0, depth.filter = 0),
    "subject.var"
  )
  expect_error(
    adapt(phyobj, cond.var = "condition", subject.var = "unique_subject",
          prev.filter = 0, depth.filter = 0),
    "repeated observations"
  )
})

# ── Stress tests ──────────────────────────────────────────────────────────────

test_that("adapt handles extreme zero inflation without crashing", {
  # 75 % of counts will be zero; models should converge or return NA, not error
  sim <- simulate_adaptme_phyloseq(
    seed = 99, n_subjects = 8, n_timepoints = 2, n_taxa = 6,
    n_da = 1, effect_log10 = 0.8, zero_inflation = 0.75
  )
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
  expect_false(inherits(output, "error"),
               info = if (inherits(output, "error")) conditionMessage(output) else "")
})

test_that("adapt works on a minimal dataset (4 subjects, 2 timepoints)", {
  sim <- simulate_adaptme_phyloseq(
    seed = 100, n_subjects = 4, n_timepoints = 2, n_taxa = 4,
    n_da = 1, effect_log10 = 0.8, zero_inflation = 0.05
  )
  output <- suppressMessages(adapt(
    input_data   = sim$phyobj,
    cond.var     = "condition",
    adj.var      = "confounder",
    subject.var  = "subject",
    prev.filter  = 0,
    depth.filter = 0
  ))
  expect_s4_class(output, "DAresult")
  expect_true(nrow(output@details) >= 1)
})

test_that("adapt handles unbalanced visits gracefully", {
  sim <- simulate_adaptme_phyloseq(
    seed = 103, n_subjects = 8, n_timepoints = 3, n_taxa = 6,
    n_da = 1, effect_log10 = 0.6, unbalanced = TRUE, drop_fraction = 0.5
  )
  output <- suppressMessages(adapt(
    input_data   = sim$phyobj,
    cond.var     = "condition",
    subject.var  = "subject",
    prev.filter  = 0,
    depth.filter = 0
  ))
  expect_s4_class(output, "DAresult")
})
