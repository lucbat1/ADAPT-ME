# ADAPT-ME

**Analysis of Microbiome Differential Abundance by Pooling Tobit Models with Mixed Effects**

`ADAPT-ME` extends the [ADAPT](https://github.com/mkbwang/ADAPT) framework for longitudinal and repeated-measures microbiome studies. It adds a subject-level random intercept to the Tobit count-ratio model, resolving the independence assumption that limits cross-sectional methods when applied to paired data.

## Key innovations over ADAPT

1. **Mixed-effects Tobit models** — Subject random intercepts are integrated out via 15-node Gauss-Hermite quadrature, removing between-subject baseline variation and improving power in paired and longitudinal designs.
2. **Stability-aware reference selection** — The BUM-guided reference selection score combines fold-change distance with within-subject standard deviation, preventing temporally unstable taxa from being selected as references.
3. **Cross-sectional ablation** — `adapt_crosssectional()` provides an isolated fixed-effects baseline (no random intercept, no `subject.var`) for method comparison and single-timepoint studies.

## Installation

```r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("lucbat1/ADAPT-ME", build_vignettes = TRUE)
```

## Quick start

```r
library(ADAPTME)

# Longitudinal / paired data
result <- adapt(
  input_data   = my_phyloseq,
  cond.var     = "Timepoint",
  base.cond    = "baseline",
  adj.var      = "Sex",
  subject.var  = "SubjectID",
  prev.filter  = 0.10,
  depth.filter = 1000,
  alpha        = 0.05
)

summary(result, select = "da")
plot(result, n.label = 5)

# Cross-sectional data (no repeated measures)
result_cs <- adapt_crosssectional(
  input_data   = my_phyloseq,
  cond.var     = "CaseStatus",
  base.cond    = "Control",
  prev.filter  = 0.10,
  depth.filter = 1000
)
```

## Input format

`adapt()` and `adapt_crosssectional()` accept a **phyloseq** object with:
- A count table of raw integer counts (not relative abundances or CLR-transformed)
- Sample metadata containing the condition variable, subject variable (for `adapt()`), and any adjustment covariates
- Taxonomy table (optional, used in `summary()` output)

## Simulation validation

`ADAPT-ME` ships with a full simulation validation harness:

```r
# Run all 6 validation scenarios (25 replicates each)
validation <- run_adaptme_simulation_validation(
  scenarios    = adaptme_validation_scenarios(),
  n_replicates = 25,
  n_subjects   = 40,
  n_timepoints = 3,
  n_taxa       = 40,
  alpha        = 0.05
)
summarize_adaptme_validation(validation)
```

From the command line:

```sh
Rscript inst/scripts/run-all-tests.R
```

Results are written to `adaptme-benchmark-results/`.

## Real data validation

Validated on the iHMP HMP2 16S stool dataset (EH1037, 42 subjects × 2 visits).
ADAPT-ME detected 13 DA taxa at FDR < 5%; the cross-sectional ablation detected none —
demonstrating the critical contribution of the subject random intercept in paired designs.

## Citation

If you use ADAPT-ME, please cite:

> Batista, L. and Kaushik, A. (2025). ADAPT-ME: Analysis of Microbiome Differential Abundance by Pooling Tobit Mixed-Effects Models. *Bioinformatics Advances*. doi: to be added upon acceptance.

And the original ADAPT method:

> Wang, M., Fontaine, S., Jiang, H., & Li, G. (2024). ADAPT: Analysis of Microbiome Differential Abundance by Pooling Tobit Models. *bioRxiv*. https://doi.org/10.1101/2024.05.14.594186
