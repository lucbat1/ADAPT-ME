#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ADAPTME))

read_env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (identical(value, "")) {
    return(default)
  }
  as.integer(value)
}

read_env_numeric <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (identical(value, "")) {
    return(default)
  }
  as.numeric(value)
}

output_dir <- Sys.getenv("ADAPTME_VALIDATION_DIR", unset = "adaptme-validation-results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_replicates <- read_env_int("ADAPTME_REPLICATES", 25)
n_subjects <- read_env_int("ADAPTME_SUBJECTS", 40)
n_timepoints <- read_env_int("ADAPTME_TIMEPOINTS", 3)
n_taxa <- read_env_int("ADAPTME_TAXA", 40)
seed <- read_env_int("ADAPTME_SEED", 1)
alpha <- read_env_numeric("ADAPTME_ALPHA", 0.05)

cat("Running ADAPT-ME simulation validation\n")
cat(sprintf("  replicates: %d\n", n_replicates))
cat(sprintf("  subjects:   %d\n", n_subjects))
cat(sprintf("  timepoints: %d\n", n_timepoints))
cat(sprintf("  taxa:       %d\n", n_taxa))
cat(sprintf("  alpha:      %.3f\n", alpha))
cat(sprintf("  output:     %s\n", output_dir))

validation <- run_adaptme_simulation_validation(
  scenarios = adaptme_validation_scenarios(),
  n_replicates = n_replicates,
  seed = seed,
  n_subjects = n_subjects,
  n_timepoints = n_timepoints,
  n_taxa = n_taxa,
  alpha = alpha
)
summary <- summarize_adaptme_validation(validation)

write.csv(validation, file.path(output_dir, "replicate-metrics.csv"), row.names = FALSE)
write.csv(summary, file.path(output_dir, "scenario-summary.csv"), row.names = FALSE)

cat("Done.\n")
print(summary)
