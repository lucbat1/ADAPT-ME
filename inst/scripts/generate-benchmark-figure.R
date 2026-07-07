#!/usr/bin/env Rscript
# generate-benchmark-figure.R
#
# Generates the benchmark performance panel figure from pre-computed results.
# Does NOT run any models — reads only from:
#   adaptme-benchmark-results/benchmark-summary.csv
#   adaptme-benchmark-results/benchmark-replicate-results.csv
#
# Output: figures/benchmark-panel.pdf  (and .png)
#
# Run from the package root:
#   Rscript inst/scripts/generate-benchmark-figure.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

output_dir <- "figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load data ─────────────────────────────────────────────────────────────────
summary_df   <- read.csv("adaptme-benchmark-results/benchmark-summary.csv",
                         stringsAsFactors = FALSE)
replicate_df <- read.csv("adaptme-benchmark-results/benchmark-replicate-results.csv",
                         stringsAsFactors = FALSE)

# ── Shared display mappings ───────────────────────────────────────────────────
method_order <- c("ADAPT-ME", "ANCOM-BC2", "MaAsLin2", "Naive LM", "ZINQ-L")
method_map   <- c(adaptme  = "ADAPT-ME",
                  ancombc2 = "ANCOM-BC2",
                  maaslin2 = "MaAsLin2",
                  naive_lm = "Naive LM",
                  zinql    = "ZINQ-L")
method_colors <- c("ADAPT-ME"   = "#E69F00",
                   "ANCOM-BC2"  = "#56B4E9",
                   "MaAsLin2"   = "#009E73",
                   "Naive LM"   = "#999999",
                   "ZINQ-L"     = "#CC79A7")

scenario_order_20 <- c("moderate_signal", "strong_signal",
                        "moderate_sparse", "sparse",
                        "unbalanced",      "confounded")
scenario_order_50 <- c("moderate_signal", "strong_signal", "sparse")

scenario_labels <- c(moderate_signal = "Moderate\nSignal",
                     strong_signal   = "Strong\nSignal",
                     moderate_sparse = "Moderate\nSparse",
                     sparse          = "Sparse",
                     unbalanced      = "Unbalanced",
                     confounded      = "Confounded")

apply_maps <- function(df, scenarios) {
  df <- df[df$scenario %in% scenarios, ]
  df$method_label   <- factor(method_map[df$method], levels = method_order)
  df$scenario_label <- factor(scenario_labels[df$scenario],
                              levels = scenario_labels[scenarios])
  df
}

# ── Panel A: Null FPR at n=20 (from replicate data) ──────────────────────────
null_rep <- replicate_df[replicate_df$scenario == "null" &
                           replicate_df$n_subjects == 20, ]
null_fpr <- aggregate(fpr ~ method, data = null_rep, FUN = mean, na.rm = TRUE)
null_fpr$method_label <- factor(method_map[null_fpr$method], levels = method_order)
null_fpr$fpr_pct      <- null_fpr$fpr * 100

p_fpr <- ggplot(null_fpr, aes(x = method_label, y = fpr_pct, fill = method_label)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 5, linetype = "dashed", colour = "red", linewidth = 0.5,
             na.rm = TRUE) +
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_y_continuous(limits = c(0, max(null_fpr$fpr_pct) * 1.25 + 0.5),
                     expand  = c(0, 0)) +
  labs(title = "Per-Taxon Type I Error Rate (Null Scenario, n = 20)",
       x = NULL, y = "Per-taxon FPR (%)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
        plot.title  = element_text(size = 10, face = "bold"))

# ── Panel B: FDR at n=20 signal scenarios ────────────────────────────────────
df20 <- apply_maps(summary_df[summary_df$n_subjects == 20, ], scenario_order_20)
df20$fdr_pct   <- df20$fdr   * 100
df20$power_pct <- df20$power * 100

p_fdr20 <- ggplot(df20, aes(x = scenario_label, y = fdr_pct, fill = method_label)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_hline(yintercept = 5, linetype = "dashed", colour = "red", linewidth = 0.5,
             na.rm = TRUE) +
  scale_fill_manual(values = method_colors, name = NULL) +
  scale_y_continuous(limits = c(0, 70), expand = c(0, 0)) +
  labs(title = "A.  FDR, Signal Scenarios (n = 20)",
       x = NULL, y = "FDR (%)") +
  theme_bw(base_size = 10) +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 8),
        axis.text.x      = element_text(size = 8),
        plot.title       = element_text(size = 10, face = "bold"))

# ── Panel C: Power at n=20 ────────────────────────────────────────────────────
p_power20 <- ggplot(df20, aes(x = scenario_label, y = power_pct, fill = method_label)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  scale_fill_manual(values = method_colors, name = NULL) +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
  labs(title = "A.  Power, Signal Scenarios (n = 20)",
       x = NULL, y = "Power (%)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x     = element_text(size = 8),
        plot.title      = element_text(size = 10, face = "bold"))

# ── Panel D: FDR at n=50 ─────────────────────────────────────────────────────
df50 <- apply_maps(summary_df[summary_df$n_subjects == 50, ], scenario_order_50)
df50$fdr_pct   <- df50$fdr   * 100
df50$power_pct <- df50$power * 100

p_fdr50 <- ggplot(df50, aes(x = scenario_label, y = fdr_pct, fill = method_label)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_hline(yintercept = 5, linetype = "dashed", colour = "red", linewidth = 0.5,
             na.rm = TRUE) +
  scale_fill_manual(values = method_colors, name = NULL) +
  scale_y_continuous(limits = c(0, 70), expand = c(0, 0)) +
  labs(title = "B.  FDR, Signal Scenarios (n = 50)",
       x = NULL, y = "FDR (%)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x     = element_text(size = 9),
        plot.title      = element_text(size = 10, face = "bold"))

# ── Panel E: Power at n=50 ───────────────────────────────────────────────────
p_power50 <- ggplot(df50, aes(x = scenario_label, y = power_pct, fill = method_label)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  scale_fill_manual(values = method_colors, name = NULL) +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
  labs(title = "B.  Power, Signal Scenarios (n = 50)",
       x = NULL, y = "Power (%)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x     = element_text(size = 9),
        plot.title      = element_text(size = 10, face = "bold"))

# ── Figure: Null FPR (standalone) ────────────────────────────────────────────
ggsave(file.path(output_dir, "fig-null-fpr.pdf"),
       p_fpr, width = 6, height = 4, device = "pdf")
ggsave(file.path(output_dir, "fig-null-fpr.png"),
       p_fpr, width = 6, height = 4, dpi = 300, device = "png")

# ── Figure: FDR (n=20 | n=50) ────────────────────────────────────────────────
fig_fdr <- (p_fdr20 | p_fdr50) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(output_dir, "fig-fdr.pdf"),
       fig_fdr, width = 13, height = 5, device = "pdf")
ggsave(file.path(output_dir, "fig-fdr.png"),
       fig_fdr, width = 13, height = 5, dpi = 300, device = "png")

# ── Figure: Power (n=20 | n=50) ──────────────────────────────────────────────
fig_power <- (p_power20 | p_power50) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(output_dir, "fig-power.pdf"),
       fig_power, width = 13, height = 5, device = "pdf")
ggsave(file.path(output_dir, "fig-power.png"),
       fig_power, width = 13, height = 5, dpi = 300, device = "png")

cat("Figures saved to", output_dir, "\n")
