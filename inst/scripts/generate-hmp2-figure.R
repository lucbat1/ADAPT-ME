#!/usr/bin/env Rscript
# generate-hmp2-figure.R
#
# Generates the HMP2 ME vs CS side-by-side volcano figure from pre-computed results.
# Does NOT re-run any models — reads only from:
#   adaptme-hmp2-results/hmp2-me-all-taxa-results.csv
#   adaptme-hmp2-results/hmp2-cs-all-taxa-results.csv
#   adaptme-hmp2-results/hmp2-phyloseq.rds  (for genus labels on DA taxa)
#
# Output: figures/hmp2-volcano-comparison.pdf  (and .png)
#
# Run from the package root:
#   Rscript inst/scripts/generate-hmp2-figure.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(phyloseq)
})

output_dir <- "figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load results ──────────────────────────────────────────────────────────────
me_df <- read.csv("adaptme-hmp2-results/hmp2-me-all-taxa-results.csv",
                  stringsAsFactors = FALSE, row.names = 1)
cs_df <- read.csv("adaptme-hmp2-results/hmp2-cs-all-taxa-results.csv",
                  stringsAsFactors = FALSE, row.names = 1)

# ── Load taxonomy from cached phyloseq ────────────────────────────────────────
phyobj  <- readRDS("adaptme-hmp2-results/hmp2-phyloseq.rds")
tax_tbl <- as.data.frame(tax_table(phyobj), stringsAsFactors = FALSE)
tax_tbl$Taxa <- rownames(tax_tbl)

# ── Annotate results with genus ───────────────────────────────────────────────
annotate_results <- function(df, alpha = 0.05) {
  df <- df[!is.na(df$pval), ]
  df$neglog10p <- -log10(df$pval)
  df$is_da     <- !is.na(df$adjusted_pval) & df$adjusted_pval < alpha

  df <- merge(df, tax_tbl[, c("Taxa", "Genus", "Phylum")],
              by = "Taxa", all.x = TRUE)

  # Build display label for DA taxa: use Genus if available, else OTU ID
  df$label <- ""
  da_idx   <- which(df$is_da)
  if (length(da_idx) > 0) {
    genus_clean <- ifelse(
      !is.na(df$Genus[da_idx]) & nchar(trimws(df$Genus[da_idx])) > 0,
      df$Genus[da_idx],
      df$Taxa[da_idx]
    )
    df$label[da_idx] <- genus_clean
  }
  df
}

me_ann <- annotate_results(me_df)
cs_ann <- annotate_results(cs_df)

n_me_da <- sum(me_ann$is_da)
n_cs_da <- sum(cs_ann$is_da)

cat(sprintf("ADAPT-ME DA taxa: %d\n", n_me_da))
cat(sprintf("ADAPT-CS DA taxa: %d\n", n_cs_da))

# ── Shared plot limits ────────────────────────────────────────────────────────
all_fc    <- c(me_ann$log10foldchange, cs_ann$log10foldchange)
all_negp  <- c(me_ann$neglog10p,      cs_ann$neglog10p)
x_lim     <- c(min(all_fc,   na.rm = TRUE) * 1.05,
                max(all_fc,   na.rm = TRUE) * 1.05)
y_lim     <- c(0, max(all_negp, na.rm = TRUE) * 1.10)

# ── Build one volcano panel ───────────────────────────────────────────────────
make_volcano <- function(df, panel_label, x_lab) {
  p <- ggplot(df, aes(x = log10foldchange, y = neglog10p)) +
    geom_point(aes(colour = is_da), alpha = 0.75, size = 1.4) +
    scale_colour_manual(
      values = c("FALSE" = "#AAAAAA", "TRUE" = "#E63946"),
      guide  = "none"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "steelblue", linewidth = 0.4) +
    coord_cartesian(xlim = x_lim, ylim = y_lim) +
    labs(title = panel_label,
         x     = x_lab,
         y     = expression(-log[10](italic(p)))) +
    theme_bw(base_size = 11) +
    theme(plot.title  = element_text(size = 11, face = "bold"),
          axis.title  = element_text(size = 10),
          axis.text   = element_text(size = 9))

  # Add labels only if there are DA taxa
  if (any(df$label != "")) {
    p <- p + geom_label_repel(
      data          = df[df$label != "", ],
      aes(label     = label),
      size          = 2.8,
      max.overlaps  = 20,
      box.padding   = 0.4,
      point.padding = 0.3,
      segment.color = "grey50",
      label.size    = 0.2,
      fill          = alpha("white", 0.8)
    )
  }
  p
}

p_me <- make_volcano(
  me_ann,
  panel_label = sprintf("ADAPT-ME: Mixed-Effects Tobit (%d DA taxa)", n_me_da),
  x_lab       = expression("Log"[10] ~ "Fold Change (Visit 2 vs Visit 1)")
)

p_cs <- make_volcano(
  cs_ann,
  panel_label = sprintf("Cross-Sectional Ablation (%d DA taxa)", n_cs_da),
  x_lab       = expression("Log"[10] ~ "Fold Change (Visit 2 vs Visit 1)")
)

# ── Save figures separately ───────────────────────────────────────────────────
ggsave(file.path(output_dir, "fig-hmp-adaptme.pdf"),
       p_me, width = 6, height = 5, device = "pdf")
ggsave(file.path(output_dir, "fig-hmp-adaptme.png"),
       p_me, width = 6, height = 5, dpi = 300, device = "png")

ggsave(file.path(output_dir, "fig-hmp-ablation.pdf"),
       p_cs, width = 6, height = 5, device = "pdf")
ggsave(file.path(output_dir, "fig-hmp-ablation.png"),
       p_cs, width = 6, height = 5, dpi = 300, device = "png")

cat(sprintf("HMP figures saved to %s\n", output_dir))
