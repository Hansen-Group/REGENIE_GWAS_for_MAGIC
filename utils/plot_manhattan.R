#!/usr/bin/env Rscript
#
# REGENIE Manhattan Plot Generator
# This script creates Manhattan plots from REGENIE GWAS results
#
# Usage: Rscript plot_manhattan.R <results_file> <output_file>
#
# Arguments:
#   results_file: Path to REGENIE results file (.regenie)
#   output_file: Path for output plot (e.g., manhattan_plot.png)
#

# Load required libraries
suppressPackageStartupMessages({
    if (!require("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
    }
    if (!require("data.table", quietly = TRUE)) {
        stop("Package 'data.table' is required. Install with: install.packages('data.table')")
    }
    if (!require("dplyr", quietly = TRUE)) {
        stop("Package 'dplyr' is required. Install with: install.packages('dplyr')")
    }
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    cat("Usage: Rscript plot_manhattan.R <results_file> <output_file>\n")
    cat("\nExample:\n")
    cat("  Rscript plot_manhattan.R results.regenie manhattan_plot.png\n")
    quit(status = 1)
}

results_file <- args[1]
output_file <- args[2]

# Check if input file exists
if (!file.exists(results_file)) {
    stop(paste("Error: File not found:", results_file))
}

cat("Reading GWAS results from:", results_file, "\n")

# Read results
df <- tryCatch({
    fread(results_file)
}, error = function(e) {
    stop(paste("Error reading file:", e$message))
})

cat("Loaded", nrow(df), "variants\n")

# Determine p-value column
if ("LOG10P" %in% colnames(df)) {
    df$P <- 10^(-df$LOG10P)
} else if ("P" %in% colnames(df)) {
    # P-value already present
} else {
    stop("Error: No p-value column found (expected 'LOG10P' or 'P')")
}

# Ensure required columns exist
required_cols <- c("CHROM", "GENPOS", "P")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
    stop(paste("Error: Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# Convert chromosome to numeric (handle X, Y, MT)
df$CHR <- as.character(df$CHROM)
df$CHR <- gsub("X", "23", df$CHR)
df$CHR <- gsub("Y", "24", df$CHR)
df$CHR <- gsub("^MT$|^M$", "26", df$CHR)
df$CHR <- as.numeric(df$CHR)

# Remove NA chromosomes
df <- df[!is.na(df$CHR), ]

# Calculate -log10(P)
df$LOG10P <- -log10(df$P)

# Remove infinite values
df <- df[is.finite(df$LOG10P), ]

cat("Processing", nrow(df), "variants after filtering\n")

# Order by chromosome and position
df <- df[order(df$CHR, df$GENPOS), ]

# Create cumulative position for x-axis
df$BPcum <- 0

chrom_lengths <- df %>%
    group_by(CHR) %>%
    summarise(max_pos = max(GENPOS, na.rm = TRUE)) %>%
    mutate(total = cumsum(as.numeric(max_pos)) - max_pos)

# Assign cumulative positions
for (chr in unique(df$CHR)) {
    chr_offset <- chrom_lengths$total[chrom_lengths$CHR == chr]
    df$BPcum[df$CHR == chr] <- df$GENPOS[df$CHR == chr] + chr_offset
}

# Calculate chromosome center positions for labels
axis_df <- df %>%
    group_by(CHR) %>%
    summarize(center = (max(BPcum) + min(BPcum)) / 2)

cat("Creating Manhattan plot...\n")

# Define significance thresholds
gwas_threshold <- -log10(5e-8)
suggestive_threshold <- -log10(1e-5)

# Create color palette (alternate colors by chromosome)
color_palette <- rep(c("#1f78b4", "#33a02c"), ceiling(max(df$CHR) / 2))

# Create Manhattan plot
p <- ggplot(df, aes(x = BPcum, y = LOG10P)) +
    geom_point(aes(color = factor(CHR)), alpha = 0.6, size = 1.2) +
    scale_color_manual(values = color_palette) +
    scale_x_continuous(
        label = axis_df$CHR,
        breaks = axis_df$center,
        expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    geom_hline(yintercept = gwas_threshold, color = "red", linetype = "dashed", linewidth = 0.5) +
    geom_hline(yintercept = suggestive_threshold, color = "blue", linetype = "dashed", linewidth = 0.5) +
    labs(
        title = "GWAS Manhattan Plot",
        subtitle = paste("REGENIE results -", nrow(df), "variants"),
        x = "Chromosome",
        y = expression(-log[10](italic(P)))
    ) +
    theme_minimal() +
    theme(
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
    )

# Save plot
cat("Saving plot to:", output_file, "\n")

# Determine output format from file extension
output_ext <- tools::file_ext(output_file)

if (output_ext == "png") {
    ggsave(output_file, plot = p, width = 14, height = 6, dpi = 300)
} else if (output_ext == "pdf") {
    ggsave(output_file, plot = p, width = 14, height = 6)
} else if (output_ext %in% c("jpg", "jpeg")) {
    ggsave(output_file, plot = p, width = 14, height = 6, dpi = 300)
} else {
    # Default to PNG
    ggsave(output_file, plot = p, width = 14, height = 6, dpi = 300)
}

cat("Manhattan plot created successfully!\n")
cat("\nLegend:\n")
cat("  Red dashed line: Genome-wide significance (P = 5e-8)\n")
cat("  Blue dashed line: Suggestive threshold (P = 1e-5)\n")

# Print summary of significant hits
n_gwas_sig <- sum(df$P < 5e-8)
n_suggestive <- sum(df$P < 1e-5)

cat("\nResults summary:\n")
cat("  Genome-wide significant variants:", n_gwas_sig, "\n")
cat("  Suggestive variants:", n_suggestive, "\n")
