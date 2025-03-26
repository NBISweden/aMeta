# This is a script for plotting KrakenUniq abundance matrix.
# Run this script as:
# Rscript plot_krakenuniq_abundance_matrix.R in_dir out_dir

args <- commandArgs(trailingOnly = TRUE)
in_dir <- as.character(args[1])
out_dir <- as.character(args[2])

library("pheatmap")

# Read abundance matrix and convert to numeric matrix
ku_abundance <- read.delim(paste0(in_dir, "/krakenuniq_abundance_matrix.txt"),
                           header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
ku_abundance <- as.matrix(ku_abundance)

# Create a version of the matrix for display, rounded for readability
ku_abundance_display <- round(ku_abundance)

# Function for tuning font size on microbial abundance heatmap
my_fontsize <- function(mat) {
  if (nrow(mat) < 50) {
    return(12)
  } else if (nrow(mat) < 100) {
    return(10)
  } else if (nrow(mat) < 150) {
    return(8)
  } else {
    return(6)
  }
}

# Absolute abundance heatmap
pdf(paste0(out_dir, "/krakenuniq_absolute_abundance_heatmap.pdf"), paper = "a4r", width = 297, height = 210)
if (nrow(ku_abundance) > 1 & ncol(ku_abundance) > 1) {
  pheatmap(ku_abundance, display_numbers = ku_abundance_display, fontsize = my_fontsize(ku_abundance),
           main = "KrakenUniq Absolute Microbial Abundance", cluster_rows = FALSE, cluster_cols = FALSE,
           number_format = "%i")
} else {
  pheatmap(ku_abundance, display_numbers = ku_abundance_display, fontsize = 8,
           main = "KrakenUniq Absolute Microbial Abundance", cluster_rows = FALSE, cluster_cols = FALSE,
           number_format = "%i", breaks = c(0, 1))
}
dev.off()

# Normalise by sequencing depth (column-wise)
for (i in 1:ncol(ku_abundance)) {
  ku_abundance[, i] <- ku_abundance[, i] / sum(ku_abundance[, i])
}
ku_abundance_display <- round(ku_abundance, digits = 3)

# Normalised abundance heatmap
pdf(paste0(out_dir, "/krakenuniq_normalized_abundance_heatmap.pdf"), paper = "a4r", width = 297, height = 210)
if (nrow(ku_abundance) > 1 & ncol(ku_abundance) > 1) {
  pheatmap(ku_abundance, display_numbers = ku_abundance_display, fontsize = my_fontsize(ku_abundance),
           main = "KrakenUniq Normalized Microbial Abundance", cluster_rows = FALSE, cluster_cols = FALSE,
           number_format = "%.3f")
} else {
  pheatmap(ku_abundance, display_numbers = ku_abundance_display, fontsize = 8,
           main = "KrakenUniq Normalized Microbial Abundance", cluster_rows = FALSE, cluster_cols = FALSE,
           number_format = "%.3f", breaks = c(0, 1))
}
dev.off()

# Top 20 species × their most enriched samples (including Homo sapiens)
top_n <- 20
species_sums <- rowSums(ku_abundance_display)

# Get top N species (including Homo sapiens)
top_species <- names(sort(species_sums, decreasing = TRUE))[1:min(top_n, length(species_sums))]

# Subset matrices
ku_top <- ku_abundance[top_species, , drop = FALSE]
ku_top_display <- ku_abundance_display[top_species, , drop = FALSE]

# Get most enriched sample per species
top_samples <- apply(ku_top, 1, function(x) names(which.max(x)))
selected_samples <- unique(top_samples)[1:min(top_n, length(unique(top_samples)))]

# Subset to selected samples
ku_top <- ku_top[, selected_samples, drop = FALSE]
ku_top_display <- ku_top_display[, selected_samples, drop = FALSE]

# Plot focused heatmap (20xN)
if (nrow(ku_top) > 0 & ncol(ku_top) > 0) {
  pdf(paste0(out_dir, "/krakenuniq_top", top_n, "x", ncol(ku_top), "_abundance_heatmap.pdf"), paper = "a4r", width = 297, height = 210)
  if (nrow(ku_top) > 1 & ncol(ku_top) > 1) {
    pheatmap(ku_top, display_numbers = ku_top_display, fontsize = my_fontsize(ku_top),
             main = paste0("Top ", nrow(ku_top), " Species × Top ", ncol(ku_top), " Samples"),
             cluster_rows = FALSE, cluster_cols = FALSE, number_format = "%i")
  } else {
    pheatmap(ku_top, display_numbers = ku_top_display, fontsize = 8,
             main = paste0("Top ", nrow(ku_top), " Species × Top ", ncol(ku_top), " Samples"),
             cluster_rows = FALSE, cluster_cols = FALSE, number_format = "%i", breaks = c(0, 1))
  }
  dev.off()
}
