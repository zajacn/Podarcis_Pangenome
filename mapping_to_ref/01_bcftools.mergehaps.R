#!/usr/bin/env Rscript

# Script to merge two haplotype VCF files into a single phased diploid VCF
# Usage: Rscript 01_bcftools.mergehaps.R hap1.vcf.gz hap2.vcf.gz output.vcf

# Load required libraries
library(VariantAnnotation)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript merge_haplotype_vcfs_v3.R hap1.vcf.gz hap2.vcf.gz output.vcf")
}

hap1_file <- args[1]
hap2_file <- args[2]
output_file <- args[3]

cat("Reading haplotype 1 VCF:", hap1_file, "\n")
cat("Reading haplotype 2 VCF:", hap2_file, "\n")

# Read VCF files
hap1_vcf <- readVcf(hap1_file, "genome")
hap2_vcf <- readVcf(hap2_file, "genome")

# Extract positions
hap1_ranges <- rowRanges(hap1_vcf)
hap2_ranges <- rowRanges(hap2_vcf)

cat("Found", length(hap1_ranges), "variants in haplotype 1\n")
cat("Found", length(hap2_ranges), "variants in haplotype 2\n")

# Create position keys for matching (VECTORIZED)
hap1_keys <- paste(seqnames(hap1_ranges), start(hap1_ranges), sep = ":")
hap2_keys <- paste(seqnames(hap2_ranges), start(hap2_ranges), sep = ":")

# Get all unique positions
all_keys <- unique(c(hap1_keys, hap2_keys))
cat("Total unique variant positions:", length(all_keys), "\n")

# VECTORIZED APPROACH - much faster!
cat("Assigning genotypes (vectorized)...\n")

# Check which keys are in each haplotype (vectorized operation)
in_hap1 <- all_keys %in% hap1_keys
in_hap2 <- all_keys %in% hap2_keys

# Assign genotypes based on presence in haplotypes (vectorized)
gt_values <- character(length(all_keys))
gt_values[in_hap1 & in_hap2] <- "1|1"    # Present in both
gt_values[in_hap1 & !in_hap2] <- "1|0"   # Present only in hap1
gt_values[!in_hap1 & in_hap2] <- "0|1"   # Present only in hap2

# Create named vector for easy lookup
names(gt_values) <- all_keys

cat("Building output VCF...\n")

# Build output VCF
# Start with all variants from hap1
output_vcf <- hap1_vcf

# Rename sample to a common name
colnames(output_vcf) <- "SAMPLE"

# Add variants that are only in hap2
hap2_only_keys <- setdiff(hap2_keys, hap1_keys)
if (length(hap2_only_keys) > 0) {
  hap2_only_idx <- which(hap2_keys %in% hap2_only_keys)
  hap2_only_vcf <- hap2_vcf[hap2_only_idx]
  
  # Rename sample to match
  colnames(hap2_only_vcf) <- "SAMPLE"
  
  # Combine
  output_vcf <- rbind(output_vcf, hap2_only_vcf)
}

cat("Sorting variants...\n")

# Sort by position
output_vcf <- output_vcf[order(as.character(seqnames(rowRanges(output_vcf))), 
                               start(rowRanges(output_vcf)))]

# Update genotypes
output_keys <- paste(seqnames(rowRanges(output_vcf)), 
                     start(rowRanges(output_vcf)), 
                     sep = ":")

# Get genotypes for output positions (vectorized lookup)
final_gt_values <- gt_values[output_keys]
gt_matrix <- matrix(final_gt_values, ncol = 1)
colnames(gt_matrix) <- "SAMPLE"

# Update the GT field
geno(output_vcf)$GT <- gt_matrix

# Write output
cat("Writing phased VCF to:", output_file, "\n")
writeVcf(output_vcf, output_file)

cat("\nDone! Phased VCF created successfully.\n")
cat("\n=== Genotype Summary ===\n")
print(table(final_gt_values))
cat("\nBreakdown:\n")
cat("1|0 = Present in haplotype 1 only:", sum(final_gt_values == "1|0"), "\n")
cat("0|1 = Present in haplotype 2 only:", sum(final_gt_values == "0|1"), "\n")
cat("1|1 = Present in both haplotypes:", sum(final_gt_values == "1|1"), "\n")