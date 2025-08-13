library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(readr)
library(pheatmap)
library(ggplot2)
library(tidyr)

# Load SNV table
snv_data <- read_tsv("C:/Users/hesam/Desktop/Mutational_Pattern/Clonal/Clonal_DeSig_Input.txt")

# Convert SNV table to per-sample VCF files
output_dir <- "C:/Users/hesam/Desktop/Mutational_Pattern/Clonal/vcf_output/"
dir.create(output_dir, showWarnings = FALSE)
samples <- unique(snv_data$Tumor_Sample_Barcode)

for (sample in samples) {
  sample_data <- snv_data %>% filter(Tumor_Sample_Barcode == sample)
  vcf_df <- data.frame(
    CHROM = sample_data$chr,
    POS = sample_data$pos,
    ID = ".",
    REF = sample_data$ref,
    ALT = sample_data$alt,
    QUAL = ".",
    FILTER = "PASS",
    INFO = "."
  )
  vcf_file <- file.path(output_dir, paste0(sample, ".vcf"))
  con <- file(vcf_file, open = "wt")
  writeLines("##fileformat=VCFv4.2", con)
  writeLines("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", con)
  write.table(vcf_df, file = con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  close(con)
}

# Load per-sample VCFs as GRanges
vcf_files <- list.files(output_dir, pattern = "\\.vcf$", full.names = TRUE)
sample_names <- tools::file_path_sans_ext(basename(vcf_files))
grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)

# Generate mutation matrix
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)

# Load COSMIC v3.3.1 signatures and fit to samples
cosmic_path <- "C:/Users/hesam/Desktop/Mutational_Pattern/COSMIC_v3.3.1_SBS_GRCh37.txt"
signatures_v31 <- read.table(cosmic_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
signatures_v31 <- as.matrix(signatures_v31)
fit_res <- fit_to_signatures(mut_mat, signatures_v31)

# Save signature contributions
write.csv(fit_res$contribution, file = "C:/Users/hesam/Desktop/Mutational_Pattern/contribution_v3.1.1_Clonal.csv")

# Load contribution table for subclonal analysis
contrib <- read.csv("Subclonal/contribution_v3.1.1_Subclonal.csv", row.names = 1, check.names = FALSE)

# Normalize contributions per sample
norm_contrib <- apply(contrib, 2, function(x) if (sum(x) > 0) x / sum(x) else x)
norm_contrib <- as.data.frame(norm_contrib)

# Filter out signatures with zero total contribution
sig_totals <- rowSums(norm_contrib)
norm_contrib_filtered <- norm_contrib[sig_totals > 0, , drop = FALSE]

# Keep only signatures with max contribution above threshold
min_sig_contrib <- 0.4
sig_max <- apply(norm_contrib_filtered, 1, max)
norm_contrib_significant <- norm_contrib_filtered[sig_max >= min_sig_contrib, , drop = FALSE]

# Plot heatmap of significant signature contributions
png("Subclonal/signature_contribution_heatmap_significant_300dpi.png", 
    width = 2400, height = 1800, res = 300)
pheatmap(norm_contrib_significant, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         main = "Signature Contribution in Subclonal Extracted with MutationalSignatures (N=277)")
dev.off()
