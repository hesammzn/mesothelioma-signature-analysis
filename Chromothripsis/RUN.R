library(data.table)
library(GenomicRanges)
library(mclust)
library(hdp)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)

# loading CNV data
cnv <- fread("Sig_CNV/Clonal/Input_Clonal_CNV.csv")
cnv <- cnv[(end - start) > 50000]
cnv <- cnv[!chromosome %in% c("X", "Y", "chrX", "chrY")]
cnv[, total_cn := nMajor + nMinor]
cnv[, seg_size := end - start]
cnv[, sample := as.character(sample)]

# removing Ig loci
cnv[, row_id := .I]
igh <- cnv[chromosome == 14 & start > 106032614 & end < 108288051]
igk <- cnv[chromosome == 2  & start > 87090568  & end < 93274235]
igl <- cnv[chromosome == 22 & start > 21080474 & end < 26065085]
cnv <- cnv[!row_id %in% rbind(igh, igk, igl)$row_id]

# loading centromere positions (hg19)
centro <- fread("CentromerePosition_hg19.txt")[, .(chrom, start = start, end = end)]
centro$chrom <- gsub("chr", "", centro$chrom)

# defining thresholds (Maclachlan et al., Table S2)
thr_MB    <- c(3, 6, 31)
thr_COUNT <- c(0, 1, 2, 3, 9)
thr_JUMP  <- c(1, 3, 8)
thr_OSCI  <- c(1, 4, 9, 38)
thr_SIZE  <- c(1.65e5, 5.02e5, 1.588e6, 7.113e6, 2.18e7, 5.396e7, 6.758e7, 1.053e8, 1.4223e8, 2.491e8)

# creating binning helper
bin_segment <- function(values, thresholds) {
  cut(values, breaks = c(-Inf, thresholds, Inf), labels = FALSE, right = TRUE)
}

# initializing outputs
count_list <- list()
band_count_list <- list()
samples <- unique(cnv$sample)

# iterating samples and computing features
for (s in samples) {
  dt <- cnv[sample == s][order(chromosome, start)]

  # binning SIZE
  dt[, SIZE_cat := bin_segment(seg_size, thr_SIZE)]

  # binning COUNT
  dt[, COUNT_cat := bin_segment(total_cn, thr_COUNT)]

  # binning JUMP
  dt[, prev := data.table::shift(total_cn)]
  dt[, jump := abs(total_cn - prev)]
  dt[is.na(jump), jump := 0]
  dt[, JUMP_cat := bin_segment(jump, thr_JUMP)]

  # counting MB (breakpoints per 10Mb)
  gr <- GRanges(dt$chromosome, IRanges(start = dt$end, end = dt$end))
  chrmax <- dt[, .(max_end = max(end)), by = chromosome]
  seqlens <- setNames(chrmax$max_end, chrmax$chromosome)
  wins <- tileGenome(seqlens, tilewidth = 1e7, cut.last.tile.in.chrom = TRUE)
  ol <- findOverlaps(gr, wins)
  dt[, MB_bin := NA_integer_]
  dt[queryHits(ol), MB_bin := subjectHits(ol)]
  mb_table <- as.data.table(table(dt$MB_bin))
  mb_table[, MB_cat := bin_segment(as.numeric(N), thr_MB)]
  mb_count <- table(mb_table$MB_cat)

  # counting BAND (breakpoints per arm)
  centro_clean <- centro[, .(chr = chrom, cent_start = start, cent_end = end)]
  common_chr <- intersect(unique(dt$chromosome), centro_clean$chr)
  dt <- dt[chromosome %in% common_chr]
  centro_clean <- centro_clean[chr %in% common_chr]
  p_arms <- centro_clean[, .(chr, start = 0, end = cent_start - 1)]
  q_arms <- centro_clean[, .(chr, start = cent_end + 1, end = 2.5e8)]
  gr_p <- GRanges(seqnames = p_arms$chr, ranges = IRanges(start = p_arms$start, end = p_arms$end))
  gr_q <- GRanges(seqnames = q_arms$chr, ranges = IRanges(start = q_arms$start, end = q_arms$end))
  gr_breaks <- GRanges(seqnames = dt$chromosome, ranges = IRanges(start = dt$end, end = dt$end))
  n_p <- length(unique(queryHits(findOverlaps(gr_breaks, gr_p))))
  n_q <- length(unique(queryHits(findOverlaps(gr_breaks, gr_q))))
  band_count_list[[s]] <- data.table(sample = s, band_count = n_p + n_q)

  # computing oscillation chain length
  cnvvec <- round(dt$total_cn)
  osc_counts <- c()
  if (length(cnvvec) >= 4) {
    cnt <- 0
    for (j in 3:length(cnvvec)) {
      if (j == length(cnvvec)) {
        osc_counts <- c(osc_counts, cnt); cnt <- 0
      } else {
        if (abs(cnvvec[j] - cnvvec[j-2]) <= 1 && cnvvec[j] != cnvvec[j-1]) {
          cnt <- cnt + 1
        } else {
          osc_counts <- c(osc_counts, cnt); cnt <- 0
        }
      }
    }
  }
  if (length(osc_counts) > 0) {
    osc_cat_vec <- bin_segment(osc_counts, thr_OSCI)
    tab_osci <- table(osc_cat_vec)
  } else {
    tab_osci <- integer()
  }

  # aggregating per-sample counts
  out <- data.table(sample = s)
  for (f in c("SIZE", "COUNT", "JUMP")) {
    tab <- table(dt[[paste0(f, "_cat")]])
    for (i in 1:max(as.numeric(names(tab)))) out[[paste0(f, "_", i)]] <- tab[as.character(i)]
  }
  for (i in 1:max(as.numeric(names(mb_count)))) out[[paste0("MB_", i)]] <- mb_count[as.character(i)]
  for (i in as.integer(names(tab_osci))) out[[paste0("OSCI_", i)]] <- as.integer(tab_osci[as.character(i)])
  for (col in names(out)) if (is.numeric(out[[col]])) out[[col]][is.na(out[[col]])] <- 0
  count_list[[s]] <- out
}

# combining samples
final <- rbindlist(count_list, fill = TRUE)
final[is.na(final)] <- 0

# clustering BAND with mclust
band_counts <- rbindlist(band_count_list)
band_counts <- band_counts[band_count > 0]
mclust_band <- Mclust(band_counts$band_count, G = 3, verbose = FALSE)
band_counts[, band_cat := mclust_band$classification]
band_counts[, band_feature := paste0("BAND_", band_cat)]
band_table <- dcast(band_counts, sample ~ band_feature, value.var = "band_count", fun.aggregate = sum)
for (col in paste0("BAND_", 1:3)) if (!col %in% colnames(band_table)) band_table[[col]] <- 0

# merging BAND into matrix
final <- merge(final, band_table, by = "sample", all.x = TRUE)
for (col in paste0("BAND_", 1:3)) final[[col]][is.na(final[[col]])] <- 0

# enforcing 28-feature set
expected_cols <- c(paste0("MB_", 1:3),
                   paste0("COUNT_", 1:5),
                   paste0("JUMP_", 1:3),
                   paste0("BAND_", 1:3),
                   paste0("OSCI_", 1:4),
                   paste0("SIZE_", 1:10))
for (col in expected_cols) if (!col %in% colnames(final)) final[[col]] <- 0
final <- final[, c("sample", expected_cols), with = FALSE]

# saving feature matrix
fwrite(final, "CN_feature28_Clonal.csv")
message("\u2705 CN signature count matrix with improved BAND saved.")

# loading feature matrix for HDP
df <- fread("CN_feature28_Clonal.csv")
rownames(df) <- df$sample
df[, sample := NULL]
mat_final <- as.matrix(apply(as.matrix(df), 2, as.numeric))
class(mat_final) <- "matrix"

# initializing HDP
n_cat <- ncol(mat_final)
hdp <- hdp_init(ppindex = 0, cpindex = 1, hh = rep(1 / n_cat, n_cat), alphaa = 1, alphab = 1)

# adding DPs per sample
hdp <- hdp_adddp(hdp, numdp = nrow(mat_final), ppindex = rep(1, nrow(mat_final)), cpindex = rep(1, nrow(mat_final)))

# setting data
hdp <- hdp_setdata(hdp, dpindex = 2:(nrow(mat_final) + 1), data = as.data.frame(mat_final))

# activating DPs
hdp <- dp_activate(hdp, dpindex = 1:(nrow(mat_final) + 1), initcc = 10, seed = 1234)

# running posterior (6 chains)
chlist <- vector("list", 6)
for (i in 1:6) {
  chlist[[i]] <- hdp_posterior(hdp, burnin = 20000, n = 100, space = 50, cpiter = 3, seed = i * 10000)
}

# combining chains
multi_chain <- hdp_multi_chain(chlist)
saveRDS(multi_chain, "mut_example_multi_CNV_Clonal.rds")

# creating output dir
if (!dir.exists("Clonal_CN_Figures")) dir.create("Clonal_CN_Figures")

# plotting CN category heatmap
hdp_final <- fread("CN_feature28_Clonal.csv")
rownames(hdp_final) <- hdp_final$sample
hdp_final2 <- hdp_final[, -1]
pheatmap(t(hdp_final2), show_colnames = FALSE)
png(file.path("Clonal_CN_Figures", "CNV_Category_Heatmap.png"), width = 1800, height = 1800, res = 300)
pheatmap(t(hdp_final2), show_colnames = FALSE)
dev.off()

# defining colors for category barplots
mat_final <- hdp_final2
n <- 30
gg_color_hue <- function(n) { hues = seq(15, 375, length = n + 1); hcl(h = hues, l = 65, c = 100)[1:n] }
cols <- gg_color_hue(n)
cnv_colors <- cols[c(1, 5, 10, 14, 21, 24)]
cnv_colors_final <- c(rep(cnv_colors[1], 3),
                      rep(cnv_colors[2], 5),
                      rep(cnv_colors[3], 3),
                      rep(cnv_colors[4], 3),
                      rep(cnv_colors[5], 4),
                      rep(cnv_colors[6], 10))
channel_names2 <- c("1","2","3","1","2","3","4","5","1","2","3","1","2","3",
                    "1","2","3","4","1","2","3","4","5","6","7","8","9","10")

# plotting CN category totals
x <- barplot(colSums(mat_final), las = 2, col = cnv_colors_final, border = NA, xaxt = "n", cex.axis = 1.15)
axis(1, at = x, label = rep("", 28), mgp = c(3, 1, 0.2))
mtext(1, at = x, text = channel_names2, col = cnv_colors_final, padj = 1.5, cex = 1.15)
png(file.path("Clonal_CN_Figures", "CNV_Category_Barplot.png"), width = 1800, height = 1000, res = 300)
x <- barplot(colSums(mat_final), las = 2, col = cnv_colors_final, border = NA, xaxt = "n", cex.axis = 1.15)
axis(1, at = x, label = rep("", 28), mgp = c(3, 1, 0.2))
mtext(1, at = x, text = channel_names2, col = cnv_colors_final, padj = 1.5, cex = 1.15)
dev.off()

# extracting components
mut_example_multi <- readRDS("mut_example_multi_CNV_Clonal.rds")
mut_example_multi_0.85_10 <- hdp_extract_components(mut_example_multi, cos.merge = 0.85, min.sample = 10)
mut_example_multi <- mut_example_multi_0.85_10

# plotting component sizes
par(mfrow = c(1,1), mar = c(5,4,4,2))
plot_comp_size(mut_example_multi, bty = "L")
png(file.path("Clonal_CN_Figures", "Signature_Component_Sizes.png"), width = 1600, height = 1200, res = 300)
par(mfrow = c(1,1), mar = c(5,4,4,2))
plot_comp_size(mut_example_multi, bty = "L")
dev.off()

# exporting signature profiles
mut_example_multi_plot <- mut_example_multi
posteriorMeans_plot <- t(comp_categ_distn(mut_example_multi_plot)[[1]])
rownames(posteriorMeans_plot) <- channel_names2
plotnames <- c("CN-SIG1", "CN-SIG2", "CN-SIG3", "CN-SIG4")

for (i in 1:4) {
  x <- barplot(posteriorMeans_plot[, i], las = 2, col = cnv_colors_final, border = NA, xaxt = "n",
               cex.axis = 1.5, main = plotnames[i], ylim = c(0, 0.25), cex.main = 2)
  axis(1, at = x, label = rep("", 28), mgp = c(3, 1, 0.2))
  mtext(1, at = x, text = channel_names2, col = cnv_colors_final, padj = 1.5, cex = 1.5)
}
rownames(posteriorMeans_plot) <- colnames(hdp_final2)
write.csv(posteriorMeans_plot, file = "CN_Signature_Profiles_Clonal.csv", row.names = TRUE)

# exporting sample-by-signature exposures
sample_signature <- comp_dp_distn(mut_example_multi)[[1]]
sample_signature <- sample_signature[-1, ]
df_ids <- fread("CN_feature28_Clonal.csv")
rownames(sample_signature) <- df_ids$sample
colnames(sample_signature) <- plotnames
write.csv(sample_signature, file = "Sample_by_CN_Signature_Exposure_Clonal.csv", row.names = TRUE)

# plotting CN-SIG4 exposure per sample
df_sig4 <- data.frame(Sample = rownames(sample_signature), CN_SIG4 = sample_signature[, "CN-SIG4"])
ggplot(df_sig4, aes(x = reorder(Sample, -CN_SIG4), y = CN_SIG4)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() + ylab("CN-SIG4 Exposure") + xlab("Sample") +
  ggtitle("CN-SIG4 Exposure Across Samples") + theme_minimal()
png(file.path("Clonal_CN_Figures", "CN_SIG4_Exposure_Across_Samples.png"), width = 2200, height = 2200, res = 300)
print(
  ggplot(df_sig4, aes(x = reorder(Sample, -CN_SIG4), y = CN_SIG4)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() + ylab("CN-SIG4 Exposure") + xlab("Sample") +
    ggtitle("CN-SIG4 Exposure Across Samples") +
    theme_minimal() +
    theme(axis.text.y = element_blank())
)
dev.off()

# plotting average signature exposure
exposures <- sample_signature
signature_means <- colMeans(exposures)
barplot(signature_means, col = rainbow(ncol(exposures)), las = 2, main = "Average Signature Exposure")
png(file.path("Clonal_CN_Figures", "Average_Signature_Exposure.png"), width = 1200, height = 900, res = 300)
barplot(signature_means, col = rainbow(ncol(exposures)), las = 2, main = "Average Signature Exposure")
dev.off()

# composing signature barplots with legend
png("CN_Signature_Barplots_clean_Clonal.png", width = 3800, height = 1450, res = 300)
par(oma = c(2, 0, 5, 0))
layout(matrix(c(1,2,3,4,5,5,5,5), nrow = 2, byrow = TRUE), heights = c(4, 1))
par(mar = c(4,4,4,1))
for (i in 1:4) {
  x <- barplot(posteriorMeans_plot[, i], las = 2, col = cnv_colors_final, border = NA, xaxt = "n",
               cex.axis = 1.2, main = plotnames[i], ylim = c(0, 0.25), cex.main = 2)
  axis(1, at = x, label = rep("", 28), mgp = c(3, 1, 0.2))
  mtext(1, at = x, text = channel_names2, col = cnv_colors_final, padj = 1, cex = 0.8)
}
mtext("De Novo CN signatures in Clonal segments", side = 3, outer = TRUE, line = 1.8, cex = 2, font = 2)
legend_colors <- cnv_colors_final[c(1,4,9,12,16,20)]
par(mar = c(0,0,0,0)); plot.new()
legend(x = 0.13, y = 0.7,
       legend = c("CN breakpoints per 10 Mb",
                  "absolute CN of segments",
                  "difference in CN between adjacent segments"),
       fill = legend_colors[1:3], border = NA, cex = 1.35, bty = "n", xpd = NA, text.font = 1)
legend(x = 0.55, y = 0.7,
       legend = c("CN breakpoints per chromosome arm",
                  "lengths of oscillating CN segment chains",
                  "size of CN segments"),
       fill = legend_colors[4:6], border = NA, cex = 1.35, bty = "n", xpd = NA, text.font = 1)
dev.off()
