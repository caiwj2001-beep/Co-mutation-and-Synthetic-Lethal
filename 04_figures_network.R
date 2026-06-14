#!/usr/bin/env Rscript
# ============================================================
# 共突变模式与合成致死 — Phase 4: Network & Figure Generation
# Co-mutation network visualization + publication figures
# ============================================================
.libPaths(c("/home/caiwj2001/R/library", .libPaths()))
library(data.table)
library(ggplot2)
library(igraph)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(ggrepel)
library(scales)

DATA_DIR <- "/home/caiwj2001/肿瘤突变特征与治疗反应/msk_impact_50k_2026"
OUT_DIR  <- "/home/caiwj2001/共突变模式与合成致死"

# ---- 1. Load data ----
cat("Loading data...\n")
mut <- fread(file.path(DATA_DIR, "data_mutations.txt"), sep = "\t", quote = "", na.strings = c("NA",""))
sam <- fread(file.path(DATA_DIR, "data_clinical_sample.txt"), sep = "\t", skip = 4, header = TRUE)
pat <- fread(file.path(DATA_DIR, "data_clinical_patient.txt"), sep = "\t", skip = 4, header = TRUE)

# Filter to coding mutations
coding_vc <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", 
               "Frame_Shift_Ins", "Splice_Site", "In_Frame_Del", "In_Frame_Ins",
               "Translation_Start_Site", "Nonstop_Mutation")
mut_filt <- mut[Variant_Classification %in% coding_vc]

# ---- 2. Cancer-type-stratified mutation frequencies ----
cat("\n--- FIGURE 1: Pan-cancer mutation landscape ---\n")

# Top genes across all cancers
gene_freq <- mut_filt[, .N, by = Hugo_Symbol][order(-N)]
top25 <- head(gene_freq, 25)$Hugo_Symbol

# Cancer-type-specific mutation frequency for top genes
top_cancer_types <- sam[, .N, by = CANCER_TYPE][order(-N)][1:15, CANCER_TYPE]

# Build matrix: gene × cancer type
# Map samples to cancer types
mut_filt[, PATIENT := gsub("-T[0-9]+-.*$", "", Tumor_Sample_Barcode)]
sam_patient <- unique(sam[, .(PATIENT_ID, CANCER_TYPE)])

ct_matrix <- matrix(0, nrow = length(top25), ncol = length(top_cancer_types))
rownames(ct_matrix) <- top25
colnames(ct_matrix) <- top_cancer_types

for (i in seq_along(top_cancer_types)) {
  ct <- top_cancer_types[i]
  ct_patients <- sam_patient[CANCER_TYPE == ct, unique(PATIENT_ID)]
  ct_mut <- mut_filt[Hugo_Symbol %in% top25 & PATIENT %in% ct_patients]
  ct_gene_freq <- ct_mut[, .(N_MUT = uniqueN(PATIENT)), by = Hugo_Symbol]
  ct_gene_freq[, FREQ := 100 * N_MUT / length(ct_patients)]
  
  for (j in seq_along(top25)) {
    g <- top25[j]
    if (g %in% ct_gene_freq$Hugo_Symbol) {
      ct_matrix[j, i] <- ct_gene_freq[Hugo_Symbol == g, FREQ]
    }
  }
}

# Figure 1: Pan-cancer mutation landscape heatmap
png(file.path(OUT_DIR, "figures", "Fig1_PanCancer_Mutation_Landscape.png"),
    width = 2000, height = 2400, res = 300)
pheatmap(ct_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(brewer.pal(9, "Reds"))(100),
         main = "Pan-cancer Mutation Landscape (MSK-IMPACT 50K)",
         fontsize_row = 10,
         fontsize_col = 10,
         display_numbers = FALSE,
         angle_col = 45,
         legend_labels = "Mutation\nFrequency (%)")
dev.off()

tiff(file.path(OUT_DIR, "figures", "Fig1_PanCancer_Mutation_Landscape.tiff"),
     width = 2000, height = 2400, res = 300, compression = "lzw")
pheatmap(ct_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(brewer.pal(9, "Reds"))(100),
         main = "Pan-cancer Mutation Landscape (MSK-IMPACT 50K)",
         fontsize_row = 10,
         fontsize_col = 10,
         display_numbers = FALSE,
         angle_col = 45)
dev.off()

cat("  [OK] Fig1: Pan-cancer mutation landscape\n")

# ---- 3. FIGURE 2: Gene-wise mutation frequency bar plot ----
cat("\n--- FIGURE 2: Top mutated genes ---\n")
gene_freq_top <- head(gene_freq, 30)
gene_freq_top[, Hugo_Symbol := factor(Hugo_Symbol, levels = rev(Hugo_Symbol))]
gene_freq_top[, Pct := 100 * N / length(unique(mut_filt$Tumor_Sample_Barcode))]

p2 <- ggplot(gene_freq_top, aes(x = N, y = Hugo_Symbol)) +
  geom_bar(stat = "identity", fill = "#2166AC", width = 0.7) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", N, Pct)), hjust = -0.1, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "Top 30 Frequently Mutated Genes in Pan-cancer",
       subtitle = sprintf("MSK-IMPACT 50K Cohort (n=%d samples)", 
                          length(unique(mut_filt$Tumor_Sample_Barcode))),
       x = "Number of mutations (coding)", y = "") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank())

ggsave(file.path(OUT_DIR, "figures", "Fig2_Top_Mutated_Genes.png"),
       p2, width = 10, height = 8, dpi = 300)
ggsave(file.path(OUT_DIR, "figures", "Fig2_Top_Mutated_Genes.tiff"),
       p2, width = 10, height = 8, dpi = 600, compression = "lzw")
cat("  [OK] Fig2: Top mutated genes\n")

# ---- 4. FIGURE 3: Co-mutation network ----
cat("\n--- FIGURE 3: Co-mutation Network ---\n")

# Build co-mutation matrix for top 30 genes (simple count-based)
top30 <- head(gene_freq, 30)$Hugo_Symbol

# Get patient×gene binary matrix
mut_target <- mut_filt[Hugo_Symbol %in% top30]
gene_sample <- mut_target[, .(MUTATED = 1), by = .(Hugo_Symbol, PATIENT = gsub("-T[0-9]+-.*$", "", Tumor_Sample_Barcode))]
gene_sample <- unique(gene_sample)

patients_all <- unique(gene_sample$PATIENT)
mat30 <- matrix(0, nrow = length(top30), ncol = length(patients_all))
rownames(mat30) <- top30
colnames(mat30) <- patients_all

for (i in seq_len(nrow(gene_sample))) {
  g <- gene_sample$Hugo_Symbol[i]
  p <- gene_sample$PATIENT[i]
  if (g %in% rownames(mat30)) mat30[g, p] <- 1
}

# Pairwise co-occurrence count matrix
co_counts <- mat30 %*% t(mat30)
diag(co_counts) <- 0

# Also compute Jaccard similarity
jaccard <- matrix(0, nrow = 30, ncol = 30)
rownames(jaccard) <- colnames(jaccard) <- top30
for (i in 1:30) {
  for (j in 1:30) {
    if (i >= j) next
    a <- co_counts[i, j]
    b <- sum(mat30[i,]) - a
    c <- sum(mat30[j,]) - a
    jaccard[i, j] <- a / (a + b + c)
    jaccard[j, i] <- jaccard[i, j]
  }
}

# Create network graph
# Use significant edges only (Jaccard > threshold)
threshold <- quantile(jaccard[upper.tri(jaccard)], 0.95)
cat(sprintf("  Jaccard similarity threshold (95th percentile): %.4f\n", threshold))

# Build adjacency matrix
adj <- jaccard
adj[jaccard < threshold] <- 0

# Create igraph
g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE)

# Node sizes proportional to mutation frequency
node_sizes <- gene_freq[match(V(g)$name, Hugo_Symbol), N]
node_sizes <- log10(node_sizes) * 5

# Node colors by functional category
func_cats <- list(
  "Tumor Suppressor" = c("TP53", "APC", "PTEN", "RB1", "CDKN2A", "NF1", "SMAD4", "VHL", "STK11"),
  "Oncogene" = c("KRAS", "PIK3CA", "BRAF", "EGFR", "CTNNB1", "NRAS", "IDH1", "ERBB2"),
  "Chromatin Remodeling" = c("ARID1A", "KMT2D", "KMT2C", "KMT2B", "KMT2A", 
                               "SMARCA4", "CREBBP", "SETD2", "ATRX", "SPEN"),
  "DNA Repair" = c("ATM", "BRCA2", "BRCA1", "MSH6", "MSH2", "BRIP1")
)

node_colors <- rep("grey70", length(V(g)$name))
names(node_colors) <- V(g)$name
for (cat_name in names(func_cats)) {
  genes_in_cat <- func_cats[[cat_name]]
  node_colors[names(node_colors) %in% genes_in_cat] <- switch(cat_name,
    "Tumor Suppressor" = "#E41A1C",
    "Oncogene" = "#377EB8",
    "Chromatin Remodeling" = "#4DAF4A",
    "DNA Repair" = "#FF7F00"
  )
}

# Edge widths
edge_widths <- E(g)$weight * 10

# Layout
set.seed(42)
layout <- layout_with_fr(g)

# Plot network
png(file.path(OUT_DIR, "figures", "Fig3_CoMutation_Network.png"),
    width = 2400, height = 2400, res = 300)
par(mar = c(0, 0, 0, 0))
plot(g,
     layout = layout,
     vertex.size = node_sizes,
     vertex.color = node_colors,
     vertex.label = V(g)$name,
     vertex.label.cex = 0.8,
     vertex.label.color = "black",
     vertex.frame.color = "grey50",
     edge.width = edge_widths,
     edge.color = adjustcolor("grey50", alpha.f = 0.5),
     main = "Pan-cancer Co-mutation Network")
legend("topleft",
       legend = names(func_cats),
       col = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00"),
       pch = 19, pt.cex = 2, bty = "n", cex = 0.9)
dev.off()

tiff(file.path(OUT_DIR, "figures", "Fig3_CoMutation_Network.tiff"),
     width = 2400, height = 2400, res = 300, compression = "lzw")
par(mar = c(0, 0, 0, 0))
plot(g,
     layout = layout,
     vertex.size = node_sizes,
     vertex.color = node_colors,
     vertex.label = V(g)$name,
     vertex.label.cex = 0.8,
     vertex.label.color = "black",
     vertex.frame.color = "grey50",
     edge.width = edge_widths,
     edge.color = adjustcolor("grey50", alpha.f = 0.5),
     main = "Pan-cancer Co-mutation Network")
legend("topleft",
       legend = names(func_cats),
       col = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00"),
       pch = 19, pt.cex = 2, bty = "n", cex = 0.9)
dev.off()

cat("  [OK] Fig3: Co-mutation network\n")

# ---- 5. FIGURE 4: Mutual exclusivity vs co-occurrence summary ----
# (This will be enhanced when co_mutation_results.RData is available)
cat("\n--- FIGURE 4: Will be generated from co_mutation_results.RData ---\n")
cat("  (Requires Phase 2b results to complete)\n")

# ---- 6. Extended Data: Cancer-type-specific co-mutation heatmaps ----
cat("\n--- Extended Data: Cancer-type-specific co-mutation ---\n")

# Focus on top 5 cancer types for detailed analysis
major_ct <- head(top_cancer_types, 5)
for (ct in major_ct) {
  ct_patients <- sam_patient[CANCER_TYPE == ct, unique(PATIENT_ID)]
  ct_mut <- mut_filt[Hugo_Symbol %in% top25 & PATIENT %in% ct_patients]
  
  if (length(ct_patients) >= 100) {
    # Build mutation matrix for this cancer type
    ct_gs <- ct_mut[, .(MUTATED = 1), by = .(Hugo_Symbol, PATIENT)]
    ct_gs <- unique(ct_gs)
    
    ct_mat <- matrix(0, nrow = length(top25), ncol = length(ct_patients))
    rownames(ct_mat) <- top25
    colnames(ct_mat) <- ct_patients
    
    for (i in seq_len(nrow(ct_gs))) {
      g <- ct_gs$Hugo_Symbol[i]
      p <- ct_gs$PATIENT[i]
      if (g %in% rownames(ct_mat)) ct_mat[g, p] <- 1
    }
    
    # Keep only genes with >3% mutation rate
    gene_rates <- rowMeans(ct_mat) * 100
    genes_keep <- names(gene_rates[gene_rates >= 3])
    
    if (length(genes_keep) >= 5) {
      ct_mat_plot <- ct_mat[genes_keep, , drop = FALSE]
      
      # Sort samples by mutation count
      sample_order <- order(-colSums(ct_mat_plot))
      ct_mat_plot <- ct_mat_plot[, sample_order, drop = FALSE]
      
      # Co-mutation heatmap (showing top 100 samples for clarity)
      n_show <- min(100, ncol(ct_mat_plot))
      
      png(file.path(OUT_DIR, "figures", 
                    sprintf("ED_CoMutation_%s.png", gsub("[ /]", "_", ct))),
          width = 1600, height = 1200, res = 200)
      pheatmap(ct_mat_plot[, 1:n_show, drop = FALSE],
               cluster_rows = TRUE, cluster_cols = FALSE,
               color = c("white", "#2166AC"),
               main = sprintf("Co-mutation Patterns in %s", ct),
               fontsize_row = 8, fontsize_col = 4,
               legend = FALSE, show_colnames = FALSE)
      dev.off()
      
      cat(sprintf("  [OK] Extended Data: %s co-mutation heatmap\n", ct))
    }
  }
}

cat("\n============= FIGURE GENERATION COMPLETE =============\n")
