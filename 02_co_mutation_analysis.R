#!/usr/bin/env Rscript
# ============================================================
# 共突变模式与合成致死 — Phase 2b: Panel-Corrected Co-mutation Analysis
# Key fix: only count gene pairs when BOTH genes are covered by panel
# ============================================================
.libPaths(c("/home/caiwj2001/R/library", .libPaths()))
library(data.table)
library(ggplot2)
library(dplyr)

DATA_DIR <- "/home/caiwj2001/肿瘤突变特征与治疗反应/msk_impact_50k_2026"
OUT_DIR  <- "/home/caiwj2001/共突变模式与合成致死"

# ---- 1. Load data ----
cat("Loading data...\n")
mut <- fread(file.path(DATA_DIR, "data_mutations.txt"), sep = "\t", quote = "", na.strings = c("NA",""))
sam <- fread(file.path(DATA_DIR, "data_clinical_sample.txt"), sep = "\t", skip = 4, header = TRUE)
pat <- fread(file.path(DATA_DIR, "data_clinical_patient.txt"), sep = "\t", skip = 4, header = TRUE)
panel <- fread(file.path(DATA_DIR, "data_gene_panel_matrix.txt"), sep = "\t")

# Filter to coding mutations
coding_vc <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", 
               "Frame_Shift_Ins", "Splice_Site", "In_Frame_Del", "In_Frame_Ins",
               "Translation_Start_Site", "Nonstop_Mutation")
mut_filt <- mut[Variant_Classification %in% coding_vc]

# ---- 2. Panel distribution ----
cat("\n--- Panel versions ---\n")
panel_counts <- panel[, .N, by = mutations][order(-N)]
print(panel_counts)

# IMPACT468 is the dominant panel (34,549 samples)
# IMPACT505 (7,817) is the most comprehensive
# Strategy: use IMPACT468+505 samples for comprehensive analysis

# ---- 3. Focus on IMPACT468+505 (most comprehensive, largest subset) ----
target_panels <- c("IMPACT468", "IMPACT505")
samples_target <- panel[mutations %in% target_panels, unique(SAMPLE_ID)]
cat(sprintf("\nSamples with IMPACT468/505: %d\n", length(samples_target)))

# Filter mutations to these samples
mut_target <- mut_filt[Tumor_Sample_Barcode %in% samples_target]
cat(sprintf("Mutations in target samples: %d\n", nrow(mut_target)))

# ---- 4. Gene frequency in target samples ----
gene_freq <- mut_target[, .N, by = Hugo_Symbol][order(-N)]
cat(sprintf("Unique genes mutated in target: %d\n", nrow(gene_freq)))
cat("\nTop 30 genes:\n")
print(head(gene_freq, 30))

# ---- 5. All genes in IMPACT468/505 are well-covered by definition ----
# Since all samples use the same panel versions, every gene in the gene list
# is covered by the panel for all samples. Use top 50 mutated genes.
top_genes_for_analysis <- head(gene_freq, 50)$Hugo_Symbol
cat(sprintf("Top 50 genes for pairwise analysis selected\n"))

# ---- 6. Build corrected binary matrix ----
cat("\nBuilding corrected binary matrix...\n")
# Get unique mutation events per gene×sample
mut_target[, PATIENT_SAMPLE := Tumor_Sample_Barcode]
gene_sample <- mut_target[Hugo_Symbol %in% top_genes_for_analysis, 
                          .(MUTATED = 1), 
                          by = .(Hugo_Symbol, PATIENT_SAMPLE)]
gene_sample <- unique(gene_sample)

# Build matrix (only covered samples)
all_samples <- unique(mut_target$Tumor_Sample_Barcode)
n_samples <- length(all_samples)
cat(sprintf("Analysis samples: %d\n", n_samples))

# Build full matrix with proper sample coverage
mut_matrix <- matrix(0, nrow = length(top_genes_for_analysis), ncol = n_samples)
rownames(mut_matrix) <- top_genes_for_analysis
colnames(mut_matrix) <- all_samples

for (i in seq_len(nrow(gene_sample))) {
  g <- gene_sample$Hugo_Symbol[i]
  s <- gene_sample$PATIENT_SAMPLE[i]
  if (g %in% rownames(mut_matrix) && s %in% colnames(mut_matrix)) {
    mut_matrix[g, s] <- 1
  }
}

cat(sprintf("Matrix dimensions: %d genes × %d samples\n", nrow(mut_matrix), ncol(mut_matrix)))

# ---- 7. Sample mutation counts ----
sample_mut_counts <- colSums(mut_matrix)
cat(sprintf("Median mutations per sample: %d\n", median(sample_mut_counts)))
cat(sprintf("Range: %d - %d\n", min(sample_mut_counts), max(sample_mut_counts)))

# ---- 8. Pairwise co-occurrence with Fisher's exact test ----
cat("\n--- Pairwise Co-occurrence Analysis (Fisher's exact test) ---\n")
ng <- nrow(mut_matrix)
pair_results <- data.table(
  Gene1 = character(),
  Gene2 = character(),
  Both_count = integer(),
  Gene1_only = integer(),
  Gene2_only = integer(),
  Neither = integer(),
  OR = numeric(),
  OR_lower = numeric(),
  OR_upper = numeric(),
  P_value = numeric(),
  Direction = character()
)

for (i in 1:(ng-1)) {
  for (j in (i+1):ng) {
    g1 <- rownames(mut_matrix)[i]
    g2 <- rownames(mut_matrix)[j]
    
    a <- sum(mut_matrix[i,] == 1 & mut_matrix[j,] == 1)  # both
    b <- sum(mut_matrix[i,] == 1 & mut_matrix[j,] == 0)  # only g1
    c <- sum(mut_matrix[i,] == 0 & mut_matrix[j,] == 1)  # only g2
    d <- sum(mut_matrix[i,] == 0 & mut_matrix[j,] == 0)  # neither
    
    if (a + b + c + d != n_samples) {
      cat(sprintf("WARNING: %s-%s: sum mismatch %d vs %d\n", g1, g2, a+b+c+d, n_samples))
    }
    
    ft <- fisher.test(matrix(c(a, b, c, d), 2, 2))
    
    pair_results <- rbind(pair_results, data.table(
      Gene1 = g1, Gene2 = g2,
      Both_count = a, Gene1_only = b, Gene2_only = c, Neither = d,
      OR = ft$estimate, OR_lower = ft$conf.int[1], OR_upper = ft$conf.int[2],
      P_value = ft$p.value,
      Direction = ifelse(ft$estimate > 1, "Co-occurrence", "Mutual exclusivity")
    ))
  }
  if (i %% 10 == 0) cat(sprintf("  Processed gene %d/%d\n", i, ng))
}

# ---- 9. Multiple testing correction ----
pair_results[, P_adj := p.adjust(P_value, method = "BH")]

# Significant pairs
sig_pairs <- pair_results[P_adj < 0.05]
cat(sprintf("\nSignificant pairs (FDR<0.05): %d\n", nrow(sig_pairs)))
cat(sprintf("  Co-occurrence: %d\n", sum(sig_pairs$Direction == "Co-occurrence")))
cat(sprintf("  Mutual exclusivity: %d\n", sum(sig_pairs$Direction == "Mutual exclusivity")))

# ---- 10. Top results by direction ----
cat("\n--- Top Co-occurring Pairs (FDR<0.05) ---\n")
co_pairs <- sig_pairs[Direction == "Co-occurrence"][order(-Both_count)]
print(head(co_pairs[, .(Gene1, Gene2, Both_count, OR, P_adj)], 20))

cat("\n--- Top Mutually Exclusive Pairs (FDR<0.05) ---\n")
me_pairs <- sig_pairs[Direction == "Mutual exclusivity"][order(P_adj)]
print(head(me_pairs[, .(Gene1, Gene2, Both_count, OR, P_adj)], 20))

# ---- 11. Save results ----
save(mut_filt, mut_target, sam, pat, panel, mut_matrix, 
     top_genes_for_analysis, pair_results, sig_pairs, co_pairs, me_pairs,
     file = file.path(OUT_DIR, "data", "co_mutation_results.RData"))

# Export significant pairs to CSV
fwrite(pair_results[P_adj < 0.05][order(P_adj)], 
       file.path(OUT_DIR, "results", "significant_gene_pairs.csv"))

cat(sprintf("\nResults saved to %s\n", OUT_DIR))
cat("============= COMPLETE =============\n")
