#!/usr/bin/env Rscript
# ============================================================
# 共突变模式与合成致死 — Phase 3: Survival Analysis
# Co-mutation × Prognosis: KM curves, Cox regression
# ============================================================
.libPaths(c("/home/caiwj2001/R/library", .libPaths()))
library(data.table)
library(ggplot2)
library(survival)
library(dplyr)

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

# ---- 2. Prepare survival data ----
# Merge sample→patient info
pat_surv <- pat[, .(PATIENT_ID, SEX, OS_STATUS, OS_MONTHS, AGE_AT_DX)]
pat_surv[, OS_MONTHS := as.numeric(OS_MONTHS)]
pat_surv[, AGE := as.numeric(AGE_AT_DX)]
pat_surv[, EVENT := fifelse(OS_STATUS == "DECEASED", 1, 0)]
pat_surv <- pat_surv[!is.na(OS_MONTHS) & OS_MONTHS > 0 & !is.na(EVENT)]

cat(sprintf("Patients with valid survival: %d\n", nrow(pat_surv)))
cat(sprintf("  Events: %d (%.1f%%)\n", sum(pat_surv$EVENT), 100*mean(pat_surv$EVENT)))

# ---- 3. Define co-mutation pairs to analyze ----
# Key pairs based on literature and biological relevance
co_mutation_pairs <- list(
  list(genes = c("TP53", "KRAS"), name = "TP53-KRAS"),
  list(genes = c("KRAS", "STK11"), name = "KRAS-STK11"),
  list(genes = c("KRAS", "KEAP1"), name = "KRAS-KEAP1"),
  list(genes = c("TP53", "PIK3CA"), name = "TP53-PIK3CA"),
  list(genes = c("ARID1A", "PIK3CA"), name = "ARID1A-PIK3CA"),
  list(genes = c("TP53", "ATM"), name = "TP53-ATM"),
  list(genes = c("APC", "KRAS"), name = "APC-KRAS"),
  list(genes = c("APC", "TP53"), name = "APC-TP53"),
  list(genes = c("TP53", "RB1"), name = "TP53-RB1"),
  list(genes = c("PIK3CA", "PTEN"), name = "PIK3CA-PTEN")
)

# ---- 4. Build mutation status per patient ----
# For each gene of interest, determine if patient has mutation
genes_of_interest <- unique(unlist(lapply(co_mutation_pairs, `[[`, "genes")))
cat(sprintf("Genes of interest: %s\n", paste(genes_of_interest, collapse=", ")))

# Get patient-level mutation status
mut_per_patient <- mut_filt[Hugo_Symbol %in% genes_of_interest, 
                             .(MUTATED = 1), 
                             by = .(PATIENT = gsub("-T[0-9]+-.*$", "", Tumor_Sample_Barcode),
                                    Hugo_Symbol)]
mut_per_patient <- unique(mut_per_patient)

# ---- 5. For each co-mutation pair, run survival analysis (pan-cancer) ----
cat("\n============= SURVIVAL ANALYSIS =============\n")

surv_results <- list()

for (pair in co_mutation_pairs) {
  g1 <- pair$genes[1]
  g2 <- pair$genes[2]
  pname <- pair$name
  
  # Get mutation status for each gene
  mut_g1 <- mut_per_patient[Hugo_Symbol == g1, unique(PATIENT)]
  mut_g2 <- mut_per_patient[Hugo_Symbol == g2, unique(PATIENT)]
  
  # Classify patients
  pat_surv[, GROUP := "Neither"]
  pat_surv[PATIENT_ID %in% setdiff(mut_g1, mut_g2), GROUP := g1]
  pat_surv[PATIENT_ID %in% setdiff(mut_g2, mut_g1), GROUP := g2]
  pat_surv[PATIENT_ID %in% intersect(mut_g1, mut_g2), GROUP := "Both"]
  
  # Counts
  grp_counts <- pat_surv[, .N, by = GROUP]
  both_n <- grp_counts[GROUP == "Both", N]
  g1_n <- grp_counts[GROUP == g1, N]
  g2_n <- grp_counts[GROUP == g2, N]
  
  cat(sprintf("\n--- %s ---\n", pname))
  cat(sprintf("  Both: %d, %s only: %d, %s only: %d, Neither: %d\n",
              both_n, g1, g1_n, g2, g2_n, 
              grp_counts[GROUP == "Neither", N]))
  
  if (both_n >= 10) {
    # Kaplan-Meier
    surv_obj <- Surv(pat_surv$OS_MONTHS, pat_surv$EVENT)
    
    # Log-rank test
    lr_test <- survdiff(surv_obj ~ GROUP, data = pat_surv)
    p_logrank <- 1 - pchisq(lr_test$chisq, df = length(unique(pat_surv$GROUP)) - 1)
    cat(sprintf("  Log-rank P = %.4f\n", p_logrank))
    
    # Cox regression: Both vs Neither (reference)
    pat_surv[, GROUP_REF := relevel(factor(GROUP), ref = "Neither")]
    cox_fit <- coxph(Surv(OS_MONTHS, EVENT) ~ GROUP_REF, data = pat_surv)
    cox_summary <- summary(cox_fit)
    
    # Both vs Neither HR
    both_coef <- grep("Both", rownames(cox_summary$coefficients), value = TRUE)
    if (length(both_coef) > 0) {
      hr_both <- cox_summary$coefficients[both_coef[1], "exp(coef)"]
      p_both <- cox_summary$coefficients[both_coef[1], "Pr(>|z|)"]
      cat(sprintf("  HR (Both vs Neither) = %.2f, P = %.4f\n", hr_both, p_both))
    }
    
    # Multivariate: adjust for age, sex
    # Only run if enough events
    if (sum(pat_surv$EVENT) >= 50) {
      cox_multi <- coxph(Surv(OS_MONTHS, EVENT) ~ GROUP_REF + AGE + SEX, data = pat_surv)
      cox_multi_sum <- summary(cox_multi)
      if (length(both_coef) > 0) {
        hr_adj <- cox_multi_sum$coefficients[both_coef[1], "exp(coef)"]
        p_adj <- cox_multi_sum$coefficients[both_coef[1], "Pr(>|z|)"]
        cat(sprintf("  Adjusted HR (Both vs Neither) = %.2f, P = %.4f\n", hr_adj, p_adj))
      }
    }
    
    # Store results
    surv_results[[pname]] <- list(
      pair = pname,
      both_n = both_n,
      logrank_p = p_logrank,
      hr = if (exists("hr_both")) hr_both else NA,
      hr_p = if (exists("p_both")) p_both else NA,
      hr_adj = if (exists("hr_adj")) hr_adj else NA,
      hr_adj_p = if (exists("p_adj")) p_adj else NA
    )
    
    # ---- Generate KM curve ----
    km_fit <- survfit(surv_obj ~ GROUP, data = pat_surv)
    
    # Custom KM plot using base ggplot2 (no survminer dependency)
    km_data <- data.frame(
      time = km_fit$time,
      surv = km_fit$surv,
      group = rep(names(km_fit$strata), km_fit$strata)
    )
    # Clean group names
    km_data$group <- gsub("GROUP=", "", km_data$group)
    
    # Color palette
    group_colors <- c("Neither" = "grey60", "Both" = "#E41A1C")
    if (g1 %in% unique(km_data$group)) group_colors[g1] <- "#377EB8"
    if (g2 %in% unique(km_data$group)) group_colors[g2] <- "#4DAF4A"
    
    p <- ggplot(km_data, aes(x = time, y = surv, color = group)) +
      geom_step(size = 1.0) +
      scale_color_manual(values = group_colors) +
      scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 12)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
      labs(
        title = sprintf("Co-mutation: %s + %s (Pan-cancer)", g1, g2),
        subtitle = sprintf("Log-rank P = %.4f, Both (n=%d) vs Neither HR=%.2f", 
                          p_logrank, both_n, if(exists("hr_both")) hr_both else NA),
        x = "Overall Survival (months)",
        y = "Survival Probability",
        color = "Mutation Status"
      ) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    ggsave(file.path(OUT_DIR, "figures", sprintf("KM_%s_pancancer.png", pname)),
           p, width = 8, height = 6, dpi = 300)
    ggsave(file.path(OUT_DIR, "figures", sprintf("KM_%s_pancancer.tiff", pname)),
           p, width = 8, height = 6, dpi = 600, compression = "lzw")
    
    cat(sprintf("  KM plot saved: KM_%s_pancancer.png/tiff\n", pname))
  } else {
    cat(sprintf("  SKIPPED: only %d patients with both mutations\n", both_n))
    
    surv_results[[pname]] <- list(
      pair = pname,
      both_n = both_n,
      logrank_p = NA,
      hr = NA, hr_p = NA, hr_adj = NA, hr_adj_p = NA
    )
  }
}

# ---- 6. Cancer-type-specific analysis for top pair ----
cat("\n============= CANCER-TYPE-SPECIFIC ANALYSIS =============\n")
# For TP53-KRAS, analyze in key cancer types

# Map sample to cancer type
sam_ct <- sam[, .(PATIENT_ID, CANCER_TYPE)]
sam_ct <- unique(sam_ct)

# Merge cancer type with survival
pat_surv_ct <- merge(pat_surv, sam_ct, by = "PATIENT_ID", all.x = TRUE)

# Focus on NSCLC
nsclc_pat <- pat_surv_ct[grepl("Non-Small Cell Lung Cancer|Lung Adenocarcinoma", CANCER_TYPE)]

if (nrow(nsclc_pat) >= 100) {
  cat(sprintf("\nNSCLC cohort: %d patients with survival\n", nrow(nsclc_pat)))
  
  # KRAS-TP53 in NSCLC
  kras_mut <- mut_per_patient[Hugo_Symbol == "KRAS", unique(PATIENT)]
  tp53_mut <- mut_per_patient[Hugo_Symbol == "TP53", unique(PATIENT)]
  
  nsclc_pat[, GROUP := "Neither"]
  nsclc_pat[PATIENT_ID %in% setdiff(tp53_mut, kras_mut), GROUP := "TP53 only"]
  nsclc_pat[PATIENT_ID %in% setdiff(kras_mut, tp53_mut), GROUP := "KRAS only"]
  nsclc_pat[PATIENT_ID %in% intersect(kras_mut, tp53_mut), GROUP := "Both"]
  
  grp_nsclc <- nsclc_pat[, .N, by = GROUP]
  cat("NSCLC KRAS-TP53:\n")
  print(grp_nsclc)
  
  both_nsclc <- grp_nsclc[GROUP == "Both", N]
  if (both_nsclc >= 10) {
    nsclc_pat[, GROUP_REF := relevel(factor(GROUP), ref = "Neither")]
    cox_nsclc <- coxph(Surv(OS_MONTHS, EVENT) ~ GROUP_REF + AGE + SEX, data = nsclc_pat)
    cat("\nNSCLC Multivariate Cox (adjusted):\n")
    print(summary(cox_nsclc)$coefficients)
    
    # KM plot for NSCLC
    km_nsclc <- survfit(Surv(OS_MONTHS, EVENT) ~ GROUP, data = nsclc_pat)
    km_nsclc_data <- data.frame(
      time = km_nsclc$time,
      surv = km_nsclc$surv,
      group = gsub("GROUP=", "", rep(names(km_nsclc$strata), km_nsclc$strata))
    )
    
    p2 <- ggplot(km_nsclc_data, aes(x = time, y = surv, color = group)) +
      geom_step(size = 1.0) +
      scale_color_manual(values = c("Neither" = "grey60", "Both" = "#E41A1C",
                                     "KRAS only" = "#377EB8", "TP53 only" = "#4DAF4A")) +
      scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 12)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
      labs(title = "KRAS-TP53 Co-mutation in NSCLC",
           subtitle = "MSK-IMPACT 50K Cohort",
           x = "Overall Survival (months)", y = "Survival Probability") +
      theme_bw(base_size = 12) + theme(legend.position = "bottom")
    
    ggsave(file.path(OUT_DIR, "figures", "KM_KRAS_TP53_NSCLC.png"),
           p2, width = 8, height = 6, dpi = 300)
    ggsave(file.path(OUT_DIR, "figures", "KM_KRAS_TP53_NSCLC.tiff"),
           p2, width = 8, height = 6, dpi = 600, compression = "lzw")
    cat("  NSCLC KM plot saved\n")
  }
}

# ---- 7. Summary table ----
cat("\n============= SURVIVAL RESULTS SUMMARY =============\n")
surv_df <- rbindlist(lapply(surv_results, as.data.table), fill = TRUE)
print(surv_df)

fwrite(surv_df, file.path(OUT_DIR, "results", "survival_results.csv"))

# ---- 8. Forest plot of co-mutation effects ----
cat("\nGenerating forest plot...\n")
forest_data <- surv_df[!is.na(hr)]

p_forest <- ggplot(forest_data, aes(x = hr, y = reorder(pair, hr))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = hr_p < 0.05), size = 3) +
  geom_errorbarh(aes(xmin = pmax(0.1, hr * 0.7), xmax = pmin(5, hr * 1.3)), height = 0.2) +
  scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey50"),
                     labels = c("TRUE" = "P<0.05", "FALSE" = "n.s.")) +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3)) +
  labs(title = "Pan-cancer Co-mutation Effects on Overall Survival",
       subtitle = "MSK-IMPACT 50K (n=44,081 with survival data)",
       x = "Hazard Ratio (Both vs Neither, log scale)",
       y = "",
       color = "") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "figures", "Forest_CoMutation_OS.png"),
       p_forest, width = 10, height = 5, dpi = 300)
ggsave(file.path(OUT_DIR, "figures", "Forest_CoMutation_OS.tiff"),
       p_forest, width = 10, height = 5, dpi = 600, compression = "lzw")

cat("============= SURVIVAL ANALYSIS COMPLETE =============\n")
