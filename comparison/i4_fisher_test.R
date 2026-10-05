### I4 vs mouse PBMC (404.csv) -- Fisher's exact test ###
### Mirrors fisher_test.R exactly, adapted for I4's per-celltype DEG table
### (no pooled p_val_adj / avg_log2FC column). Base R only.
### Input: the I4 PBMC DEG table from Kim et al., Nat Commun 2024 (Inspiration4),
###   Supplementary Data 15, sheet "DEGs" (immediate post-flight R+1 vs pre-flight,
###   per-celltype log2FC and padj for 9 PBMC cell types), saved as CSV with the
###   column header row (human_gene, CD4_T, ..., other, CD4_T.padj, ..., other.padj)
###   as the first line. Extra columns in the sheet are ignored.
### Mouse input: ~/404.csv (PBMC DEG table with human_ortholog column, from
###   convert_orthologs.py), same as fisher_test.R / overlap_heatmap.R.

### 1. Set file paths ###
path_to_i4_csv    <- "~/I4_DEGs.csv"   # I4 per-celltype DEG table (Supplementary Data 15)
path_to_mouse_csv <- "~/404.csv"       # Mouse PBMC DEG table (has human_ortholog col)

celltypes <- c("CD4_T","CD8_T","other_T","B","NK","CD14_Mono","CD16_Mono","DC","other")

### 2. Read raw data ###
i4_raw    <- read.csv(path_to_i4_csv)
mice_deg  <- read.csv(path_to_mouse_csv)

# The I4 table contains a handful of exact duplicate rows (same gene, identical
# values in every celltype column) -- collapse those before anything else.
i4_key <- do.call(paste, c(i4_raw[c("human_gene", celltypes, paste0(celltypes, ".padj"))], sep = "\r"))
i4_raw <- i4_raw[!duplicated(i4_key), ]

cat("I4 unique genes (raw table)     :", length(unique(i4_raw$human_gene)), "\n")
cat("Mouse rows count                :", nrow(mice_deg), "\n")

### 3. Define the I4 "DEG set" (human_sig analog) ###
# A gene counts as an I4 DEG if it passes padj < 0.05 AND |log2FC| > 0.1 in EVERY
# one of the 9 annotated PBMC celltypes (same thresholds as the rest of the
# cross-species analysis, applied per celltype). Yields 152 genes.
padj_cols <- paste0(celltypes, ".padj")
sig_in_all <- apply(i4_raw[padj_cols] < 0.05, 1, all) &
              apply(abs(i4_raw[celltypes]) > 0.1, 1, all)
i4_sig_genes <- unique(i4_raw$human_gene[sig_in_all & !is.na(sig_in_all)])
cat("I4 genes significant in ALL 9 celltypes (I4 DEG set):", length(i4_sig_genes), "\n\n")

### 4. Fisher's exact test for a given mouse age group ###
calc_fisher_for_age <- function(i4_all_genes, i4_sig_genes, mice_deg, age_filter) {
  cat("========================\n")
  cat("Processing age group:", age_filter, "\n")
  cat("========================\n")

  # 4.1 Filter mouse genes by age and 'Overall' celltype (identical to fisher_test.R)
  mice_age_all <- mice_deg[mice_deg$celltype == "Overall" & mice_deg$age == age_filter, ]

  if (nrow(mice_age_all) == 0) {
    warning(paste0("No rows found in Mouse data for age = ", age_filter, " and celltype = 'Overall'"))
    return(NULL)
  }

  mouse_all_genes <- unique(mice_age_all$human_ortholog)

  cat("I4 gene count (background)        :", length(i4_all_genes), "\n")
  cat("Mouse human_ortholog count         :", length(mouse_all_genes), "\n")

  # 4.2 Define background (intersection of genes present in both tables)
  matched_genes <- intersect(i4_all_genes, mouse_all_genes)
  cat("Background matched genes count     :", length(matched_genes), "\n")

  if (length(matched_genes) == 0) {
    warning("No matched genes found. Cannot perform Fisher test.")
    return(NULL)
  }

  # 4.3 I4-significant genes within matched background
  i4_sig <- intersect(i4_sig_genes, matched_genes)

  # 4.4 Mouse-significant genes within matched background
  #     Row-level filter THEN unique() -- if a human_ortholog has multiple mouse
  #     paralog rows (common in this table) with disagreeing directions, the
  #     ortholog still counts as significant as long as >=1 row passes the
  #     threshold. Direction is intentionally ignored here (matches fisher_test.R).
  mouse_sig_rows <- mice_age_all[
    mice_age_all$human_ortholog %in% matched_genes &
    mice_age_all$pvals_adj < 0.05 &
    abs(mice_age_all$logfoldchanges) > 0.1,
  ]
  mouse_sig <- unique(mouse_sig_rows$human_ortholog)

  cat("Significant I4 genes (in background)   :", length(i4_sig), "\n")
  cat("Significant Mouse genes (in background):", length(mouse_sig), "\n")

  # 4.5 Counts for contingency table
  overlap_genes    <- intersect(i4_sig, mouse_sig)
  only_i4_genes    <- setdiff(i4_sig, mouse_sig)
  only_mouse_genes <- setdiff(mouse_sig, i4_sig)
  neither_count    <- length(matched_genes) -
    length(overlap_genes) - length(only_i4_genes) - length(only_mouse_genes)
  if (neither_count < 0) {
    warning("neither_count < 0. Check matched_genes definition. Setting to 0.")
    neither_count <- 0
  }

  cat("\nCounts for Fisher (limited to matched_genes background):\n")
  cat("Overlap (Both significant) :", length(overlap_genes), "\n")
  cat("Only I4 significant        :", length(only_i4_genes), "\n")
  cat("Only Mouse significant     :", length(only_mouse_genes), "\n")
  cat("Neither significant        :", neither_count, "\n\n")

  # 4.6 Construct 2x2 contingency table
  cont_tbl <- matrix(
    c(length(overlap_genes), length(only_i4_genes),
      length(only_mouse_genes), neither_count),
    nrow = 2, byrow = TRUE
  )
  rownames(cont_tbl) <- c("I4_sig", "I4_not_sig")
  colnames(cont_tbl) <- c("Mouse_sig", "Mouse_not_sig")
  cat("2x2 Contingency Table:\n")
  print(cont_tbl)

  # 4.7 Fisher's exact test (one-sided, enrichment)
  fisher_res <- fisher.test(cont_tbl, alternative = "greater")
  cat("\nFisher's exact test p value (age =", age_filter, "):",
      signif(fisher_res$p.value, 5), "\n\n")

  invisible(list(
    age = age_filter, background = matched_genes, cont_table = cont_tbl,
    overlap_genes = overlap_genes, only_i4 = length(only_i4_genes),
    only_mouse = length(only_mouse_genes), fisher = fisher_res
  ))
}

### 5. Run for old and young mouse groups ###
i4_all_genes <- unique(i4_raw$human_gene)

res_old   <- calc_fisher_for_age(i4_all_genes, i4_sig_genes, mice_deg, "old")
res_young <- calc_fisher_for_age(i4_all_genes, i4_sig_genes, mice_deg, "young")

cat("--- All Fisher tests completed ---\n")

### 6. Save overlap gene lists for downstream heatmap / reporting ###
if (!is.null(res_old))   writeLines(res_old$overlap_genes,   "~/i4_pbmc_overlap_old.txt")
if (!is.null(res_young)) writeLines(res_young$overlap_genes, "~/i4_pbmc_overlap_young.txt")

### 7. Venn-diagram display counts ###
# Convention used for the Venn diagrams (Fig 6, Fig S5): circles show the full
# (unrestricted) significant gene sets on each side, while the p-value printed
# under each Venn comes from the background-restricted Fisher test above.
for (a in c("old", "young")) {
  ma <- mice_deg[mice_deg$celltype == "Overall" & mice_deg$age == a, ]
  msig <- unique(ma$human_ortholog[ma$pvals_adj < 0.05 & abs(ma$logfoldchanges) > 0.1])
  cat(sprintf("Venn display (%s): I4-only = %d | overlap = %d | mouse-only = %d\n",
              a, length(setdiff(i4_sig_genes, msig)), length(intersect(i4_sig_genes, msig)),
              length(setdiff(msig, i4_sig_genes))))
}
