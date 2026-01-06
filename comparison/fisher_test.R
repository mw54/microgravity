library(dplyr)

### 1. Set file paths ###
path_to_1_csv <- "~/41467_2023_42013_MOESM4_ESM.csv"  # Human DEG table
path_to_2_csv <- "~/402.csv"                          # Mouse genome table (converted) ##403.csv 404.csv 405.csv

### 2. Read raw data ###
human_deg <- read.csv(path_to_1_csv)
mice_deg  <- read.csv(path_to_2_csv)

cat("Human genes count:", length(unique(human_deg$gene)), "\n")
cat("Mouse rows count :", nrow(mice_deg), "\n")

### 3. Fisher's exact test for specific age group ###
calc_fisher_for_age <- function(human_deg, mice_deg, age_filter) {
  cat("\n========================\n")
  cat("Processing age group:", age_filter, "\n")
  cat("========================\n")
  
  # 3.1 Filter mouse genes by age and 'Overall' celltype
  mice_age_all <- mice_deg %>%
    filter(
      celltype == "Overall",
      age == age_filter
    )
  
  if (nrow(mice_age_all) == 0) {
    warning(paste0("No rows found in Mouse data for age = ", age_filter, " and celltype = 'Overall'"))
    return(NULL)
  }
  
  human_all_genes <- unique(human_deg$gene)
  mouse_all_genes <- unique(mice_age_all$human_ortholog)
  
  cat("Human gene count        :", length(human_all_genes), "\n")
  cat("Mouse human_ortholog count:", length(mouse_all_genes), "\n")
  
  # 3.2 Define background (intersection of genes present in both)
  matched_genes <- intersect(human_all_genes, mouse_all_genes)
  cat("Background matched genes count :", length(matched_genes), "\n")
  
  if (length(matched_genes) == 0) {
    warning("No matched genes found. Cannot perform Fisher test.")
    return(NULL)
  }
  
  # 3.3 Identify significant human genes (within matched background)
  human_sig <- human_deg %>%
    filter(
      gene %in% matched_genes,
      p_val_adj < 0.05,
      abs(avg_log2FC) > 0.1
    ) %>%
    pull(gene) %>%
    unique()
  
  # 3.4 Identify significant mouse genes (within matched background)
  mouse_sig <- mice_age_all %>%
    filter(
      human_ortholog %in% matched_genes,
      pvals_adj < 0.05,
      abs(logfoldchanges) > 0.1
    ) %>%
    pull(human_ortholog) %>%
    unique()
  
  cat("Significant Human genes (in background) :", length(human_sig), "\n")
  cat("Significant Mouse genes (in background) :", length(mouse_sig), "\n")
  
  # 3.5 Calculate counts for contingency table
  overlap_genes    <- intersect(human_sig, mouse_sig)
  only_human_genes <- setdiff(human_sig, mouse_sig)
  only_mouse_genes <- setdiff(mouse_sig, human_sig)
  neither_count    <- length(matched_genes) -
    length(overlap_genes) -
    length(only_human_genes) -
    length(only_mouse_genes)
  
  if (neither_count < 0) {
    warning("neither_count < 0. Check matched_genes definition. Setting to 0.")
    neither_count <- 0
  }
  
  cat("\nCounts for Fisher (limited to matched_genes background):\n")
  cat("Overlap (Both significant) :", length(overlap_genes),    "\n")
  cat("Only Human significant     :", length(only_human_genes), "\n")
  cat("Only Mouse significant     :", length(only_mouse_genes), "\n")
  cat("Neither significant        :", neither_count,           "\n\n")
  
  # 3.6 Construct 2x2 contingency table
  cont_tbl <- matrix(
    c(
      length(overlap_genes),
      length(only_human_genes),
      length(only_mouse_genes),
      neither_count
    ),
    nrow = 2,
    byrow = TRUE
  )
  rownames(cont_tbl) <- c("Human_sig",   "Human_not_sig")
  colnames(cont_tbl) <- c("Mouse_sig",   "Mouse_not_sig")
  
  cat("2x2 Contingency Table:\n")
  print(cont_tbl)
  
  # 3.7 Perform Fisher's exact test (greater = enrichment)
  fisher_res <- fisher.test(cont_tbl, alternative = "greater")
  cat("\nFisher's exact test p value (age =", age_filter, "):",
      signif(fisher_res$p.value, 5), "\n")
  
  invisible(list(
    age           = age_filter,
    background    = matched_genes,
    cont_table    = cont_tbl,
    overlap_genes = overlap_genes,
    fisher        = fisher_res
  ))
}

### 4. Run Fisher test for old and young groups ###
res_old <- calc_fisher_for_age(
  human_deg  = human_deg,
  mice_deg   = mice_deg,
  age_filter = "old"
)

res_young <- calc_fisher_for_age(
  human_deg  = human_deg,
  mice_deg   = mice_deg,
  age_filter = "young"
)

cat("\n--- All Fisher tests completed ---\n")
