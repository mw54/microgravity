### I4 vs mouse PBMC (404.csv) -- overlap heatmap ###
### Mirrors overlap_heatmap.R exactly, adapted for I4's per-celltype DEG table.
### Base R + pheatmap only.
### Input: the I4 PBMC DEG table from Kim et al., Nat Commun 2024 (Inspiration4),
###   Supplementary Data 15, sheet "DEGs" (immediate post-flight R+1 vs pre-flight,
###   per-celltype log2FC and padj for 9 PBMC cell types), saved as CSV with the
###   column header row (human_gene, CD4_T, ..., other, CD4_T.padj, ..., other.padj)
###   as the first line. Extra columns in the sheet are ignored.
### Mouse input: ~/404.csv (PBMC DEG table with human_ortholog column, from
###   convert_orthologs.py), same as fisher_test.R / overlap_heatmap.R.

library(pheatmap)

### 1. Set file paths ###
path_to_i4_csv    <- "~/I4_DEGs.csv"
path_to_mouse_csv <- "~/404.csv"
output_file_base  <- "~/I4_PBMC/"

celltypes <- c("CD4_T","CD8_T","other_T","B","NK","CD14_Mono","CD16_Mono","DC","other")
padj_cols <- paste0(celltypes, ".padj")

output_file_old   <- paste0(output_file_base, "I4_PBMC_old_concordant.pdf")
output_file_young <- paste0(output_file_base, "I4_PBMC_young_concordant.pdf")

if (!dir.exists(output_file_base)) dir.create(output_file_base, recursive = TRUE)

### 2. Load & dedupe I4 table (see i4_fisher_test.R for rationale) ###
i4_raw <- read.csv(path_to_i4_csv)
i4_key <- do.call(paste, c(i4_raw[c("human_gene", celltypes, padj_cols)], sep = "\r"))
i4_raw <- i4_raw[!duplicated(i4_key), ]

### 3. Build the I4 "human_filtered" analog ###
# Gene must pass padj<0.05 AND |log2FC|>0.1 in ALL 9 celltypes (152 genes;
# same DEG-set definition as the Fisher test / Venn diagrams).
# log2FC representative value = mean across the 9 celltypes -- justified
# because 147/152 genes have a fully concordant sign (the other 5: 8/9) across celltypes
# (checked separately); this mirrors how overlap_heatmap.R collapses
# duplicate mouse paralog rows with mean(logfoldchanges).
sig_in_all <- apply(i4_raw[padj_cols] < 0.05, 1, all) &
              apply(abs(i4_raw[celltypes]) > 0.1, 1, all)
i4_sig <- i4_raw[sig_in_all & !is.na(sig_in_all), ]
i4_filtered <- data.frame(
  gene      = i4_sig$human_gene,
  avg_log2FC = rowMeans(i4_sig[celltypes], na.rm = TRUE)
)
cat("I4 filtered gene count:", nrow(i4_filtered), "\n")

# Direction-consistency sanity check (reported, not filtered on here --
# concordance with the MOUSE side is what determines inclusion in the heatmap)
signs <- sign(as.matrix(i4_sig[celltypes]))
consistent <- rowSums(signs == signs[, 1]) == ncol(signs)
cat("I4 genes with fully consistent sign across all 9 celltypes:",
    sum(consistent), "/", nrow(i4_sig), "\n\n")

### 4. Load mouse data ###
mice_deg <- read.csv(path_to_mouse_csv)

### 5. Plotting function (same structure as overlap_heatmap.R) ###
generate_overlap_heatmap <- function(i4_data, all_mice_data, age_filter, output_filename) {

  message(paste("\n--- Processing:", age_filter, "group ---"))

  mice_filtered_rows <- all_mice_data[
    all_mice_data$celltype == "Overall" &
    all_mice_data$age == age_filter &
    all_mice_data$pvals_adj < 0.05 &
    abs(all_mice_data$logfoldchanges) > 0.1,
  ]
  # Collapse duplicate human_ortholog rows (mouse paralogs) via mean --
  # identical behaviour to group_by(human_ortholog) %>% summarise(mean(...))
  mice_filtered <- aggregate(
    logfoldchanges ~ human_ortholog, data = mice_filtered_rows, FUN = mean, na.rm = TRUE
  )

  message(paste("Mice:", age_filter, "filtered significant genes:", nrow(mice_filtered)))
  if (nrow(mice_filtered) == 0) { warning("No significant mouse genes found. Skipping."); return(invisible(NULL)) }

  overlap_data <- merge(i4_data, mice_filtered, by.x = "gene", by.y = "human_ortholog")
  message(paste("Found", nrow(overlap_data), "overlapping significant genes."))
  if (nrow(overlap_data) == 0) return(invisible(NULL))

  concordant_data <- overlap_data[sign(overlap_data$avg_log2FC) == sign(overlap_data$logfoldchanges), ]
  message(paste("Concordant genes:", nrow(concordant_data)))
  if (nrow(concordant_data) == 0) return(invisible(NULL))

  plot_data <- concordant_data
  colnames(plot_data)[colnames(plot_data) == "avg_log2FC"]     <- "Human"
  colnames(plot_data)[colnames(plot_data) == "logfoldchanges"] <- "Mouse"

  plot_matrix <- as.matrix(plot_data[, c("Human", "Mouse")])
  rownames(plot_matrix) <- plot_data$gene

  my_colors <- colorRampPalette(c("blue", "white", "red"))(100)
  max_lfc <- max(abs(plot_matrix), na.rm = TRUE)
  if (max_lfc == 0) max_lfc <- 1
  my_breaks <- seq(-max_lfc, max_lfc, length.out = length(my_colors) + 1)

  message(paste("Plotting heatmap, saving to:", output_filename))

  pheatmap(
    plot_matrix,
    color = my_colors,
    breaks = my_breaks,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    fontsize_row = 8,
    main = paste("I4 vs Mouse(PBMC) -", age_filter),
    legend = TRUE,
    legend_breaks = c(-max_lfc, 0, max_lfc),
    legend_labels = c(paste0("-", round(max_lfc, 2)), "0", paste0("+", round(max_lfc, 2))),
    filename = output_filename,
    width = 4,                                  # same as overlap_heatmap.R
    height = 2 + 0.42 * nrow(plot_matrix)       # scaled to gene count (few genes here)
  )

  message("Heatmap generated.")
  invisible(concordant_data)
}

### 6. Run for both age groups ###
res_old <- tryCatch({
  generate_overlap_heatmap(i4_filtered, mice_deg, "old", output_file_old)
}, error = function(e) { message("Error (old): ", e$message); NULL })

res_young <- tryCatch({
  generate_overlap_heatmap(i4_filtered, mice_deg, "young", output_file_young)
}, error = function(e) { message("Error (young): ", e$message); NULL })

if (!is.null(res_old))   write.csv(res_old,   paste0(output_file_base, "I4_PBMC_old_concordant_genes.csv"), row.names = FALSE)
if (!is.null(res_young)) write.csv(res_young, paste0(output_file_base, "I4_PBMC_young_concordant_genes.csv"), row.names = FALSE)

message("\n--- All processing completed ---")
