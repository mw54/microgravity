library(dplyr)
library(pheatmap)
library(tibble)

### 1. Set file paths ###
path_to_1_csv <- "~/41467_2023_42013_MOESM4_ESM.csv"
path_to_2_csv <- "~/402.csv" ##403.csv 404.csv 405.csv

output_file_base <- "~/260104/"
output_file_old <- paste0(output_file_base, "402_old_concordant.pdf")
output_file_young <- paste0(output_file_base, "402_young_concordant.pdf")

if (!dir.exists(dirname(output_file_base))) {
  dir.create(dirname(output_file_base), recursive = TRUE)
}

### 2. Define plotting function ###
generate_overlap_heatmap <- function(human_data, all_mice_data, age_filter, output_filename) {
  
  message(paste("\n--- Processing:", age_filter, "group ---"))
  
  mice_filtered <- all_mice_data %>%
    filter(celltype == "Overall", age == age_filter, pvals_adj < 0.05, abs(logfoldchanges) > 0.1) %>%
    group_by(human_ortholog) %>%
    summarise(logfoldchanges = mean(logfoldchanges, na.rm = TRUE))
  
  message(paste("Mice:", age_filter, "filtered significant genes:", nrow(mice_filtered)))
  
  if (nrow(mice_filtered) == 0) {
    warning("No significant mouse genes found. Skipping.")
    return()
  }
  
  overlap_data <- inner_join(
    human_data,
    mice_filtered,
    by = c("gene" = "human_ortholog")
  )
  
  message(paste("Found", nrow(overlap_data), "overlapping significant genes."))
  
  if (nrow(overlap_data) == 0) return()
  
  concordant_data <- overlap_data %>%
    filter(sign(avg_log2FC) == sign(logfoldchanges))
  
  message(paste("Concordant genes:", nrow(concordant_data)))
  
  if (nrow(concordant_data) == 0) return()
  
  plot_data <- concordant_data %>%
    rename(Human = avg_log2FC, Mouse = logfoldchanges)
  
  plot_matrix <- plot_data %>%
    column_to_rownames(var = "gene") %>%
    as.matrix()
  
  my_colors <- colorRampPalette(c("blue", "white", "red"))(100)
  
  max_lfc <- max(abs(plot_matrix), na.rm = TRUE)
  if (max_lfc == 0) max_lfc <- 1
  
  my_breaks <- seq(-max_lfc, max_lfc, length.out = length(my_colors) + 1)
  
  message(paste("Plotting heatmap, saving to:", output_filename))
  
  # --- Plotting ---
  pheatmap(
    plot_matrix,
    color = my_colors,
    breaks = my_breaks,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    fontsize_row = 8,
    main = paste("Human vs Mouse(Spleen) -", age_filter),
    legend = TRUE,
    legend_breaks = c(-max_lfc, 0, max_lfc),
    legend_labels = c(paste0("-", round(max_lfc, 2)),
                      "0",
                      paste0("+", round(max_lfc, 2))),
    filename = output_filename,
    width = 4,   # Width adjusted to 4
    height = 14
  )
  # --- End Plotting ---
  
  message("✔ Heatmap generated.")
}

### 3. Main execution flow ###
tryCatch({
  human_deg <- read.csv(path_to_1_csv)
  mice_deg <- read.csv(path_to_2_csv)
  message("✔ Data loaded successfully.")
  
  human_filtered <- human_deg %>%
    filter(abs(avg_log2FC) > 0.1, p_val_adj < 0.05) %>%
    group_by(gene) %>%
    summarise(avg_log2FC = mean(avg_log2FC, na.rm = TRUE))
  
  message(paste("Human significant genes:", nrow(human_filtered)))
  
  generate_overlap_heatmap(
    human_data = human_filtered,
    all_mice_data = mice_deg,
    age_filter = "old",
    output_filename = output_file_old
  )
  
  generate_overlap_heatmap(
    human_data = human_filtered,
    all_mice_data = mice_deg,
    age_filter = "young",
    output_filename = output_file_young
  )
  
  message("\n--- All processing completed ---")
  
}, error = function(e) {
  message("❌ Error occurred: ", e$message)
})
