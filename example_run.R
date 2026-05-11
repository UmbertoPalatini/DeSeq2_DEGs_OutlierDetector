# example_run.R
# ------------------------------------------------------------------------------
# End-to-end demo of the outlier annotation on the bundled synthetic dataset.
# From the repository root run:
#
#   Rscript example_run.R
#
# It will:
#   1. load the annotation script and the heatmap helper
#   2. annotate the example DEGs
#   3. print a summary to the console
#   4. write an annotated results table and a heatmap PNG under example_data/
# ------------------------------------------------------------------------------

this_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  f    <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f) > 0 && nzchar(f[1])) dirname(normalizePath(f[1])) else getwd()
}, error = function(e) getwd())
setwd(this_dir)

source("DEGs_outlier_detection_postDeSeq.R")
source("plot_outlier_heatmap.R")

res <- analyze_differentialabundance_outliers(
  results_file           = "example_data/example_results.tsv",
  normalized_counts_file = "example_data/example_normalized_counts.tsv",
  sample_sheet_file      = "example_data/example_samplesheet.tsv"
)

cat("\n--- Summary ---\n")
print(res$summary)

cat("\n--- First 10 annotated rows ---\n")
show_cols <- c("baseMean", "log2FoldChange", "padj",
               "outlier_flag", "outlier_type", "group_pattern")
show_cols <- intersect(show_cols, colnames(res$comprehensive_results))
print(head(res$comprehensive_results[, show_cols, drop = FALSE], 10))

write.table(
  res$comprehensive_results,
  "example_data/example_annotated_results.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("\nAnnotated results written to example_data/example_annotated_results.tsv\n")

plot_outlier_heatmap(
  annotation             = res,
  normalized_counts_file = "example_data/example_normalized_counts.tsv",
  sample_sheet_file      = "example_data/example_samplesheet.tsv",
  output_file            = "example_data/example_heatmap.png"
)
