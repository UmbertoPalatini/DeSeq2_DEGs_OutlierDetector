# plot_outlier_heatmap.R
# ------------------------------------------------------------------------------
# Visual diagnostic for the post-DESeq2 outlier annotation. Given the output of
# analyze_differentialabundance_outliers() plus the same normalized counts and
# sample sheet used for the annotation, it plots a row-scaled heatmap of the
# annotated DEGs with:
#
#   * column annotation  : contrast group (from the sample sheet)
#   * row annotation     : outlier_flag (Clean / Outliers_Present / Artifact)
#                          and, optionally, outlier_type
#
# The colour of the row-annotation bar immediately shows which genes have been
# flagged and which are retained, on the same axes as the expression pattern
# that drove the flag.
#
# Usage:
#   source("DEGs_outlier_detection_postDeSeq.R")
#   source("plot_outlier_heatmap.R")
#
#   res <- analyze_differentialabundance_outliers(...)
#   plot_outlier_heatmap(
#     annotation             = res,
#     normalized_counts_file = "example_data/example_normalized_counts.tsv",
#     sample_sheet_file      = "example_data/example_samplesheet.tsv",
#     output_file            = "example_data/example_heatmap.png"
#   )
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(pheatmap)
})

plot_outlier_heatmap <- function(annotation,
                                 normalized_counts_file,
                                 sample_sheet_file,
                                 output_file = NULL,
                                 top_n       = NULL,
                                 width       = 9,
                                 height      = 8,
                                 dpi         = 200) {

  comp <- annotation$comprehensive_results
  if (!is.data.frame(comp)) {
    stop("annotation$comprehensive_results not found - pass the full return ",
         "value of analyze_differentialabundance_outliers().")
  }

  counts <- read.table(normalized_counts_file, sep = "\t", header = TRUE,
                       row.names = 1, check.names = FALSE)
  sheet  <- read.table(sample_sheet_file, sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE)

  if (!"sample_id" %in% colnames(sheet)) {
    if ("sample" %in% colnames(sheet)) sheet$sample_id <- sheet$sample
    else colnames(sheet)[1] <- "sample_id"
  }
  common <- intersect(colnames(counts), sheet$sample_id)
  counts <- counts[, common, drop = FALSE]
  sheet  <- sheet[match(common, sheet$sample_id), , drop = FALSE]

  group_col <- intersect(c("condition", "group", "treatment",
                           "Condition", "Group", "Treatment"),
                         colnames(sheet))[1]
  if (is.na(group_col)) group_col <- colnames(sheet)[2]

  genes_plot <- intersect(rownames(comp), rownames(counts))
  if (!is.null(top_n) && length(genes_plot) > top_n) {
    genes_plot <- genes_plot[seq_len(top_n)]
  }

  mat     <- as.matrix(counts[genes_plot, , drop = FALSE])
  mat_log <- log2(mat + 1)

  row_ann <- data.frame(
    outlier_flag = factor(comp[genes_plot, "outlier_flag"],
                          levels = c("Clean", "Outliers_Present", "Artifact")),
    outlier_type = ifelse(is.na(comp[genes_plot, "outlier_type"]),
                          "none", comp[genes_plot, "outlier_type"]),
    row.names = genes_plot
  )

  col_ann <- data.frame(
    group = sheet[[group_col]],
    row.names = sheet$sample_id
  )

  flag_colors <- c(Clean            = "#4C9F70",
                   Outliers_Present = "#F0B94B",
                   Artifact         = "#D9534F")
  type_levels <- unique(row_ann$outlier_type)
  type_palette <- c(none              = "#DDDDDD",
                    high_expression   = "#8C6BB1",
                    low_expression    = "#6BAED6",
                    moderate          = "#BDBDBD",
                    sparse_expression = "#FD8D3C")
  type_palette <- type_palette[names(type_palette) %in% type_levels]

  group_levels <- unique(col_ann$group)
  group_palette <- setNames(
    c("#5B8DB8", "#C05555", "#7BAF7B", "#B08ECA")[seq_along(group_levels)],
    group_levels
  )

  ann_colors <- list(
    outlier_flag = flag_colors,
    outlier_type = type_palette,
    group        = group_palette
  )

  gaps <- NULL
  if (is.factor(row_ann$outlier_flag) || is.character(row_ann$outlier_flag)) {
    ord     <- order(row_ann$outlier_flag)
    mat_log <- mat_log[ord, , drop = FALSE]
    row_ann <- row_ann[ord, , drop = FALSE]
    flag_vec <- as.character(row_ann$outlier_flag)
    gaps <- cumsum(rle(flag_vec)$lengths)
    gaps <- gaps[-length(gaps)]
  }

  ph <- pheatmap(
    mat_log,
    cluster_rows      = FALSE,
    cluster_cols      = TRUE,
    scale             = "row",
    annotation_row    = row_ann,
    annotation_col    = col_ann,
    annotation_colors = ann_colors,
    gaps_row          = gaps,
    show_rownames     = nrow(mat_log) <= 60,
    show_colnames     = TRUE,
    fontsize_row      = 7,
    fontsize_col      = 8,
    main              = "DEG expression, annotated by outlier status",
    silent            = !is.null(output_file)
  )

  if (!is.null(output_file)) {
    ext <- tolower(tools::file_ext(output_file))
    if (ext == "pdf") {
      pdf(output_file, width = width, height = height)
    } else {
      png(output_file, width = width * dpi, height = height * dpi, res = dpi)
    }
    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
    dev.off()
    cat("Heatmap written to:", output_file, "\n")
  }

  invisible(ph)
}
