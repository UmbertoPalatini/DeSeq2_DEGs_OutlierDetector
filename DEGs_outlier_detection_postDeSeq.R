# DEGs_outlier_detection_postDeSeq.R
# ------------------------------------------------------------------------------
# Post-hoc outlier / artifact filter for DESeq2 (or edgeR) differential
# expression results, tailored to small-n RNA-seq designs where a handful of
# samples with aberrant expression can dominate a GLM fit and produce
# statistically significant but biologically implausible DEGs.
#
# Given (i) a results table, (ii) a library-size-normalized counts matrix, and
# (iii) a sample sheet assigning each sample to a contrast group, the script
# annotates each significant DEG with one of three flags:
#
#   Clean            - no per-sample irregularities detected
#   Outliers_Present - outlier sample(s) detected but the group-level pattern
#                      remains consistent; signal is retained
#   Artifact         - the DEG call is attributable to a single- or few-sample
#                      effect, or expression confined to a vanishingly small
#                      number of samples across the contrast
#
# The script never re-runs differential expression. It only annotates. The
# caller decides whether to drop Artifact-flagged genes from downstream
# analysis.
#
# Entry point:
#   analyze_differentialabundance_outliers(
#     results_file, normalized_counts_file, sample_sheet_file, ...)
#
# See README.md for scientific rationale, algorithm details, and examples.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
})

analyze_differentialabundance_outliers <- function(results_file,
                                                   normalized_counts_file,
                                                   sample_sheet_file,
                                                   padj_threshold          = 0.05,
                                                   min_samples_consistent  = 2,
                                                   min_expressing_samples  = 3) {

  cat("=== POST-DESEQ2 OUTLIER / ARTIFACT ANNOTATION ===\n")

  cat("Loading input files...\n")
  results_table     <- read.table(results_file,           sep = "\t", header = TRUE, row.names = 1)
  normalized_counts <- read.table(normalized_counts_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  sample_data       <- read.table(sample_sheet_file,      sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  cat("- Results table     :", nrow(results_table),     "genes\n")
  cat("- Normalized counts :", nrow(normalized_counts), "genes x",
                               ncol(normalized_counts), "samples\n")
  cat("- Sample sheet      :", nrow(sample_data),       "samples\n\n")

  if (!"sample_id" %in% colnames(sample_data)) {
    if ("sample" %in% colnames(sample_data)) {
      sample_data$sample_id <- sample_data$sample
    } else {
      colnames(sample_data)[1] <- "sample_id"
    }
  }

  common_samples <- intersect(colnames(normalized_counts), sample_data$sample_id)
  if (length(common_samples) == 0) {
    stop("No matching samples between normalized counts and sample sheet.")
  }
  cat("Matching samples between counts and sample sheet:", length(common_samples), "\n")

  normalized_counts     <- normalized_counts[, common_samples, drop = FALSE]
  sample_data_filtered  <- sample_data[sample_data$sample_id %in% common_samples, ]
  rownames(sample_data_filtered) <- sample_data_filtered$sample_id
  sample_data_filtered  <- sample_data_filtered[colnames(normalized_counts), ]

  group_col <- detect_group_column(sample_data_filtered)
  cat("Group column used for contrast:", group_col, "\n\n")

  detect_outlier_driven_degs_from_data(
    results_table, normalized_counts, sample_data_filtered,
    group_col, min_samples_consistent, padj_threshold, min_expressing_samples
  )
}

detect_group_column <- function(sample_data) {
  candidates <- c("condition", "group", "treatment",
                  "Condition", "Group", "Treatment")
  for (col in candidates) {
    if (col %in% colnames(sample_data)) return(col)
  }
  if (ncol(sample_data) >= 2) return(colnames(sample_data)[2])
  NULL
}

# ------------------------------------------------------------------------------
# Rule 1 - single-sample outlier detection on the normalized counts of each
# significant DEG. An outlier is a sample whose value is an extreme departure
# from the distribution of the remaining samples, judged by complementary
# criteria so that both moderate spikes and very sparse-then-spike patterns
# are caught.
# ------------------------------------------------------------------------------
identify_outliers <- function(counts) {
  median_val <- median(counts)
  mean_val   <- mean(counts)

  # Ratio relative to the central tendency of the full vector.
  ratio_hits <- which(counts > (median_val * 20))

  # Fallback when the median is zero (many non-expressing samples).
  if (median_val == 0 && mean_val > 0) {
    ratio_hits <- which(counts > (mean_val * 10))
  }

  # Extreme departures relative to the weakest non-zero signal -
  # catches the pattern "most samples at detection floor, one sample very high".
  non_zero <- counts[counts > 0]
  if (length(non_zero) > 1) {
    extreme_hits <- which(counts > (min(non_zero) * 1000))
    ratio_hits   <- unique(c(ratio_hits, extreme_hits))
  }

  # IQR-based criterion tuned to be conservative (3*IQR above Q3).
  q1  <- quantile(counts, 0.25)
  q3  <- quantile(counts, 0.75)
  iqr <- q3 - q1
  if (iqr > 0) {
    iqr_hits   <- which(counts > (q3 + 3 * iqr))
    ratio_hits <- unique(c(ratio_hits, iqr_hits))
  }

  if (length(ratio_hits) == 0) return(list(has_outliers = FALSE))

  outlier_samples    <- names(counts)[ratio_hits]
  outlier_values     <- counts[ratio_hits]
  non_outlier_values <- counts[-ratio_hits]

  list(
    has_outliers       = TRUE,
    outlier_samples    = outlier_samples,
    outlier_values     = outlier_values,
    non_outlier_values = non_outlier_values,
    outlier_type       = classify_outlier_type(outlier_values, non_outlier_values)
  )
}

classify_outlier_type <- function(outlier_values, non_outlier_values) {
  om <- median(outlier_values)
  nm <- median(non_outlier_values)
  if (om > nm * 2)       "high_expression"
  else if (om < nm / 2)  "low_expression"
  else                   "moderate"
}

# Decide whether detected outliers are likely driving the DEG call. A DEG
# is flagged as Artifact when the outlier sample(s) sit in one contrast
# group AND the departure is large in absolute or relative terms, leaving
# the other samples in that group indistinguishable from the opposite group.
assess_outlier_dependence <- function(counts, sample_data, outlier_info,
                                      min_samples, group_col) {
  outlier_samples    <- outlier_info$outlier_samples
  outlier_values     <- outlier_info$outlier_values
  non_outlier_values <- outlier_info$non_outlier_values

  outlier_median     <- median(outlier_values)
  non_outlier_median <- median(non_outlier_values)
  non_outlier_max    <- max(non_outlier_values)

  median_ratio <- if (non_outlier_median > 0) outlier_median / non_outlier_median else Inf
  max_ratio    <- if (non_outlier_max    > 0) outlier_median / non_outlier_max    else Inf

  if (!is.null(group_col) && group_col %in% colnames(sample_data)) {
    groups <- sample_data[[group_col]]
    names(groups) <- rownames(sample_data)

    outlier_groups       <- groups[outlier_samples]
    outlier_group_table  <- table(outlier_groups)
    max_per_group        <- max(outlier_group_table)

    non_outlier_samples  <- setdiff(names(counts), outlier_samples)
    group_sizes          <- table(groups[non_outlier_samples])
    min_group_size       <- if (length(group_sizes) > 0) min(group_sizes) else 0

    likely_artifact <- (
      length(outlier_samples) <= 2 &&
      min_group_size >= min_samples &&
      (
        outlier_median > 10000 ||
        median_ratio   > 100   ||
        (max_per_group == length(outlier_samples) &&
          (median_ratio > 50 || (max_ratio > 10 && outlier_median > 1000)))
      )
    )

    group_pattern <- paste0(
      "Outliers in: ",
      paste(names(outlier_group_table), "(", outlier_group_table, ")", collapse = ", ")
    )
  } else {
    # Without group information, apply a stricter magnitude-only criterion.
    likely_artifact <- (length(outlier_samples) <= 1 &&
                        (median_ratio > 100 || outlier_median > 5000))
    group_pattern   <- "No group information available"
  }

  list(likely_artifact = likely_artifact, group_pattern = group_pattern)
}

# ------------------------------------------------------------------------------
# Main loop. For each significant DEG: check rule 1 (single-sample outliers),
# then overlay rule 2 (sparse_expression).
# ------------------------------------------------------------------------------
detect_outlier_driven_degs_from_data <- function(results_table, normalized_counts,
                                                 sample_data, group_col,
                                                 min_samples_consistent,
                                                 padj_threshold,
                                                 min_expressing_samples) {

  cat("=== ANNOTATION PASS ===\n")

  padj_col <- pick_padj_column(results_table)
  lfc_col  <- pick_lfc_column(results_table)

  sig_idx    <- which(results_table[[padj_col]] < padj_threshold &
                      !is.na(results_table[[padj_col]]))
  sig_genes  <- rownames(results_table)[sig_idx]
  cat("Significant DEGs at padj <", padj_threshold, ":", length(sig_genes), "\n")

  common_genes <- intersect(sig_genes, rownames(normalized_counts))
  if (length(common_genes) == 0) {
    stop("No common genes between results table and normalized counts.")
  }
  cat("Genes present in both results and counts :", length(common_genes), "\n\n")

  comprehensive_results <- results_table[common_genes, , drop = FALSE]
  comprehensive_results <- cbind(gene_id = rownames(comprehensive_results),
                                 comprehensive_results)
  comprehensive_results$has_outliers      <- FALSE
  comprehensive_results$outlier_type      <- NA_character_
  comprehensive_results$outlier_samples   <- NA_character_
  comprehensive_results$outlier_values    <- NA_character_
  comprehensive_results$non_outlier_range <- NA_character_
  comprehensive_results$group_pattern     <- NA_character_
  comprehensive_results$likely_artifact   <- FALSE
  comprehensive_results$outlier_flag      <- "Clean"

  progress_interval <- max(1, floor(length(common_genes) / 20))

  # ---- Rule 1 : per-gene single-sample outlier ----
  for (i in seq_along(common_genes)) {
    gene <- common_genes[i]
    if (i %% progress_interval == 0) {
      cat("Processing gene", i, "of", length(common_genes), "\n")
    }

    gene_counts <- as.numeric(normalized_counts[gene, ])
    names(gene_counts) <- colnames(normalized_counts)
    if (var(gene_counts) == 0) next

    outlier_info <- identify_outliers(gene_counts)
    if (!outlier_info$has_outliers) next

    assessment <- assess_outlier_dependence(
      gene_counts, sample_data, outlier_info,
      min_samples_consistent, group_col
    )

    comprehensive_results[gene, "has_outliers"]      <- TRUE
    comprehensive_results[gene, "outlier_type"]      <- outlier_info$outlier_type
    comprehensive_results[gene, "outlier_samples"]   <- paste(outlier_info$outlier_samples, collapse = ";")
    comprehensive_results[gene, "outlier_values"]    <- paste(round(outlier_info$outlier_values, 2), collapse = ";")
    comprehensive_results[gene, "non_outlier_range"] <- paste(round(range(outlier_info$non_outlier_values), 2), collapse = "-")
    comprehensive_results[gene, "group_pattern"]     <- assessment$group_pattern
    comprehensive_results[gene, "likely_artifact"]   <- assessment$likely_artifact
    comprehensive_results[gene, "outlier_flag"]      <- if (assessment$likely_artifact) "Artifact" else "Outliers_Present"
  }

  # ---- Rule 2 : sparse_expression ----
  # A DEG whose entire non-zero signal comes from a very small number of
  # samples across the full contrast cannot inform a group-level inference,
  # because within-group variance cannot be meaningfully estimated.
  for (gene in rownames(comprehensive_results)) {
    gene_counts  <- as.numeric(normalized_counts[gene, ])
    n_expressing <- sum(gene_counts > 0, na.rm = TRUE)
    if (n_expressing < min_expressing_samples) {
      expressing <- colnames(normalized_counts)[gene_counts > 0]
      comprehensive_results[gene, "has_outliers"]     <- TRUE
      comprehensive_results[gene, "outlier_type"]     <- "sparse_expression"
      comprehensive_results[gene, "outlier_samples"]  <- paste(expressing, collapse = ";")
      comprehensive_results[gene, "outlier_values"]   <- paste(round(gene_counts[gene_counts > 0], 2), collapse = ";")
      comprehensive_results[gene, "group_pattern"]    <- paste0(
        "Expressed in ", n_expressing, "/", length(gene_counts), " samples"
      )
      comprehensive_results[gene, "likely_artifact"]  <- TRUE
      comprehensive_results[gene, "outlier_flag"]     <- "Artifact"
    }
  }

  comprehensive_results <- comprehensive_results[order(
    comprehensive_results$outlier_flag != "Clean",
    comprehensive_results$likely_artifact,
    comprehensive_results[[padj_col]]
  ), ]

  outlier_flags <- comprehensive_results[comprehensive_results$has_outliers,
                                         c("outlier_type", "outlier_samples",
                                           "outlier_values", "non_outlier_range",
                                           "group_pattern", "likely_artifact",
                                           lfc_col, padj_col), drop = FALSE]
  outlier_flags$gene <- rownames(outlier_flags)
  outlier_flags <- outlier_flags[, c("gene", setdiff(colnames(outlier_flags), "gene"))]
  outlier_flags <- outlier_flags[order(-outlier_flags$likely_artifact,
                                       outlier_flags[[padj_col]]), ]

  cat("\n=== SUMMARY ===\n")
  print(table(comprehensive_results$outlier_flag))

  list(
    outlier_flags         = outlier_flags,
    comprehensive_results = comprehensive_results,
    summary = list(
      total_degs       = length(common_genes),
      degs_with_flags  = nrow(outlier_flags),
      likely_artifacts = sum(comprehensive_results$outlier_flag == "Artifact"),
      outliers_present = sum(comprehensive_results$outlier_flag == "Outliers_Present"),
      clean_degs       = sum(comprehensive_results$outlier_flag == "Clean")
    )
  )
}

# ------------------------------------------------------------------------------
# Column-name helpers - accept DESeq2 and edgeR/limma conventions.
# ------------------------------------------------------------------------------
pick_padj_column <- function(results_table) {
  for (col in c("padj", "adj.P.Val")) {
    if (col %in% colnames(results_table)) return(col)
  }
  hits <- grep("adj|padj|fdr|BH", colnames(results_table), ignore.case = TRUE, value = TRUE)
  if (length(hits) > 0) return(hits[1])
  stop("Cannot find an adjusted p-value column in results table.")
}

pick_lfc_column <- function(results_table) {
  for (col in c("log2FoldChange", "logFC")) {
    if (col %in% colnames(results_table)) return(col)
  }
  NA_character_
}

cat("Script loaded. Entry point:\n")
cat("  res <- analyze_differentialabundance_outliers(\n")
cat("           results_file           = 'deseq2.results.tsv',\n")
cat("           normalized_counts_file = 'deseq2.normalised_counts.tsv',\n")
cat("           sample_sheet_file      = 'samplesheet.tsv')\n")
