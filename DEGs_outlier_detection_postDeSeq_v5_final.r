# Fixed assess function with correct signature and better artifact logic
library(DESeq2)
library(dplyr)

# Function to work with nf-core/differentialabundance output files
analyze_differentialabundance_outliers <- function(results_file, 
                                                   normalized_counts_file,
                                                   sample_sheet_file,
                                                   outlier_threshold = 3, 
                                                   min_samples_consistent = 2,
                                                   padj_threshold = 0.05) {
  
  cat("=== ANALYZING NF-CORE/DIFFERENTIALABUNDANCE OUTPUT FOR OUTLIERS ===\n")
  
  # Load the files
  cat("Loading files...\n")
  results_table <- read.table(results_file, sep = "\t", header = TRUE, row.names = 1)
  normalized_counts <- read.table(normalized_counts_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  sample_data <- read.table(sample_sheet_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  cat("Files loaded successfully:\n")
  cat("- Results table:", nrow(results_table), "genes\n")
  cat("- Normalized counts:", nrow(normalized_counts), "genes,", ncol(normalized_counts), "samples\n")
  cat("- Sample data:", nrow(sample_data), "samples\n\n")
  
  # Ensure sample names match between files
  if(!"sample_id" %in% colnames(sample_data)) {
    if("sample" %in% colnames(sample_data)) {
      sample_data$sample_id <- sample_data$sample
    } else if(ncol(sample_data) >= 1) {
      colnames(sample_data)[1] <- "sample_id"
    }
  }
  
  # Match samples between normalized counts and sample data
  common_samples <- intersect(colnames(normalized_counts), sample_data$sample_id)
  cat("Matching samples found:", length(common_samples), "\n")
  
  if(length(common_samples) == 0) {
    stop("No matching samples between normalized counts and sample data. Check sample names.")
  }
  
  # Filter to common samples
  normalized_counts <- normalized_counts[, common_samples, drop = FALSE]
  sample_data_filtered <- sample_data[sample_data$sample_id %in% common_samples, ]
  rownames(sample_data_filtered) <- sample_data_filtered$sample_id
  sample_data_filtered <- sample_data_filtered[colnames(normalized_counts), ]
  
  # Detect group/condition column
  group_col <- detect_group_column(sample_data_filtered)
  cat("Using group column:", group_col, "\n\n")
  
  # Run outlier detection
  return(detect_outlier_driven_degs_from_data(
    results_table, normalized_counts, sample_data_filtered, 
    group_col, outlier_threshold, min_samples_consistent, padj_threshold
  ))
}

# Helper function to detect the group/condition column
detect_group_column <- function(sample_data) {
  possible_cols <- c("condition", "group", "treatment", "Group", "Condition", "Treatment")
  
  for(col in possible_cols) {
    if(col %in% colnames(sample_data)) {
      return(col)
    }
  }
  
  if(ncol(sample_data) >= 2) {
    return(colnames(sample_data)[2])
  }
  
  return(NULL)
}

# Core outlier detection function
detect_outlier_driven_degs_from_data <- function(results_table, normalized_counts, 
                                                 sample_data, group_col,
                                                 outlier_threshold = 3, 
                                                 min_samples_consistent = 2,
                                                 padj_threshold = 0.05) {
  
  cat("=== OUTLIER-DRIVEN DEG DETECTION ===\n")
  cat("Analyzing DESeq2 results for single-outlier artifacts\n\n")
  
  # Filter to significant DEGs
  if("padj" %in% colnames(results_table)) {
    sig_genes <- rownames(results_table)[which(results_table$padj < padj_threshold & !is.na(results_table$padj))]
  } else if("adj.P.Val" %in% colnames(results_table)) {
    sig_genes <- rownames(results_table)[which(results_table$adj.P.Val < padj_threshold & !is.na(results_table$adj.P.Val))]
  } else {
    pval_cols <- grep("adj|padj|fdr|BH", colnames(results_table), ignore.case = TRUE)
    if(length(pval_cols) > 0) {
      pval_col <- colnames(results_table)[pval_cols[1]]
      sig_genes <- rownames(results_table)[which(results_table[, pval_col] < padj_threshold & !is.na(results_table[, pval_col]))]
      cat("Using p-value column:", pval_col, "\n")
    } else {
      stop("Cannot find adjusted p-value column in results table")
    }
  }
  
  cat("Analyzing", length(sig_genes), "significant DEGs\n")
  
  # Ensure we have matching genes between results and counts
  common_genes <- intersect(sig_genes, rownames(normalized_counts))
  cat("Common genes between results and counts:", length(common_genes), "\n\n")
  
  if(length(common_genes) == 0) {
    stop("No common genes between results table and normalized counts")
  }
  
  # Initialize results
  outlier_flags <- data.frame(
    gene = character(),
    outlier_type = character(),
    outlier_samples = character(),
    outlier_values = character(),
    non_outlier_range = character(),
    group_pattern = character(),
    likely_artifact = logical(),
    log2FoldChange = numeric(),
    padj = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Progress tracking
  progress_interval <- max(1, floor(length(common_genes) / 20))
  
  for(i in seq_along(common_genes)) {
    gene <- common_genes[i]
    
    if(i %% progress_interval == 0) {
      cat("Processing gene", i, "of", length(common_genes), "\n")
    }
    
    gene_counts <- as.numeric(normalized_counts[gene, ])
    names(gene_counts) <- colnames(normalized_counts)
    
    # Skip if no variance
    if(var(gene_counts) == 0) next
    
    # Identify outliers using multiple methods
    outlier_info <- identify_outliers(gene_counts, sample_data, outlier_threshold)
    
    # DEBUG: Print details for genes we know have outliers
    if(gene %in% c("AAEL020754", "AAEL010235", "AAEL000793")) {
      cat("DEBUG - Gene:", gene, "\n")
      cat("Counts:", paste(round(gene_counts, 2), collapse = ", "), "\n")
      
      # Test simple outlier detection manually
      max_val <- max(gene_counts)
      median_val <- median(gene_counts)
      ratio <- max_val / median_val
      cat("Max value:", max_val, "Median:", median_val, "Ratio:", round(ratio, 1), "\n")
      
      outlier_info <- identify_outliers(gene_counts, sample_data, outlier_threshold)
      cat("Has outliers detected:", outlier_info$has_outliers, "\n")
      if(outlier_info$has_outliers) {
        cat("Outlier samples:", paste(outlier_info$outlier_samples, collapse = ", "), "\n")
        cat("Outlier values:", paste(round(outlier_info$outlier_values, 2), collapse = ", "), "\n")
      }
      cat("---\n")
    }
    
    if(outlier_info$has_outliers) {
      # Check if DEG signal depends on outliers
      artifact_assessment <- assess_outlier_dependence_fixed(
        gene_counts, sample_data, outlier_info, min_samples_consistent, group_col
      )
      
      # Get result statistics
      gene_result <- results_table[gene, ]
      if("log2FoldChange" %in% colnames(results_table)) {
        lfc <- gene_result$log2FoldChange
      } else if("logFC" %in% colnames(results_table)) {
        lfc <- gene_result$logFC
      } else {
        lfc <- NA
      }
      
      if("padj" %in% colnames(results_table)) {
        padj_val <- gene_result$padj
      } else if("adj.P.Val" %in% colnames(results_table)) {
        padj_val <- gene_result$adj.P.Val
      } else {
        padj_val <- NA
      }
      
      outlier_flags <- rbind(outlier_flags, data.frame(
        gene = gene,
        outlier_type = outlier_info$outlier_type,
        outlier_samples = paste(outlier_info$outlier_samples, collapse = ";"),
        outlier_values = paste(round(outlier_info$outlier_values, 2), collapse = ";"),
        non_outlier_range = paste(round(range(outlier_info$non_outlier_values), 2), collapse = "-"),
        group_pattern = artifact_assessment$group_pattern,
        likely_artifact = artifact_assessment$likely_artifact,
        log2FoldChange = lfc,
        padj = padj_val,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  cat("\nSUMMARY:\n")
  cat("- Total significant DEGs analyzed:", length(common_genes), "\n")
  cat("- DEGs with outliers detected:", nrow(outlier_flags), "\n")
  cat("- Likely artifacts (single outlier driven):", sum(outlier_flags$likely_artifact), "\n")
  cat("- Robust DEGs (consistent pattern):", sum(!outlier_flags$likely_artifact), "\n\n")
  
  # Create comprehensive results table with outlier analysis
  cat("Creating comprehensive results table...\n")
  
  # Start with the original results table for common genes
  comprehensive_results <- results_table[common_genes, , drop = FALSE]
  
  # Add outlier analysis columns
  comprehensive_results$has_outliers <- FALSE
  comprehensive_results$outlier_type <- NA
  comprehensive_results$outlier_samples <- NA
  comprehensive_results$outlier_values <- NA
  comprehensive_results$non_outlier_range <- NA
  comprehensive_results$group_pattern <- NA
  comprehensive_results$likely_artifact <- FALSE
  comprehensive_results$outlier_flag <- "Clean"
  
  # Fill in outlier information for genes that have outliers
  if(nrow(outlier_flags) > 0) {
    for(i in 1:nrow(outlier_flags)) {
      gene <- outlier_flags$gene[i]
      if(gene %in% rownames(comprehensive_results)) {
        comprehensive_results[gene, "has_outliers"] <- TRUE
        comprehensive_results[gene, "outlier_type"] <- outlier_flags$outlier_type[i]
        comprehensive_results[gene, "outlier_samples"] <- outlier_flags$outlier_samples[i]
        comprehensive_results[gene, "outlier_values"] <- outlier_flags$outlier_values[i]
        comprehensive_results[gene, "non_outlier_range"] <- outlier_flags$non_outlier_range[i]
        comprehensive_results[gene, "group_pattern"] <- outlier_flags$group_pattern[i]
        comprehensive_results[gene, "likely_artifact"] <- outlier_flags$likely_artifact[i]
        comprehensive_results[gene, "outlier_flag"] <- ifelse(outlier_flags$likely_artifact[i], "Artifact", "Outliers_Present")
      }
    }
  }
  
  # Sort by outlier flag (artifacts first) then by padj
  if("padj" %in% colnames(comprehensive_results)) {
    comprehensive_results <- comprehensive_results[order(comprehensive_results$outlier_flag != "Clean", 
                                                         comprehensive_results$likely_artifact, 
                                                         comprehensive_results$padj), ]
  } else if("adj.P.Val" %in% colnames(comprehensive_results)) {
    comprehensive_results <- comprehensive_results[order(comprehensive_results$outlier_flag != "Clean", 
                                                         comprehensive_results$likely_artifact, 
                                                         comprehensive_results$adj.P.Val), ]
  }
  
  # Print summary table
  cat("\nOutlier Analysis Summary by Gene:\n")
  summary_table <- table(comprehensive_results$outlier_flag)
  print(summary_table)
  
  cat("\nTop 10 genes (artifacts listed first):\n")
  if("padj" %in% colnames(comprehensive_results)) {
    print_cols <- c("outlier_flag", "likely_artifact", "log2FoldChange", "padj", "outlier_samples")
  } else if("adj.P.Val" %in% colnames(comprehensive_results)) {
    print_cols <- c("outlier_flag", "likely_artifact", "logFC", "adj.P.Val", "outlier_samples")
  } else {
    print_cols <- c("outlier_flag", "likely_artifact", "outlier_samples")
  }
  
  available_cols <- intersect(print_cols, colnames(comprehensive_results))
  print(head(comprehensive_results[, available_cols, drop = FALSE], 10))
  
  # Sort by likelihood of being artifact
  outlier_flags <- outlier_flags[order(-outlier_flags$likely_artifact, outlier_flags$padj), ]
  
  return(list(
    outlier_flags = outlier_flags,
    comprehensive_results = comprehensive_results,
    summary = list(
      total_degs = length(common_genes),
      degs_with_outliers = nrow(outlier_flags),
      likely_artifacts = sum(outlier_flags$likely_artifact),
      robust_degs = sum(!outlier_flags$likely_artifact),
      clean_degs = sum(comprehensive_results$outlier_flag == "Clean")
    )
  ))
}

# Fixed assess function with correct signature and better artifact logic
assess_outlier_dependence_fixed <- function(counts, sample_data, outlier_info, min_samples, group_col) {
  outlier_samples <- outlier_info$outlier_samples
  outlier_values <- outlier_info$outlier_values
  non_outlier_values <- outlier_info$non_outlier_values
  
  # Calculate outlier magnitude
  outlier_median <- median(outlier_values)
  non_outlier_median <- median(non_outlier_values)
  non_outlier_max <- max(non_outlier_values)
  
  # Calculate ratios to assess severity
  if(non_outlier_median > 0) {
    median_ratio <- outlier_median / non_outlier_median
  } else {
    median_ratio <- Inf
  }
  
  if(non_outlier_max > 0) {
    max_ratio <- outlier_median / non_outlier_max
  } else {
    max_ratio <- Inf
  }
  
  if(!is.null(group_col) && group_col %in% colnames(sample_data)) {
    groups <- sample_data[[group_col]]
    names(groups) <- rownames(sample_data)
    
    # Check group distribution of outliers
    outlier_groups <- groups[outlier_samples]
    outlier_group_table <- table(outlier_groups)
    
    # Check if outliers are concentrated in one group
    max_outliers_per_group <- max(outlier_group_table)
    
    # Check if non-outlier samples show consistent pattern
    non_outlier_samples <- setdiff(names(counts), outlier_samples)
    non_outlier_groups <- groups[non_outlier_samples]
    
    # Calculate group sizes without outliers
    group_sizes <- table(non_outlier_groups)
    min_group_size <- ifelse(length(group_sizes) > 0, min(group_sizes), 0)
    
    # IMPROVED artifact criteria - consider magnitude AND pattern
    likely_artifact <- (
      length(outlier_samples) <= 2 &&  # Very few outliers
      min_group_size >= min_samples &&  # Other groups have enough samples
      (
        # Extreme magnitude outliers (regardless of group distribution)
        outlier_median > 10000 ||  # Very extreme absolute values
        median_ratio > 100 ||       # Extreme ratios
        # OR traditional pattern-based detection for moderate outliers
        (max_outliers_per_group == length(outlier_samples) && 
         (median_ratio > 50 || (max_ratio > 10 && outlier_median > 1000)))
      )
    )
    
    group_pattern <- paste("Outliers in:", paste(names(outlier_group_table), "(", outlier_group_table, ")", collapse=", "))
    
  } else {
    # No group information available - be more conservative
    likely_artifact <- (
      length(outlier_samples) <= 1 &&
      (median_ratio > 100 || outlier_median > 5000)  # Very high thresholds without group info
    )
    group_pattern <- "No group information available"
  }
  
  return(list(
    likely_artifact = likely_artifact,
    group_pattern = group_pattern
  ))
}

identify_outliers <- function(counts, sample_data, threshold) {
  # Balanced approach: catch real technical artifacts while preserving biological differences
  
  # Calculate simple statistics
  median_val <- median(counts)
  mean_val <- mean(counts)
  max_val <- max(counts)
  
  # Method 1: Moderate ratio test - >20x median (between original 5x and ultra-conservative 1000x)
  simple_outliers <- which(counts > (median_val * 20))
  
  # Method 2: If median is 0 (many zeros), use >10x mean
  if(median_val == 0 && mean_val > 0) {
    simple_outliers <- which(counts > (mean_val * 10))
  }
  
  # Method 3: Absolute threshold - anything >1000x the minimum non-zero value
  non_zero_vals <- counts[counts > 0]
  if(length(non_zero_vals) > 1) {
    min_nonzero <- min(non_zero_vals)
    extreme_outliers <- which(counts > (min_nonzero * 1000))
    simple_outliers <- unique(c(simple_outliers, extreme_outliers))
  }
  
  # Method 4: More conservative IQR method - use 3*IQR instead of 1.5*IQR
  q1 <- quantile(counts, 0.25)
  q3 <- quantile(counts, 0.75)
  iqr <- q3 - q1
  if(iqr > 0) {
    iqr_outliers <- which(counts > (q3 + 3 * iqr))
    simple_outliers <- unique(c(simple_outliers, iqr_outliers))
  }
  
  has_outliers <- length(simple_outliers) > 0
  
  if(has_outliers) {
    outlier_samples <- names(counts)[simple_outliers]
    outlier_values <- counts[simple_outliers]
    non_outlier_values <- counts[-simple_outliers]
    
    # Determine outlier type
    outlier_median <- median(outlier_values)
    non_outlier_median <- median(non_outlier_values)
    
    if(outlier_median > non_outlier_median * 2) {
      outlier_type <- "high_expression"
    } else if(outlier_median < non_outlier_median / 2) {
      outlier_type <- "low_expression"
    } else {
      outlier_type <- "moderate"
    }
    
    return(list(
      has_outliers = TRUE,
      outlier_samples = outlier_samples,
      outlier_values = outlier_values,
      non_outlier_values = non_outlier_values,
      outlier_type = outlier_type
    ))
  } else {
    return(list(has_outliers = FALSE))
  }
}

cat("Script loaded successfully. Use:\n")
cat("outlier_analysis <- analyze_differentialabundance_outliers(\n")
cat("  results_file = 'your_results.tsv',\n")
cat("  normalized_counts_file = 'your_counts.tsv',\n")
cat("  sample_sheet_file = 'your_samplesheet.tsv'\n")
cat(")\n")