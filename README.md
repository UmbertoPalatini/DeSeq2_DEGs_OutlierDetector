# Post-DESeq2 Outlier-Driven DEG Detection

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A post-hoc analysis tool that identifies differentially expressed genes (DEGs) whose statistical significance may be driven by single-sample outliers, regardless of whether these outliers are techical errors or real biological "chaos", rather than consistent biological differences. This approach helps distinguish between technical artifacts and legitimate biological variation in RNA-seq differential expression results.

## Background

While DESeq2 and edgeR are powerful tools for differential expression analysis, they can produce inflated false positive rates when outlier samples violate distributional assumptions. Rather than applying aggressive pre-filtering that might remove legitimate biological variation, this script performs targeted post-analysis to identify genes where significance depends on extreme outliers.
Some references:
- PMID: 33709073
- PMID: 40926263
- PMID: 25516281
- PMID: 39478636

## Key Features

- **Post-hoc analysis**: Works with existing DESeq2/edgeR results without re-running differential expression
- **Balanced outlier detection**: Uses multiple methods tuned to catch technical artifacts while preserving biological variation
- **Group-aware assessment**: Considers experimental design when evaluating outlier patterns
- **Comprehensive output**: Provides both flagged outliers and complete annotated results table
- **nf-core compatible**: Designed to work seamlessly with nf-core/differentialabundance pipeline output

## Installation

### Dependencies

```r
# Required packages
install.packages(c("dplyr"))

# DESeq2 for data handling (optional, only needed for column names)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
```

### Download

```bash
# Download the script
wget https://raw.githubusercontent.com/yourusername/repo/main/DEGs_outlier_detection_postDeSeq_v5_final.r

# Or clone the repository
git clone https://github.com/yourusername/outlier-deg-detection.git
```

## Quick Start

```r
# Load the script
source("DEGs_outlier_detection_postDeSeq_v5_final.r")

# Run analysis on your data
outlier_analysis <- analyze_differentialabundance_outliers(
  results_file = "deseq2.results.tsv",
  normalized_counts_file = "deseq2.normalised_counts.tsv",
  sample_sheet_file = "samplesheet.tsv"
)

# View results
View(outlier_analysis$comprehensive_results)

# Save annotated results
write.table(outlier_analysis$comprehensive_results, 
           "deg_results_with_outlier_analysis.tsv", 
           sep = "\t", quote = FALSE, row.names = TRUE)
```

## Input Requirements

### Required Files

| File | Description | Expected Format |
|------|-------------|-----------------|
| **Results Table** | DESeq2/edgeR differential expression results | TSV with `padj` or `adj.P.Val` column |
| **Normalized Counts** | Library-size normalized expression matrix | TSV, genes as rows, samples as columns |
| **Sample Sheet** | Sample metadata with group assignments | TSV with sample IDs and condition/group column |

### File Formats

**Results Table** (from DESeq2/edgeR):
```
gene_id         baseMean    log2FoldChange  padj
GENE001         1234.5      2.1            0.001
GENE002         567.8       -1.8           0.05
```

**Normalized Counts**:
```
gene_id         Sample1     Sample2     Sample3     Sample4
GENE001         120.5       89.2        156.1       203.7
GENE002         45.2        67.1        23.8        89.4
```

**Sample Sheet**:
```
sample_id       condition   batch
Sample1         control     1
Sample2         control     1
Sample3         treatment   2
Sample4         treatment   2
```

## Algorithm Overview

### Outlier Detection Methods

The script uses multiple complementary approaches to identify outliers:

1. **Ratio-based Detection**: Flags values >20x median (balanced threshold)
2. **Zero-heavy Data Handling**: Uses 10x mean for datasets with many zeros
3. **Extreme Value Detection**: Flags values >1000x minimum non-zero value
4. **Conservative IQR**: Uses 3×IQR threshold (instead of standard 1.5×IQR)

**Threshold Philosophy**: The detection methods are calibrated to catch extreme technical artifacts (e.g., single samples with 100-1000x higher expression than background) while preserving legitimate biological differences between experimental groups (typically 2-20x fold changes).

### Artifact Assessment

Genes are classified as likely artifacts when **ALL** of the following conditions are met:

**Required conditions:**
- ≤2 outlier samples total
- Sufficient remaining samples in each group (≥ `min_samples_consistent`)
- All outliers concentrated in a single experimental group

**AND at least ONE magnitude criterion:**
- Extreme absolute values (>10,000 normalized counts)
- Extreme ratio (>100x median of all samples)
- High ratio (>50x median) OR (>10x max non-outlier value AND >1000 absolute)

## Output Structure

### Main Results: `comprehensive_results`

Enhanced version of your original results table with added outlier analysis columns:

| Column | Description |
|--------|-------------|
| `outlier_flag` | **"Clean"**, **"Outliers_Present"**, or **"Artifact"** |
| `likely_artifact` | Boolean flag for probable false positives |
| `has_outliers` | Whether outliers were detected |
| `outlier_type` | Type of outlier: "high_expression", "low_expression", or "moderate" |
| `outlier_samples` | Which samples are outliers (semicolon-separated) |
| `outlier_values` | The actual outlier expression values (semicolon-separated) |
| `non_outlier_range` | Min-max range of non-outlier expression values |
| `group_pattern` | Distribution of outliers across experimental groups |

### Detailed Analysis: `outlier_flags`

Subset containing only genes with detected outliers, with full diagnostic information.

### Summary Statistics: `summary`

Count of genes in each category for quick assessment.

## Interpretation Guide

### Gene Categories

**Clean Genes** (`outlier_flag = "Clean"`)
- No outliers detected
- High confidence in differential expression call
- Safe to include in downstream analysis

**Outliers Present** (`outlier_flag = "Outliers_Present"`)
- Outliers detected but likely biological variation
- Examples: Moderate fold-changes (5-20x) with consistent patterns
- Generally safe to include, consider manual review

**Likely Artifacts** (`outlier_flag = "Artifact"`)
- Extreme outliers driving statistical significance
- Examples: Single samples with 100x+ expression vs background
- **Recommend exclusion** from downstream analysis

### Example Cases

**Clear Artifact**:
```
Gene: EXAMPLE001
Values: [0, 0, 0, 0, 8844, 0, 0, 0, 0, 1, 0, 0]
Classification: Artifact (single 8844x outlier vs ~0.5 background)
```

**Biological Variation**:
```
Gene: EXAMPLE002  
Values: [15210, 7255, 7864, 3865, 3054, 2163, 2779, 2284, 3531, 1849, 518, 2263, 1628]
Classification: Clean (5.5x ratio - legitimate group differences)
```

**Robust Signal**:
```
Gene: EXAMPLE003
Values: [120, 110, 135, 125, 450, 478, 523, 467, 445, 489, 501, 456]
Classification: Clean (consistent group differences)
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `padj_threshold` | 0.05 | Significance threshold for DEG inclusion |
| `outlier_threshold` | 3 | Legacy parameter (passed to functions but not used in current implementation) |
| `min_samples_consistent` | 2 | Minimum non-outlier samples required per group |

## Usage Examples

### Basic Analysis

```r
# Standard analysis
results <- analyze_differentialabundance_outliers(
  results_file = "my_deseq_results.tsv",
  normalized_counts_file = "my_normalized_counts.tsv",
  sample_sheet_file = "my_samples.tsv"
)

# Check summary
print(results$summary)
#> $total_degs: 150
#> $clean_degs: 120
#> $likely_artifacts: 8
#> $robust_degs: 22
```

### Filter Results for Clean Gene List

```r
# Get high-confidence DEGs (exclude artifacts)
clean_degs <- results$comprehensive_results[
  results$comprehensive_results$outlier_flag != "Artifact", 
]

# Or be more conservative (only clean genes)
very_clean_degs <- results$comprehensive_results[
  results$comprehensive_results$outlier_flag == "Clean", 
]

# Save filtered results
write.table(clean_degs, "filtered_deg_results.tsv", 
           sep = "\t", quote = FALSE, row.names = TRUE)
```

### Manual Review of Flagged Genes

```r
# View only flagged artifacts
artifacts <- results$outlier_flags[results$outlier_flags$likely_artifact, ]
View(artifacts)

# Export for manual review
write.table(artifacts, "genes_flagged_as_artifacts.tsv", 
           sep = "\t", quote = FALSE, row.names = FALSE)
```

### Test Outlier Detection on Custom Data

```r
# Test the outlier detection function directly
test_counts <- c(15210.78, 7255.144, 7863.838, 3865.132, 3054.153, 
                 2162.5, 2779.326, 2284.098, 3530.547, 1848.999, 
                 517.8544, 2262.534, 1627.743)
names(test_counts) <- paste0("Sample", 1:13)

result <- identify_outliers(test_counts, NULL, 3)
print(result$has_outliers)  # Should be FALSE for biological variation
```

## Workflow Integration

### With nf-core/differentialabundance

```bash
# After running nf-core/differentialabundance
nextflow run nf-core/differentialabundance --input samplesheet.csv --contrasts contrasts.csv
```

```r
# Point to nf-core outputs
results <- analyze_differentialabundance_outliers(
  results_file = "results/deseq2/condition_treatment_vs_control/deseq2.condition_treatment_vs_control.results.tsv",
  normalized_counts_file = "results/deseq2/normalised_counts/deseq2.normalised_counts.tsv",
  sample_sheet_file = "samplesheet.csv"  # Your original input
)
```

### Integration with Standard DESeq2 Workflow

```r
# After standard DESeq2 analysis
dds <- DESeq(dds)
res <- results(dds)

# Save intermediate files for outlier analysis
write.table(res, "deseq2_results.tsv", sep = "\t", quote = FALSE)
write.table(counts(dds, normalized = TRUE), "normalized_counts.tsv", sep = "\t", quote = FALSE)

# Run outlier detection
outlier_results <- analyze_differentialabundance_outliers(
  results_file = "deseq2_results.tsv",
  normalized_counts_file = "normalized_counts.tsv", 
  sample_sheet_file = "sample_metadata.tsv"
)
```

## Tuning and Validation

### Threshold Sensitivity

The current thresholds are designed to be balanced:
- **20x median**: Catches extreme single-sample spikes while preserving group differences
- **3×IQR**: Conservative IQR detection to avoid false positives
- **1000x minimum**: Extreme ratio detection for very sparse data

### Validation Approach

1. **Visual inspection**: Examine flagged vs clean genes manually
2. **Biological validation**: Check if flagged genes make biological sense
3. **Threshold adjustment**: Modify parameters based on your data characteristics

### Custom Threshold Example

```r
# For more sensitive detection (catches more outliers)
# Modify the identify_outliers function thresholds:
# - Change 20x to 15x median
# - Change 3*IQR to 2*IQR

# For more conservative detection (catches fewer outliers)  
# - Change 20x to 50x median
# - Change 3*IQR to 5*IQR
```

## Limitations

- **Small sample sizes**: Less effective with <4 samples per group
- **Legitimate extreme biology**: May flag genuine biological extremes (rare but possible)
- **Count data optimized**: Designed specifically for RNA-seq count matrices
- **Manual review recommended**: Automated flags should be verified by domain experts
- **Threshold dependency**: Performance depends on appropriate threshold selection for your data type

## Performance

- **Runtime**: ~1-2 seconds per 1000 genes analyzed
- **Memory**: Minimal additional memory beyond input data
- **Scalability**: Tested with datasets up to 50,000 genes and 200 samples

## Troubleshooting

### Common Issues

**"No matching samples" error**:
- Check sample name consistency between files
- Ensure no extra spaces or special characters
- Verify column headers match expected format

**"Cannot find adjusted p-value column"**:
- Ensure results file contains `padj` or `adj.P.Val` column
- Check that results file is properly formatted

**Too many/few outliers detected**:
- Algorithm thresholds may need adjustment for your data
- Consider manual inspection of borderline cases
- Review threshold sensitivity section above

### Getting Help

1. **Check file formats**: Ensure input files match expected structure
2. **Verify sample matching**: Run diagnostic checks on sample names
3. **Test with known cases**: Use the custom data testing example above
4. **Manual validation**: Examine a few flagged cases to verify behavior

## Contributing

Issues and pull requests welcome! Please ensure:

- Test with sample datasets
- Document any parameter changes
- Include examples for new features
- Update README for significant changes

## Citation

If you use this script in your research, please cite this repository and consider citing relevant differential expression methodology papers that inform the approach.

## License

MIT License - see LICENSE file for details.