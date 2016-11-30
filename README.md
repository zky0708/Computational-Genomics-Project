# Systematic Analysis of Tumor Purity
This project presents a systematic analysis of four tumor-purity-estimating methods on about 7500 patients across 15 cancer types from The Cancer Genome Atlas (TCGA). 

1. purity_estimate.R -- Preprocess input data; get tumor estimates of each sample using four different methods.
2. comparison.R -- Compare the four estimates by computing Pearson correlations; plot figures used in report.
3. clinical_associations.R -- Examine associations of purity levels with clinical features; survival analysis; evaluate performance of each method; plot figures used in report.
