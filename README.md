# Systematic Analysis of Tumor Purity
This is a course project for COMS W4761-- Computational Genomics, systematically analyzing and comparing four popular tumor-purity-estimating methods using genomic datasets of around 7500 patients across 15 cancer types from The Cancer Genome Atlas (TCGA). 

1. purity_estimate.R -- Preprocess input data; get tumor estimates of each sample using four different methods.
2. comparison.R -- Compare the four estimates by computing Pearson correlations; plot figures used in report.
3. clinical_associations.R -- Examine associations of purity levels with clinical features; survival analysis; evaluate performance of each method; plot figures used in report.
