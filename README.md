# INMA Cancer Proteomics

This repository contains advanced technical work in bioinformatics and differential expression analysis developed in 2024 in collaboration with Maria Sancho Albero within the NAIVE project (Catedra Samca).

## Project scope

The analysis focuses on preprocessing, quality control, statistical modeling, and visualization of high-throughput proteomics and nanoparticle-protein interaction data using R.

## Current structure

- Codes: active scripts with the latest corrected/improved versions.
- src: reusable functions for IO, preprocessing, statistics, and visualization.
- config: pipeline parameters and contrasts.
- scripts: production pipeline scripts.
- scripts/legacy: archived historical scripts preserved for traceability.
- notebooks/exploratory: ad hoc analysis notebooks.
- tests/testthat: basic unit tests.
- Inputs and Outputs: data and generated artifacts from the analysis.

## Pipeline order

Main execution order in Codes:

1. Codes/001_Data_formatting.R
2. Codes/002_Setting_data.R
3. Codes/003_PCA_analysis.R
4. Codes/004_DE_Analyses.R
5. Codes/005_Expression_levels.R
6. Codes/006_Volcano_plots.R
7. Codes/007_Data_formatting_new_data.R
8. Codes/008_Setting_new_data.R
9. Codes/009_PCA_analysis.R
10. Codes/010_DE_Analyses_added.R
11. Codes/011_Venn_Diagrams_added.R
12. Codes/012_Bar_plots_per_contrasts.R
13. Codes/013_DE_Analyses_reg_out.R
14. Codes/014_PCA_analysis_after_reg_out.R

## Notes

- Differential expression output names now use mean_gt_X to avoid invalid Windows file names.
- Codes contains the current maintained versions, while scripts/legacy keeps historical snapshots.
