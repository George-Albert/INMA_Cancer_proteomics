# INMA Cancer Proteomics

This repository contains advanced technical work in bioinformatics and differential expression analysis developed in 2024 in collaboration with Maria Sancho Albero within the NAIVE project (Catedra Samca).

## Project scope

The analysis focuses on preprocessing, quality control, statistical modeling, and visualization of high-throughput proteomics and nanoparticle-protein interaction data using R.

## Current structure

- src: reusable functions for IO, preprocessing, statistics, and visualization.
- config: pipeline parameters and contrasts.
- scripts: production pipeline scripts.
- scripts/legacy: original scripts preserved for traceability.
- notebooks/exploratory: ad hoc analysis notebooks.
- tests/testthat: basic unit tests.
- Inputs and Outputs: data and generated artifacts from the analysis.

## Pipeline order

1. scripts/001_data_formatting.R
2. scripts/002_setting_data.R
3. scripts/003_pca_analysis.R
4. scripts/004_de_analyses.R

## Notes

- Differential expression output names now use mean_gt_X to avoid invalid Windows file names.
- The legacy scripts are kept unchanged under scripts/legacy for reproducibility and comparison.
