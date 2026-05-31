# INMA Cancer Proteomics

![Language](https://img.shields.io/badge/language-R-blue)
![Platform](https://img.shields.io/badge/platform-Windows-lightgrey)
![Tests](https://img.shields.io/badge/tests-testthat%20passing-brightgreen)

Comprehensive R-based pipeline for proteomics data curation, quality control, dimensionality reduction, differential expression analysis, and downstream biological interpretation.

This repository consolidates technical work developed in 2024 and later refined to improve robustness, reproducibility, and maintainability for nanoparticle-protein interaction and cancer-related omics analyses.

## Table of Contents

- [Project goals](#project-goals)
- [Repository organization](#repository-organization)
- [Main pipeline in Codes](#main-pipeline-in-codes)
- [Recent technical improvements](#recent-technical-improvements)
- [Reproducibility](#reproducibility)
- [Scientific and collaboration context](#scientific-and-collaboration-context)
- [Selected publications relevant to this domain](#selected-publications-relevant-to-this-domain)
- [External references](#external-references)
- [How to cite this repository](#how-to-cite-this-repository)
- [License](#license)

## Project goals

- Standardize heterogeneous input tables into analysis-ready matrices.
- Build reproducible metadata and feature annotations.
- Perform PCA and differential expression analyses across multiple contrasts.
- Extend analysis to added datasets and regression-adjusted variants.
- Generate interpretable outputs (volcano plots, expression summaries, Venn diagrams, contrast bar plots).
- Maintain reusable statistical and utility modules in a shared source layer.

## Repository organization

- `Codes`: Active analysis scripts. This is the main maintained workflow.
- `src`: Reusable functions for IO, preprocessing, statistics, utilities, and visualization helpers.
- `config`: YAML configuration files for paths, contrasts, and thresholds.
- `scripts`: Additional pipeline scripts and historical variants.
- `scripts/legacy`: Archived snapshots retained for traceability only.
- `Inputs`: Source and intermediate data files used by the analyses.
- `Outputs`: Generated figures and processed analysis outputs.
- `tests`: Unit tests implemented with testthat.

## Main pipeline in Codes

Recommended execution order:

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

## Recent technical improvements

- Centralized shared helpers in src/utils.R and src/statistics.R to reduce duplicated code.
- Standardized output naming with mean_gt_X patterns for better cross-platform compatibility (including Windows).
- Added stricter alignment checks between metadata and read matrices in extended analyses.
- Cleaned inconsistent imports and removed duplicated local helper definitions.
- Hardened file matching patterns and directory handling in plotting/report scripts.
- Fixed test path resolution so testthat executes correctly from project root.

## Reproducibility

### Environment

- OS tested: Windows.
- R installation tested through R 4.6.0.
- Prefer running with `Rscript` available in PATH.
- If `Rscript` is not in PATH, configure your shell environment first instead of hardcoding user-specific absolute paths in commands.

Quick check:

```powershell
Rscript --version
```

### Run scripts

Run any step with:

```powershell
Rscript --vanilla Codes/001_Data_formatting.R
```

Replace the script path according to the pipeline order.

### Run tests

From project root:

```powershell
Rscript --vanilla -e "testthat::test_dir('tests/testthat')"
```

### Reproducibility checklist

- Run commands from repository root.
- Use the ordered pipeline in `Codes/`.
- Keep parameterization in `config/paths.yml`, `config/contrasts.yml`, and `config/thresholds.yml`.
- Keep generated artifacts under `Outputs/` and preserve input snapshots under `Inputs/`.
- Run tests after changes and record results in commit messages or PR notes.

### Traceability notes

- Active, maintained workflow: `Codes/` + `src/` + `config/`.
- Historical snapshots: `scripts/legacy/` (kept for provenance; may contain old assumptions and machine-specific paths).
- For reproducible analyses, do not execute legacy scripts unless you are explicitly reproducing historical behavior.

## Scientific and collaboration context

This repository is linked to research activity in nanomedicine and extracellular vesicle-enabled therapeutics associated with the University of Zaragoza ecosystem (INMA, CIBER-BBN, IIS Aragon) and related collaborators.

Public institutional and researcher-profile sources indicate:

- Maria Sancho-Albero is affiliated with INMA (CSIC-Universidad de Zaragoza) and UNIZAR structures.
- Her work focuses on nanomedicine, extracellular vesicles, targeted delivery, and cancer-oriented bioorthogonal or nanocarrier strategies.
- Institutional news from NANBIOSIS reports the award of an ERC Starting Grant (project SEVEN) for metastasis-targeted nanotherapy research.

Note: This repository is an analysis/codebase resource. Scientific claims and career milestones should always be interpreted from primary institutional and publication sources listed below.

## Selected publications relevant to this domain

Examples of peer-reviewed papers associated with the researcher profile and topic area:

- Cancer-derived exosomes loaded with ultrathin palladium nanosheets for targeted bioorthogonal catalysis. Nature Catalysis (2019). DOI: 10.1038/s41929-019-0333-4
- Exosome origin determines cell targeting and the transfer of therapeutic nanoparticles towards target cells. Journal of Nanobiotechnology (2019). DOI: 10.1186/s12951-018-0437-z
- Isolation of exosomes from whole blood by a new microfluidic device: proof of concept application in the diagnosis and monitoring of pancreatic cancer. Journal of Nanobiotechnology (2020). DOI: 10.1186/s12951-020-00701-7
- Transfer of photothermal nanoparticles using stem cell derived small extracellular vesicles for in vivo treatment of primary and multinodular tumours. Journal of Extracellular Vesicles (2022). DOI: 10.1002/jev2.12193
- Exosomes loaded with ultrasmall Pt nanoparticles: a novel low-toxicity alternative to cisplatin. Journal of Nanobiotechnology (2022). DOI: 10.1186/s12951-022-01675-4
- Supramolecular Nucleic Acid-Based Organosilica Nanoparticles Responsive to Physical and Biological Inputs. Journal of the American Chemical Society (2023). DOI: 10.1021/jacs.3c04345

For a current and complete list, check ORCID and bibliographic indexes.

## External references

- INMA profile (Sancho Albero, Maria): https://inma.unizar-csic.es/investigadores/sancho-albero-maria/
- NFP Research Group profile: https://nfp.unizar.es/post_teams/maria-sancho-albero/
- ORCID profile: https://orcid.org/0000-0001-8762-5457
- NANBIOSIS ERC announcement: https://www.nanbiosis.es/maria-sancho-nanbiosis-researcher-and-ana-serrano-awarded-erc-starting-grants/
- OpenAlex author entry: https://api.openalex.org/authors?filter=orcid:0000-0001-8762-5457
- OpenAlex works list: https://api.openalex.org/works?filter=authorships.author.orcid:0000-0001-8762-5457&sort=cited_by_count:desc&per-page=15

## How to cite this repository

If you use this repository, please cite it in your methods section and include the repository URL plus access date.

Suggested citation format:

```text
INMA Cancer Proteomics Repository. Comprehensive R pipeline for proteomics and differential expression analyses.
GitHub repository: https://github.com/<owner>/<repo>
Accessed: YYYY-MM-DD
```

If needed, this repository can also include a `CITATION.cff` file in a future update.

## License

See LICENSE for repository licensing terms.
