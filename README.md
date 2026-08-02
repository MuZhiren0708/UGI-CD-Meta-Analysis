# Prevalence of Upper Gastrointestinal Involvement in Crohn's Disease: A Systematic Review and Meta-Analysis

## Overview

This repository contains the study-level dataset and analysis code for a systematic review and machine learning–enhanced meta-analysis estimating the global pooled prevalence of upper gastrointestinal (UGI) involvement in Crohn's disease (CD).

## Study Details

- **Design:** Systematic review and meta-analysis (PRISMA 2020 compliant)
- **Registration:** PROSPERO (CRD420261339554)
- **Studies included:** 55 independent cohorts
- **Total patients:** 22,465
- **Countries:** 25
- **Pooled prevalence:** 11.6% (95% CI: 9.4%–14.2%)

## Repository Structure

```
UGI-CD-Meta-Analysis/
├── data/
│   └── extracted_data.csv      # Study-level extracted data (55 cohorts)
├── code/
│   └── UGI_CD_Meta_Analysis.R  # Complete R analysis code (traditional meta + ML enhanced meta)
├── README.md
└── LICENSE
```

## Data Description

The dataset (`data/extracted_data.csv`) contains 55 rows (one per cohort) with the following columns:

| Column | Description |
|--------|-------------|
| Study | First author and year |
| Country / Region | Country or region of study |
| Study design | Retrospective cohort, Prospective cohort, Cross-sectional, or Case-control |
| Study Period | Years of data collection |
| Mean age | Mean age of participants (years) |
| Male | Male proportion (%) |
| Diagnostic Modality for UGI | Method used to detect UGI involvement |
| CD Patients | Total number of Crohn's disease patients |
| UGI Cases | Number of patients with UGI involvement |
| Prevalence of UGI Phenotype | UGI prevalence (%) |
| NOS Score | Newcastle-Ottawa Scale quality score (0–9) |

## Analysis Components

### 1. Traditional Meta-Analysis
- Random-effects model (DerSimonian–Laird estimator)
- Subgroup analyses (age group, geographic region, diagnostic modality, study design, study quality)
- REML multivariate meta-regression
- Publication bias assessment (funnel plot, Egger's test, trim-and-fill)
- Leave-one-out sensitivity analysis

### 2. Machine Learning Enhanced Meta-Analysis
- Four models compared: LASSO regression, Random Forest, XGBoost, conventional meta-regression
- 10-fold cross-validation with grid search hyperparameter optimization
- Independent hold-out validation (20% of studies)
- SHAP (SHapley Additive exPlanations) analysis for heterogeneity attribution
- K-means clustering for unsupervised risk stratification

## Requirements

### R (All Analysis)
- R >= 4.0
- Packages: `ggplot2`, `patchwork`, `dplyr`, `svglite`, `ggbeeswarm`, `tidyverse`, `xgboost`, `caret`, `countrycode`, `shapviz`

## Citation

If you use this data or code, please cite:

> [Citation to be added upon publication]

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions regarding the data or analysis, please open an issue in this repository.
