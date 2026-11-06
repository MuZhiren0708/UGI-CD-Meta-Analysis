# Prevalence of Montreal L4 Involvement in Crohn's Disease

Study-level data and reproducible R code for **Prevalence of Montreal L4 Involvement in Crohn's Disease: A Meta-Analysis with Exploratory Machine Learning**.

This revised release supersedes the earlier extraction and analysis files in this repository. The input contains **50 source reports and 55 original cohort records**. After outcome reconciliation and removal or replacement of repeated/incompatible estimates, **41 estimates from 37 reports** enter the primary synthesis: **17,812 assessed patients and 2,921 cases**. These records form **33 conservative source-family validation blocks**. The 55 displayed records are not 55 independent cohorts.

The outcome is broad, study-reported Montreal L4 involvement, which may include proximal small-bowel disease. It is not uniformly verified strict esophagogastroduodenal involvement. Available anatomical subsets and diagnostic-evidence restrictions are analysed separately. All shared records are aggregate study/cohort data, not individual patient records.

## Main results

- Inverse-variance, logit-transformed DerSimonian-Laird pooled prevalence: **12.78% (95% CI 10.30–15.74%)**.
- ML uses source-reported age summary, male percentage and log CD denominator, with **10 repetitions of nested grouped five-fold cross-validation**. Related publications/source families remain in the same fold. Imputation and scaling are fitted within training folds.
- Predictive R²: random forest **0.191**, XGBoost **−0.059**, LASSO **−0.056**, REML meta-regression **0.234**, and the L3-extended LASSO **−0.145**. Negative values are retained. RMSE and MAE are in percentage points. These are continuous cohort-level predictions; no clinical classification threshold, individual lifetime trajectory or external validation is claimed.
- The supplied NOS item allocations for Tarrant (2008) sum to **8**. That total now enters the NOS ≥7 subgroup: **19 estimates; 12,739 patients; 2,051 cases; 12.51% (95% CI 9.31–16.61%)**. This reconciles the supplied star allocations; it is not a new independent quality assessment.

## Files

| Path | Contents |
|---|---|
| `data/table1_55_rows.csv`, `data/all_source_rows.csv` | All original cohort records, inclusion roles, source locators and reconciled totals |
| `data/extracted_data.csv` | Alias of the current 55-row dataset, replacing the obsolete extraction at this path |
| `data/analysis_data.csv` | The 41 estimates used in the primary synthesis and ML analysis |
| `data/data_dictionary.csv` | Definitions, units and missing-value interpretation for all 55 input fields |
| `data/source_audit_50.csv`, `data/source_membership.csv` | Report-level audit and inclusion/replacement decisions |
| `data/NOS_item_scores.csv`, `data/NOS_arithmetic_audit.json` | Supplied item stars and verified arithmetic totals for all 55 rows |
| `data/anatomic_phenotype_outcomes.csv` | Recoverable isolated, concomitant, EGD and proximal small-bowel outcomes |
| `analysis_results/` | Actual saved OOF predictions, fold assignments, hyperparameters, fitted models, meta-analysis results and independent metric audit |
| `ML_demographic_sensitivity/` | Saved analysis restricted by demographic-denominator alignment |
| `code/` | Both original-named entry scripts, analysis helpers, revised figure styling and the reproduction entry point |
| `figure_data/` | Saved inputs for the revised layouts of Figures 7–11 |
| `figures/SVG/` | R-generated main Figures 2–11 and Supplemental Figures 1–12 (22 SVGs) |
| `docs/Figure_style_references.md` | Sources used to inform the revised visual layouts |
| `validation/` | Release checks, the small data-change log and file checksums |
| `online_resources/` | Updated Online Resource 2 (data/results) and Online Resource 3 (scripts/models) ZIPs |

Empty values, `NA` and `NR` mean unavailable information; they do not imply zero. Original diagnostic labels are retained as descriptive labels, rather than interpreted as proof of mutually exclusive or uniformly objective tests. PDF identifiers, page locators and hashes are supplied, but source-paper PDFs are not redistributed.

## Reproduce from the saved analysis

Use R 4.5 or later with these packages installed:

```r
install.packages(c("meta", "metafor", "glmnet", "ranger", "xgboost",
                  "ggplot2", "patchwork", "svglite", "jsonlite", "cluster", "scales"))
```

Run commands from the repository root:

```text
Rscript code/run_analysis.R reproduced_results saved
```

This copies the saved ML outputs, recomputes the meta-analysis and subgroup results from the corrected inputs, regenerates all 22 SVGs, applies the newer layouts to Figures 7–11, and independently checks saved OOF metrics and fold separation. Results are written to `reproduced_results/`; committed snapshots remain in their original folders. Times New Roman should be installed to match the published artwork. The standard entry point retains editable SVG text and requires no Python or Node.js.

To refit the complete frozen ML pipeline, including the demographic-denominator sensitivity analysis:

```text
Rscript code/run_analysis.R refitted_results full
```

The seed and tuning grids are in `code/UGI_CD_Meta_Analysis.R`; saved R package versions are in `analysis_results/ML_sessionInfo.txt`. Exact refitting requires compatible package versions. Model parameters and fitted-model RDS files are included. An optional extra R library path can be set with `UGI_R_LIBRARY`.

The original entry names remain usable, with one shared output directory:

```text
Rscript code/UGI_CD_Meta_Analysis.R refitted_results all
Rscript code/Untitled1.R refitted_results
```

The wrapper is the recommended route because it also applies the latest Figure 7–11 styling. Historical Figure 1 search-flow records cannot be reconstructed from these aggregate inputs and are not generated by the scripts.

## This release was checked

The `saved` workflow completed successfully on 2026-11-06 using R 4.5.3. It regenerated 22 SVGs, reproduced the corrected NOS subgroup, retained 55 displayed forest-plot records, and independently recomputed OOF metrics to numerical precision. The saved ML models were not refitted solely for this repository synchronization: outcome counts, modelling predictors, fold assignments and cohort membership did not change. The four NOS input fields and the affected quality-subgroup outputs were synchronized with the later supplied supplemental star allocations.

This release preserves remaining source limitations, including heterogeneous anatomical definitions, incomplete diagnostic detail, some ambiguous coexistence definitions, possible source overlap and differing demographic assessment times. It does not resolve those limitations by inventing observations.

The repository retains its existing MIT licence. Please cite the accompanying manuscript when using these data or scripts; the final publication citation will be added when available.
