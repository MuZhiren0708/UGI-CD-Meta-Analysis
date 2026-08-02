# =============================================================================
# Traditional Meta-Analysis: UGI Involvement in Crohn's Disease
# =============================================================================
# This script performs the traditional meta-analysis components:
#   - Pooled prevalence (random-effects, DerSimonian-Laird)
#   - Subgroup analyses (age, region, modality, design, quality)
#   - Multivariate meta-regression (REML)
#   - Publication bias (funnel plot, Egger's test, trim-and-fill)
#   - Leave-one-out sensitivity analysis
# =============================================================================

# --- Load libraries ---
library(meta)
library(metafor)
library(dplyr)

# --- Load data ---
data <- read.csv("data/extracted_data.csv", stringsAsFactors = FALSE)

# --- Data preparation ---
# Calculate logit-transformed prevalence and standard error
data$prevalence_prop <- data$UGI.Cases / data$CD.Patients
data <- escalc(
  measure = "PLO",
  xi = UGI.Cases,
  ni = CD.Patients,
  data = data
)

# --- 1. Pooled Prevalence (Random-Effects Model) ---
meta_result <- rma(
  yi,
  vi,
  data = data,
  method = "DL",
  slab = Study
)

# Back-transform to prevalence
pooled_prev <- predict(meta_result, transf = transf.ilogit)
cat("Pooled Prevalence:", round(pooled_prev$pred * 100, 2), "%\n")
cat("95% CI:", round(pooled_prev$ci.lb * 100, 2), "% -", round(pooled_prev$ci.ub * 100, 2), "%\n")
cat("I²:", round(meta_result$I2 * 100, 1), "%\n")
cat("Q =", round(meta_result$QE, 2), ", p =", meta_result$QEp, "\n")

# --- 2. Subgroup Analyses ---
# Age group subgroup
# (Define age groups based on mean age: pediatric <=16, adult 17-59, elderly >=60)
data$age_group <- cut(data$Mean.age,
  breaks = c(-Inf, 16, 59, Inf),
  labels = c("Paediatric", "Adult", "Elderly")
)

subgroup_age <- rma(yi, vi, data = data, method = "DL", subset = age_group == "Paediatric")
cat("\n--- Paediatric ---\n")
print(predict(subgroup_age, transf = transf.ilogit))

subgroup_adult <- rma(yi, vi, data = data, method = "DL", subset = age_group == "Adult")
cat("\n--- Adult ---\n")
print(predict(subgroup_adult, transf = transf.ilogit))

# --- 3. Multivariate Meta-Regression (REML) ---
# reg_result <- rma(yi, vi,
#   mods = ~ Mean.age + Male + NOS.Score + factor(Diagnostic.Modality.for.UGI),
#   data = data,
#   method = "REML"
# )
# summary(reg_result)

# --- 4. Publication Bias ---
# Egger's test
egger_test <- regtest(meta_result)
cat("\nEgger's test: t =", round(egger_test$zval, 2), ", p =", egger_test$pval, "\n")

# Trim-and-fill
tf <- trimfill(meta_result)
cat("Trim-and-fill imputed studies:", tf$k0, "\n")

# --- 5. Leave-One-Out Sensitivity Analysis ---
loo <- leave1out(meta_result, transf = transf.ilogit)
print(loo)

# --- Forest Plot ---
# png("figures/forest_plot.png", width = 1200, height = 2000, res = 150)
# forest(meta_result, transf = transf.ilogit,
#        slab = data$Study, xlab = "Prevalence (%)",
#        main = "UGI Involvement in CD - Pooled Prevalence")
# dev.off()

# --- Funnel Plot ---
# png("figures/funnel_plot.png", width = 800, height = 600, res = 150)
# funnel(meta_result, main = "Funnel Plot")
# dev.off()
