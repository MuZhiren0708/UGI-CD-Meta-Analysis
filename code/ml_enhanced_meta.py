# =============================================================================
# Machine Learning Enhanced Meta-Analysis: UGI Involvement in Crohn's Disease
# =============================================================================
# This script performs the ML-enhanced meta-analysis components:
#   - Four ML models: LASSO, Random Forest, XGBoost, conventional meta-regression
#   - 10-fold cross-validation with grid search
#   - Independent hold-out validation (20%)
#   - SHAP analysis for heterogeneity attribution
#   - K-means clustering for risk stratification
# =============================================================================

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.linear_model import Lasso
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_absolute_error
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import xgboost as xgb
import shap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Load data ---
data = pd.read_csv("data/extracted_data.csv")

# --- Feature engineering ---
# Parse mean age (handle entries like "35 (22.5)")
def parse_age(val):
    if pd.isna(val):
        return np.nan
    s = str(val).split("(")[0].strip()
    try:
        return float(s)
    except ValueError:
        return np.nan

data["Mean_age_num"] = data["Mean age"].apply(parse_age)

# Encode categorical variables
data["Study_design_enc"] = pd.Categorical(data["Study design"]).codes
data["Modality_enc"] = pd.Categorical(data["Diagnostic Modality for UGI"]).codes
data["Country_enc"] = pd.Categorical(data["Country / Region"]).codes

# Target variable
y = data["Prevalence of UGI Phenotype"].values

# Feature set
feature_cols = ["Mean_age_num", "Male", "CD Patients", "NOS Score",
                "Study_design_enc", "Modality_enc"]
X = data[feature_cols].copy()
X = X.fillna(X.median())

# --- Train-test split (20% hold-out) ---
np.random.seed(42)
n_studies = len(data)
holdout_idx = np.random.choice(n_studies, size=int(0.2 * n_studies), replace=False)
train_mask = np.array([i not in holdout_idx for i in range(n_studies)])

X_train, X_test = X.values[train_mask], X.values[~train_mask]
y_train, y_test = y[train_mask], y[~train_mask]

scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# --- Model definitions ---
models = {
    "LASSO": {
        "model": Lasso(max_iter=10000),
        "params": {"alpha": [0.001, 0.01, 0.1, 1.0]},
        "scaled": True
    },
    "Random Forest": {
        "model": RandomForestRegressor(random_state=42),
        "params": {"n_estimators": [50, 100, 200],
                   "max_depth": [3, 5, None]},
        "scaled": False
    },
    "XGBoost": {
        "model": xgb.XGBRegressor(random_state=42, verbosity=0),
        "params": {"n_estimators": [50, 100, 200],
                   "max_depth": [3, 5, 7],
                   "learning_rate": [0.01, 0.1, 0.2]},
        "scaled": False
    }
}

# --- Cross-validation and model comparison ---
cv = KFold(n_splits=10, shuffle=True, random_state=42)
results = {}

print("=" * 60)
print("Model Performance Comparison (10-Fold CV)")
print("=" * 60)

for name, cfg in models.items():
    X_tr = X_train_scaled if cfg["scaled"] else X_train
    X_te = X_test_scaled if cfg["scaled"] else X_test

    gs = GridSearchCV(cfg["model"], cfg["params"], cv=cv,
                      scoring="r2", n_jobs=-1)
    gs.fit(X_tr, y_train)

    y_pred_test = gs.predict(X_te)
    test_r2 = r2_score(y_test, y_pred_test)
    test_mae = mean_absolute_error(y_test, y_pred_test)

    results[name] = {
        "best_model": gs.best_estimator_,
        "cv_r2": gs.best_score_,
        "test_r2": test_r2,
        "test_mae": test_mae
    }

    print(f"\n{name}:")
    print(f"  CV R² = {gs.best_score_:.4f}")
    print(f"  Hold-out R² = {test_r2:.4f}")
    print(f"  Hold-out MAE = {test_mae:.4f}")

# Conventional meta-regression (linear regression baseline)
from sklearn.linear_model import LinearRegression
lr = LinearRegression()
lr.fit(X_train_scaled, y_train)
lr_pred = lr.predict(X_test_scaled)
print(f"\nConventional Meta-Regression:")
print(f"  Hold-out R² = {r2_score(y_test, lr_pred):.4f}")
print(f"  Hold-out MAE = {mean_absolute_error(y_test, lr_pred):.4f}")

# --- SHAP Analysis (best model: XGBoost) ---
print("\n" + "=" * 60)
print("SHAP Analysis (XGBoost)")
print("=" * 60)

best_model = results["XGBoost"]["best_model"]
explainer = shap.TreeExplainer(best_model)
shap_values = explainer.shap_values(X_train)

# Summary plot
plt.figure(figsize=(10, 6))
shap.summary_plot(shap_values, X_train, feature_names=feature_cols, show=False)
plt.tight_layout()
plt.savefig("figures/shap_summary.png", dpi=300, bbox_inches="tight")
plt.close()

# Feature importance (mean |SHAP|)
mean_shap = np.abs(shap_values).mean(axis=0)
for feat, val in sorted(zip(feature_cols, mean_shap), key=lambda x: -x[1]):
    print(f"  {feat}: {val:.4f}")

# --- K-means Clustering ---
print("\n" + "=" * 60)
print("K-means Clustering")
print("=" * 60)

cluster_features = ["Mean_age_num", "Male", "CD Patients", "NOS Score",
                    "Prevalence of UGI Phenotype"]
X_cluster = data[cluster_features].copy()
X_cluster = X_cluster.fillna(X_cluster.median())
X_cluster_scaled = StandardScaler().fit_transform(X_cluster)

# Determine optimal k (elbow + silhouette)
from sklearn.metrics import silhouette_score
inertias = []
silhouettes = []
k_range = range(2, 8)
for k in k_range:
    km = KMeans(n_clusters=k, random_state=42, n_init=10)
    km.fit(X_cluster_scaled)
    inertias.append(km.inertia_)
    silhouettes.append(silhouette_score(X_cluster_scaled, km.labels_))

optimal_k = k_range[np.argmax(silhouettes)]
print(f"Optimal k = {optimal_k} (silhouette = {max(silhouettes):.4f})")

# Final clustering
km_final = KMeans(n_clusters=optimal_k, random_state=42, n_init=10)
data["Cluster"] = km_final.fit_predict(X_cluster_scaled)

# PCA visualization
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_cluster_scaled)

plt.figure(figsize=(8, 6))
for c in range(optimal_k):
    mask = data["Cluster"] == c
    plt.scatter(X_pca[mask, 0], X_pca[mask, 1], label=f"Cluster {c}", s=60)
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
plt.title("K-means Clustering (PCA Visualization)")
plt.legend()
plt.tight_layout()
plt.savefig("figures/kmeans_pca.png", dpi=300, bbox_inches="tight")
plt.close()

print("\nAnalysis complete. Figures saved to figures/")
