# 🛒 ShopCo Customer Analytics & Retention Strategy Framework

[![R](https://img.shields.io/badge/R-4.0+-198CE7.svg)](https://www.r-project.org/)
[![Libraries](https://img.shields.io/badge/mice%20|%20survival%20|%20pROC-brightgreen.svg)]()
[![Domain](https://img.shields.io/badge/Retail%20Analytics-Survival%20Analysis-orange.svg)]()

This project is dedicated to building an end-to-end, enterprise-grade customer data analytics and retention strategy framework for ShopCo, a large national retail chain. Moving beyond traditional monolithic analysis, the project introduces advanced missing value imputation, dimensionality reduction, unsupervised clustering, machine learning classification, and biostatistical Survival Analysis. This approach successfully tackles core challenges in retail data, such as dimensional redundancy, survivorship bias, and churn metric miscalculation.

---

## 💡 Market Background & Core Pain Points

When profiling customers and designing retention interventions, data scientists and business analysts typically face the following fatal challenges:
* **Complete-Case Bias:** The traditional approach of dropping rows with missing values causes severe statistical bias. For example, older demographics often omit specific survey fields; directly deleting these records artificially underestimates the value of senior high-net-worth cohorts.
* **Metric Ambiguity & Collinearity:** Customer satisfaction surveys contain massive amounts of highly correlated, redundant metrics. Feeding these directly into predictive models without dimensionality reduction triggers multicollinearity, rendering business attribution uninterpretable.
* **Misleading Binary Churn Metrics:** Traditional "churn (0/1)" metrics completely ignore the varying lengths of the customer lifecycle (statistical right-censoring), causing standard retention rate calculations to fail entirely.

---

## 🏗️ Algorithm Architecture: Three-Tier Evolution Approach

To combat these business traps, this project developed a **progressive three-phase analytical architecture**, upgrading from foundational data engineering to advanced survival modeling:

### 📍 Phase 1: Data Robustness & Orthogonal Dimensionality Reduction
Built a solid data pipeline and feature engineering framework to resolve underlying data quality issues.
* **Multiple Imputation:** Engineered a Multiple Imputation by Chained Equations (MICE) pipeline using the `mice` package. Applied Predictive Mean Matching (PMM) for continuous variables and Proportional Odds Logistic Regression (POLR) for ordered categorical data, scientifically imputing MCAR/MAR gaps while rigorously preserving original data distributions.
* **Orthogonal Feature Extraction:** Utilized Principal Component Analysis (PCA) and Factor Analysis (with Varimax rotation) to perfectly decouple 6 noisy satisfaction metrics into two independent, highly actionable business dimensions: "Value" and "Experience".

### 📍 Phase 2: Segmentation & Behavioral Attribution
Accurately profiled customer value and spending drivers via unsupervised learning and regression analysis.
* **High-Tier Clustering:** Combined Ward's Hierarchical Clustering with K-Means algorithms to engineer a 3-tier customer segmentation model (Tech Enthusiasts, Household Shoppers, Balanced Shoppers), validating inter-cluster significance via ANOVA.
* **Non-linear Spending Drivers:** Discarded poor-fitting linear age assumptions in favor of a Polynomial Regression model. This accurately captured the parabolic spending curve across the customer lifecycle, pinpointing the absolute peak spending age at exactly 43.3 years.

### 📍 Phase 3: Conversion Prediction & Survival Analysis
Introduced probabilistic classification and advanced time-series survival analysis to optimize marketing budgets.
* **Threshold-Optimized Classifier:** Deployed Logistic and Probit link function models to predict premium membership upgrade probabilities. By plotting the ROC curve and calculating Youden's Index, the model identified the optimal business threshold balancing sensitivity and specificity.
* **Cox Proportional Hazards Model:** Adapted medical survival analysis techniques, utilizing Kaplan-Meier curves instead of naive binary churn rates to properly account for customer "follow-up time" (censoring). Validated the constant hazard assumption via Schoenfeld residuals.

---

## 📊 Core Optimization Achievements

Through extensive pooled testing across imputed datasets (Rubin's Rules), the models generated highly actionable commercial insights:

* 🚀 **Conversion Threshold Optimization:** In premium upgrade prediction, by dropping the classification threshold from the default 0.5 down to the statistically optimal **0.378**, the model skyrocketed the capture rate (Sensitivity) of potential upgraders from 57.0% to **75.0%**—without sacrificing overall accuracy (maintained at ~70%). This significantly minimized subscription revenue loss caused by "false negatives".
* 🛡️ **Scientific Quantification of Retention ROI:** The Cox survival model provided highly compelling evidence: customers randomly assigned to the retention program demonstrated a Hazard Ratio of 0.541. Controlling for age, the retention initiative successfully and significantly reduced the customer churn rate by **46%**.

---

## 📂 Project File Structure

```text
.
├── main.R                                 # Core execution script (MICE, PCA, K-Means, Logit/Cox models)
├── ShopCo_Customer_Analytics_Report.pdf   # Executive management summary and detailed statistical report
└── data/
    └── project_data_AA6PNR.csv            # Raw customer sample dataset (features & response variables)
