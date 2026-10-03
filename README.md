# ShopCo Customer Analytics & Retention Strategy

## Overview
This repository contains a comprehensive end-to-end customer analytics project for a fictitious retail chain, ShopCo. The project analyzes customer behavior, segmentation, and retention using advanced statistical and machine learning techniques in R. 

## Key Features & Methodology
* **Missing Data Imputation:** Handled MCAR and MAR data gaps using Multiple Imputation by Chained Equations (`mice`), effectively eliminating complete-case bias.
* **Dimensionality Reduction:** Applied Principal Component Analysis (PCA) and Factor Analysis (Varimax rotation) to deconstruct survey data into orthogonal "Value" and "Experience" dimensions.
* **Customer Segmentation:** Engineered a 3-tier customer segmentation model (Tech Enthusiasts, Household Shoppers, Balanced Shoppers) using a hybrid approach of K-Means and Ward’s Hierarchical Clustering.
* **Predictive Modeling:** Developed pooled polynomial regression (to predict peak spending age) and logistic regression models (to classify premium upgrade likelihood), optimizing decision thresholds via Youden's Index.
* **Survival Analysis:** Quantified the ROI of a randomized customer retention program using Kaplan-Meier curves and Cox Proportional Hazards modeling, proving a 46% reduction in churn rate.

## Repository Structure
* `main.R`: The complete reproducible R script covering data processing, modeling, and visualization.
* `ShopCo_Customer_Analytics_Report.pdf`: The detailed management summary and technical report.
* `data/`: Contains the customer dataset used for this analysis.

## Tech Stack
**Language:** R
**Key Libraries:** `mice`, `survival`, `pROC`, `cluster`, `ca`
