# Task 1: Data quality and missing values
library(mice)

d <- read.csv("project_data_AA6PNR.csv", stringsAsFactors = TRUE)
d$education_level <- factor(d$education_level,
                            ordered = TRUE,
                            levels = c("high_school", "bachelor", "graduate"))

missing_summary <- data.frame(
  Missing_Count = colSums(is.na(d)),
  Missing_Percentage = round(colMeans(is.na(d)) * 100, 2)
)
cat("\n--- Missing Data Summary ---\n")
print(subset(missing_summary, Missing_Count > 0))

vars_with_na <- names(d)[colSums(is.na(d)) > 0]
other_vars   <- setdiff(names(d), "customer_id")   # ID carries no information

miss_tests <- do.call(rbind, lapply(vars_with_na, function(v) {
  R <- factor(is.na(d[[v]]), levels = c(FALSE, TRUE),
              labels = c("Observed", "Missing"))
  do.call(rbind, lapply(setdiff(other_vars, v), function(x) {
    ok <- !is.na(d[[x]])
    if (is.numeric(d[[x]]) && length(unique(d[[x]])) > 2) {
      p <- t.test(d[[x]][ok] ~ R[ok])$p.value
      test <- "t-test"
    } else {
      p <- suppressWarnings(chisq.test(table(R[ok], d[[x]][ok]))$p.value)
      test <- "chi-square"
    }
    data.frame(missing_var = v, compared_with = x, test = test,
               p_value = round(p, 4))
  }))
}))
cat("\n--- Missingness tests: all comparisons ---\n")
print(miss_tests, row.names = FALSE)
cat("\n--- Missingness tests: significant associations (p < 0.05) ---\n")
print(subset(miss_tests, p_value < 0.05))

cat("\n--- Missing education_level by region (row %) ---\n")
print(round(prop.table(table(d$region, is.na(d$education_level)), 1) * 100, 1))

cat("\n--- Mean age by missingness in spend_electronics ---\n")
print(tapply(d$age, is.na(d$spend_electronics), mean))


boxplot(age ~ is.na(spend_electronics), data = d,
        main = "Age vs Missingness in spend_electronics",
        names = c("Observed", "Missing"), xlab = "",
        ylab = "Age", col = c("lightblue", "salmon"))
dev.copy(png, "boxplot_age.png", width = 800, height = 600, res = 120)
dev.off()


init <- mice(d, maxit = 0)
print(init$method)


pred <- init$predictorMatrix
pred[, "customer_id"] <- 0

imp <- mice(d, m = 20, maxit = 10, seed = 2026,
            predictorMatrix = pred, printFlag = FALSE)


print(plot(imp))
dev.copy(png, "mice_convergence.png", width = 900, height = 900, res = 120)
dev.off()

print(stripplot(imp, spend_electronics ~ .imp, pch = 20, cex = 1.2))
dev.copy(png, "mice_stripplot.png", width = 900, height = 700, res = 120)
dev.off()


cat("\nFigures saved in:", getwd(), "\n")
print(list.files(pattern = "\\.png$"))

mean_complete <- mean(d$spend_electronics, na.rm = TRUE)

fit_mean <- with(imp, lm(spend_electronics ~ 1))
pooled_results <- summary(pool(fit_mean), conf.int = TRUE)
mean_pooled <- pooled_results$estimate  
cat(sprintf("\nMean of spend_electronics (Complete cases): %.2f\n", mean_complete)) 
cat(sprintf("Mean of spend_electronics (Pooled imputed): %.2f  [95%% CI %.2f, %.2f]\n",             
            mean_pooled, pooled_results$`2.5 %`, pooled_results$`97.5 %`))

cat("\n--- Observed spend_electronics by age group ---\n")
print(round(tapply(d$spend_electronics,
                   cut(d$age, c(17, 35, 45, 55, 65, 100)),
                   mean, na.rm = TRUE), 2))

# Task 2: What does "satisfaction" mean?

d_complete <- mice::complete(imp, 1)

sat_items  <- c("sat_price", "sat_value", "sat_promotions",
                "sat_staff", "sat_store", "sat_checkout")
spend_vars <- c("spend_grocery", "spend_electronics", "annual_spend")
sat_data   <- d_complete[, sat_items]

pca_result <- prcomp(sat_data, scale. = TRUE)

cat("\n--- PCA: variance explained ---\n")
print(summary(pca_result))

cat("\n--- PCA: eigenvalues ---\n")
print(round(pca_result$sdev^2, 3))

cat("\n--- PCA loadings (first two components) ---\n")
print(round(pca_result$rotation[, 1:2], 3))


plot(pca_result, type = "lines",
     main = "Scree Plot of Satisfaction PCA",
     col = "blue", lwd = 2, pch = 19)
abline(h = 1, lty = 2, col = "grey40")
dev.copy(png, "pca_scree.png", width = 800, height = 600, res = 120)
dev.off()

fa_result <- factanal(sat_data, factors = 2, rotation = "varimax",
                      scores = "regression")
cat("\n--- Factor analysis (varimax) ---\n")
print(fa_result, cutoff = 0.3)   # loadings, variance explained, fit test


L <- unclass(fa_result$loadings)
f_value <- which.max(abs(L["sat_price", ]))
f_exp   <- setdiff(1:2, f_value)

d_complete$Sat_Value      <- fa_result$scores[, f_value] * sign(sum(L[, f_value]))
d_complete$Sat_Experience <- fa_result$scores[, f_exp]   * sign(sum(L[, f_exp]))

cat("\nFactor", f_value, "= Sat_Value, Factor", f_exp, "= Sat_Experience\n")

sat_matrix   <- scale(as.matrix(sat_data))
spend_matrix <- scale(as.matrix(d_complete[, spend_vars]))

cca_result <- cancor(sat_matrix, spend_matrix)

cat("\n--- Canonical correlations ---\n")
print(round(cca_result$cor, 4))


n <- nrow(sat_matrix); p <- ncol(sat_matrix); q <- ncol(spend_matrix)
r <- cca_result$cor
bartlett <- t(sapply(seq_along(r), function(k) {
  chi <- -(n - 1 - (p + q + 1) / 2) * sum(log(1 - r[k:length(r)]^2))
  df  <- (p - k + 1) * (q - k + 1)
  c(from_k = k, chi_sq = chi, df = df,
    p_value = pchisq(chi, df, lower.tail = FALSE))
}))
cat("\n--- Bartlett test (H0: correlations k, k+1, ... are all zero) ---\n")
print(round(bartlett, 4))

U1 <- sat_matrix   %*% cca_result$xcoef[, 1]
V1 <- spend_matrix %*% cca_result$ycoef[, 1]
s  <- sign(cor(V1, spend_matrix[, "annual_spend"]))[1]
U1 <- U1 * s
V1 <- V1 * s

cat("\n--- Standardised canonical weights, first pair ---\n")
print(round(cca_result$xcoef[, 1] * s * sqrt(n - 1), 3))
print(round(cca_result$ycoef[, 1] * s * sqrt(n - 1), 3))

cat("\n--- Structure correlations (variable vs. its own canonical variate) ---\n")
print(round(cor(sat_matrix, U1), 3))
print(round(cor(spend_matrix, V1), 3))


stab <- t(sapply(1:imp$m, function(i) {
  di <- mice::complete(imp, i)
  pc <- prcomp(di[, sat_items], scale. = TRUE)
  c(cum_var_PC1_PC2 = sum(pc$sdev[1:2]^2) / length(sat_items),
    first_can_cor   = cancor(scale(di[, sat_items]),
                             scale(di[, spend_vars]))$cor[1])
}))
cat("\n--- Range across the 20 imputations (min / max) ---\n")
print(round(apply(stab, 2, range), 4))

# Task 3: Channels and customer segments

library(ca)
library(cluster)   # silhouette()


channel_table <- table(d_complete$region, d_complete$preferred_channel)
cat("\n--- Cross-tabulation: Region vs Preferred Channel ---\n")
print(channel_table)

cat("\n--- Row profiles (% of each region using each channel) ---\n")
print(round(prop.table(channel_table, 1) * 100, 1))

chi_test <- chisq.test(channel_table)
cat("\n--- Chi-square Test ---\n")
print(chi_test)

cat("\n--- Minimum expected count ---\n")
print(round(min(chi_test$expected), 1))

cat("\n--- Standardised residuals ---\n")
print(round(chi_test$stdres, 2))

ca_result <- ca(channel_table)
cat("\n--- Correspondence Analysis Summary ---\n")
print(summary(ca_result))

plot(ca_result, main = "Correspondence Analysis: Region vs Channel")
dev.copy(png, "ca_plot.png", width = 900, height = 800, res = 120)
dev.off()

cluster_data <- d_complete[, c("spend_grocery", "spend_electronics")]
cluster_data_scaled <- scale(cluster_data)

cat("\n--- Correlation between the two clustering variables ---\n")
print(round(cor(cluster_data$spend_grocery, cluster_data$spend_electronics), 3))

dist_matrix <- dist(cluster_data_scaled, method = "euclidean")
hc_result <- hclust(dist_matrix, method = "ward.D2")

plot(hc_result, labels = FALSE, hang = -1,
     main = "Dendrogram for Hierarchical Clustering",
     xlab = "", sub = "")
rect.hclust(hc_result, k = 3, border = c("red", "blue", "darkgreen"))
dev.copy(png, "dendrogram.png", width = 900, height = 700, res = 120)
dev.off()

wss <- numeric(10)
for (k in 1:10) {
  set.seed(2026)
  wss[k] <- kmeans(cluster_data_scaled, centers = k, nstart = 25)$tot.withinss
}
cat("\n--- Total within-cluster sum of squares, k = 1..10 ---\n")
print(round(wss, 1))

plot(1:10, wss, type = "b", pch = 19, col = "red", lwd = 2,
     main = "Elbow Method for Optimal k",
     xlab = "Number of Clusters (k)",
     ylab = "Total Within-Cluster Sum of Squares")
dev.copy(png, "elbow_curve.png", width = 800, height = 600, res = 120)
dev.off()

sil_width <- sapply(2:6, function(k) {
  set.seed(2026)
  km <- kmeans(cluster_data_scaled, centers = k, nstart = 25)
  mean(silhouette(km$cluster, dist_matrix)[, "sil_width"])
})
names(sil_width) <- 2:6
cat("\n--- Average silhouette width by k ---\n")
print(round(sil_width, 3))

k_opt <- 3
hc_clusters <- cutree(hc_result, k = k_opt)

set.seed(2026)
kmeans_result <- kmeans(cluster_data_scaled, centers = k_opt, nstart = 25)
kmeans_clusters <- kmeans_result$cluster


agree_table <- table(HC = hc_clusters, KMeans = kmeans_clusters)
cat("\n--- Agreement between Hierarchical and K-Means ---\n")
print(agree_table)
cat(sprintf("Agreement: %d of %d customers (%.1f%%)\n",
            sum(apply(agree_table, 1, max)), sum(agree_table),
            100 * sum(apply(agree_table, 1, max)) / sum(agree_table)))

centres <- aggregate(cluster_data, by = list(cl = kmeans_clusters), FUN = mean)
seg_names <- character(k_opt)
seg_names[centres$cl[which.max(centres$spend_electronics)]] <- "Tech"
seg_names[centres$cl[which.max(centres$spend_grocery)]]     <- "Household"
seg_names[seg_names == ""] <- "Balanced"

d_complete$Segment <- factor(seg_names[kmeans_clusters],
                             levels = c("Tech", "Household", "Balanced"))

profile_vars <- c("spend_grocery", "spend_electronics", "annual_spend",
                  "age", "Sat_Value", "Sat_Experience")
cluster_profiles <- aggregate(d_complete[, profile_vars],
                              by = list(Segment = d_complete$Segment),
                              FUN = mean)
cluster_profiles$Size <- as.numeric(table(d_complete$Segment))
cluster_profiles[, -1] <- round(cluster_profiles[, -1], 2)

cat("\n--- Segment Profiles ---\n")
print(cluster_profiles)


test_vars <- c("annual_spend", "age", "Sat_Value", "Sat_Experience")
anova_p <- sapply(test_vars, function(v) {
  summary(aov(d_complete[[v]] ~ d_complete$Segment))[[1]][["Pr(>F)"]][1]
})
cat("\n--- One-way ANOVA p-values: differences between segments ---\n")
print(round(anova_p, 4))

# Task 4: What drives annual spend?

lm_age_linear <- lm(annual_spend ~ age, data = d)
cat("\n--- Linear age model ---\n")
print(summary(lm_age_linear))

plot(d$age, resid(lm_age_linear),
     main = "Residuals vs Age (Linear Model)",
     xlab = "Age", ylab = "Residuals",
     pch = 20, col = "darkgray")
abline(h = 0, col = "red", lwd = 2)
lines(loess.smooth(d$age, resid(lm_age_linear)), col = "blue", lwd = 2)
dev.copy(png, "residuals_age.png", width = 800, height = 600, res = 120)
dev.off()

lm_age_quad <- lm(annual_spend ~ age + I(age^2), data = d)
cat("\n--- Quadratic age model ---\n")
print(summary(lm_age_quad))

cat("\n--- F-test: linear vs quadratic ---\n")
print(anova(lm_age_linear, lm_age_quad))
cat(sprintf("AIC Linear: %.2f | AIC Quadratic: %.2f\n",
            AIC(lm_age_linear), AIC(lm_age_quad)))

b <- coef(lm_age_quad)["age"]
a <- coef(lm_age_quad)["I(age^2)"]
peak_age <- -b / (2 * a)
cat(sprintf("\nPredicted spending is highest at age: %.1f\n", peak_age))

pred_ages <- data.frame(age = c(20, 30, round(peak_age), 55, 70))
cat("\n--- Predicted annual spend (EUR) by age, quadratic model ---\n")
print(round(cbind(pred_ages, predicted = predict(lm_age_quad, pred_ages)), 2))


lm_final <- lm(annual_spend ~ age + I(age^2) + Sat_Value + Sat_Experience + Segment,
               data = d_complete)
cat("\n--- Final Regression Model Summary ---\n")
print(summary(lm_final))

cat("\n--- 95% confidence intervals (EUR) ---\n")
print(round(confint(lm_final), 2))

lm_no_age <- update(lm_final, . ~ . - age - I(age^2))
lm_no_seg <- update(lm_final, . ~ . - Segment)
cat("\n--- F-test: age terms jointly ---\n")
print(anova(lm_no_age, lm_final))
cat("\n--- F-test: Segment jointly ---\n")
print(anova(lm_no_seg, lm_final))

cat("\n--- SD of satisfaction scores ---\n")
print(round(sapply(d_complete[, c("Sat_Value", "Sat_Experience")], sd), 3))

bf <- coef(lm_final)
cat(sprintf("\nPeak age in the final model: %.1f\n",
            -bf["age"] / (2 * bf["I(age^2)"])))

imp_ext4 <- cbind(imp, data.frame(Sat_Value      = d_complete$Sat_Value,
                                  Sat_Experience = d_complete$Sat_Experience,
                                  Segment        = d_complete$Segment))
fit_final_pooled <- with(imp_ext4,
                         lm(annual_spend ~ age + I(age^2) + Sat_Value +
                              Sat_Experience + Segment))
cat("\n--- Pooled final model (Rubin's rules) ---\n")
print(summary(pool(fit_final_pooled), conf.int = TRUE))


# Task 5: Who upgrades to premium?

library(pROC)

imp_extended <- cbind(imp, data.frame(Sat_Value      = d_complete$Sat_Value,
                                      Sat_Experience = d_complete$Sat_Experience))

old_contr <- options(contrasts = c("contr.treatment", "contr.treatment"))

fit_logit <- with(imp_extended,
                  glm(premium_upgrade ~ education_level + spend_grocery +
                        Sat_Value + Sat_Experience,
                      family = binomial(link = "logit")))
pooled_summary_logit <- summary(pool(fit_logit), conf.int = TRUE)

odds_ratios <- data.frame(
  Predictor  = pooled_summary_logit$term,
  Odds_Ratio = round(exp(pooled_summary_logit$estimate), 3),   CI_Lower   = round(exp(pooled_summary_logit$`2.5 %`), 3),
  CI_Upper   = round(exp(pooled_summary_logit$`97.5 %`), 3),
  P_value    = round(pooled_summary_logit$p.value, 4)
)
cat("\n--- Odds Ratios with 95% CI (pooled logistic regression) ---\n")
print(odds_ratios)

g <- pooled_summary_logit[pooled_summary_logit$term == "spend_grocery", ]
cat(sprintf("\nOR per 10 EUR grocery spend: %.3f [95%% CI %.3f, %.3f]\n",
            exp(10 * g$estimate), exp(10 * g$`2.5 %`), exp(10 * g$`97.5 %`)))

fit_probit <- with(imp_extended,
                   glm(premium_upgrade ~ education_level + spend_grocery +
                         Sat_Value + Sat_Experience,
                       family = binomial(link = "probit")))
pooled_summary_probit <- summary(pool(fit_probit))

coef_compare <- data.frame(
  Predictor   = pooled_summary_logit$term,
  Logit_Coef  = round(pooled_summary_logit$estimate, 4),
  Probit_Coef = round(pooled_summary_probit$estimate, 4),
  Ratio       = round(pooled_summary_logit$estimate /
                        pooled_summary_probit$estimate, 2),
  Logit_p     = round(pooled_summary_logit$p.value, 4),
  Probit_p    = round(pooled_summary_probit$p.value, 4)
)
cat("\n--- Coefficient comparison (logit vs probit) ---\n")
print(coef_compare)

options(old_contr)   # restore default contrasts


prob_logit  <- rowMeans(sapply(fit_logit$analyses,  fitted))
prob_probit <- rowMeans(sapply(fit_probit$analyses, fitted))

cat("\n--- Predicted probabilities: logit vs probit ---\n")
cat(sprintf("Correlation: %.4f\n", cor(prob_logit, prob_probit)))
cat(sprintf("Largest absolute difference: %.4f\n",
            max(abs(prob_logit - prob_probit))))
cat(sprintf("Customers classified differently at 0.5: %d\n",
            sum((prob_logit >= 0.5) != (prob_probit >= 0.5))))

actual_y <- d$premium_upgrade          # no missing values

cat(sprintf("\nShare of premium upgraders: %.3f\n", mean(actual_y)))

roc_curve <- roc(actual_y, prob_logit, quiet = TRUE)
cat(sprintf("AUC: %.4f [95%% CI %.4f, %.4f]\n",
            auc(roc_curve), ci.auc(roc_curve)[1], ci.auc(roc_curve)[3]))

plot(roc_curve, main = "ROC Curve for Premium Upgrade (Logit)",
     col = "darkorange", lwd = 2, print.auc = TRUE)
dev.copy(png, "roc_curve.png", width = 800, height = 800, res = 120)
dev.off()


youden <- coords(roc_curve, "best", best.method = "youden",
                 ret = c("threshold", "sensitivity", "specificity"))
cat("\n--- Youden-optimal threshold ---\n")
print(round(youden, 3))

conf_stats <- function(threshold) {
  pred <- factor(as.numeric(prob_logit >= threshold), levels = c(0, 1))
  tab  <- table(Predicted = pred, Actual = factor(actual_y, levels = c(0, 1)))
  cat(sprintf("\n--- Confusion matrix (threshold = %.3f) ---\n", threshold))
  print(tab)
  cat(sprintf("Accuracy: %.3f | Sensitivity: %.3f | Specificity: %.3f\n",
              sum(diag(tab)) / sum(tab),
              tab["1", "1"] / sum(tab[, "1"]),
              tab["0", "0"] / sum(tab[, "0"])))
}
conf_stats(0.5)
conf_stats(youden$threshold[1])

# Task 6: Survival analysis - did the retention program work?

library(survival)


surv_obj <- Surv(time = d$months_observed, event = d$churn_event)

km_fit <- survfit(surv_obj ~ loyalty_program, data = d)

cat("\n--- Kaplan-Meier: overview by group ---\n")
print(km_fit)

cat("\n--- Estimated share still a customer at 3, 6 and 10 months ---\n")
print(summary(km_fit, times = c(3, 6, 10))[c("strata", "time", "surv")])

plot(km_fit, col = c("red", "blue"), lty = 1:2, lwd = 2,
     main = "Kaplan-Meier Survival Curves by Loyalty Program",
     xlab = "Time (Months)", ylab = "Probability of Still Being a Customer")
legend("bottomleft", legend = c("Not enrolled (0)", "Enrolled (1)"),
       col = c("red", "blue"), lty = 1:2, lwd = 2)
dev.copy(png, "km_curve.png", width = 900, height = 700, res = 120)
dev.off()

cat("\n--- Log-Rank Test ---\n")
print(survdiff(surv_obj ~ loyalty_program, data = d))


cox_model <- coxph(surv_obj ~ loyalty_program + age, data = d)
cat("\n--- Cox Proportional Hazards Model Summary ---\n")
print(summary(cox_model))

cat(sprintf("\nHazard ratio for 10 years of age: %.3f\n",
            exp(10 * coef(cox_model)["age"])))

cat("\n--- Proportional Hazards Assumption Test (cox.zph) ---\n")
print(cox.zph(cox_model))

cat("\n--- Naive comparison: share with churn_event = 1 by group ---\n")
print(round(prop.table(table(Program = d$loyalty_program,
                             Churned = d$churn_event), 1), 3))

cat("\n--- Follow-up time by group (months) ---\n")
print(round(do.call(rbind, tapply(d$months_observed, d$loyalty_program, summary)), 2))

cat("\n--- Customers censored before month 6 ---\n")
print(sum(d$churn_event == 0 & d$months_observed < 6))


d_6m <- subset(d, !(churn_event == 0 & months_observed < 6))
d_6m$churn_6m <- ifelse(d_6m$churn_event == 1 & d_6m$months_observed <= 6, 1, 0)

cat(sprintf("\nCustomers in the 6-month analysis: %d (excluded: %d)\n",
            nrow(d_6m), nrow(d) - nrow(d_6m)))

cat("\n--- Share leaving within 6 months by group ---\n")
print(round(tapply(d_6m$churn_6m, d_6m$loyalty_program, mean), 3))

logit_6m  <- glm(churn_6m ~ loyalty_program + age, data = d_6m,
                 family = binomial(link = "logit"))
probit_6m <- glm(churn_6m ~ loyalty_program + age, data = d_6m,
                 family = binomial(link = "probit"))

cat("\n--- Logit (left within 6 months) ---\n")
print(summary(logit_6m))
cat("\n--- Probit (left within 6 months) ---\n")
print(summary(probit_6m))

cat("\n--- Logit odds ratios with 95% CI ---\n")
print(round(exp(cbind(OR = coef(logit_6m), confint.default(logit_6m))), 3))

cat("\n--- Coefficient ratio logit / probit ---\n")
print(round(coef(logit_6m) / coef(probit_6m), 2))