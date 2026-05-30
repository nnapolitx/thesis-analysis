library(dplyr)
library(biotools)
library(car)
library(knitr)
library(ggplot2)
library(emmeans)
library(MVN)
library(sandwich)
library(lmtest)
select <- dplyr::select

# ---- Load and factor ----
# Must have run and have saved objects from steps 1-3 for this  
# script to work correctly. Also be sure that 'select' is set to 
# dplyr as the 'car' package masks it.

res_analysis <- read.csv("data-clean/res_analysis.csv")

res_analysis <- res_analysis <- res_analysis %>%
  mutate(
    grade = factor(grade),
    condition = factor(condition),
    cond_grade = factor(cond_grade),
    status = factor(status),
    condition2 = factor(condition2),
    cond_grade2 = factor(cond_grade2),
    # add in HNF accuracy gain by points (currently only exists by %)
    hnf_score_gain = hnf_total_post - hnf_total_pre
  )

# Set reference levels
res_analysis$condition <- relevel(res_analysis$condition, 
                                  ref = "control_1")
res_analysis$cond_grade <- relevel(res_analysis$cond_grade, 
                                   ref = "control_1_fifth")
res_analysis$condition2 <- relevel(res_analysis$condition2, 
                                   ref = "control")
res_analysis$cond_grade2 <- relevel(res_analysis$cond_grade2, 
                                    ref = "control_fifth")
res_analysis$status <- relevel(res_analysis$status,
                               ref = "participant")

# Verify collapsing
kable(table(res_analysis$condition, res_analysis$condition2))
kable(table(res_analysis$cond_grade, res_analysis$cond_grade2))

# ---- Define DVs and analysis sample ----
dvs <- c("gain_corsi", "gain_hnf_overall", "gain_dccs",
               "gain_precalc", "gain_tejas")

dv_names <- c("Corsi (Working Memory)",
               "H&F Overall (Inhibitory Control)",
               "DCCS (Cognitive Flexibility)",
               "Woodcock-Munoz (Math)",
               "Tejas LEE (Literacy)")

pretest_vars <- c("total_points_pre_c",
                  "nih_score_overall_pre_c",
                  "dccs_tv_score_pre_c",
                  "precalc_puntaje_total_pre_c",
                  "tejas_puntaje_total_pre_c")

res_mancova <- res_analysis %>%
  filter(status == "participant")

# ---- Center pretest scores ----
res_mancova <- res_mancova %>%
  mutate(
    total_points_pre_c = scale(total_points_pre, center = TRUE, 
                               scale = FALSE),
    nih_score_overall_pre_c = scale(nih_score_overall_pre, 
                                    center = TRUE, scale = FALSE),
    nih_score_overall_57_pre_c = scale(nih_score_overall_57_pre, 
                                       center = TRUE, scale = FALSE),
    ies_overall_pre_c = scale(ies_overall_pre, center = TRUE, 
                              scale = FALSE),
    hnf_total_pre_c = scale(hnf_total_pre, center = TRUE, 
                            scale = FALSE),
    dccs_tv_score_pre_c = scale(dccs_tv_score_pre, center = TRUE, 
                                scale = FALSE),
    precalc_puntaje_total_pre_c = scale(precalc_puntaje_total_pre, 
                                        center = TRUE, scale = FALSE),
    tejas_puntaje_total_pre_c = scale(tejas_puntaje_total_pre, 
                                      center = TRUE, scale = FALSE)
  )

# ---- MVN checks: full sample ----
dv_matrix <- as.matrix(res_mancova[, dvs])

mardia_result <- mvn(data = dv_matrix, mvn_test = "mardia",
                     univariate_test = "SW", desc = FALSE)
summary(mardia_result)
# Does not comply with multivar normality

hz_result <- mvn(data = dv_matrix, mvn_test = "hz", desc = FALSE)
summary(hz_result)
# Supports mardia test

# ---- MVN checks: by grade ----
for (g in c("kinder", "fifth")) {
  cat("\nGrade:", g, "\n")
  grade_matrix <- res_mancova %>%
    filter(grade == g) %>%
    dplyr::select(all_of(dvs)) %>%
    as.matrix()
  mardia_grade <- mvn(data = grade_matrix, mvn_test = "mardia",
                      desc = FALSE)
  summary(mardia_grade)
}

# Kindergarten does comply with multivariate normality, while fifth
# grade does not.

# ---- MVN checks: by condition2 ----
for (cond in levels(res_mancova$condition2)) {
  cat("\nCondition:", cond, "\n")
  cond_matrix <- res_mancova %>%
    filter(condition2 == cond) %>%
    dplyr::select(all_of(dvs)) %>%
    as.matrix()
  mardia_cond <- mvn(data = cond_matrix, mvn_test = "mardia",
                     desc = FALSE)
  summary(mardia_cond)
}

# None of the conditions comply with multivar norm checks, but 
# considering the developmental differences by age, this is 
# expected.

# ---- MVN checks by condition*grade ----
for (condgrade in levels(res_mancova$cond_grade2)) {
  cat("\nCondition*grade:", condgrade, "\n")
  condgrade_matrix <- res_mancova %>%
    filter(cond_grade2 == condgrade) %>%
    dplyr::select(all_of(dvs)) %>%
    as.matrix()
  mardia_condgrade <- mvn(data = condgrade_matrix, mvn_test = "mardia",
                     desc = FALSE)
  summary(mardia_condgrade)
}
# By classroom, only the digital fifth grade and paper kindergarten
# classrooms do not full comply with multivariate normality checks.

# ---- Multivariate outliers (Mahalanobis distance) ----
mah_dist <- mahalanobis(dv_matrix,
                        colMeans(dv_matrix, na.rm = TRUE),
                        cov(dv_matrix, use = "complete.obs"))

# Critical value: chi-square with df = number of DVs at p < .001
crit_val <- qchisq(0.999, df = length(dvs))
cat("Critical value (chi2, df =", length(dvs), ", p = .001):",
    round(crit_val, 3), "\n")

outliers_mv <- res_mancova %>%
  mutate(mahal_d = mah_dist) %>%
  filter(mahal_d > crit_val) %>%
  select(id, cond_grade2, all_of(dvs), mahal_d)

if (nrow(outliers_mv) > 0) {
  kable(outliers_mv %>% mutate(mahal_d = round(mahal_d, 3)),
        caption = "Multivariate outliers on gain scores")
}

# There are three outliers who had gains above 3SD: one in fifth 
# digital on the HnF; one in mixed Kinder on woodcock; and one in 
# paper kinder on woodcock. These outliers are likely affecting
# the multivariate normality assumption given they correspond to 
# the classrooms that did not comply with the assumption.

# ---- Cook's distance and leverage for flagged outlier (dig504) ----
# Using cond_grade2 consistent with updated variable structure
fit_lm <- lm(gain_hnf_overall ~ cond_grade2, data = res_mancova)
fit_no <- lm(gain_hnf_overall ~ cond_grade2,
              data = subset(res_mancova, id != "dig504"))

cooks <- cooks.distance(fit_lm)
cooks_threshold <- 4 / nrow(model.frame(fit_lm))
lev_threshold <- 2 * length(coef(fit_lm)) / nobs(fit_lm)

cat("\nMost influential case (Cook's distance):\n")
cat("ID:", res_mancova$id[which.max(cooks)], "\n")
cat("Cook's D:", round(cooks[res_mancova$id == "dig504"], 5), "\n")
cat("Threshold (4/n):", round(cooks_threshold, 5), "\n")

lev <- hatvalues(fit_lm)
cat("\nLeverage for dig504:",
    round(lev[res_mancova$id == "dig504"], 5), "\n")
cat("Leverage threshold (2k/n):", round(lev_threshold, 5), "\n")

# Compare coefficients with and without dig504
cat("\nCoefficient comparison (full vs. without dig504):\n")
coef_compare <- data.frame(
  full    = round(coef(fit_lm), 4),
  no_dig504 = round(coef(fit_no), 4)
) %>%
  mutate(difference = round(full - no_dig504, 4))
kable(coef_compare, caption = "Cook's distance sensitivity: dig504")

# Formal model summaries
cat("\nModel with dig504:\n")
print(summary(fit_lm))
cat("\nModel without dig504:\n")
print(summary(fit_no))

# Decision: retain outlier
cat("\nDecision: dig504 retained. Cook's D exceeds 4/n threshold",
    "but leverage is within acceptable limits and coefficient",
    "differences do not produce significant changes to the model's",
    "interpretation.\n")

# ---- Control group homogeneity checks ----
# Note: control_2 is kindergarten only, so all comparisons
# are restricted to kindergarten participants.
res_control_kinder <- res_mancova %>%
  filter(condition %in% c("control_1", "control_2"),
         grade == "kinder")

cat("Control 1 n:", sum(res_control_kinder$condition == "control_1"), "\n")
cat("Control 2 n:", sum(res_control_kinder$condition == "control_2"), "\n")

# ---- Wilcoxon tests: baseline scores ----
baseline_vars <- c("total_points_pre", "nih_score_overall_pre", 
                   "dccs_tv_score_pre", "precalc_puntaje_total_pre", 
                   "tejas_puntaje_total_pre")
baseline_names <- c("Corsi pretest", "Hnf pretest", "DCCS pretest", 
                   "Woodcock pretest", "TejasLee pretest")

cat("\n=== Wilcoxon tests: baseline scores (kindergarten only) ===\n")
for (i in seq_along(baseline_vars)) {
  cat("\n", baseline_names[i], "\n")
  print(wilcox.test(get(baseline_vars[i]) ~ condition,
                    data = res_control_kinder))
}

# ---- Wilcoxon tests: gain scores ----
cat("\n=== Wilcoxon tests: gain scores (kindergarten only) ===\n")
for (i in seq_along(dvs)) {
  cat("\n", dv_names[i], "\n")
  print(wilcox.test(get(dvs[i]) ~ condition,
                    data = res_control_kinder))
}

# ---- Levene's tests: variance homogeneity ----
cat("\n=== Levene's tests: gain scores (kindergarten only) ===\n")
for (i in seq_along(dvs)) {
  cat("\n", dv_names[i], "\n")
  print(leveneTest(get(dvs[i]) ~ condition,
                   data = res_control_kinder))
}

cat("\n=== Levene's tests: baseline scores (kindergarten only) ===\n")
for (i in seq_along(baseline_vars)) {
  cat("\n", baseline_names[i], "\n")
  print(leveneTest(get(baseline_vars[i]) ~ condition,
                   data = res_control_kinder))
}

# ---- Box's M test ----
# Full sample with condition2 (primary)
cat("=== Box's M: full sample, condition2 ===\n")
print(boxM(dv_matrix, res_mancova$condition2))

# Full sample with original condition (for comparison/transparency)
cat("\n=== Box's M: full sample, original condition ===\n")
print(boxM(dv_matrix, res_mancova$condition))

# By grade using condition2
cat("\n=== Box's M by grade (condition2) ===\n")
for (g in c("kinder", "fifth")) {
  cat("\nGrade:", g, "\n")
  grade_data   <- res_mancova %>% filter(grade == g)
  grade_matrix <- as.matrix(grade_data[, dvs])
  
  # Drop unused condition levels for fifth grade
  grade_cond <- droplevels(grade_data$condition2)
  print(boxM(grade_matrix, grade_cond))
}

# by grade using all conditions:
for (g in c("kinder", "fifth")) {
  cat("\nGrade:", g, "\n")
  grade_data   <- res_mancova %>% filter(grade == g)
  grade_matrix <- as.matrix(grade_data[, dvs])
  
  # Drop unused condition levels for fifth grade
  grade_cond <- droplevels(grade_data$condition)
  print(boxM(grade_matrix, grade_cond))
}

# ---- Homogeneity of regression slopes ----

# Test condition2 x pretest interaction for each matched DV
cat("=== Homogeneity of regression slopes: condition2 x pretest ===\n")
slope_results <- data.frame()

for (i in seq_along(dvs)) {
  formula_int <- as.formula(
    paste(dvs[i], "~ condition2 *", baseline_vars[i]))
  
  model_int   <- lm(formula_int, data = res_mancova)
  anova_int   <- Anova(model_int, type = "III")
  int_row     <- anova_int[grep("condition2:", rownames(anova_int)), ]
  
  slope_results <- bind_rows(slope_results, data.frame(
    DV        = dv_names[i],
    Covariate = baseline_names[i],
    F_value   = round(int_row[["F value"]], 3),
    df1       = int_row[["Df"]],
    df2       = model_int$df.residual,
    p_value   = round(int_row[["Pr(>F)"]], 4)
  ))
}

kable(slope_results,
      col.names = c("DV", "Covariate", "F", "df1", "df2", "p"),
      caption   = "Homogeneity of regression slopes: condition2 x pretest")

# Test condition2 x ed_level interaction for each DV
cat("\n=== Homogeneity of regression slopes: condition2 x ed_level ===\n")
ed_slope_results <- data.frame()

for (i in seq_along(dvs)) {
  formula_ed <- as.formula(
    paste(dvs[i], "~ condition2 * ed_level"))
  
  model_ed  <- lm(formula_ed, data = res_mancova)
  anova_ed  <- Anova(model_ed, type = "III")
  int_row   <- anova_ed[grep("condition2:", rownames(anova_ed)), ]
  
  ed_slope_results <- bind_rows(ed_slope_results, data.frame(
    DV        = dv_names[i],
    Covariate = "ed_level",
    F_value   = round(int_row[["F value"]], 3),
    df1       = int_row[["Df"]],
    df2       = model_ed$df.residual,
    p_value   = round(int_row[["Pr(>F)"]], 4)
  ))
}

kable(ed_slope_results,
      col.names = c("DV", "Covariate", "F", "df1", "df2", "p"),
      caption   = "Homogeneity of regression slopes: condition2 x ed_level")

# ---- MANCOVA FULL SAMPLE ----
dv_matrix_full <- cbind(
  res_mancova$gain_corsi,
  res_mancova$gain_hnf_overall,
  res_mancova$gain_dccs,
  res_mancova$gain_precalc,
  res_mancova$gain_tejas
)
colnames(dv_matrix_full) <- c("gain_corsi", "gain_hnf_overall",
                              "gain_dccs", "gain_precalc", 
                              "gain_tejas")

mancova_full <- lm(dv_matrix_full ~ condition2 * grade +
                     total_points_pre_c +
                     nih_score_overall_pre_c +
                     dccs_tv_score_pre_c +
                     precalc_puntaje_total_pre_c +
                     tejas_puntaje_total_pre_c,
                   data = res_mancova)

qr(model.matrix(mancova_full))$rank
ncol(model.matrix(mancova_full))

summary(Manova(mancova_full, type = "III"), multivariate = TRUE)
# Condition is significant, suggesting that different conditions
# regardless of grade have significantly different gains. Grade is
# not signitificant, however grade x condition is significant, 
# which would suggest that different conditions had significantly
# different gains between the two grade levels.

# Adding in ed_level just for sensitivity
mancova_full_s <- lm(dv_matrix_full ~ condition2 * grade +
                     total_points_pre_c +
                     nih_score_overall_pre_c +
                     dccs_tv_score_pre_c +
                     precalc_puntaje_total_pre_c +
                     tejas_puntaje_total_pre_c + ed_level,
                   data = res_mancova)

qr(model.matrix(mancova_full_s))$rank
ncol(model.matrix(mancova_full_s))

summary(Manova(mancova_full_s, type = "III"), multivariate = TRUE)
# Differences NS. Ed_level also NS.

# ---- ANCOVAs FULL SAMPLE ----
ancova_results_full <- data.frame()

for (i in seq_along(dvs)) {
  dv      <- dvs[i]
  pretest <- pretest_vars[i]
  
  formula_ancova <- as.formula(
    paste(dv, "~ condition2 * grade +", pretest))
  
  model_ancova <- lm(formula_ancova, data = res_mancova)
  anova_out    <- Anova(model_ancova, type = "III")
  
  # Extract rows for condition2, grade, interaction, and pretest
  cond_row  <- anova_out["condition2", ]
  grade_row <- anova_out["grade", ]
  int_row   <- anova_out["condition2:grade", ]
  pre_row   <- anova_out[pretest, ]
  
  # Partial eta squared for each term
  ss_resid   <- anova_out["Residuals", "Sum Sq"]
  eta2_cond  <- cond_row[["Sum Sq"]]  / (cond_row[["Sum Sq"]]  + ss_resid)
  eta2_grade <- grade_row[["Sum Sq"]] / (grade_row[["Sum Sq"]] + ss_resid)
  eta2_int   <- int_row[["Sum Sq"]]   / (int_row[["Sum Sq"]]   + ss_resid)
  
  ancova_results_full <- bind_rows(ancova_results_full, data.frame(
    DV           = dv_names[i],
    F_cond       = round(cond_row[["F value"]], 3),
    p_cond       = round(cond_row[["Pr(>F)"]], 4),
    eta2_cond    = round(eta2_cond, 3),
    F_grade      = round(grade_row[["F value"]], 3),
    p_grade      = round(grade_row[["Pr(>F)"]], 4),
    eta2_grade   = round(eta2_grade, 3),
    F_int        = round(int_row[["F value"]], 3),
    p_int        = round(int_row[["Pr(>F)"]], 4),
    eta2_int     = round(eta2_int, 3)
  ))
}

# Apply FDR and Bonferroni corrections separately for each term
ancova_results_full <- ancova_results_full %>%
  mutate(
    p_cond_fdr   = round(p.adjust(p_cond,  method = "fdr"), 4),
    p_cond_bonf  = round(p.adjust(p_cond,  method = "bonferroni"), 4),
    p_int_fdr    = round(p.adjust(p_int,   method = "fdr"), 4),
    p_int_bonf   = round(p.adjust(p_int,   method = "bonferroni"), 4)
  )

# Table: condition2 main effect
kable(ancova_results_full %>%
        select(DV, F_cond, p_cond, eta2_cond, 
               p_cond_fdr, p_cond_bonf),
      col.names = c("Outcome", "F", "p", "partial eta2",
                    "p (FDR)", "p (Bonferroni)"),
      caption   = "Full sample ANCOVAs: condition2 main effect")

# Table: condition2 x grade interaction
kable(ancova_results_full %>%
        select(DV, F_int, p_int, eta2_int,
               p_int_fdr, p_int_bonf),
      col.names = c("Outcome", "F", "p", "partial eta2",
                    "p (FDR)", "p (Bonferroni)"),
      caption   = "Full sample ANCOVAs: condition2 x grade interaction")

# Full ANCOVA summaries to see all terms
cat("\n--- Full ANCOVA summaries ---\n")
for (i in seq_along(dvs)) {
  dv      <- dvs[i]
  pretest <- pretest_vars[i]
  
  cat("\n", dv_names[i], "\n")
  
  formula_ancova <- as.formula(
    paste(dv, "~ condition2 * grade +", pretest))
  
  model_ancova <- lm(formula_ancova, data = res_mancova)
  print(Anova(model_ancova, type = "III"))
}

# ---- Pairwise comparisons: Full sample ----
# H&F - significant condition effect, no significant interaction
model_hnf_full <- lm(gain_hnf_overall ~ condition2 * grade + 
                       nih_score_overall_pre_c,
                     data = res_mancova)
emm_hnf_full   <- emmeans(model_hnf_full, ~ condition2)
pairs_hnf_full <- pairs(emm_hnf_full, adjust = "holm")
cat("=== H&F: Pairwise comparisons (full sample) ===\n")
print(summary(pairs_hnf_full))

# CBTT - significant effect, no interaction, no bonf correction
model_corsi_full <- lm(gain_corsi ~ condition2 * grade + 
                         total_points_pre_c,
                       data = res_mancova)
emm_corsi_full   <- emmeans(model_corsi_full, ~ condition2)
pairs_corsi_full <- pairs(emm_corsi_full, adjust = "holm")
cat("\n=== CBTT: Pairwise comparisons (full sample) ===\n")
print(summary(pairs_corsi_full))

# ---- MANCOVA KINDERGARTEN SUBGROUP ----
res_kinder <- res_mancova %>%
  filter(grade == "kinder") %>%
  mutate(condition2 = droplevels(condition2))

cat("Kindergarten n:", nrow(res_kinder), "\n")
cat("Condition distribution:\n")
print(table(res_kinder$condition2))

dv_matrix_kinder <- as.matrix(res_kinder[, dvs])
colnames(dv_matrix_kinder) <- dvs

# Confirm rank
mm_kinder <- model.matrix(~ condition2 +
                            total_points_pre_c +
                            nih_score_overall_pre_c +
                            dccs_tv_score_pre_c +
                            precalc_puntaje_total_pre_c +
                            tejas_puntaje_total_pre_c,
                          data = res_kinder)
cat("\nModel matrix rank:", qr(mm_kinder)$rank, "\n")
cat("Model matrix columns:", ncol(mm_kinder), "\n\n")

# Primary model
mancova_kinder <- lm(dv_matrix_kinder ~ condition2 +
                       total_points_pre_c +
                       nih_score_overall_pre_c +
                       dccs_tv_score_pre_c +
                       precalc_puntaje_total_pre_c +
                       tejas_puntaje_total_pre_c,
                     data = res_kinder)

summary(Manova(mancova_kinder, type = "III"), multivariate = TRUE)

# Sensitivity model with ed_level
mancova_kinder_s <- lm(dv_matrix_kinder ~ condition2 +
                         total_points_pre_c +
                         nih_score_overall_pre_c +
                         dccs_tv_score_pre_c +
                         precalc_puntaje_total_pre_c +
                         tejas_puntaje_total_pre_c +
                         ed_level,
                       data = res_kinder)

summary(Manova(mancova_kinder_s, type = "III"), multivariate = TRUE)


# ---- ANCOVAs KINDERGARTEN SUBGROUP ----
ancova_results_kinder <- data.frame()

for (i in seq_along(dvs)) {
  dv      <- dvs[i]
  pretest <- pretest_vars[i]
  
  formula_ancova <- as.formula(
    paste(dv, "~ condition2 +", pretest))
  
  model_ancova <- lm(formula_ancova, data = res_kinder)
  anova_out <- Anova(model_ancova, type = "III")
  
  # Extract rows for condition2, grade, interaction, and pretest
  cond_row <- anova_out["condition2", ]
  pre_row <- anova_out[pretest, ]
  
  # Partial eta squared for each term
  ss_resid <- anova_out["Residuals", "Sum Sq"]
  eta2_cond  <- cond_row[["Sum Sq"]]  / (cond_row[["Sum Sq"]]  + ss_resid)
  
  ancova_results_kinder <- bind_rows(ancova_results_kinder, data.frame(
    DV = dv_names[i],
    F_cond = round(cond_row[["F value"]], 3),
    p_cond = round(cond_row[["Pr(>F)"]], 4),
    eta2_cond = round(eta2_cond, 3)
  ))
}

# Apply FDR and Bonferroni corrections separately for each term
ancova_results_kinder <- ancova_results_kinder %>%
  mutate(
    p_cond_fdr   = round(p.adjust(p_cond,  method = "fdr"), 4),
    p_cond_bonf  = round(p.adjust(p_cond,  method = "bonferroni"), 4)
  )

# Table: condition2 main effect
kable(ancova_results_kinder %>%
        select(DV, F_cond, p_cond, eta2_cond, 
               p_cond_fdr, p_cond_bonf),
      col.names = c("Outcome", "F", "p", "partial eta2",
                    "p (FDR)", "p (Bonferroni)"),
      caption   = "Kindergarten subgroup ANCOVAs: condition2 main effect")

# ANCOVA summaries to see all terms
cat("\n--- ANCOVA summaries ---\n")
for (i in seq_along(dvs)) {
  dv      <- dvs[i]
  pretest <- pretest_vars[i]
  
  cat("\n", dv_names[i], "\n")
  
  formula_ancova <- as.formula(
    paste(dv, "~ condition2 +", pretest))
  
  model_ancova <- lm(formula_ancova, data = res_kinder)
  print(Anova(model_ancova, type = "III"))
}

# ---- Pairwise comparisons: Kindergarten ----
# DCCS (Cognitive Flexibility) - primary, survived correction
model_dccs_k <- lm(gain_dccs ~ condition2 + dccs_tv_score_pre_c,
                   data = res_kinder)
emm_dccs_k   <- emmeans(model_dccs_k, ~ condition2)
pairs_dccs_k <- pairs(emm_dccs_k, adjust = "holm")
cat("=== DCCS: Pairwise comparisons (kindergarten) ===\n")
print(summary(pairs_dccs_k))

# Corsi (Working Memory) - secondary, survived FDR but not Bonferroni
model_corsi_k <- lm(gain_corsi ~ condition2 + total_points_pre_c,
                    data = res_kinder)
emm_corsi_k   <- emmeans(model_corsi_k, ~ condition2)
pairs_corsi_k <- pairs(emm_corsi_k, adjust = "holm")
cat("\n=== Corsi: Pairwise comparisons (kindergarten, secondary) ===\n")
print(summary(pairs_corsi_k))

# ---- Combined summary table ----
all_pairs_kinder <- bind_rows(
  as.data.frame(summary(pairs_dccs_k)) %>%
    mutate(Outcome = "DCCS", Status = "Primary"),
  as.data.frame(summary(pairs_corsi_k)) %>%
    mutate(Outcome = "Corsi", Status = "Secondary")
) %>%
  select(Outcome, Status, contrast, estimate, SE, df, t.ratio, p.value) %>%
  mutate(
    estimate = round(estimate, 3),
    SE       = round(SE, 3),
    t.ratio  = round(t.ratio, 3),
    p.value  = round(p.value, 4)
  )

kable(all_pairs_kinder,
      col.names = c("Outcome", "Status", "Contrast",
                    "Estimate", "SE", "df", "t", "p (Holm)"),
      caption   = "Pairwise comparisons for significant kindergarten ANCOVAs")

# ---- MANCOVA FIFTH SUBGROUP ----
res_fifth <- res_mancova %>%
  filter(grade == "fifth") %>%
  mutate(condition2 = droplevels(condition2))

cat("Fifth n:", nrow(res_fifth), "\n")
cat("Condition distribution:\n")
print(table(res_fifth$condition2))

dv_matrix_fifth <- as.matrix(res_fifth[, dvs])
colnames(dv_matrix_fifth) <- dvs

# Confirm rank
mm_fifth <- model.matrix(~ condition2 +
                            total_points_pre_c +
                            nih_score_overall_pre_c +
                            dccs_tv_score_pre_c +
                            precalc_puntaje_total_pre_c +
                            tejas_puntaje_total_pre_c,
                          data = res_fifth)
cat("\nModel matrix rank:", qr(mm_fifth)$rank, "\n")
cat("Model matrix columns:", ncol(mm_fifth), "\n\n")

# Primary model
mancova_fifth <- lm(dv_matrix_fifth ~ condition2 +
                       total_points_pre_c +
                       nih_score_overall_pre_c +
                       dccs_tv_score_pre_c +
                       precalc_puntaje_total_pre_c +
                       tejas_puntaje_total_pre_c,
                     data = res_fifth)

summary(Manova(mancova_fifth, type = "III"), multivariate = TRUE)

# Sensitivity model with ed_level
mancova_fifth_s <- lm(dv_matrix_fifth ~ condition2 +
                         total_points_pre_c +
                         nih_score_overall_pre_c +
                         dccs_tv_score_pre_c +
                         precalc_puntaje_total_pre_c +
                         tejas_puntaje_total_pre_c +
                         ed_level,
                       data = res_fifth)

summary(Manova(mancova_fifth_s, type = "III"), multivariate = TRUE)

# ---- ANCOVAs FIFTH SUBGROUP ----
ancova_results_fifth <- data.frame()

for (i in seq_along(dvs)) {
  dv      <- dvs[i]
  pretest <- pretest_vars[i]
  
  formula_ancova <- as.formula(
    paste(dv, "~ condition2 +", pretest))
  
  model_ancova <- lm(formula_ancova, data = res_fifth)
  anova_out <- Anova(model_ancova, type = "III")
  
  # Extract rows for condition2, grade, interaction, and pretest
  cond_row <- anova_out["condition2", ]
  pre_row <- anova_out[pretest, ]
  
  # Partial eta squared for each term
  ss_resid <- anova_out["Residuals", "Sum Sq"]
  eta2_cond  <- cond_row[["Sum Sq"]]  / (cond_row[["Sum Sq"]]  + ss_resid)
  
  ancova_results_fifth <- bind_rows(ancova_results_fifth, data.frame(
    DV = dv_names[i],
    F_cond = round(cond_row[["F value"]], 3),
    p_cond = round(cond_row[["Pr(>F)"]], 4),
    eta2_cond = round(eta2_cond, 3)
  ))
}

# Apply FDR and Bonferroni corrections separately for each term
ancova_results_fifth <- ancova_results_fifth %>%
  mutate(
    p_cond_fdr   = round(p.adjust(p_cond,  method = "fdr"), 4),
    p_cond_bonf  = round(p.adjust(p_cond,  method = "bonferroni"), 4)
  )

# Table: condition2 main effect
kable(ancova_results_fifth %>%
        select(DV, F_cond, p_cond, eta2_cond, 
               p_cond_fdr, p_cond_bonf),
      col.names = c("Outcome", "F", "p", "partial eta2",
                    "p (FDR)", "p (Bonferroni)"),
      caption = "Fifth grade sample ANCOVAs: condition2 main effect")

# Full ANCOVA summaries to see all terms
cat("\n--- Full ANCOVA summaries ---\n")
for (i in seq_along(dvs)) {
  dv      <- dvs[i]
  pretest <- pretest_vars[i]
  
  cat("\n", dv_names[i], "\n")
  
  formula_ancova <- as.formula(
    paste(dv, "~ condition2 +", pretest))
  
  model_ancova <- lm(formula_ancova, data = res_fifth)
  print(Anova(model_ancova, type = "III"))
}


# ---- Pairwise comparisons: Fifth grade ----

# H&F Overall (Inhibitory Control) - primary, survived both corrections
model_hnf_5 <- lm(gain_hnf_overall ~ condition2 + nih_score_overall_pre_c,
                  data = res_fifth)
emm_hnf_5   <- emmeans(model_hnf_5, ~ condition2)
pairs_hnf_5 <- pairs(emm_hnf_5, adjust = "holm")
cat("=== H&F Overall: Pairwise comparisons (fifth grade) ===\n")
print(summary(pairs_hnf_5))

# DCCS (Cognitive Flexibility) - primary, survived both corrections
model_dccs_5 <- lm(gain_dccs ~ condition2 + dccs_tv_score_pre_c,
                   data = res_fifth)
emm_dccs_5   <- emmeans(model_dccs_5, ~ condition2)
pairs_dccs_5 <- pairs(emm_dccs_5, adjust = "holm")
cat("\n=== DCCS: Pairwise comparisons (fifth grade) ===\n")
print(summary(pairs_dccs_5))

# ---- Combined summary table ----
all_pairs_fifth <- bind_rows(
  as.data.frame(summary(pairs_hnf_5)) %>%
    mutate(Outcome = "H&F Overall", Status = "Primary"),
  as.data.frame(summary(pairs_dccs_5)) %>%
    mutate(Outcome = "DCCS", Status = "Primary")
) %>%
  select(Outcome, Status, contrast, estimate, SE, df, t.ratio, p.value) %>%
  mutate(
    estimate = round(estimate, 3),
    SE       = round(SE, 3),
    t.ratio  = round(t.ratio, 3),
    p.value  = round(p.value, 4)
  )

kable(all_pairs_fifth,
      col.names = c("Outcome", "Status", "Contrast",
                    "Estimate", "SE", "df", "t", "p (Holm)"),
      caption   = "Pairwise comparisons for significant fifth grade ANCOVAs")
# ---- SAVE rdata ----
save.image(file = "rdata_files/MANCOVA_res_29_5.RData")
# ---- Plots ----
library(ggplot2)
library(emmeans)
library(dplyr)
# 
# ---- Shared theme ----
theme_publication <- theme_bw() +
  theme(
    text              = element_text(family = "Arial", size = 11),
    plot.title        = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle     = element_text(size = 10, hjust = 0, color = "grey40"),
    axis.title        = element_text(size = 11),
    axis.text         = element_text(size = 10),
    strip.text        = element_text(face = "bold", size = 11),
    strip.background  = element_rect(fill = "grey92", color = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "none"
  )

condition_colors <- c(
  "control - digital" = "#2166AC",
  "control - mixed"   = "#4DAC26",
  "control - paper"   = "#D01C8B",
  "digital - mixed"   = "#F1B6DA",
  "digital - paper"   = "#B8E186",
  "mixed - paper"     = "#92C5DE"
)

# ---- Helper: significance label ----
sig_label <- function(p) {
  case_when(
    p < .001 ~ "***",
    p < .01  ~ "**",
    p < .05  ~ "*",
    p < .10  ~ "\u2020",
    TRUE     ~ "ns"
  )
}

# PLOT 1: Forest plot - Kindergarten pairwise comparisons
# Re-run models to ensure objects are available
model_dccs_k  <- lm(gain_dccs  ~ condition2 + dccs_tv_score_pre_c,
                    data = res_kinder)
model_corsi_k <- lm(gain_corsi ~ condition2 + total_points_pre_c,
                    data = res_kinder)

pairs_dccs_k  <- pairs(emmeans(model_dccs_k,  ~ condition2), adjust = "holm")
pairs_corsi_k <- pairs(emmeans(model_corsi_k, ~ condition2), adjust = "holm")

kinder_pairs <- bind_rows(
  as.data.frame(summary(pairs_dccs_k))  %>% mutate(Outcome = "DCCS\n(Cognitive Flexibility)",
                                                   Status = "Primary"),
  as.data.frame(summary(pairs_corsi_k)) %>% mutate(Outcome = "Corsi\n(Working Memory)",
                                                   Status = "Secondary")
) %>%
  mutate(
    contrast  = factor(contrast, levels = rev(unique(contrast))),
    sig       = sig_label(p.value),
    lower     = estimate - SE,
    upper     = estimate + SE,
    significant = p.value < .05
  )

p_kinder <- ggplot(kinder_pairs,
                   aes(x = estimate, y = contrast,
                       color = contrast, alpha = significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.25, linewidth = 0.8) +
  geom_point(size = 3) +
  geom_text(aes(label = sig, x = upper + 0.05),
            hjust = 0, size = 3.5, color = "black") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4)) +
  scale_color_manual(values = condition_colors) +
  facet_wrap(~ paste0(Status, ": ", Outcome), scales = "free_x", ncol = 2) +
  labs(
    title    = "Kindergarten: Pairwise Contrasts by Outcome",
    subtitle = "Estimates \u00b1 SE (Holm-corrected). \u2020p < .10, *p < .05, **p < .01, ***p < .001, ns = not significant",
    x        = "Estimated Difference in Gain Score",
    y        = NULL
  ) +
  theme_publication

ggsave("plots/forest_kinder.png",
       plot   = p_kinder,
       width  = 10,
       height = 5.5,
       dpi    = 600,
       bg     = "white")




# PLOT 2: Forest plot - Fifth grade pairwise comparisons
model_hnf_5  <- lm(gain_hnf_overall ~ condition2 + nih_score_overall_pre_c,
                   data = res_fifth)
model_dccs_5 <- lm(gain_dccs ~ condition2 + dccs_tv_score_pre_c,
                   data = res_fifth)

pairs_hnf_5  <- pairs(emmeans(model_hnf_5,  ~ condition2), adjust = "holm")
pairs_dccs_5 <- pairs(emmeans(model_dccs_5, ~ condition2), adjust = "holm")

fifth_pairs <- bind_rows(
  as.data.frame(summary(pairs_hnf_5))  %>% mutate(Outcome = "H&F\n(Inhibitory Control)",
                                                  Status = "Primary"),
  as.data.frame(summary(pairs_dccs_5)) %>% mutate(Outcome = "DCCS\n(Cognitive Flexibility)",
                                                  Status = "Primary")
) %>%
  mutate(
    contrast    = factor(contrast, levels = rev(unique(contrast))),
    sig         = sig_label(p.value),
    lower       = estimate - SE,
    upper       = estimate + SE,
    significant = p.value < .05
  )

p_fifth <- ggplot(fifth_pairs,
                  aes(x = estimate, y = contrast,
                      color = contrast, alpha = significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.25, linewidth = 0.8) +
  geom_point(size = 3) +
  geom_text(aes(label = sig, x = upper + 0.03),
            hjust = 0, size = 3.5, color = "black") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4)) +
  scale_color_manual(values = condition_colors) +
  facet_wrap(~ paste0(Status, ": ", Outcome), scales = "free_x", ncol = 2) +
  labs(
    title    = "Fifth Grade: Pairwise Contrasts by Outcome",
    subtitle = "Estimates \u00b1 SE (Holm-corrected). \u2020p < .10, *p < .05, **p < .01, ***p < .001, ns = not significant",
    x        = "Estimated Difference in Gain Score",
    y        = NULL
  ) +
  theme_publication

ggsave("plots/forest_fifth.png",
       plot   = p_fifth,
       width  = 10,
       height = 5.5,
       dpi    = 600,
       bg     = "white")


# PLOT 3: Interaction plot - DCCS by condition and grade
# Full sample model with interaction for emmeans
model_dccs_full <- lm(gain_dccs ~ condition2 * grade + dccs_tv_score_pre_c,
                      data = res_mancova)

emm_dccs_int <- emmeans(model_dccs_full, ~ condition2 | grade)
dccs_int_df  <- as.data.frame(emm_dccs_int) %>%
  rename(estimate = emmean, SE = SE) %>%
  mutate(
    grade      = factor(grade,
                        levels = c("kinder", "fifth"),
                        labels = c("Kindergarten", "Fifth Grade")),
    condition2 = factor(condition2,
                        levels = c("control", "digital", "paper", "mixed"),
                        labels = c("Control", "Digital", "Paper", "Mixed")),
    lower      = estimate - SE,
    upper      = estimate + SE
  )

p_interaction <- ggplot(dccs_int_df,
                        aes(x = condition2, y = estimate,
                            group = grade, color = grade, shape = grade)) +
  geom_line(linewidth = 0.9, position = position_dodge(0.15)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.12, linewidth = 0.7,
                position = position_dodge(0.15)) +
  geom_point(size = 3.5, position = position_dodge(0.15)) +
  scale_color_manual(
    values = c("Kindergarten" = "#2166AC", "Fifth Grade" = "#D01C8B"),
    name   = "Grade"
  ) +
  scale_shape_manual(
    values = c("Kindergarten" = 16, "Fifth Grade" = 17),
    name   = "Grade"
  ) +
  labs(
    title    = "DCCS Cognitive Flexibility: Condition \u00d7 Grade Interaction",
    subtitle = "Estimated marginal means \u00b1 SE, adjusted for DCCS pretest",
    x        = "Condition",
    y        = "Adjusted Gain Score (EMM)"
  ) +
  theme_publication +
  theme(legend.position = "right")

ggsave("plots/interaction_dccs.png",
       plot   = p_interaction,
       width  = 7,
       height = 5,
       dpi    = 600,
       bg     = "white")
