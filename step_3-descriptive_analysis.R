library(dplyr)
library(ggplot2)
library(psych)
library(knitr)
library(kableExtra)
library(tidyr)
library(pwr)
library(effectsize)
select <- dplyr::select

# ---- Load full analysis data frames -----
res_analysis <- read.csv("data-clean/res_analysis_outlr.csv")
str(res_analysis)

res_atypical <- read.csv("data-clean/res_atypical_outlr.csv")
str(res_atypical)

res_full <- read.csv("data-clean/res_full_outlr.csv")
str(res_full)

# For reliability analyses, need to load item-level data frames
corsi_pre <- read.csv("data-clean/corsi_pre_c.csv")
str(corsi_pre)

precalc_pre <- read.csv("data-clean/precalc_pre_c.csv")
str(precalc_pre)

tejas_pre <- read.csv("data-clean/tejas_pre_c.csv")
str(tejas_pre)

hnf_pre <- read.csv("data-clean/hnf_pre_c.csv")
str(hnf_pre)

dccs_pre <- read.csv("data-clean/dccs_pre_c.csv")
str(dccs_pre)

# ---- Factor variables ----
factor_grouping_vars <- function(df) {
  df <- df %>%
    mutate(
      sex       = factor(sex),
      grade     = factor(grade),
      school    = factor(school),
      class_id  = factor(class_id),
      condition = factor(condition,
                         levels = c("control_1", "control_2",
                                    "digital", "paper", "mixed")),
      status    = factor(status,
                         levels = c("participant", "drop-out",
                                    "atypical", "rejected"))
    )
  
  # Add a column for a single kindergarten control group
  df <- df %>%
    mutate(
      condition2 = factor(
        case_when(
          condition %in% c("control_1", "control_2") ~ "control",
          TRUE ~ as.character(condition)))
    )
  
  df$condition <- relevel(df$condition, ref = "control_1")
  df$status <- relevel(df$status, ref = "participant")
  if ("control" %in% levels(df$condition2)) {
    df$condition2 <- relevel(df$condition2, ref = "control")
  }
  
  return(df)
}

res_analysis <- factor_grouping_vars(res_analysis)
res_full     <- factor_grouping_vars(res_full)
res_atypical <- factor_grouping_vars(res_atypical)
str(res_analysis)
str(res_full)
str(res_atypical)


# ---- demographic descriptives ----
age_by_grade <- res_analysis %>%
  filter(!is.na(age_mths)) %>%
  group_by(grade) %>%
  summarise(
    n = n(),
    mean = round(mean(age_mths), 2),
    sd = round(sd(age_mths),   2),
    min = round(min(age_mths), 2),
    max = round(max(age_mths), 2)
  )
kable(age_by_grade)

age_by_group <- res_analysis %>%
  filter(!is.na(age_mths)) %>%
  mutate(cond_grade = paste(condition, grade, sep = "_")) %>%
  group_by(cond_grade) %>%
  summarise(
    n = n(),
    mean = round(mean(age_mths), 2),
    sd = round(sd(age_mths), 2),
    min = round(min(age_mths), 2),
    max = round(max(age_mths), 2)
  )

kable(age_by_group)

# Check for significant differences:
res_analysis <- res_analysis %>%
  mutate(cond_grade = paste(condition, grade, sep = "_"),
         cond_grade2 = factor(
           case_when(
             cond_grade %in% c("control_1_fifth", 
                               "control_2_fifth") ~ "control_fifth",
             cond_grade %in% c("control_1_kinder", 
                               "control_2_kinder") ~ "control_kinder",
             TRUE ~ as.character(cond_grade))))

age_anova_kinder <- aov(age_mths ~ condition, 
                       data = res_analysis %>% 
                         filter(grade == "kinder", !is.na(age_mths)))
summary(age_anova_kinder)

age_anova_fifth <- aov(age_mths ~ condition,
                       data = res_analysis %>% 
                         filter(grade == "fifth", !is.na(age_mths)))
summary(age_anova_fifth)

# ---- Sex counts by condition x grade ----
sex_table <- res_analysis %>%
  filter(status == "participant") %>%
  mutate(cond_grade = paste(condition, grade, sep = "_")) %>%
  group_by(cond_grade, sex) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = sex, values_from = n, values_fill = 0) %>%
  mutate(total = F + M,
         pct_F = round(F / total * 100, 1),
         pct_M = round(M / total * 100, 1))

kable(sex_table)
sum(sex_table$total)

# ---- Chi-square test: sex distribution across condition x grade ----
sex_chisq_matrix <- res_analysis %>%
  filter(status == "participant", !is.na(sex)) %>%
  mutate(cond_grade = paste(condition, grade, sep = "_")) %>%
  group_by(cond_grade, sex) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = sex, values_from = n, values_fill = 0)

chisq_sex <- chisq.test(
  as.matrix(sex_chisq_matrix[, c("F", "M")])
)
print(chisq_sex)

# ---- Baseline checks for covariates ----
# Parent education level (SES indicator)
ed_table <- res_analysis %>%
  filter(status == "participant", !is.na(ed_level)) %>%
  group_by(condition, ed_level) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = ed_level, 
              values_from = n, 
              values_fill = 0)

kable(ed_table)

# Chi-square test for ed_level distribution across conditions
ed_chisq_matrix <- res_analysis %>%
  filter(status == "participant", !is.na(ed_level)) %>%
  group_by(condition, ed_level) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = ed_level, 
              values_from = n, 
              values_fill = 0)

chisq_ed <- chisq.test(as.matrix(ed_chisq_matrix[, -1]))
print(chisq_ed)

# Also check by grade
ed_table_grade <- res_analysis %>%
  filter(status == "participant", !is.na(ed_level)) %>%
  group_by(condition, grade, ed_level) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = ed_level,
              values_from = n,
              values_fill = 0)

kable(ed_table_grade)

# Chi-square test for ed_level distribution across conditions
ed_chisq_matrix2 <- res_analysis %>%
  filter(status == "participant", !is.na(ed_level)) %>%
  group_by(cond_grade, ed_level) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = ed_level, 
              values_from = n, 
              values_fill = 0)
chisq_ed2 <- chisq.test(as.matrix(ed_chisq_matrix2[, -1]))
print(chisq_ed2)

fisher_ed <- fisher.test(
  as.matrix(ed_chisq_matrix[, -1]),
  simulate.p.value = TRUE,
  B = 10000
)
print(fisher_ed)

# ---- Kindergarten only ----
cat("=== Kindergarten: ed_level distribution by condition ===\n")

ed_kinder <- res_analysis %>%
  filter(status == "participant", !is.na(ed_level),
         grade == "kinder") %>%
  group_by(condition, ed_level) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = ed_level, values_from = n,
              values_fill = 0)
kable(ed_kinder)

fisher_kinder <- fisher.test(
  as.matrix(ed_kinder[, -1]),
  simulate.p.value = TRUE,
  B = 10000
)
cat("Kindergarten Fisher's exact test:\n")
print(fisher_kinder)

# ---- Fifth grade only ----
cat("\n=== Fifth grade: ed_level distribution by condition ===\n")

ed_fifth <- res_analysis %>%
  filter(status == "participant", !is.na(ed_level),
         grade == "fifth") %>%
  group_by(condition, ed_level) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = ed_level, values_from = n,
              values_fill = 0)
kable(ed_fifth)

fisher_fifth <- fisher.test(
  as.matrix(ed_fifth[, -1]),
  simulate.p.value = TRUE,
  B = 10000
)
cat("Fifth grade Fisher's exact test:\n")
print(fisher_fifth)

# ---- Pairwise Fisher's exact tests for ed_level by condition, by grade ----

run_pairwise_fisher_ed <- function(grade_filter) {
  
  dat <- res_analysis %>%
    filter(status == "participant",
           !is.na(ed_level),
           grade == grade_filter)
  
  conditions <- unique(as.character(dat$condition))
  
  pairwise_results <- combn(conditions, 2, function(pair) {
    pair_dat <- dat %>%
      filter(condition %in% pair) %>%
      group_by(condition, ed_level) %>%
      summarise(n = n(), .groups = "drop") %>%
      pivot_wider(names_from  = ed_level,
                  values_from = n,
                  values_fill = 0)
    
    result <- fisher.test(
      as.matrix(pair_dat[, -1]),
      simulate.p.value = TRUE,
      B = 10000
    )
    
    data.frame(
      comparison = paste(pair, collapse = " vs "),
      p_value    = round(result$p.value, 4)
    )
  }, simplify = FALSE)
  
  bind_rows(pairwise_results) %>%
    mutate(p_adj = round(p.adjust(p_value, method = "holm"), 4),
           grade = grade_filter)
}

# Run for each grade
pairwise_ed_kinder <- run_pairwise_fisher_ed("kinder")
pairwise_ed_fifth  <- run_pairwise_fisher_ed("fifth")

# Print results
cat("--- Kindergarten: pairwise ed_level comparisons ---\n")
kable(pairwise_ed_kinder %>% select(comparison, p_value, p_adj),
      col.names = c("Comparison", "p-value", "p-adjusted (Holm)"),
      caption = "Pairwise Fisher's exact tests: parent education by condition (kindergarten)")
# Pairwise shows no sig diff between any group.

cat("\n--- Fifth grade: pairwise ed_level comparisons ---\n")
kable(pairwise_ed_fifth %>% select(comparison, p_value, p_adj),
      col.names = c("Comparison", "p-value", "p-adjusted (Holm)"),
      caption = "Pairwise Fisher's exact tests: parent education by condition (fifth grade)")

# ---- Descriptive: mean ed_level by condition and grade ----
cat("\n=== Mean ed_level by condition and grade ===\n")
ed_means <- res_analysis %>%
  filter(status == "participant", !is.na(ed_level)) %>%
  mutate(ed_level = as.numeric(as.character(ed_level))) %>%
  group_by(condition, grade) %>%
  summarise(
    n       = n(),
    mean_ed = round(mean(ed_level), 2),
    sd_ed   = round(sd(ed_level),   2),
    .groups = "drop"
  )
kable(ed_means,
      caption = "Mean parent education level by condition and grade")

# ---- Baseline check for attendance (intervention groups only) ----
# Descriptives and filter out control
att_desc <- res_analysis %>%
  filter(status == "participant",
         !is.na(attendence),
         !condition %in% c("control_1", "control_2")) %>%
  group_by(cond_grade) %>%
  summarise(
    n      = n(),
    mean   = round(mean(attendence), 2),
    median = round(median(attendence), 2),
    sd     = round(sd(attendence), 2),
    min    = min(attendence),
    max    = max(attendence)
  )

kable(att_desc)

# ANOVA: does attendance differ across intervention conditions?
att_anova <- aov(attendence ~ cond_grade,
                 data = res_analysis %>%
                   filter(status == "participant",
                          !is.na(attendence),
                          !condition %in% c("control_1", "control_2")))

summary(att_anova)

# ---- Reliability Analysis ----
# Filter to analysis sample participants only
analysis_ids <- res_analysis %>%
  filter(status == "participant") %>%
  pull(id)

# ---- Corsi Test-retest reliability for Corsi using control groups only ----
corsi_controls <- res_analysis %>%
  filter(condition %in% c("control_1", "control_2"))

corsi_testretest <- cor.test(
  corsi_controls$total_points_pre,
  corsi_controls$total_points_post,
  method = "pearson"
)
print(corsi_testretest)

cor.test(corsi_controls$total_points_pre[corsi_controls$grade == "fifth"],
         corsi_controls$total_points_post[corsi_controls$grade == "fifth"])

cor.test(corsi_controls$total_points_pre[corsi_controls$grade == "kinder"],
         corsi_controls$total_points_post[corsi_controls$grade == "kinder"])

# ---- H&F: block-level scores (hearts, flowers, mixed) ----
hnf_block_items <- hnf_pre %>%
  filter(id %in% analysis_ids) %>%
  select(nih_score_hearts, nih_score_flowers, nih_score_mixed)

cat("\n--- H&F Block-Level Reliability ---\n")
hnf_alpha <- psych::alpha(hnf_block_items, na.rm = TRUE)
print(hnf_alpha$total)

hnf_omega <- psych::omega(hnf_block_items, nfactors = 1, plot = FALSE)
cat("McDonald's Omega (H&F):", round(hnf_omega$omega.tot, 3), "\n")

# ---- DCCS: trial-level accuracy from dccs_pre ----
dccs_items <- dccs_pre %>%
  filter(id %in% analysis_ids) %>%
  select(starts_with("respuesta_pregunta_"))

cat("\n--- DCCS Trial-Level Reliability ---\n")
dccs_alpha <- psych::alpha(dccs_items, na.rm = TRUE)
print(dccs_alpha$total)

dccs_omega <- psych::omega(dccs_items, nfactors = 1, plot = FALSE)
cat("McDonald's Omega (DCCS):", round(dccs_omega$omega.tot, 3), "\n")

# ---- Woodcock-munoz: item-level point columns from precalc_pre ----
precalc_items <- precalc_pre %>%
  filter(id %in% analysis_ids) %>%
  select(starts_with("puntaje_")) %>%
  select(-matches("ambito")) %>%
  select(-puntaje_total, -puntaje_25)

cat("\n--- Precalc Item-Level Reliability ---\n")
precalc_alpha <- psych::alpha(precalc_items, na.rm = TRUE)
print(precalc_alpha$total)

precalc_omega <- psych::omega(precalc_items, nfactors = 1, plot = FALSE)
cat("McDonald's Omega (Precalc):", round(precalc_omega$omega.tot, 3), "\n")

# ---- TejasLee: item-level point columns from tejas_pre ----
tejas_items <- tejas_pre %>%
  filter(id %in% analysis_ids) %>%
  select(starts_with("puntaje_")) %>%
  select(-matches("ambito"))%>%
  select(-puntaje_total)

cat("\n--- Tejas Item-Level Reliability ---\n")
tejas_alpha <- psych::alpha(tejas_items, na.rm = TRUE)
print(tejas_alpha$total)

tejas_omega <- psych::omega(tejas_items, nfactors = 1, plot = FALSE)
cat("McDonald's Omega (Tejas):", round(tejas_omega$omega.tot, 3), "\n")

# ---- reliability summary ----
reliability_summary <- data.frame(
  Instrument = c("Hearts & Flowers", "DCCS", 
                 "Woodcock-Munoz", "Tejas LEE"),
  Measure = c("3 block scores (hearts, flowers, mixed)",
              "24 trial accuracy scores",
              "72 item scores",
              "72 item scores"),
  Std_Alpha = c(round(hnf_alpha$total$std.alpha, 3),
                round(dccs_alpha$total$std.alpha, 3),
                round(precalc_alpha$total$std.alpha, 3),
                round(tejas_alpha$total$std.alpha, 3)),
  Omega = c(round(hnf_omega$omega.tot, 3),
            round(dccs_omega$omega.tot, 3),
            round(precalc_omega$omega.tot, 3),
            round(tejas_omega$omega.tot, 3))
)

kable(reliability_summary)


# ---- Pretest outcome descriptives ----
# Helper function: mean, median, SD, skewness, kurtosis
desc_stats <- function(x) {
  x <- x[!is.na(x)]
  data.frame(
    n        = length(x),
    mean     = round(mean(x), 2),
    median   = round(median(x), 2),
    sd       = round(sd(x), 2),
    skewness = round(psych::skew(x), 2),
    kurtosis = round(psych::kurtosi(x), 2)
  )
}

# ---- Function to build descriptives table for one variable ----
build_desc_table <- function(df, var, label) {
  df %>%
    filter(!is.na(.data[[var]])) %>%
    group_by(cond_grade) %>%
    summarise(
      n        = n(),
      mean     = round(mean(.data[[var]]),          2),
      median   = round(median(.data[[var]]),        2),
      sd       = round(sd(.data[[var]]),            2),
      skewness = round(psych::skew(.data[[var]]),   2),
      kurtosis = round(psych::kurtosi(.data[[var]]),2),
      .groups  = "drop"
    ) %>%
    mutate(outcome = label) %>%
    select(outcome, everything())
}

# ---- Corsi descriptives ----
corsi_desc <- bind_rows(
  build_desc_table(res_analysis, "total_points_pre", "Corsi Pre"),
  build_desc_table(res_analysis, "total_points_post", "Corsi Post"),
  build_desc_table(res_analysis, "gain_corsi", "Corsi Gain")
)

kable(corsi_desc)

# test for differences at baseline
corsipre_anova_kinder <- aov(total_points_pre ~ cond_grade, 
                             data = res_analysis%>% 
                               filter(grade == "kinder"))
summary(corsipre_anova_kinder)

# Not significant, but app. 0.09 is worth looking at pairwise
pair_k <- res_analysis %>% 
  filter(grade == "kinder")
pair_5 <- res_analysis %>% 
  filter(grade == "fifth")

pairwise.t.test(pair_k$total_points_pre, pair_k$condition)

corsipre_anova_fifth <- aov(total_points_pre ~ cond_grade, 
                             data = res_analysis%>% 
                               filter(grade == "fifth"))
summary(corsipre_anova_fifth)
# NS

# ---- HnF descriptives ----
# need to add a gain score for accuracy
res_analysis <- res_analysis %>%
  mutate(
    acc_prop_overall_gain = (hnf_total_post/57) - (hnf_total_pre/57)
  ) %>%
  relocate(acc_prop_overall_gain, .after = acc_prop_overall_post)

hnf_desc <- bind_rows(
  # Primary TV score (block average)
  build_desc_table(res_analysis, "nih_score_overall_pre",  "H&F TV Score Pre"),
  build_desc_table(res_analysis, "nih_score_overall_post", "H&F TV Score Post"),
  build_desc_table(res_analysis, "gain_hnf_overall", "H&F TV Score Gain"),
  # Sensitivity: 57-trial TV score
  build_desc_table(res_analysis, "nih_score_overall_57_pre", "H&F 57-Trial TV Pre"),
  build_desc_table(res_analysis, "nih_score_overall_57_post", "H&F 57-Trial TV Post"),
  build_desc_table(res_analysis, "gain_hnf_57", "H&F 57-Trial TV Gain"),
  # Sensitivity: IES score
  build_desc_table(res_analysis, "ies_overall_pre", "H&F IES Pre"),
  build_desc_table(res_analysis, "ies_overall_post", "H&F IES Post"),
  build_desc_table(res_analysis, "gain_hnf_ies", "H&F IES Gain"),
  # Block accuracy
  build_desc_table(res_analysis, "acc_prop_hearts_pre", "Acc Hearts Pre"),
  build_desc_table(res_analysis, "acc_prop_hearts_post", "Acc Hearts Post"),
  build_desc_table(res_analysis, "acc_prop_flowers_pre", "Acc Flowers Pre"),
  build_desc_table(res_analysis, "acc_prop_flowers_post", "Acc Flowers Post"),
  build_desc_table(res_analysis, "acc_prop_mixed_pre", "Acc Mixed Pre"),
  build_desc_table(res_analysis, "acc_prop_mixed_post", "Acc Mixed Post"),
  build_desc_table(res_analysis, "acc_prop_overall_pre", "Acc Overall Pre"),
  build_desc_table(res_analysis, "acc_prop_overall_post", "Acc Overall Post"),
  build_desc_table(res_analysis, "acc_prop_overall_gain", "Acc Overall Gain"),
  # Block RT means
  build_desc_table(res_analysis, "rt_mean_hearts_pre",    "RT Hearts Pre"),
  build_desc_table(res_analysis, "rt_mean_hearts_post",   "RT Hearts Post"),
  build_desc_table(res_analysis, "rt_mean_flowers_pre",   "RT Flowers Pre"),
  build_desc_table(res_analysis, "rt_mean_flowers_post",  "RT Flowers Post"),
  build_desc_table(res_analysis, "rt_mean_mixed_pre",     "RT Mixed Pre"),
  build_desc_table(res_analysis, "rt_mean_mixed_post",    "RT Mixed Post"),
  build_desc_table(res_analysis, "rt_mean_overall_pre",   "RT Overall Pre"),
  build_desc_table(res_analysis, "rt_mean_overall_post",  "RT Overall Post")
)

kable(hnf_desc)

# Significance checks on overall on composite score
hnfpre_anova_kinder <- aov(nih_score_overall_pre ~ cond_grade, 
                             data = res_analysis%>% 
                               filter(grade == "kinder"))
summary(hnfpre_anova_kinder)
# Not significant, but app. 0.15 is worth looking at pairwise
pairwise.t.test(pair_k$nih_score_overall_pre, pair_k$condition)
# Mixed group had lower scores but NS compared to any other group.

hnfpre_anova_fifth <- aov(nih_score_overall_pre ~ cond_grade, 
                            data = res_analysis%>% 
                              filter(grade == "fifth"))
summary(hnfpre_anova_fifth)

# ANOVAS by accuracy then RT
# Accuracy by block at pretest
anova_acc_overall_kinder <- aov(acc_prop_overall_pre ~ condition,
                                data = res_analysis %>% 
                                  filter(grade == "kinder"))
anova_acc_overall_fifth  <- aov(acc_prop_overall_pre ~ condition,
                                data = res_analysis %>% 
                                  filter(grade == "fifth"))

# RT overall at pretest
anova_rt_overall_kinder  <- aov(rt_mean_overall_pre ~ condition,
                                data = res_analysis %>% 
                                  filter(grade == "kinder"))
anova_rt_overall_fifth   <- aov(rt_mean_overall_pre ~ condition,
                                data = res_analysis %>% 
                                  filter(grade == "fifth"))

summary(anova_acc_overall_kinder)
# Very close, may need to check pairwise comparisons.
summary(anova_acc_overall_fifth)
# Again, close
summary(anova_rt_overall_kinder)
# statistically the same
summary(anova_rt_overall_fifth)

# pairwise accuracy kinder:
pairwise.t.test(pair_k$acc_prop_overall_pre, pair_k$condition)
# pairwise accuracy fifth
pairwise.t.test(pair_5$acc_prop_overall_pre, pair_5$condition)

# check mean scores
kable(hnf_desc)

# ---- DCCS descriptives ----
dccs_desc <- bind_rows(
  build_desc_table(res_analysis, "dccs_tv_score_pre", "DCCS Pre"),
  build_desc_table(res_analysis, "dccs_tv_score_post", "DCCS Post"),
  build_desc_table(res_analysis, "gain_dccs", "DCCS Gain"),
  build_desc_table(res_analysis, "dccs_mean_acc_pre", "DCCS Acc Pre"),
  build_desc_table(res_analysis, "dccs_mean_acc_post", "DCCS Acc Post"),
  build_desc_table(res_analysis, "dccs_mean_rt_pre", "DCCS RT Pre"),
  build_desc_table(res_analysis, "dccs_mean_rt_post", "DCCS RT Post")
)

kable(dccs_desc)

# Significance checks on overall on composite score
dccspre_anova_kinder <- aov(dccs_tv_score_pre ~ cond_grade, 
                           data = res_analysis%>% 
                             filter(grade == "kinder"))
summary(dccspre_anova_kinder)
# not significant

pairwise.t.test(pair_k$dccs_tv_score_pre, pair_k$condition)
# Pairwise shows no significant difference

dccspre_anova_fifth <- aov(dccs_tv_score_pre ~ cond_grade, 
                          data = res_analysis%>% 
                            filter(grade == "fifth"))
summary(dccspre_anova_fifth)
# NS

# ANOVAS by accuracy then RT
# Accuracy by block at pretest
anova_acc_dccs_kinder <- aov(dccs_mean_acc_pre ~ condition,
                                data = res_analysis %>% 
                                  filter(grade == "kinder"))
summary(anova_acc_dccs_kinder)

pairwise.t.test(pair_k$dccs_mean_acc_pre, pair_k$condition)

anova_acc_dccs_fifth  <- aov(dccs_mean_acc_pre ~ condition,
                                data = res_analysis %>% 
                                  filter(grade == "fifth"))
summary(anova_acc_dccs_fifth)

# RT overall at pretest
anova_rt_dccs_kinder  <- aov(dccs_mean_rt_pre ~ condition,
                                data = res_analysis %>% 
                                  filter(grade == "kinder"))
summary(anova_rt_dccs_kinder)

anova_rt_dccs_fifth   <- aov(dccs_mean_rt_pre ~ condition,
                                data = res_analysis %>% 
                                  filter(grade == "fifth"))
summary(anova_rt_dccs_fifth)


# ---- Precalc descriptives ----
precalc_desc <- bind_rows(
  build_desc_table(res_analysis, "precalc_puntaje_total_pre", "Precalc Pre"),
  build_desc_table(res_analysis, "precalc_puntaje_total_post", "Precalc Post"),
  build_desc_table(res_analysis, "gain_precalc", "Precalc Gain")
)

kable(precalc_desc)

# Significance checks on overall on composite score
precalcpre_anova_kinder <- aov(precalc_puntaje_total_pre ~ cond_grade, 
                            data = res_analysis%>% 
                              filter(grade == "kinder"))
summary(precalcpre_anova_kinder)
# Not significant, will check pairwise, just because the p<.2

precalcpre_anova_fifth <- aov(precalc_puntaje_total_pre ~ cond_grade, 
                           data = res_analysis%>% 
                             filter(grade == "fifth"))
summary(precalcpre_anova_fifth)
# NS

pairwise.t.test(pair_k$precalc_puntaje_total_pre, pair_k$condition)

# ---- Tejas descriptives ----
tejas_desc <- bind_rows(
  build_desc_table(res_analysis, "tejas_puntaje_total_pre", "Tejas Pre"),
  build_desc_table(res_analysis, "tejas_puntaje_total_post", "Tejas Post"),
  build_desc_table(res_analysis, "gain_tejas", "Tejas Gain")
)

kable(tejas_desc)

# Significance checks on overall score
tejaspre_anova_kinder <- aov(tejas_puntaje_total_pre ~ cond_grade, 
                               data = res_analysis%>% 
                                 filter(grade == "kinder"))
summary(tejaspre_anova_kinder)
# NS

tejaspre_anova_fifth <- aov(tejas_puntaje_total_pre ~ cond_grade, 
                              data = res_analysis%>% 
                                filter(grade == "fifth"))
summary(tejaspre_anova_fifth)
# NS


# ---- Complete Baseline Equivalence Summary Table ----
baseline_summary <- data.frame(
  Variable = c(
    # Demographic
    "Age in months (kindergarten)",
    "Age in months (fifth grade)",
    "Sex distribution",
    "Parent education level",
    "Attendance (intervention only)",
    # EF outcomes
    "Corsi Block Test (kindergarten)",
    "Corsi Block Test (fifth grade)",
    "H&F Overall composite (kindergarten)",
    "H&F Overall composite (fifth grade)",
    "H&F Overall accuracy (kindergarten)",
    "H&F Overall accuracy (fifth grade)",
    "H&F Overall RT (kindergarten)",
    "H&F Overall RT (fifth grade)",
    "DCCS composite (kindergarten)",
    "DCCS composite (fifth grade)",
    "DCCS accuracy (kindergarten)",
    "DCCS accuracy (fifth grade)",
    "DCCS RT (kindergarten)",
    "DCCS RT (fifth grade)",
    # Academic outcomes
    "Woodcock-Munoz total (kindergarten)",
    "Woodcock-Munoz total (fifth grade)",
    "Tejas LEE total (kindergarten)",
    "Tejas LEE total (fifth grade)"
  ),
  Test = c(
    "One-way ANOVA", "One-way ANOVA", "Chi-square", 
    "Fisher's exact", "One-way ANOVA",
    "One-way ANOVA", "One-way ANOVA",
    "One-way ANOVA", "One-way ANOVA",
    "One-way ANOVA", "One-way ANOVA",
    "One-way ANOVA", "One-way ANOVA",
    "One-way ANOVA", "One-way ANOVA",
    "One-way ANOVA", "One-way ANOVA",
    "One-way ANOVA", "One-way ANOVA",
    "One-way ANOVA", "One-way ANOVA",
    "One-way ANOVA", "One-way ANOVA"
  ),
  Statistic = c(
    "F(4, 85) = 1.12", "F(3, 90) = 1.47",
    "X²(8) = 8.97",
    "p < .001",
    "F(5, 127) = 1.04",
    "F(4, 86) = 1.92", "F(3, 92) = 1.46",
    "F(4, 86) = 1.83", "F(3, 92) = 1.58",
    "F(4, 86) = 2.26", "F(3, 92) = 2.36",
    "F(4, 86) = 1.04", "F(3, 92) = 0.26",
    "F(4, 86) = 1.17", "F(3, 92) = 0.64",
    "F(4, 86) = 0.96", "F(3, 92) = 1.48",
    "F(4, 86) = 0.81", "F(3, 92) = 0.18",
    "F(4, 86) = 1.66", "F(3, 92) = 1.32",
    "F(4, 86) = 1.49", "F(3, 92) = 0.57"
  ),
  p_value = c(
    ".351", ".227",
    ".344",
    "< .001",
    ".395",
    ".114", ".230",
    ".130", ".199",
    ".069", ".077",
    ".392", ".856",
    ".330", ".594",
    ".436", ".226",
    ".523", ".907",
    ".167", ".272",
    ".212", ".634"
  ),
  Equivalent = c(
    "Yes", "Yes",
    "Yes",
    "No*",
    "Yes",
    "Yes", "Yes",
    "Yes", "Yes",
    "Yes", "Yes",
    "Yes", "Yes",
    "Yes", "Yes",
    "Yes", "Yes",
    "Yes", "Yes",
    "Yes", "Yes",
    "Yes", "Yes"
  )
)

kable(baseline_summary,
      col.names = c("Variable", "Test", "Statistic", 
                    "p-value", "Equivalent"),
      caption = "Table 6. Baseline equivalence summary across conditions")

# ---- Save ----
# Save updated res_analysis with cond_grade variable
write.csv(res_analysis, "data-clean/res_analysis.csv", row.names = FALSE)

# Save global environment checkpoint
save.image(file = "rdata_files/may_28_descriptives_complete.RData")

# The rest of the code is just building tables for descriptive
# statistics for the manuscript and supplementary materials.
# These have been left here for now if researchers wish to modify
# the code and create new tables.

# ---- descriptives tables for manuscript ----
library(tidyverse)
library(knitr)
library(psych)
library(effsize)

# ---- Variable list ----
vars_pre <- c(
  "nih_score_overall_pre",
  "ies_overall_pre",
  "acc_prop_overall_pre",
  "rt_mean_overall_pre",
  "dccs_tv_score_pre",
  "dccs_mean_acc_pre",
  "dccs_mean_rt_pre",
  "total_points_pre",
  "precalc_puntaje_total_pre",
  "tejas_puntaje_total_pre"
)

vars_post <- c(
  "nih_score_overall_post",
  "ies_overall_post",
  "acc_prop_overall_post",
  "rt_mean_overall_post",
  "dccs_tv_score_post",
  "dccs_mean_acc_post",
  "dccs_mean_rt_post",
  "total_points_post",
  "precalc_puntaje_total_post",
  "tejas_puntaje_total_post"
)

var_labels <- c(
  "H&F TV Score",
  "H&F IES",
  "H&F Accuracy",
  "H&F RT",
  "DCCS TV Score",
  "DCCS Accuracy",
  "DCCS RT",
  "Corsi Total Points",
  "Woodcock-Munoz",
  "TejasLee"
)

# ---- Helper: descriptives + paired t-test + Cohen's d ----
compute_desc_table <- function(data, vars_pre, vars_post, var_labels) {
  
  results <- map_dfr(seq_along(vars_pre), function(i) {
    
    vp  <- vars_pre[i]
    vpo <- vars_post[i]
    lbl <- var_labels[i]
    
    x_pre  <- data[[vp]]
    x_post <- data[[vpo]]
    gain   <- x_post - x_pre
    
    # Remove pairs with any NA
    complete <- !is.na(x_pre) & !is.na(x_post)
    x_pre_c  <- x_pre[complete]
    x_post_c <- x_post[complete]
    gain_c   <- gain[complete]
    
    n <- sum(complete)
    
    # Skew and kurtosis (using psych)
    sk_pre  <- round(psych::skew(x_pre_c),  2)
    ku_pre  <- round(psych::kurtosi(x_pre_c), 2)
    sk_post <- round(psych::skew(x_post_c), 2)
    ku_post <- round(psych::kurtosi(x_post_c), 2)
    
    # Paired t-test
    tt <- t.test(x_post_c, x_pre_c, paired = TRUE)
    
    # Cohen's d (paired)
    d <- effsize::cohen.d(x_post_c, x_pre_c, paired = TRUE)
    
    data.frame(
      Variable      = lbl,
      n             = n,
      Pre_M         = round(mean(x_pre_c), 2),
      Pre_SD        = round(sd(x_pre_c), 2),
      Pre_Mdn       = round(median(x_pre_c), 2),
      Pre_Skew      = sk_pre,
      Pre_Kurt      = ku_pre,
      Post_M        = round(mean(x_post_c), 2),
      Post_SD       = round(sd(x_post_c), 2),
      Post_Mdn      = round(median(x_post_c), 2),
      Post_Skew     = sk_post,
      Post_Kurt     = ku_post,
      Gain_M        = round(mean(gain_c), 2),
      Gain_SD       = round(sd(gain_c), 2),
      t             = round(tt$statistic, 3),
      df            = tt$parameter,
      p             = round(tt$p.value, 4),
      d             = round(d$estimate, 3)
    )
  })
  
  results
}

# ---- Normality check helper for baseline equivalence ----
check_normality <- function(data, var, group_var = "condition") {
  data %>%
    filter(!is.na(.data[[var]])) %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      n         = n(),
      W         = tryCatch(shapiro.test(.data[[var]])$statistic, error = function(e) NA),
      p_shapiro = tryCatch(shapiro.test(.data[[var]])$p.value,   error = function(e) NA),
      .groups   = "drop"
    ) %>%
    mutate(variable = var, normal = p_shapiro > .05)
}

# ---- Baseline equivalence helper ----
# Uses Welch ANOVA if normal, Kruskal-Wallis if not
baseline_equiv <- function(data, vars_pre, var_labels) {
  
  map_dfr(seq_along(vars_pre), function(i) {
    
    var <- vars_pre[i]
    lbl <- var_labels[i]
    
    d <- data %>% filter(!is.na(.data[[var]]))
    
    # Check normality per group
    norm_check <- d %>%
      group_by(condition) %>%
      summarise(p = tryCatch(shapiro.test(.data[[var]])$p.value,
                             error = function(e) NA),
                .groups = "drop")
    
    all_normal <- all(norm_check$p > .05, na.rm = TRUE)
    
    if (all_normal) {
      # Welch ANOVA
      fit <- oneway.test(as.formula(paste(var, "~ condition")),
                         data = d, var.equal = FALSE)
      data.frame(
        Variable = lbl,
        Test     = "Welch ANOVA",
        stat     = round(fit$statistic, 3),
        df1      = round(fit$parameter[1], 1),
        df2      = round(fit$parameter[2], 1),
        p        = round(fit$p.value, 4)
      )
    } else {
      # Kruskal-Wallis
      fit <- kruskal.test(as.formula(paste(var, "~ condition")), data = d)
      data.frame(
        Variable = lbl,
        Test     = "Kruskal-Wallis",
        stat     = round(fit$statistic, 3),
        df1      = fit$parameter,
        df2      = NA,
        p        = round(fit$p.value, 4)
      )
    }
  })
}

# ---- Run for kindergarten ----
kinder_data <- res_analysis %>%
  filter(status == "participant", grade == "kinder")

cat("=== KINDERGARTEN: Descriptives + Pre-Post Change ===\n")
kinder_desc <- compute_desc_table(kinder_data, vars_pre, vars_post, var_labels)
kable(kinder_desc)

cat("\n=== KINDERGARTEN: Baseline Equivalence ===\n")
kinder_equiv <- baseline_equiv(kinder_data, vars_pre, var_labels)
kable(kinder_equiv)

# ---- Run for fifth grade ----
fifth_data <- res_analysis %>%
  filter(status == "participant", grade == "fifth")

cat("\n=== FIFTH GRADE: Descriptives + Pre-Post Change ===\n")
fifth_desc <- compute_desc_table(fifth_data, vars_pre, vars_post, var_labels)
kable(fifth_desc)

cat("\n=== FIFTH GRADE: Baseline Equivalence ===\n")
fifth_equiv <- baseline_equiv(fifth_data, vars_pre, var_labels)
kable(fifth_equiv)

# ---- Modified function with group label ----
compute_desc_table_by_group <- function(data, vars_pre, vars_post, 
                                        var_labels, group_label) {
  results <- compute_desc_table(data, vars_pre, vars_post, var_labels)
  results$Group <- group_label
  results <- results %>% select(Group, everything())
  return(results)
}

# ---- Kindergarten: all conditions ----
kinder_groups <- list(
  "control_1" = "Control 1",
  "control_2" = "Control 2",
  "digital"   = "Digital",
  "paper"     = "Paper",
  "mixed"     = "Mixed"
)

kinder_all <- map_dfr(names(kinder_groups), function(cond) {
  d <- kinder_data %>% filter(condition == cond)
  compute_desc_table_by_group(d, vars_pre, vars_post, 
                              var_labels, kinder_groups[[cond]])
})

kable(kinder_all)

# ---- Fifth grade: all conditions ----
fifth_groups <- list(
  "control_1" = "Control 1",
  "digital"   = "Digital",
  "paper"     = "Paper",
  "mixed"     = "Mixed"
)

fifth_all <- map_dfr(names(fifth_groups), function(cond) {
  d <- fifth_data %>% filter(condition == cond)
  compute_desc_table_by_group(d, vars_pre, vars_post,
                              var_labels, fifth_groups[[cond]])
})

kable(fifth_all)
