library(dplyr)
library(ggplot2)
library(psych)
library(knitr)
library(kableExtra)
library(tidyr)

res <- read.csv("data-clean/res.csv")
str(res)
table(res$condition)
# ---- Factor sex, grade, class_id, condition and status ----
res <- res %>%
  mutate(
    sex = factor(sex),
    grade = factor(grade),
    school = factor(school),
    class_id = factor(class_id),
    condition = factor(condition, 
                       levels = c("control_1", "control_2", 
                                  "digital", "paper", "mixed")),
    status = factor(status, levels = c("participant", "drop-out", 
                                  "atypical", "rejected")),
    ed_level = factor(ed_level)
  )

# Relevel condition so control_1 is reference
res$condition <- relevel(res$condition, ref = "control_1")
res$status    <- relevel(res$status,    ref = "participant")

# Define analysis samples
# Full sample: excludes only "rejected" (declined participation)
res_full <- res %>% filter(status != "rejected")

# Primary analysis sample: excludes rejected AND atypical
res_analysis <- res %>% filter(status %in% c("participant", 
                                             "drop-out"))

# Atypical sample (for supplementary analyses)
res_atypical <- res %>% filter(status == "atypical")

# Quick check
table(res$status)
table(res_full$status)
table(res_analysis$status)

# ---- Get full sample summary (invited vs. participant) ----
attrition_table <- res %>%
  group_by(class_id) %>%
  summarise(
    invited     = n(),
    participant = sum(status == "participant"),
    rejected    = sum(status == "rejected"),
    atypical    = sum(status == "atypical"),
    dropout     = sum(status == "drop-out")
  ) %>%
  ungroup() %>%
  bind_rows(
    summarise(.,
              class_id    = "TOTAL",
              invited     = sum(invited),
              participant = sum(participant),
              rejected    = sum(rejected),
              atypical    = sum(atypical),
              dropout     = sum(dropout))
  )

kable(attrition_table)

# Chi-square test: are dropout/attrition rates equal across conditions?
# Collapse to condition level for chi-square
attrition_by_condition <- res %>%
  group_by(condition) %>%
  summarise(
    participant = sum(status == "participant"),
    dropout     = sum(status == "drop-out"),
    atypical    = sum(status == "atypical"),
    rejected    = sum(status == "rejected")
  )

# Chi-square on participant vs dropout counts by condition
chisq_dropout <- chisq.test(
  rbind(attrition_by_condition$participant,
        attrition_by_condition$dropout)
)
print(chisq_dropout)
# too small n for chi-square
fisher_dropout <- fisher.test(
  rbind(attrition_by_condition$participant,
        attrition_by_condition$dropout),
  simulate.p.value = TRUE
)
print(fisher_dropout)

# ---- Check missing patterns for DV's ----
primary_dvs <- c("gain_corsi", "gain_hnf_overall", 
                 "gain_dccs", "gain_precalc", "gain_tejas")

# Overall missing data summary
missing_summary <- res_analysis %>%
  select(all_of(primary_dvs)) %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  tidyr::pivot_longer(everything(), 
                      names_to  = "variable", 
                      values_to = "n_missing")

kable(missing_summary)

# Check if missingness clusters by condition
missing_by_condition <- res_analysis %>%
  group_by(condition) %>%
  summarise(across(all_of(primary_dvs), ~ sum(is.na(.))))

kable(missing_by_condition)

# Check if missingness clusters by school
missing_by_school <- res_analysis %>%
  group_by(school) %>%
  summarise(across(all_of(primary_dvs), ~ sum(is.na(.))))

kable(missing_by_school)

# Check if missingness clusters by grade
missing_by_grade <- res_analysis %>%
  group_by(grade) %>%
  summarise(across(all_of(primary_dvs), ~ sum(is.na(.))))

kable(missing_by_grade)

# Analysis sample cell sizes
cell_sizes <- res_analysis %>%
  filter(status == "participant") %>%
  group_by(condition, grade) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = grade, values_from = n) %>%
  mutate(total = fifth + kinder)

# Add atypical counts per cell
cell_atypical <- res_atypical %>%
  group_by(condition, grade) %>%
  summarise(n_atypical = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from  = grade, 
                     values_from = n_atypical,
                     names_prefix = "atypical_")

cell_summary <- left_join(cell_sizes, cell_atypical, by = "condition") %>%
  mutate(across(starts_with("atypical_"), ~ replace_na(., 0)))

kable(cell_summary)

# ---- outlier detection on DVs ----
primary_dvs <- c("gain_corsi", "gain_hnf_overall", "gain_hnf_57",
                 "gain_dccs", "gain_precalc", "gain_tejas")

# ---- Compute mean and SD for each DV in analysis sample ----
dv_stats <- res_analysis %>%
  summarise(across(all_of(primary_dvs),
                   list(mean = ~ mean(., na.rm = TRUE),
                        sd   = ~ sd(.,   na.rm = TRUE))))

# Print for reference
print(dv_stats)

# ---- Flag cases beyond 3 SDs ----
flag_outliers <- function(df, dvs, stats) {
  df_flagged <- df %>%
    mutate(across(all_of(dvs), 
                  ~ abs(. - get(paste0(cur_column(), "_mean"), 
                                envir = as.environment(as.list(stats)))) /
                    get(paste0(cur_column(), "_sd"),
                        envir = as.environment(as.list(stats))),
                  .names = "z_{.col}"))
  return(df_flagged)
}

# Compute z-scores for each gain score manually for clarity
res_analysis <- res_analysis %>%
  mutate(
    z_gain_corsi = (gain_corsi - mean(gain_corsi, na.rm = TRUE)) /
      sd(gain_corsi, na.rm = TRUE),
    z_gain_hnf_overall = (gain_hnf_overall - mean(gain_hnf_overall, na.rm = TRUE)) / 
      sd(gain_hnf_overall, na.rm = TRUE),
    z_gain_hnf_57 = (gain_hnf_57 - mean(gain_hnf_57, na.rm = TRUE)) /
      sd(gain_hnf_57, na.rm = TRUE),
    z_gain_hnf_ies = (gain_hnf_ies - mean(gain_hnf_ies, na.rm = TRUE)) /
      sd(gain_hnf_ies, na.rm = TRUE),
    z_gain_dccs = (gain_dccs - mean(gain_dccs, na.rm = TRUE)) /
      sd(gain_dccs, na.rm = TRUE),
    z_gain_precalc = (gain_precalc - mean(gain_precalc, na.rm = TRUE)) /
      sd(gain_precalc, na.rm = TRUE),
    z_gain_tejas = (gain_tejas - mean(gain_tejas, na.rm = TRUE)) /
      sd(gain_tejas, na.rm = TRUE)
  )
# ---- Identify flagged cases ----
outlier_flags <- res_analysis %>%
  filter(abs(z_gain_corsi) > 3 |
           abs(z_gain_hnf_overall) > 3 |
           abs(z_gain_dccs) > 3 |
           abs(z_gain_precalc) > 3 |
           abs(z_gain_tejas) > 3) %>%
  select(id, status, condition, grade, school,
         gain_corsi, gain_hnf_overall, gain_dccs, 
         gain_precalc, gain_tejas,
         z_gain_corsi, z_gain_hnf_overall, z_gain_dccs,
         z_gain_precalc, z_gain_tejas)

cat("Number of flagged cases:", nrow(outlier_flags), "\n\n")
print(outlier_flags)
kable(outlier_flags)

# ---- Cross-reference with atypical status ----
cat("\nStatus breakdown of flagged cases:\n")
print(table(outlier_flags$status))

cat("\nCondition breakdown of flagged cases:\n")
print(table(outlier_flags$condition))

cat("\nGrade breakdown of flagged cases:\n")
print(table(outlier_flags$grade))

# ---- Summary table of flagged cases ----
kable(outlier_flags %>% 
        select(id, status, condition, grade,
               z_gain_corsi, z_gain_hnf_overall, z_gain_dccs,
               z_gain_precalc, z_gain_tejas) %>%
        mutate(across(starts_with("z_"), ~ round(., 2))))

# ---- Code control group as zero attendence ----
res_analysis <- res_analysis %>%
  mutate(attendence = ifelse(is.na(attendence) & condition %in% 
                               c("control_1", "control_2"), 0, 
                                 attendence))

# ---- Save new data frames (suffix: _outlr) and workspace ----
write.csv(res_analysis, "data-clean/res_analysis_outlr.csv", row.names = FALSE)
write.csv(res_atypical, "data-clean/res_atypical_outlr.csv", row.names = FALSE)
write.csv(res_full, "data-clean/res_full_outlr.csv", row.names = FALSE)

# Save global environment as final checkpoint
save.image(file = "rdata_files/may_28_outlr.RData")

