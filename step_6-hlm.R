library(dplyr)
library(knitr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(performance)
library(emmeans)
library(clubSandwich)
library(tidyr)
select <- dplyr::select

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

res_hlm <- res_analysis %>%
  filter(status == "participant")

# ---- Reshape to long DF for analysis ----
hlm_long <- res_hlm %>%
  select(id:school, ed_level, attendence, condition2, cond_grade2,
         total_points_pre, total_points_post,
         nih_score_overall_pre, nih_score_overall_post,
         dccs_tv_score_pre, dccs_tv_score_post,
         precalc_puntaje_total_pre, precalc_puntaje_total_post,
         tejas_puntaje_total_pre, tejas_puntaje_total_post) %>%
  pivot_longer(
    cols = c(total_points_pre, total_points_post,
             nih_score_overall_pre, nih_score_overall_post,
             dccs_tv_score_pre, dccs_tv_score_post,
             precalc_puntaje_total_pre, precalc_puntaje_total_post,
             tejas_puntaje_total_pre, tejas_puntaje_total_post),
    names_to  = c(".value", "time"),
    names_pattern = "(.+)_(pre|post)"
  )

hlm_long <- hlm_long %>%
  mutate(
    time = ifelse(time == "pre", 0, 1),
    time_lab = factor(time, levels = c(0, 1), 
                      labels = c("pre", "post"))
    )

hlm_long <- hlm_long %>%
  rename(corsi = total_points,
         hnf = nih_score_overall,
         dccs = dccs_tv_score,
         woodcock = precalc_puntaje_total,
         tejas = tejas_puntaje_total)

str(hlm_long)

# ---- HLM for Corsi Measure ----
null_corsi <- lmer(corsi ~ 1 + (1 | id), 
                   data = hlm_long, REML = TRUE)
summary(null_corsi)

icc(null_corsi)

# Cross-classified by classroom
xcnull_corsi <- lmer(corsi ~ 1 + (1 | id) + (1 | cond_grade2),
                     data = hlm_long, REML = TRUE)
summary(xcnull_corsi)
icc(xcnull_corsi)
icc(xcnull_corsi, by_group = TRUE)

# Model specification
# Unconditional growth model:
m1 <- lmer(corsi ~ time + (1 | id), data = hlm_long, REML = TRUE)
summary(m1)

m1.2 <- lmer(corsi ~ time + (time | id), data = hlm_long, REML = TRUE)
# Not identifiable

m2 <- lmer(corsi ~ time * condition2 * grade + (1 | id),
           data = hlm_long, REML = TRUE)

m3 <- lmer(corsi ~ time * condition2 * grade + (1 | id) + (1 | cond_grade2), 
           data = hlm_long, REML = TRUE)
summary(m1)
summary(m2)
summary(m3) # Negative Eigenvalue

m2s <- lmer(corsi ~ time * condition2 * grade + (1 | id) + 
              ed_level, data = hlm_long, REML = TRUE)
summary(m2s)

anova(uncond_cbtt, interact_cbtt, sensitiv_cbtt)

coef_test(m2, vcov = "CR1", cluster = hlm_long$cond_grade2)
r2(m2)

# Proportion of between-person variance explained by fixed effects
(8.512 - 3.561) / 8.512

# Kindergarten as a reference:
m2_k <- lmer(corsi ~ time * condition2 * grade + (1 | id),
             data = hlm_long %>% mutate(grade = relevel(grade, ref = "kinder")),
             REML = TRUE)
summary(m2_k)

# ---- HLM for H&F ----



# ---- HLM for DCCS ----



# ---- HLM for Woodcock ----



# ---- HLM for TejasLee ----




# ---- Plots ----



# ---- Save data/data frames ----
save.image(file = "rdata_files/hlm_31_5.RData")
write.csv(hlm_long, file = "data-clean/hlm_long.csv", row.names = F)
