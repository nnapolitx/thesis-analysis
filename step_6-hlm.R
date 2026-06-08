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

# To run this analysis, it is necessary to have run scripts from 
# step_1 to step_3-descriptive_analysis, and have saved the objects
# to the .csv format. Additionally, be certain to have run the line
# select <- dplyr::select as, if you have loaded the car package in 
# the previous script (step_5), the select method will be masked.

# ---- Load and factor data ----
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

# Center education level variable to make intercept more interpretable
hlm_long <- hlm_long %>%
  mutate(ed_level_c = ed_level - 1)

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
m1_corsi <- lmer(corsi ~ time + (1 | id), data = hlm_long, REML = TRUE)
summary(m1_corsi)

m2_corsi <- lmer(corsi ~ time * condition2 * grade + (1 | id),
           data = hlm_long, REML = TRUE)

summary(m1_corsi)
summary(m2_corsi)

m3_corsi <- lmer(corsi ~ time * condition2 * grade + (1 | id) + (1 | cond_grade2),
                 data = hlm_long, REML = TRUE)
# Eigenvalue basically zero; not interpretable. Probably has to do 
# with collinearity problem between condition and classroom.

m2s_corsi <- lmer(corsi ~ time * condition2 * grade + (1 | id) + 
                  ed_level_c + time:ed_level_c, data = hlm_long, 
                  REML = TRUE)
summary(m2s_corsi)

anova(m1_corsi, m2_corsi, m2s_corsi)

r2(m2_corsi)

# Proportion of between-person variance explained by fixed effects
var_null <- as.numeric(VarCorr(null_corsi)$id)
var_m2   <- as.numeric(VarCorr(m2_corsi)$id)
(var_null - var_m2) / var_null

# Kindergarten as a reference:
m2_k <- lmer(corsi ~ time * condition2 * grade + (1 | id),
             data = hlm_long %>% mutate(grade = relevel(grade, ref = "kinder")),
             REML = TRUE)
summary(m2_k)

# ---- HLM for H&F ----
null_hnf <- lmer(hnf ~ 1 + (1 | id), data = hlm_long, REML = TRUE)
summary(null_hnf)

icc(null_hnf)

# Cross-classified by classroom
xcnull_hnf <- lmer(hnf ~ 1 + (1 | id) + (1 | cond_grade2),
                     data = hlm_long, REML = TRUE)
summary(xcnull_hnf)
icc(xcnull_hnf)
icc(xcnull_hnf, by_group = TRUE)

# Model specification
# Unconditional growth model:
m1_hnf <- lmer(hnf ~ time + (1 | id), data = hlm_long, REML = TRUE)
summary(m1_hnf)

m2_hnf <- lmer(hnf ~ time * condition2 * grade + (1 | id),
           data = hlm_long, REML = TRUE)
summary(m2_hnf)

# m3_hnf <- lmer(hnf ~ time * condition2 * grade + (1 | id) + 
#                 (1 | cond_grade2), data = hlm_long, REML = TRUE)
# summary(m3_hnf) # Negative Eigenvalue

m2s_hnf_main <- lmer(hnf ~ time * condition2 * grade + (1 | id) + 
                      ed_level_c, data = hlm_long, 
                    REML = TRUE)
summary(m2s_hnf_main)

m2s_hnf_mod <- lmer(hnf ~ time * condition2 * grade + (1 | id) + 
                ed_level_c + time:ed_level_c, data = hlm_long, 
                REML = TRUE)
summary(m2s_hnf_mod)

anova(m1_hnf, m2_hnf, m2s_hnf_main, m2s_hnf_mod)

# coef_test(m2_hnf, vcov = "CR1", cluster = hlm_long$cond_grade2)
r2(m2_hnf)

# Proportion of between-person variance explained by fixed effects
var_null <- as.numeric(VarCorr(null_hnf)$id)
var_m2   <- as.numeric(VarCorr(m2_hnf)$id)
(var_null - var_m2) / var_null

# Kindergarten as a reference:
m2_hnfk <- lmer(hnf ~ time * condition2 * grade + (1 | id),
             data = hlm_long %>% mutate(grade = relevel(grade, ref = "kinder")),
             REML = TRUE)
summary(m2_hnfk)

# ---- HLM for DCCS ----
null_dccs <- lmer(dccs ~ 1 + (1 | id), 
                 data = hlm_long, REML = TRUE)
summary(null_dccs)

icc(null_dccs)

# Cross-classified by classroom
xcnull_dccs <- lmer(dccs ~ 1 + (1 | id) + (1 | cond_grade2),
                   data = hlm_long, REML = TRUE)
summary(xcnull_dccs)
icc(xcnull_dccs)
icc(xcnull_dccs, by_group = TRUE)

# Model specification
# Unconditional growth model:
m1_dccs <- lmer(dccs ~ time + (1 | id), data = hlm_long, REML = TRUE)
summary(m1_dccs)

m2_dccs <- lmer(dccs ~ time * condition2 * grade + (1 | id),
               data = hlm_long, REML = TRUE)
summary(m2_dccs)

m3_dccs <- lmer(dccs ~ time * condition2 * grade + (1 | id) + 
                 (1 | cond_grade2), data = hlm_long, REML = TRUE)
summary(m3_dccs) # Negative Eigenvalue

m2s_dccs_main <- lmer(dccs ~ time * condition2 * grade + (1 | id) + 
                  ed_level, data = hlm_long, REML = TRUE)
summary(m2s_dccs_main)

m2s_dccs_mod <- lmer(dccs ~ time * condition2 * grade + (1 | id) + 
                      ed_level_c + time:ed_level_c, 
                      data = hlm_long, REML = TRUE)
summary(m2s_dccs_mod)

anova(m1_dccs, m2_dccs, m2s_dccs_main, m2s_dccs_mod)

coef_test(m2_dccs, vcov = "CR1", cluster = hlm_long$cond_grade2)
r2(m2_dccs)

# Proportion of between-person variance explained by fixed effects
var_null <- as.numeric(VarCorr(null_dccs)$id)
var_m2   <- as.numeric(VarCorr(m2_dccs)$id)
(var_null - var_m2) / var_null

# Kindergarten as a reference:
m2_dccsk <- lmer(dccs ~ time * condition2 * grade + (1 | id),
                data = hlm_long %>% mutate(grade = relevel(grade, ref = "kinder")),
                REML = TRUE)
summary(m2_dccsk)

# ---- HLM for Woodcock ----
null_wdck <- lmer(woodcock ~ 1 + (1 | id), 
                  data = hlm_long, REML = TRUE)
summary(null_wdck)

icc(null_wdck)

# Cross-classified by classroom
xcnull_wdck <- lmer(woodcock ~ 1 + (1 | id) + (1 | cond_grade2),
                    data = hlm_long, REML = TRUE)
summary(xcnull_wdck)
icc(xcnull_wdck)
icc(xcnull_wdck, by_group = TRUE)

# Model specification
# Unconditional growth model:
m1_wdck <- lmer(woodcock ~ time + (1 | id), data = hlm_long, REML = TRUE)
summary(m1_wdck)

m2_wdck <- lmer(woodcock ~ time * condition2 * grade + (1 | id),
                data = hlm_long, REML = TRUE)
summary(m2_wdck)

m3_wdck <- lmer(woodcock ~ time * condition2 * grade + (1 | id) + 
                  (1 | cond_grade2), data = hlm_long, REML = TRUE)
# Negative Eigenvalue

m2s_wdck_main <- lmer(woodcock ~ time * condition2 * grade + (1 | id) + 
                      ed_level_c, data = hlm_long, REML = TRUE)
summary(m2s_wdck_main)

m2s_wdck_mod <- lmer(woodcock ~ time * condition2 * grade + (1 | id) + 
                     ed_level_c + time:ed_level_c, 
                     data = hlm_long, REML = TRUE)
summary(m2s_wdck_mod)

anova(m1_wdck, m2_wdck, m2s_wdck_main, m2s_wdck_mod)
r2(m2_wdck)

# Kindergarten as a reference:
m2_wdckk <- lmer(woodcock ~ time * condition2 * grade + (1 | id),
                 data = hlm_long %>% mutate(grade = relevel(grade, ref = "kinder")),
                 REML = TRUE)
summary(m2_wdckk)

# Proportion of between-person variance explained by fixed effects
var_null <- as.numeric(VarCorr(null_wdck)$id)
var_m2   <- as.numeric(VarCorr(m2_wdck)$id)
(var_null - var_m2) / var_null

# ---- HLM for TejasLee ----
null_tejas <- lmer(tejas ~ 1 + (1 | id), 
                  data = hlm_long, REML = TRUE)
summary(null_tejas)

icc(null_tejas)

# Cross-classified by classroom
xcnull_tejas <- lmer(tejas ~ 1 + (1 | id) + (1 | cond_grade2),
                    data = hlm_long, REML = TRUE)
summary(xcnull_tejas)
icc(xcnull_tejas)
icc(xcnull_tejas, by_group = TRUE)

# Model specification
# Unconditional growth model:
m1_tejas <- lmer(tejas ~ time + (1 | id), data = hlm_long, REML = TRUE)
summary(m1_tejas)

m2_tejas <- lmer(tejas ~ time * condition2 * grade + (1 | id),
                data = hlm_long, REML = TRUE)
summary(m2_tejas)

m3_tejas <- lmer(tejas ~ time * condition2 * grade + (1 | id) + 
                  (1 | cond_grade2), data = hlm_long, REML = TRUE)
# Negative Eigenvalue

m2s_tejas_main <- lmer(tejas ~ time * condition2 * grade + (1 | id) + 
                   ed_level_c, data = hlm_long, REML = TRUE)
summary(m2s_tejas_main)

m2s_tejas_mod <- lmer(tejas ~ time * condition2 * grade + (1 | id) + 
                       ed_level_c + time:ed_level_c, 
                     data = hlm_long, REML = TRUE)
summary(m2s_tejas_mod)


anova(m1_tejas, m2_tejas, m2s_tejas_main, m2s_tejas_mod)

r2(m2_tejas)

# Proportion of between-person variance explained by fixed effects
var_null <- as.numeric(VarCorr(null_tejas)$id)
var_m2   <- as.numeric(VarCorr(m2_tejas)$id)
(var_null - var_m2) / var_null

# Kindergarten as a reference:
m2_tejask <- lmer(tejas ~ time * condition2 * grade + (1 | id),
                 data = hlm_long %>% mutate(grade = relevel(grade, ref = "kinder")),
                 REML = TRUE)
summary(m2_tejask)

# ---- Plots ----
# ICC decomposition bar chart
icc_data <- data.frame(
  outcome = rep(c("CBTT", "H&F", "DCCS", "Woodcock-Munoz", "Tejas LEE"), 3),
  level   = rep(c("Null model", "Participant", "Classroom"), each = 5),
  icc     = c(
    # Null model
    51.3, 76.2, 49.3, 78.0, 82.9,
    # Participant
    13.4, 16.3, 17.9, 23.0, 21.8,
    # Classroom
    39.5, 61.0, 32.4, 56.3, 62.1
  )
)

icc_data$outcome <- factor(icc_data$outcome, 
                           levels = c("CBTT", "H&F", "DCCS", 
                                      "Woodcock-Munoz", "Tejas LEE"))
icc_data$level <- factor(icc_data$level, 
                         levels = c("Null model", "Participant", "Classroom"))

ggplot(icc_data, aes(x = outcome, y = icc, fill = level)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75),
           width = 0.65, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Null model"  = "#185FA5",
                               "Participant" = "#85B7EB",
                               "Classroom"   = "#1D9E75")) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "ICC (%)", fill = NULL,
       caption = "Null model ICC reflects total between-person variance.\nParticipant and classroom ICCs are from the cross-classified null model.") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "top",
    legend.key.size    = unit(0.4, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 11),
    plot.caption       = element_text(size = 9, color = "gray50",
                                      hjust = 0)
  )

# ---- Parental education x baseline score scatter plots ----
# Filter to pretest observations only
res_pre <- hlm_long %>%
  filter(time == 0)

# Common theme
ed_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "top",
    legend.key.size    = unit(0.4, "cm"),
    axis.text.x        = element_text(size = 11),
    plot.caption       = element_text(size = 9, color = "gray50", hjust = 0)
  )

# x-axis labels
ed_labels <- c("0" = "Primary", "1" = "Secondary", 
               "2" = "Vocational", "3" = "University")

# ---- H&F baseline ----
p_hnf <- ggplot(res_pre, aes(x = ed_level_c, y = hnf, 
                             color = grade, shape = grade)) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  scale_x_continuous(breaks = 0:3, labels = ed_labels) +
  scale_color_manual(values = c("fifth"  = "#185FA5",
                                "kinder" = "#1D9E75"),
                     labels = c("fifth"  = "Fifth grade",
                                "kinder" = "Kindergarten")) +
  scale_shape_manual(values = c("fifth"  = 16,
                                "kinder" = 17),
                     labels = c("fifth"  = "Fifth grade",
                                "kinder" = "Kindergarten")) +
  labs(x = "Parental education level", 
       y = "H&F baseline score",
       color = NULL, shape = NULL,
       caption = "β = 0.21, SE = 0.09, p = .029. Jitter added to reduce overplotting.") +
  ed_theme

# ---- Woodcock baseline ----
p_wdck <- ggplot(res_pre, aes(x = ed_level_c, y = woodcock,
                              color = grade, shape = grade)) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  scale_x_continuous(breaks = 0:3, labels = ed_labels) +
  scale_color_manual(values = c("fifth"  = "#185FA5",
                                "kinder" = "#1D9E75"),
                     labels = c("fifth"  = "Fifth grade",
                                "kinder" = "Kindergarten")) +
  scale_shape_manual(values = c("fifth"  = 16,
                                "kinder" = 17),
                     labels = c("fifth"  = "Fifth grade",
                                "kinder" = "Kindergarten")) +
  labs(x = "Parental education level",
       y = "Woodcock-Munoz baseline score",
       color = NULL, shape = NULL,
       caption = "β = 1.92, SE = 0.75, p = .012. Jitter added to reduce overplotting.") +
  ed_theme

# ---- Tejas baseline ----
p_tejas <- ggplot(res_pre, aes(x = ed_level_c, y = tejas,
                               color = grade, shape = grade)) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  scale_x_continuous(breaks = 0:3, labels = ed_labels) +
  scale_color_manual(values = c("fifth"  = "#185FA5",
                                "kinder" = "#1D9E75"),
                     labels = c("fifth"  = "Fifth grade",
                                "kinder" = "Kindergarten")) +
  scale_shape_manual(values = c("fifth"  = 17,
                                "kinder" = 17),
                     labels = c("fifth"  = "Fifth grade",
                                "kinder" = "Kindergarten")) +
  labs(x = "Parental education level",
       y = "Tejas LEE baseline score",
       color = NULL, shape = NULL,
       caption = "β = 2.72, SE = 0.87, p = .002. Jitter added to reduce overplotting.") +
  ed_theme

# ---- Print individually or combine ----
print(p_hnf)
print(p_wdck)
print(p_tejas)

# Optional: combine into one figure using patchwork
# library(patchwork)
# p_hnf / p_wdck / p_tejas +
#   plot_annotation(
#     title = "Parental education and baseline outcome scores",
#     caption = "Lines show linear fits with 95% CI by grade level."
#   )

# ---- Save data/data frames ----
# This may be commented out to avoid overwriting.
# save.image(file = "rdata_files/hlm_1_6.RData")
# write.csv(hlm_long, file = "data-clean/hlm_long.csv", row.names = F)

