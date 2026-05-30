library(dplyr)
library(ggplot2)
library(psych)
library(readr)
library(stringr)
library(janitor)
library(readxl)
select <- dplyr::select

# set wd just in case
setwd("C:/Users/nnapo/Documents/PhD Classes/Tesis 1-6/analysis/thesis-analysis")


# ----Load data frames from raw data files (anonymized)----
# Corsi pre and posttest
corsi_pre <- read.csv("data-raw/corsi-pretest.csv")
corsi_post <- read.csv("data-raw/corsi-post.csv")

# hearts and flowers pre and posttest
hnf_pre <- read.csv("data-raw/hnf-pretest.csv")
hnf_post <- read.csv("data-raw/hnf-post.csv")

# dimensional change card sort pre and posttest
dccs_pre <- read.csv("data-raw/dccs-pretest.csv")
dccs_post <- read.csv("data-raw/dccs-post.csv")

# Woodcock-munoz pre and posttest
precalc_pre <- read.csv("data-raw/precalc-pretest.csv")
precalc_post <- read.csv("data-raw/precalc-post.csv")

# TejasLee pre and posttest
tejas_pre <- read.csv("data-raw/tejas-pretest.csv")
tejas_post <- read.csv("data-raw/tejas-post.csv")

tejas_pre$pregunta_9 <- as.integer(tejas_pre$pregunta_9)
tejas_post$pregunta_16 <- as.integer(tejas_post$pregunta_16)
tejas_post$pregunta_17 <- as.integer(tejas_post$pregunta_17)

# Socio-demografic data and attendence
dmogrph <- read_excel("data-raw/demo.xlsx")

# ----Factoring categorical vars----
dfnames <- c("corsi_pre", "corsi_post", "dccs_pre", "dccs_post", "dmogrph", "hnf_pre", 
             "hnf_post", "precalc_pre", "precalc_post", "tejas_pre", "tejas_post")

# Function to factor the relevant variables
factor_vars <- function(df) {
  df$sex       <- factor(df$sex)
  df$grade     <- factor(df$grade)
  df$school    <- factor(df$school)
  df$condition <- relevel(factor(df$condition), ref = "control_1")
  df$status    <- relevel(factor(df$status), ref = "participant")
  return(df)
}

# Apply to all dataframes
for (name in dfnames) {
  assign(name, factor_vars(get(name)))
}

# ----Factor any other cat vars ----
dmogrph$ed_level <- relevel(factor(dmogrph$ed_level), 
                            ref = 1) 
dmogrph$esp_fl <- factor(dmogrph$esp_fl)
dmogrph$foreign_born <- factor(dmogrph$foreign_born)
dmogrph$time <- as.integer(dmogrph$time)
dmogrph$attendence <- as.integer(dmogrph$attendence)

# ----Add a classroom id----
add_class_id <- function(df) {
  df <- df %>%
    mutate(
      class_id = paste(school, grade, sep = "_"),
      .after = school
    )
  df$class_id <- as.factor(df$class_id)
  return(df)
}

for (name in dfnames) {
  assign(name, add_class_id(get(name)))
}

# check all str to see if it worked
str(corsi_pre)
str(corsi_post)
str(hnf_pre)
str(hnf_post)
str(dccs_pre)
str(dccs_post)
str(precalc_pre)
str(precalc_post)
str(tejas_pre)
str(tejas_post)
str(dmogrph)

# ----Calc composite scores for HNF and DCCS (possibly use two)----
# Remove practice blocks from hnf dfs:
hnf_pre2 <- hnf_pre %>%
  select(id:total_time, pregunta_7:pregunta_18_time, 
         pregunta_25:pregunta_69_time) %>%
  mutate(across(ends_with("_time"), ~ . * 1000))

hnf_post2 <- hnf_post %>%
  select(id:total_time, pregunta_7:pregunta_18_time, 
         pregunta_25:pregunta_69_time) %>%
  mutate(across(ends_with("_time"), ~ . * 1000))

# Calculate NIH composite score for HNF

# ---- Scoring Function for H&F ----
# Primary TV Score: Block-level average
compute_hnf_tv_zelazo <- function(df) {
  
  hearts_acc   <- paste0("pregunta_", 7:18)
  hearts_time  <- paste0("pregunta_", 7:18, "_time")
  flowers_acc  <- paste0("pregunta_", 25:36)
  flowers_time <- paste0("pregunta_", 25:36, "_time")
  mixed_acc    <- paste0("pregunta_", 37:69)
  mixed_time   <- paste0("pregunta_", 37:69, "_time")
  
  compute_block_tv <- function(acc_cols, time_cols, n_trials) {
    n_correct      <- rowSums(df[, acc_cols], na.rm = FALSE)
    accuracy_score <- n_correct * (5 / n_trials)
    accuracy_pct   <- n_correct / n_trials
    
    acc_mat  <- as.matrix(df[, acc_cols])
    time_mat <- as.matrix(df[, time_cols])
    
    time_mat[acc_mat != 1] <- NA
    time_mat[time_mat < 100 | time_mat > 10000] <- NA
    time_mat[!is.na(time_mat) & time_mat < 500]  <- 500
    time_mat[!is.na(time_mat) & time_mat > 3000] <- 3000
    
    median_rt <- apply(time_mat, 1, function(x) {
      if (all(is.na(x))) NA else median(x, na.rm = TRUE)
    })
    
    log_rt   <- log10(median_rt)
    log_min  <- log10(500)
    log_max  <- log10(3000)
    rt_score <- 5 * (log_max - log_rt) / (log_max - log_min)
    rt_score <- pmax(0, pmin(5, rt_score))
    
    tv_score <- ifelse(
      is.na(accuracy_pct), NA,
      ifelse(
        accuracy_pct > 0.80 & !is.na(rt_score),
        accuracy_score + rt_score,
        accuracy_score
      )
    )
    return(tv_score)
  }
  
  df$nih_score_hearts  <- compute_block_tv(hearts_acc, hearts_time,  12)
  df$nih_score_flowers <- compute_block_tv(flowers_acc, flowers_time, 12)
  df$nih_score_mixed   <- compute_block_tv(mixed_acc, mixed_time,   33)
  df$nih_score_overall <- (df$nih_score_hearts +
                             df$nih_score_flowers +
                             df$nih_score_mixed) / 3
  return(df)
}

# ---- Sensitivity 1: Single TV score across all 57 trials ----
compute_hnf_tv_57 <- function(df) {
  
  all_acc  <- c(paste0("pregunta_", 7:18),
                paste0("pregunta_", 25:36),
                paste0("pregunta_", 37:69))
  all_time <- c(paste0("pregunta_", 7:18,  "_time"),
                paste0("pregunta_", 25:36, "_time"),
                paste0("pregunta_", 37:69, "_time"))
  
  acc_mat        <- as.matrix(df[, all_acc])
  n_correct      <- rowSums(acc_mat, na.rm = FALSE)
  accuracy_score <- n_correct * (5 / 57)
  accuracy_pct   <- n_correct / 57
  
  time_mat <- as.matrix(df[, all_time])
  time_mat[acc_mat != 1] <- NA
  time_mat[time_mat < 100 | time_mat > 10000] <- NA
  time_mat[!is.na(time_mat) & time_mat < 500]  <- 500
  time_mat[!is.na(time_mat) & time_mat > 3000] <- 3000
  
  median_rt <- apply(time_mat, 1, function(x) {
    if (all(is.na(x))) NA else median(x, na.rm = TRUE)
  })
  
  log_rt   <- log10(median_rt)
  log_min  <- log10(500)
  log_max  <- log10(3000)
  rt_score <- 5 * (log_max - log_rt) / (log_max - log_min)
  rt_score <- pmax(0, pmin(5, rt_score))
  
  df$nih_score_overall_57 <- ifelse(
    is.na(accuracy_pct), NA,
    ifelse(
      accuracy_pct >= 0.80 & !is.na(rt_score),
      accuracy_score + rt_score,
      accuracy_score
    )
  )
  return(df)
}

# Sensitivity 2: IES score
compute_hnf_ies <- function(df) {
  
  hearts_acc   <- paste0("pregunta_", 7:18)
  hearts_time  <- paste0("pregunta_", 7:18, "_time")
  flowers_acc  <- paste0("pregunta_", 25:36)
  flowers_time <- paste0("pregunta_", 25:36, "_time")
  mixed_acc    <- paste0("pregunta_", 37:69)
  mixed_time   <- paste0("pregunta_", 37:69, "_time")
  
  compute_block_ies <- function(acc_cols, time_cols, n_trials) {
    acc_mat      <- as.matrix(df[, acc_cols])
    time_mat     <- as.matrix(df[, time_cols])
    n_correct    <- rowSums(acc_mat, na.rm = FALSE)
    accuracy_pct <- n_correct / n_trials
    
    time_mat[acc_mat != 1] <- NA
    time_mat[time_mat < 100 | time_mat > 10000] <- NA
    
    mean_rt <- rowMeans(time_mat, na.rm = TRUE)
    mean_rt[is.nan(mean_rt)] <- NA
    
    ies_score <- ifelse(
      is.na(accuracy_pct) | accuracy_pct == 0,
      NA,
      mean_rt / accuracy_pct
    )
    return(ies_score)
  }
  
  df$ies_hearts  <- compute_block_ies(hearts_acc,  hearts_time,  12)
  df$ies_flowers <- compute_block_ies(flowers_acc, flowers_time, 12)
  df$ies_mixed   <- compute_block_ies(mixed_acc,   mixed_time,   33)
  df$ies_overall <- rowMeans(
    cbind(df$ies_hearts, df$ies_flowers, df$ies_mixed),
    na.rm = FALSE
  )
  return(df)
}

# ---- Apply all scoring functions ----
hnf_pre2  <- compute_hnf_tv_zelazo(hnf_pre2)
hnf_post2 <- compute_hnf_tv_zelazo(hnf_post2)

hnf_pre2  <- compute_hnf_tv_57(hnf_pre2)
hnf_post2 <- compute_hnf_tv_57(hnf_post2)

hnf_pre2  <- compute_hnf_ies(hnf_pre2)
hnf_post2 <- compute_hnf_ies(hnf_post2)

str(hnf_pre)
str(hnf_post)

# ---- Quick sanity checks ----
cat("=== Primary TV Score (Block Average) ===\n")
summary(hnf_pre2[,  c("nih_score_hearts", "nih_score_flowers",
                      "nih_score_mixed",  "nih_score_overall")])
summary(hnf_post2[, c("nih_score_hearts", "nih_score_flowers",
                      "nih_score_mixed",  "nih_score_overall")])

cat("\n=== Sensitivity: 57-Trial TV Score ===\n")
summary(hnf_pre2[,  "nih_score_overall_57"])
summary(hnf_post2[, "nih_score_overall_57"])

cat("\n=== Sensitivity: IES Score ===\n")
summary(hnf_pre2[,  c("ies_hearts", "ies_flowers", "ies_mixed", "ies_overall")])
summary(hnf_post2[, c("ies_hearts", "ies_flowers", "ies_mixed", "ies_overall")])

# ---- DCCS TV Score (unchanged) ----
compute_dccs_tv <- function(df) {
  acc_cols  <- paste0("respuesta_pregunta_", 1:24)
  time_cols <- paste0("tiempo_pregunta_", 1:24)
  
  acc_matrix     <- as.matrix(df[, acc_cols])
  n_correct      <- rowSums(acc_matrix, na.rm = FALSE)
  accuracy_score <- n_correct * (5 / 24)
  accuracy_pct   <- n_correct / 24
  
  rt_matrix <- as.matrix(df[, time_cols])
  rt_matrix[acc_matrix != 1] <- NA
  rt_matrix[rt_matrix < 100 | rt_matrix > 10000] <- NA
  rt_matrix[rt_matrix < 500  & !is.na(rt_matrix)] <- 500
  rt_matrix[rt_matrix > 3000 & !is.na(rt_matrix)] <- 3000
  
  median_rt <- apply(rt_matrix, 1, function(x) {
    if (all(is.na(x))) NA else median(x, na.rm = TRUE)
  })
  
  log_rt   <- log10(median_rt)
  log_min  <- log10(500)
  log_max  <- log10(3000)
  rt_score <- 5 * (log_max - log_rt) / (log_max - log_min)
  rt_score <- pmax(0, pmin(5, rt_score))
  
  # Note: DCCS uses >= 0.80 consistent with Zelazo original
  tv_score <- ifelse(
    is.na(accuracy_pct), NA,
    ifelse(
      accuracy_pct >= 0.80 & !is.na(rt_score),
      accuracy_score + rt_score,
      accuracy_score
    )
  )
  
  df$dccs_tv_score <- tv_score
  return(df)
}

compute_dccs_descriptives <- function(df) {
  acc_cols  <- paste0("respuesta_pregunta_", 1:24)
  time_cols <- paste0("tiempo_pregunta_", 1:24)
  
  acc_matrix <- as.matrix(df[, acc_cols])
  rt_matrix  <- as.matrix(df[, time_cols])
  
  df$dccs_mean_acc <- rowMeans(acc_matrix, na.rm = FALSE)
  
  rt_matrix[acc_matrix != 1] <- NA
  rt_matrix[rt_matrix < 100 | rt_matrix > 10000] <- NA
  df$dccs_mean_rt <- apply(rt_matrix, 1, function(x) {
    if (all(is.na(x))) NA else mean(x, na.rm = TRUE)
  })
  return(df)
}

dccs_pre  <- compute_dccs_tv(dccs_pre)
dccs_post <- compute_dccs_tv(dccs_post)
dccs_pre  <- compute_dccs_descriptives(dccs_pre)
dccs_post <- compute_dccs_descriptives(dccs_post)

cat("\n=== DCCS TV Score ===\n")
summary(dccs_pre[,  "dccs_tv_score"])
summary(dccs_post[, "dccs_tv_score"])

# ---- H&F Descriptives (unchanged) ----
compute_hnf_descriptives <- function(df) {
  df <- df %>%
    mutate(
      time_seconds_hearts = time_seconds_hearts * 1000,
      time_seconds_flowers = time_seconds_flowers * 1000,
      time_seconds_heart_flowers = time_seconds_heart_flowers * 1000
    )
  
  df <- df %>%
    mutate(
      acc_prop_hearts  = score_hearts        / 12,
      acc_prop_flowers = score_flowers       / 12,
      acc_prop_mixed   = score_heart_flowers / 33,
      acc_prop_overall = hnf_total           / 57
    )
  
  hearts_time  <- paste0("pregunta_", 7:18,  "_time")
  flowers_time <- paste0("pregunta_", 25:36, "_time")
  mixed_time   <- paste0("pregunta_", 37:69, "_time")
  all_time     <- c(hearts_time, flowers_time, mixed_time)
  
  df$rt_mean_hearts <- rowMeans(as.matrix(df[, hearts_time]), 
                                na.rm = TRUE)
  df$rt_mean_flowers <- rowMeans(as.matrix(df[, flowers_time]), 
                                 na.rm = TRUE)
  df$rt_mean_mixed <- rowMeans(as.matrix(df[, mixed_time]), 
                               na.rm = TRUE)
  df$rt_mean_overall <- rowMeans(as.matrix(df[, all_time]), 
                                 na.rm = TRUE)
  
  df$rt_sd_hearts <- apply(as.matrix(df[, hearts_time]), 1, sd, 
                           na.rm = TRUE)
  df$rt_sd_flowers <- apply(as.matrix(df[, flowers_time]), 1, sd, 
                            na.rm = TRUE)
  df$rt_sd_mixed <- apply(as.matrix(df[, mixed_time]), 1, sd, 
                          na.rm = TRUE)
  df$rt_sd_overall <- apply(as.matrix(df[, all_time]), 1, sd, 
                            na.rm = TRUE)
  
  return(df)
}

hnf_pre2  <- compute_hnf_descriptives(hnf_pre2)
hnf_post2 <- compute_hnf_descriptives(hnf_post2)

# ---- Tejas letter score (unchanged) ----
tejas_pre$puntaje_ambito_letra  <- tejas_pre$puntaje_ambito_3  +
  tejas_pre$puntaje_ambito_5
tejas_post$puntaje_ambito_letra <- tejas_post$puntaje_ambito_3 +
  tejas_post$puntaje_ambito_5

# ---- Variable selection ----
corsi_pre_sel <- corsi_pre %>%
  select(id, total_points, fwd_att, bwd_att) %>%
  rename_with(~ paste0(., "_pre"), -id)

corsi_post_sel <- corsi_post %>%
  select(id, total_points, fwd_att, bwd_att) %>%
  rename_with(~ paste0(., "_post"), -id)

hnf_pre_sel <- hnf_pre2 %>%
  select(id,
         nih_score_overall, nih_score_hearts, nih_score_flowers,
         nih_score_mixed, nih_score_overall_57,
         ies_hearts, ies_flowers, ies_mixed, ies_overall,
         hnf_total, score_hearts, score_flowers, score_heart_flowers,
         time_seconds_hearts, time_seconds_flowers,
         time_seconds_heart_flowers, total_time,
         acc_prop_hearts, acc_prop_flowers, acc_prop_mixed, acc_prop_overall,
         rt_mean_hearts, rt_mean_flowers, rt_mean_mixed, rt_mean_overall,
         rt_sd_hearts, rt_sd_flowers, rt_sd_mixed, rt_sd_overall) %>%
  rename_with(~ paste0(., "_pre"), -id)

hnf_post_sel <- hnf_post2 %>%
  select(id,
         nih_score_overall, nih_score_hearts, nih_score_flowers,
         nih_score_mixed, nih_score_overall_57,
         ies_hearts, ies_flowers, ies_mixed, ies_overall,
         hnf_total, score_hearts, score_flowers, score_heart_flowers,
         time_seconds_hearts, time_seconds_flowers,
         time_seconds_heart_flowers, total_time,
         acc_prop_hearts, acc_prop_flowers, acc_prop_mixed, acc_prop_overall,
         rt_mean_hearts, rt_mean_flowers, rt_mean_mixed, rt_mean_overall,
         rt_sd_hearts, rt_sd_flowers, rt_sd_mixed, rt_sd_overall) %>%
  rename_with(~ paste0(., "_post"), -id)

dccs_pre_sel <- dccs_pre %>%
  select(id, dccs_tv_score, dccs_mean_acc, dccs_mean_rt) %>%
  rename_with(~ paste0(., "_pre"), -id)

dccs_post_sel <- dccs_post %>%
  select(id, dccs_tv_score, dccs_mean_acc, dccs_mean_rt) %>%
  rename_with(~ paste0(., "_post"), -id)

precalc_pre_sel <- precalc_pre %>%
  select(id, puntaje_total, puntaje_ambito_1, puntaje_ambito_2,
         puntaje_ambito_3) %>%
  rename_with(~ paste0("precalc_", ., "_pre"), -id)

precalc_post_sel <- precalc_post %>%
  select(id, puntaje_total, puntaje_ambito_1, puntaje_ambito_2,
         puntaje_ambito_3) %>%
  rename_with(~ paste0("precalc_", ., "_post"), -id)

tejas_pre_sel <- tejas_pre %>%
  select(id, puntaje_total, puntaje_ambito_1, puntaje_ambito_2,
         puntaje_ambito_3, puntaje_ambito_4, puntaje_ambito_5,
         puntaje_ambito_6, puntaje_ambito_7, puntaje_ambito_letra) %>%
  rename_with(~ paste0("tejas_", ., "_pre"), -id)

tejas_post_sel <- tejas_post %>%
  select(id, puntaje_total, puntaje_ambito_1, puntaje_ambito_2,
         puntaje_ambito_3, puntaje_ambito_4, puntaje_ambito_5,
         puntaje_ambito_6, puntaje_ambito_7, puntaje_ambito_letra) %>%
  rename_with(~ paste0("tejas_", ., "_post"), -id)

dmogrph_sel <- dmogrph %>%
  select(id, sex, grade, school, class_id, condition, status,
         age_mths, ed_level, esp_fl, foreign_born,
         activity, name, time, attendence)

# ---- Merge ----
res <- dmogrph_sel %>%
  left_join(corsi_pre_sel,    by = "id") %>%
  left_join(corsi_post_sel,   by = "id") %>%
  left_join(hnf_pre_sel,      by = "id") %>%
  left_join(hnf_post_sel,     by = "id") %>%
  left_join(dccs_pre_sel,     by = "id") %>%
  left_join(dccs_post_sel,    by = "id") %>%
  left_join(precalc_pre_sel,  by = "id") %>%
  left_join(precalc_post_sel, by = "id") %>%
  left_join(tejas_pre_sel,    by = "id") %>%
  left_join(tejas_post_sel,   by = "id")

# ---- Gain scores ----
res <- res %>%
  mutate(
    # Corsi
    gain_corsi = total_points_post - total_points_pre,
    
    # H&F primary
    gain_hnf_overall = nih_score_overall_post - nih_score_overall_pre,
    gain_hnf_hearts = nih_score_hearts_post - nih_score_hearts_pre,
    gain_hnf_flowers = nih_score_flowers_post - nih_score_flowers_pre,
    gain_hnf_mixed = nih_score_mixed_post - nih_score_mixed_pre,
    
    # H&F sensitivity
    gain_hnf_57  = nih_score_overall_57_post - nih_score_overall_57_pre,
    gain_hnf_ies = ies_overall_post - ies_overall_pre,
    
    # DCCS
    gain_dccs = dccs_tv_score_post - dccs_tv_score_pre,
    
    # Precalc
    gain_precalc = precalc_puntaje_total_post - precalc_puntaje_total_pre,
    
    # Tejas
    gain_tejas = tejas_puntaje_total_post - tejas_puntaje_total_pre
  )

# ---- Final checks ----
cat("\n=== Gain Score Summary ===\n")
summary(res[, c("gain_corsi", "gain_hnf_overall", "gain_hnf_57",
                "gain_hnf_ies", "gain_dccs",
                "gain_precalc", "gain_tejas")])

cat("\nDimensions of merged dataframe:", dim(res), "\n")
cat("Column names:\n")
print(names(res))
# ----Save dfs and workspace-----
# may be commented out to avoid overwrite
write.csv(corsi_pre, "data-clean/corsi_pre_c.csv", row.names = FALSE)
write.csv(corsi_post, "data-clean/corsi_post_c.csv", row.names = FALSE)
write.csv(hnf_pre2, "data-clean/hnf_pre_c.csv", row.names = FALSE)
write.csv(hnf_post2, "data-clean/hnf_post_c.csv", row.names = FALSE)
write.csv(dccs_pre, "data-clean/dccs_pre_c.csv", row.names = FALSE)
write.csv(dccs_post, "data-clean/dccs_post_c.csv", row.names = FALSE)
write.csv(precalc_pre, "data-clean/precalc_pre_c.csv",row.names = FALSE)
write.csv(precalc_post, "data-clean/precalc_post_c.csv", row.names = FALSE)
write.csv(tejas_pre, "data-clean/tejas_pre_c.csv", row.names = FALSE)
write.csv(tejas_post, "data-clean/tejas_post_c.csv", row.names = FALSE)
write.csv(dmogrph, "data-clean/dmo_c.csv", row.names = FALSE)
write.csv(res, "data-clean/res.csv", row.names = FALSE)

save.image(file = "data-clean/may_28.RData")
