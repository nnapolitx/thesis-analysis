library(dplyr)
library(pwr)
library(knitr)

# It is not necessary to save any objects from this script for the
# rest of analysis scripts to correctly.

# --- Load and factor data ----
res_analysis <- read.csv("data-clean/res_analysis.csv")
res_analysis <- res_analysis %>%
  mutate(
    grade = factor(grade),
    condition = factor(condition),
    cond_grade = factor(cond_grade),
    status = factor(status),
    condition2 = factor(condition2),
    cond_grade2 = factor(cond_grade2)
  )
res_analysis$condition <- relevel(res_analysis$condition, 
                                  ref = "control_1")
res_analysis$cond_grade <- relevel(res_analysis$cond_grade, 
                                   ref = "control_1_fifth")
res_analysis$condition2 <- relevel(res_analysis$condition2, 
                                  ref = "control")
res_analysis$cond_grade2 <- relevel(res_analysis$cond_grade2, 
                                   ref = "control_fifth")
# ---- Sample sizes ----
n_total <- sum(res_analysis$status == "participant")
n_kinder <- sum(res_analysis$status == "participant" & 
                      res_analysis$grade  == "kinder")
n_fifth <- sum(res_analysis$status == "participant" & 
                      res_analysis$grade  == "fifth")

cat("Total participants:", n_total, "\n")
cat("Kindergarten n:", n_kinder, "\n")
cat("Fifth grade n:", n_fifth, "\n\n")

# ---- Model parameters ----
# Full model: condition (4 levels) x grade (2 levels) + 5 pretest covariates
# Numerator df for condition main effect = k - 1 = 3
# Numerator df for grade main effect = 2 - 1 = 1
# Numerator df for condition * grade interaction = (4-1)*(2-1) = 3
# Total predictors in model = 4 (condition) + 1 (grade level) +
# 3 (interactions of grade*condition) + 5 covariates (pretests) = 12
# Denominator df = n - total predictors (12) - 1

n_pred_full    <- 12   # 3 conditions + 1 grade + 3 interaction + 5 covariates
n_pred_subgroup <- 8   # 3 condition + 5 covariates

# ---- STEP 15: Minimum detectable effect size at 80% power ---- 
# Full sample model
# Using condition main effect as the focal test (df = 3)
power_full <- pwr.f2.test(
  u         = 3,
  v         = n_total - n_pred_full - 1,
  sig.level = 0.05,
  power     = 0.80
)

cat("=== Full sample (n =", n_total, ") ===\n")
cat("Denominator df:", n_total - n_pred_full - 1, "\n")
cat("Minimum detectable f2:", round(power_full$f2, 4), "\n")
cat("Equivalent partial eta2:", 
    round(power_full$f2 / (1 + power_full$f2), 4), "\n")
cat("Approximate Cohen's d:", 
    round(2 * sqrt(power_full$f2), 4), "\n\n")

# ---- Kindergarten subgroup ----
# 4 conditions present in kinder, condition df = 3
power_kinder <- pwr.f2.test(
  u         = 3,
  v         = n_kinder - n_pred_subgroup - 1,
  sig.level = 0.05,
  power     = 0.80
)

cat("=== Kindergarten subgroup (n =", n_kinder, ") ===\n")
cat("Denominator df:", n_kinder - n_pred_subgroup - 1, "\n")
cat("Minimum detectable f2:", round(power_kinder$f2, 4), "\n")
cat("Equivalent partial eta2:", 
    round(power_kinder$f2 / (1 + power_kinder$f2), 4), "\n")
cat("Approximate Cohen's d:", 
    round(2 * sqrt(power_kinder$f2), 4), "\n\n")

# ---- Fifth grade subgroup ----
# 4 conditions in fifth grade, condition df = 3
power_fifth <- pwr.f2.test(
  u         = 3,
  v         = n_fifth - n_pred_subgroup - 1,
  sig.level = 0.05,
  power     = 0.80
)

cat("=== Fifth grade subgroup (n =", n_fifth, ") ===\n")
cat("Denominator df:", n_fifth - n_pred_subgroup - 1, "\n")
cat("Minimum detectable f2:", round(power_fifth$f2, 4), "\n")
cat("Equivalent partial eta2:", 
    round(power_fifth$f2 / (1 + power_fifth$f2), 4), "\n")
cat("Approximate Cohen's d:", 
    round(2 * sqrt(power_fifth$f2), 4), "\n\n")

# ---- STEP 16: Sensitivity analysis ----
# What power do we have for small, medium, and large effects?
# Using Cohen's conventions: f2 = 0.02 (small), 0.15 (medium), 
# 0.35 (large)

effect_sizes <- c(0.02, 0.15, 0.35)
effect_labels <- c("Small (f2 = 0.02)", 
                   "Medium (f2 = 0.15)", 
                   "Large (f2 = 0.35)")

sensitivity_results <- data.frame(
  Effect_Size  = effect_labels,
  Power_Full   = sapply(effect_sizes, function(f2) {
    round(pwr.f2.test(u = 3,
                      v = n_total - n_pred_full - 1,
                      f2 = f2,
                      sig.level = 0.05)$power, 3)
  }),
  Power_Kinder = sapply(effect_sizes, function(f2) {
    round(pwr.f2.test(u = 3,
                      v = n_kinder - n_pred_subgroup - 1,
                      f2 = f2,
                      sig.level = 0.05)$power, 3)
  }),
  Power_Fifth  = sapply(effect_sizes, function(f2) {
    round(pwr.f2.test(u = 3,
                      v = n_fifth - n_pred_subgroup - 1,
                      f2 = f2,
                      sig.level = 0.05)$power, 3)
  })
)

kable(sensitivity_results,
      col.names = c("Effect Size", 
                    paste0("Full Sample (n = ", n_total, ")"),
                    paste0("Kindergarten (n = ", n_kinder, ")"),
                    paste0("Fifth Grade (n = ", n_fifth, ")")),
      caption = "Table 7. Sensitivity analysis: observed power by effect size and sample")

# ---- Summary table: minimum detectable effects ----

mde_summary <- data.frame(
  Sample = c(paste0("Full sample (n = ", n_total, ")"),
             paste0("Kindergarten (n = ", n_kinder, ")"),
             paste0("Fifth grade (n = ", n_fifth, ")")),
  df_num = c(3, 3, 3),
  df_den = c(n_total  - n_pred_full    - 1,
             n_kinder - n_pred_subgroup - 1,
             n_fifth  - n_pred_subgroup - 1),
  min_f2    = c(round(power_full$f2,    4),
                round(power_kinder$f2,  4),
                round(power_fifth$f2,   4)),
  min_eta2  = c(round(power_full$f2    / (1 + power_full$f2),   4),
                round(power_kinder$f2  / (1 + power_kinder$f2), 4),
                round(power_fifth$f2   / (1 + power_fifth$f2),  4)),
  min_d     = c(round(2 * sqrt(power_full$f2),   4),
                round(2 * sqrt(power_kinder$f2), 4),
                round(2 * sqrt(power_fifth$f2),  4))
)

kable(mde_summary,
      col.names = c("Sample", "df (numerator)", "df (denominator)",
                    "Min. f2", "Min. partial eta2", 
                    "Approx. Cohen's d"),
      caption = "Table 2. Minimum detectable effect sizes at 80% power (alpha = .05)")


# ---- Updated minimum detectable effects summary ----
mde_summary_v2 <- data.frame(
  Sample = c(
    paste0("Full sample (n = ", n_total, ")"),
    paste0("Kindergarten full (n = ", n_kinder, ")"),
    paste0("Fifth grade (n = ", n_fifth, ")")
  ),
  df_num = c(4, 4, 3, 3),
  df_den = c(
    n_total        - n_pred_full     - 1,
    n_kinder       - n_pred_subgroup - 1,
    n_fifth        - n_pred_subgroup - 1
  ),
  min_f2 = c(
    round(power_full$f2,          4),
    round(power_kinder$f2,        4),
    round(power_fifth$f2,         4)
  ),
  min_eta2 = c(
    round(power_full$f2          / (1 + power_full$f2),          4),
    round(power_kinder$f2        / (1 + power_kinder$f2),        4),
    round(power_fifth$f2         / (1 + power_fifth$f2),         4)
  ),
  min_d = c(
    round(2 * sqrt(power_full$f2),          4),
    round(2 * sqrt(power_kinder$f2),        4),
    round(2 * sqrt(power_fifth$f2),         4)
  )
)

kable(mde_summary_v2,
      col.names = c("Sample", "df (numerator)", "df (denominator)",
                    "Min. f2", "Min. partial eta2",
                    "Approx. Cohen's d"),
      caption = "Table 4. Minimum detectable effect sizes at 80% power (alpha = .05)")
# ---- Consort style Flow Diagram ----
# Figure 1: Participant Flow Diagram
# CONSORT-style recruitment and attrition flowchart
# Saved to plots/fig1_flowdiagram.png at 600 dpi

library(ggplot2)
library(grid)

create_flow_diagram <- function() {
  
  png("plots/fig1_flowdiagram.png",
      width  = 18,
      height = 24,
      units  = "cm",
      res    = 600,
      bg     = "white")
  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 1)))
  pushViewport(viewport(x = 0.5, y = 0.5, 
                        width = 0.95, height = 0.95))
  
  draw_box <- function(x, y, w, h, label, 
                       fill = "white", 
                       border = "black",
                       fontsize = 8,
                       bold = FALSE) {
    grid.rect(x = x, y = y, width = w, height = h,
              gp = gpar(fill = fill, col = border, lwd = 1.2),
              default.units = "npc")
    grid.text(label,
              x = x, y = y,
              gp = gpar(fontsize = fontsize,
                        fontface = if (bold) "bold" else "plain",
                        col = "black"),
              default.units = "npc",
              just = "center")
  }
  
  draw_arrow <- function(x1, y1, x2, y2) {
    grid.lines(x = c(x1, x2), y = c(y1, y2),
               gp = gpar(col = "black", lwd = 1.2),
               arrow = arrow(length = unit(0.15, "cm"),
                             type = "closed",
                             angle = 20),
               default.units = "npc")
  }
  
  draw_line <- function(x1, y1, x2, y2) {
    grid.lines(x = c(x1, x2), y = c(y1, y2),
               gp = gpar(col = "black", lwd = 1.0),
               default.units = "npc")
  }
  
  # Layout parameters
  bw <- 0.55
  bh <- 0.065
  sw <- 0.30
  cx <- 0.50
  rx <- 0.82
  
  # Y positions (top to bottom) — y2 removed, remaining boxes shifted up evenly
  y1 <- 0.90   # Assessed for eligibility (was y2 content, now first box)
  y3 <- 0.72   # Consented
  y4 <- 0.54   # Intervention period
  y5 <- 0.38   # Posttest assessment
  y6 <- 0.20   # Final analysis sample
  
  # Exclusion box y positions — midpoints between adjacent main boxes
  ye1 <- (y1 - bh/2 + y3 + bh/2) / 2   # between box 1 and box 3
  ye2 <- (y3 - bh/2 + y4 + bh/2) / 2   # between box 3 and box 4
  ye3 <- (y4 - bh/2 + y5 + bh/2) / 2   # between box 4 and box 5
  
  # Draw main boxes
  
  # Box 1: Assessed for eligibility (replaces old y1 + y2)
  draw_box(cx, y1, bw, bh,
           "Assessed for eligibility and invited\n to participate: n = 222",
           fill = "#D6EAF8", bold = TRUE, fontsize = 12)
  
  # Box 3: Consented
  draw_box(cx, y3, bw, bh,
           "Consented and began pretest assessment\nn = 194  (kindergarten: 98, fifth grade: 96)",
           fill = "white", fontsize = 11)
  
  # Box 4: Intervention period
  draw_box(cx, y4, bw, bh,
           "Participated in intervention period\n(or control condition)\nn = 188",
           fill = "white", fontsize = 11)
  
  # Box 5: Posttest
  draw_box(cx, y5, bw, bh,
           "Completed posttest assessment\nn = 185",
           fill = "white", fontsize = 11)
  
  # Box 6: Final analysis sample
  draw_box(cx, y6, bw, bh + 0.04,
           paste0("Final analysis sample: n = 185\n",
                  "Kindergarten: n = 91  |  Fifth grade: n = 94\n",
                  "control 1: 37  |  control 2: 14  |  digital: 45\n",
                  "paper: 46  |  mixed: 43"),
           fill = "#D5F5E3", bold = TRUE, fontsize = 11)
  
  # Draw exclusion boxes
  draw_box(rx, ye1, sw, bh + 0.01,
           "Declined participation\n(child or parent refusal)\nn = 28",
           fill = "#FADBD8", fontsize = 10)
  
  draw_box(rx, ye2, sw, bh + 0.01,
           "Excluded: IEP diagnosis\n(neuroatypical development)\nn = 6",
           fill = "#FADBD8", fontsize = 10)
  
  draw_box(rx, ye3, sw, bh + 0.01,
           "Dropped out during study\n(MCAR confirmed, p = .116)\nn = 3",
           fill = "#FADBD8", fontsize = 10)
  
  # Main vertical arrows
  draw_arrow(cx, y1 - bh/2, cx, y3 + bh/2)
  draw_arrow(cx, y3 - bh/2, cx, y4 + bh/2)
  draw_arrow(cx, y4 - bh/2, cx, y5 + bh/2)
  draw_arrow(cx, y5 - bh/2, cx, y6 + (bh + 0.04)/2)
  
  # Exclusion connectors — horizontal lines from vertical arrow midpoint to red box left edge
  draw_line(cx, ye1, rx - sw/2, ye1)
  draw_line(cx, ye2, rx - sw/2, ye2)
  draw_line(cx, ye3, rx - sw/2, ye3)
  
  # Title
  grid.text("Figure 1. Participant Flow Diagram",
            x = 0.5, y = 0.99,
            gp = gpar(fontsize = 12, fontface = "bold"),
            default.units = "npc",
            just = "center")
  
  # Footer note
  grid.text(
    paste0("Note. MCAR = missing completely at random. ",
           "IEP = Individualized Education Plan. ",
           "Atypically\ndeveloping participants (n = 6) are excluded from primary ",
           "analyses. ",
           "Control 2 condition was\nassigned to kindergarten only ",
           "due to low enrollment at this grade level."),
    x = 0.5, y = 0.06,
    gp = gpar(fontsize = 11, fontface = "italic", col = "grey30"),
    default.units = "npc",
    just = "center")
  
  dev.off()
  cat("Flow diagram saved to plots/fig1_flowdiagram.png\n")
}

create_flow_diagram()

# ---- SAVE ----
save.image(file = "rdata_files/may_28_pwr_complete.RData")

