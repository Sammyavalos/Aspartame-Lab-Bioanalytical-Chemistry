# Aspartame calibration
# Stock: 0.0589 g aspartame in 20 mL water (10 mM).
# Serial dilutions:
# 14 mL of 10 mM + 3.5 mL water = 8 mM
# 12 mL of 8 mM + 4 mL water = 6 mM
# 10 mL of 6 mM + 5 mL water = 4 mM
# 8 mL of 4 mM + 8 mL water = 2 mM
# 6 mL of 2 mM + 6 mL water = 1 mM
# 4 mL of 1 mM + 4 mL water = 0.5 mM

concentration <- c(0.5, 1, 2, 4, 6, 8, 10)  # mM
absorbance <- c(0.121, 0.188, 0.342, 0.494, 0.677, 0.832, 0.920)

cal_model <- lm(absorbance ~ concentration)
cal_summary <- summary(cal_model)

# Bars show residual SE, not replicate SD (one reading per standard).
s_y <- cal_summary$sigma

plot(concentration, absorbance,
     type = "p",
     pch = 19,
     col = "darkblue",
     ylim = c(min(absorbance - s_y), max(absorbance + s_y)),
     main = "Aspartame Calibration Curve",
     xlab = "Aspartame concentration (mM)",
     ylab = "Absorbance",
     las = 1)

coefficients(cal_model)

abline(cal_model, col = "red", lwd = 2, lty = 1)

arrows(x0 = concentration, y0 = absorbance - s_y,
       x1 = concentration, y1 = absorbance + s_y,
       angle = 90, code = 3, length = 0.05, col = "darkblue")


legend("topleft",
       legend = c("Calibration standards", "Linear fit"),
       col = c("darkblue", "red"),
       lty = c(NA, 1),
       pch = c(19, NA),
       lwd = c(1, 2),
       bty = "n")

legend("bottomright",
       legend = c(
         sprintf("A = %.5f C + %.5f",
                 coef(cal_model)[2], coef(cal_model)[1]),
         sprintf("R² = %.4f", cal_summary$r.squared),
         sprintf("Slope relative SE = %.2f%%",
                 cal_summary$coefficients["concentration", "Std. Error"] /
                   coef(cal_model)[["concentration"]] * 100)
       ),
       bty = "n", cex = 0.85)

cal_summary

# Molar absorptivity, L/(mol*cm)
# Multiply the slope by 1000 to convert from mM to M.
path_length_cm <- 1  # Assumes a 1 cm cuvette.

epsilon <- coef(cal_model)[["concentration"]] *
  1000 / path_length_cm

epsilon_SE <- cal_summary$coefficients[
  "concentration", "Std. Error"
] * 1000 / path_length_cm

regression_SE <- cal_summary$sigma
R_squared <- cal_summary$r.squared

coef(cal_model)

c(
  molar_absorptivity = epsilon,
  molar_absorptivity_SE = epsilon_SE,
  regression_SE = regression_SE,
  R_squared = R_squared
)

slope_relative_SE_pct <- epsilon_SE / epsilon * 100
slope_relative_SE_pct

# Crystal Light Zero Sugar Lemonade
# Stock: 1.00 g powder in 20 mL (50 mg/mL).
# Stock volumes diluted to 5 mL:
# 5 mL stock       = 50 mg/mL
# 2.5 mL stock     = 25 mg/mL
# 1.25 mL stock    = 12.5 mg/mL
# 0.625 mL stock   = 6.25 mg/mL
# 0.3125 mL stock  = 3.125 mg/mL
# 0.15625 mL stock = 1.5625 mg/mL

CL_packet_mass_g <- c(4.2969, 5.0558, 4.0565)
CL_mass_mean <- mean(CL_packet_mass_g)
CL_mass_sd <- sd(CL_packet_mass_g)
CL_mass_rsd <- CL_mass_sd / CL_mass_mean * 100

CL_mass_mean  # g
CL_mass_sd    # g
CL_mass_rsd   # %

CL_stock_mass_g <- 1.00
CL_stock_vol_mL <- 20
CL_stock_mg_per_mL <- CL_stock_mass_g * 1000 / CL_stock_vol_mL
CL_stock_mg_per_mL

# Packet 1 was measured at all six concentrations.
# Packets 2 and 3 were measured at 3.125 mg/mL.
CL_data <- data.frame(
  packet = c(1, 1, 1, 1, 1, 1, 2, 3),
  sweet_mg_mL = c(50, 25, 12.5, 6.25, 3.125, 1.5625, 3.125, 3.125),
  absorbance = c(0.918, 0.937, 0.874, 0.952, 0.467, 0.425, 0.452, 0.458)
)

b0 <- coef(cal_model)[["(Intercept)"]]
b1 <- coef(cal_model)[["concentration"]]
CL_data$asp_mM <- (CL_data$absorbance - b0) / b1

MW_asp <- 294.30
CL_data$asp_M <- CL_data$asp_mM / 1000
CL_data$asp_mg_mL <- CL_data$asp_M * MW_asp

CL_data$asp_mg_in_5mL <- CL_data$asp_mg_mL * 5

CL_data$mg_per_g <- CL_data$asp_mg_mL / CL_data$sweet_mg_mL * 1000

CL_data$in_cal_range <- CL_data$absorbance >= min(absorbance) &
  CL_data$absorbance <= max(absorbance)

CL_data

# Compare packets at 3.125 mg/mL.
CL_reps <- CL_data[CL_data$sweet_mg_mL == 3.125, ]
CL_reps

n_CL <- nrow(CL_reps)
CL_mean <- mean(CL_reps$mg_per_g)
CL_sd <- sd(CL_reps$mg_per_g)

# 95% CI across the three packets
t_crit_CL <- qt(0.975, df = n_CL - 1)
CL_CI_half <- t_crit_CL * CL_sd / sqrt(n_CL)

CL_CI <- c(
  lower = CL_mean - CL_CI_half,
  upper = CL_mean + CL_CI_half
)

c(n = n_CL, mean_mg_per_g = CL_mean, SD = CL_sd, t_crit = t_crit_CL)
CL_CI

# t.test(CL_reps$mg_per_g, mu = REPORTED)

# Shared Bumble Bee calibration
# Class data from the same spectrophotometer.
# Assumed all groups used our standard concentrations;
shared_conc <- c(0.5, 1, 2, 4, 6, 8, 10)
shared_abs <- rbind(
  Tue_Greenlee_Parra  = c(0.149, 0.256, 0.260, 0.540, 1.118, 1.344, 1.713),
  Wed_Zuinga_Avalos   = c(0.121, 0.188, 0.342, 0.494, 0.677, 0.832, 0.920),
  Thu_Nikolas_Malia   = c(0.136, 0.140, 0.348, 0.495, 0.666, 0.847, 0.928),
  Wed_Boothby_Ciccoli = c(0.077, 0.133, 0.236, 0.473, 0.675, 0.828, 0.940)
)

shared <- data.frame(
  group = rep(rownames(shared_abs), each = length(shared_conc)),
  conc  = rep(shared_conc, times = nrow(shared_abs)),
  A     = as.vector(t(shared_abs))
)

grp <- factor(shared$group)
plot(shared$conc, shared$A, pch = 19, col = as.integer(grp),
     main = "Bumble Bee Calibration: Individual Readings",
     xlab = "Aspartame concentration (mM)", ylab = "Absorbance", las = 1)
legend("topleft", legend = levels(grp), col = seq_along(levels(grp)),
       pch = 19, bty = "n", cex = 0.8)

# Iterative Grubbs test on pooled residuals (alpha = 0.05)
shared_kept <- shared
repeat {
  fit_tmp <- lm(A ~ conc, data = shared_kept)
  res_tmp <- residuals(fit_tmp)
  n_tmp   <- length(res_tmp)
  G       <- max(abs(res_tmp - mean(res_tmp))) / sd(res_tmp)
  t_g     <- qt(1 - 0.05 / (2 * n_tmp), df = n_tmp - 2)
  G_crit  <- (n_tmp - 1) / sqrt(n_tmp) * sqrt(t_g^2 / (n_tmp - 2 + t_g^2))
  print(c(n = n_tmp, G = G, G_crit = G_crit))
  if (G <= G_crit) break
  worst <- which.max(abs(res_tmp - mean(res_tmp)))
  print(shared_kept[worst, ])          # reading removed as an outlier
  shared_kept <- shared_kept[-worst, ]
}
# Retained readings: mean and SD at each concentration
shared_mean <- tapply(shared_kept$A, shared_kept$conc, mean)
shared_sd   <- tapply(shared_kept$A, shared_kept$conc, sd)

shared_model   <- lm(A ~ conc, data = shared_kept)
shared_summary <- summary(shared_model)

shared_epsilon    <- coef(shared_model)[["conc"]] * 1000 / path_length_cm
shared_epsilon_SE <- shared_summary$coefficients["conc", "Std. Error"] *
  1000 / path_length_cm
c(epsilon = shared_epsilon, epsilon_SE = shared_epsilon_SE,
  regression_SE = shared_summary$sigma, R_squared = shared_summary$r.squared)

plot(shared_conc, shared_mean, pch = 19, col = "darkblue",
     ylim = c(0, max(shared_mean + shared_sd)),
     main = "Aspartame Calibration Curve (Bumble Bee Group Data)",
     xlab = "Aspartame concentration (mM)", ylab = "Absorbance", las = 1)
arrows(shared_conc, shared_mean - shared_sd,
       shared_conc, shared_mean + shared_sd,
       angle = 90, code = 3, length = 0.05, col = "darkblue")
abline(shared_model, col = "red", lwd = 2)
legend("topleft", legend = c("Mean of groups (± SD)", "Linear fit"),
       col = c("darkblue", "red"), lty = c(NA, 1), pch = c(19, NA),
       lwd = c(1, 2), bty = "n")
legend("bottomright", bty = "n", cex = 0.85,
       legend = c(sprintf("A = %.5f C + %.5f",
                          coef(shared_model)[2], coef(shared_model)[1]),
                  sprintf("R² = %.4f", shared_summary$r.squared)))

# Crystal Light results using the shared calibration
CL_shared_mM <- (CL_reps$absorbance - coef(shared_model)[["(Intercept)"]]) /
  coef(shared_model)[["conc"]]
CL_shared_mg_per_g <- CL_shared_mM / 1000 * MW_asp / CL_reps$sweet_mg_mL * 1000

CL_shared_mean <- mean(CL_shared_mg_per_g)
CL_shared_sd   <- sd(CL_shared_mg_per_g)
CL_shared_half <- qt(0.975, df = n_CL - 1) * CL_shared_sd / sqrt(n_CL)
c(mean = CL_shared_mean, SD = CL_shared_sd,
  lower = CL_shared_mean - CL_shared_half,
  upper = CL_shared_mean + CL_shared_half)

c(CL_shared_mean, CL_shared_sd)

# No manufacturer aspartame value was available, so a t-test
#was not performed; the 95% CI is reported instead.