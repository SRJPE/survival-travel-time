# ============================================================================
# SIMPLE JOINT ACOUSTIC TELEMETRY + CWT MODEL
# R + JAGS, self-contained simulation
#
# Goal:
#   Estimate ONE survival probability and ONE travel-time distribution
#   using two observation systems:
#
#      acoustic telemetry --> survival x acoustic detection
#      coded-wire tags    --> survival x trawl/capture probability
#
# Deliberately omitted:
#   routing, reaches, flow, fish size, year effects, random effects,
#   time-varying parameters, and tag failure.
#
# Packages needed:
#   install.packages(c("rjags", "coda"))
#
# System requirement:
#   JAGS must also be installed.
# ============================================================================

rm(list = ls())
set.seed(20260930)

library(rjags)
library(coda)

# ----------------------------------------------------------------------------
# OUTPUT DIRECTORY
# ----------------------------------------------------------------------------

out_dir <- "C:/Users/owner/Documents/20260924_flora_josh_mortality/sims_joint_draft_20260930/jags"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# 1. TRUTH AND DESIGN
# ----------------------------------------------------------------------------

truth <- list(
  S = 0.50,             # survival release -> endpoint
  median_T = 12,        # median travel time, days
  sigma_logT = 0.40,    # SD of log travel time
  p_AT = 0.90,          # KNOWN acoustic detection probability at endpoint
  p_CWT = 0.003         # KNOWN CWT capture probability at endpoint
)

N_AT  <- 250
N_CWT <- 60000

mu_logT <- log(truth$median_T)

# ----------------------------------------------------------------------------
# 2. SIMULATE ACOUSTIC TELEMETRY DATA
# ----------------------------------------------------------------------------
# A fish must survive AND be detected to contribute an endpoint detection.
# Travel time is generated from the same biological distribution for all fish.

alive_AT <- rbinom(N_AT, 1, truth$S)
T_AT_all <- rlnorm(N_AT, mu_logT, truth$sigma_logT)
det_AT   <- rbinom(N_AT, 1, truth$p_AT)

seen_AT <- alive_AT == 1 & det_AT == 1
T_AT    <- T_AT_all[seen_AT]
n_AT_obs <- length(T_AT)

# ----------------------------------------------------------------------------
# 3. SIMULATE CWT DATA
# ----------------------------------------------------------------------------
# A fish must survive AND be captured in the endpoint trawl to be recovered.

alive_CWT <- rbinom(N_CWT, 1, truth$S)
T_CWT_all <- rlnorm(N_CWT, mu_logT, truth$sigma_logT)
cap_CWT   <- rbinom(N_CWT, 1, truth$p_CWT)

seen_CWT <- alive_CWT == 1 & cap_CWT == 1
T_CWT    <- T_CWT_all[seen_CWT]
n_CWT_obs <- length(T_CWT)

cat("\nSIMULATED DATA\n")
cat("AT detected:", n_AT_obs, "of", N_AT, "\n")
cat("CWT recovered:", n_CWT_obs, "of", N_CWT, "\n")
cat("AT median observed travel time:", round(median(T_AT), 2), "days\n")
cat("CWT median observed travel time:", round(median(T_CWT), 2), "days\n\n")

if (n_AT_obs < 2 || n_CWT_obs < 2) {
  stop("Too few observations. Increase sample size or detection/capture probability.")
}

# ----------------------------------------------------------------------------
# 4. WRITE THE SHORT JAGS MODEL
# ----------------------------------------------------------------------------

model_string <- "
model {

  # Common survival
  S ~ dbeta(1, 1)

  # Common travel-time distribution
  # log(T) ~ Normal(mu_logT, sigma_logT^2)
  mu_logT ~ dnorm(log(10), 0.1)
  sigma_logT ~ dunif(0.05, 1.5)
  tau_logT <- pow(sigma_logT, -2)

  # Acoustic observation model
  n_AT_obs ~ dbin(S * p_AT, N_AT)

  # CWT observation model
  n_CWT_obs ~ dbin(S * p_CWT, N_CWT)

  # Observed travel times from both tag types inform the SAME distribution
  for (i in 1:n_AT_obs) {
    logT_AT[i] ~ dnorm(mu_logT, tau_logT)
  }

  for (i in 1:n_CWT_obs) {
    logT_CWT[i] ~ dnorm(mu_logT, tau_logT)
  }

  # Derived quantities
  median_T <- exp(mu_logT)
  mean_T <- exp(mu_logT + 0.5 * pow(sigma_logT, 2))
}
"

model_file <- file.path(out_dir, "simple_joint_AT_CWT.jags")
writeLines(model_string, model_file)

# ----------------------------------------------------------------------------
# 5. FIT THE JOINT MODEL
# ----------------------------------------------------------------------------

jags_data <- list(
  N_AT = N_AT,
  n_AT_obs = n_AT_obs,
  logT_AT = log(T_AT),
  p_AT = truth$p_AT,

  N_CWT = N_CWT,
  n_CWT_obs = n_CWT_obs,
  logT_CWT = log(T_CWT),
  p_CWT = truth$p_CWT
)

inits <- function() {
  list(
    S = runif(1, 0.35, 0.65),
    mu_logT = rnorm(1, log(12), 0.15),
    sigma_logT = runif(1, 0.25, 0.60)
  )
}

jm <- jags.model(
  file = model_file,
  data = jags_data,
  inits = inits,
  n.chains = 3,
  n.adapt = 1000
)

update(jm, 2000)

samp <- coda.samples(
  model = jm,
  variable.names = c("S", "mu_logT", "sigma_logT", "median_T", "mean_T"),
  n.iter = 6000,
  thin = 2
)

print(summary(samp))
print(gelman.diag(samp, multivariate = FALSE))

draws <- as.matrix(samp)

# ----------------------------------------------------------------------------
# 6. SMALL PARAMETER-RECOVERY TABLE
# ----------------------------------------------------------------------------

q <- function(x) quantile(x, c(0.025, 0.50, 0.975))

recovery <- rbind(
  S = c(
    truth = truth$S,
    mean = mean(draws[, "S"]),
    q(draws[, "S"])
  ),
  median_T = c(
    truth = truth$median_T,
    mean = mean(draws[, "median_T"]),
    q(draws[, "median_T"])
  ),
  sigma_logT = c(
    truth = truth$sigma_logT,
    mean = mean(draws[, "sigma_logT"]),
    q(draws[, "sigma_logT"])
  )
)

colnames(recovery) <- c("truth", "post_mean", "q2.5", "q50", "q97.5")

cat("\nPARAMETER RECOVERY\n")
print(round(recovery, 3))

write.csv(recovery, file.path(out_dir, "simple_joint_parameter_recovery.csv"))

# ----------------------------------------------------------------------------
# 7. FIGURE 1: OBSERVED TRAVEL TIMES + TRUE/POSTERIOR DISTRIBUTIONS
# ----------------------------------------------------------------------------

post_mu  <- mean(draws[, "mu_logT"])
post_sig <- mean(draws[, "sigma_logT"])

xmax <- max(c(T_AT, T_CWT, qlnorm(0.995, post_mu, post_sig)))
xx <- seq(0.05, xmax, length.out = 500)

png(file.path(out_dir, "figure1_travel_time.png"), width = 1600, height = 1000, res = 180)

plot(
  density(T_AT, from = 0),
  lwd = 3,
  xlim = c(0, xmax),
  ylim = c(0, max(density(T_AT)$y, density(T_CWT)$y) * 1.15),
  xlab = "Travel time to endpoint (days)",
  ylab = "Density",
  main = "Observed travel times from two tag types"
)

lines(density(T_CWT, from = 0), lwd = 3, lty = 2)
lines(
  xx,
  dlnorm(xx, log(truth$median_T), truth$sigma_logT),
  lwd = 2,
  lty = 3
)
lines(
  xx,
  dlnorm(xx, post_mu, post_sig),
  lwd = 2
)

legend(
  "topright",
  legend = c(
    "Acoustic detections",
    "CWT recoveries",
    "True travel-time distribution",
    "Joint posterior mean"
  ),
  lty = c(1, 2, 3, 1),
  lwd = c(3, 3, 2, 2),
  bty = "n"
)

dev.off()

# ----------------------------------------------------------------------------
# 8. FIGURE 2: POSTERIOR RECOVERY OF SURVIVAL AND MEDIAN TRAVEL TIME
# ----------------------------------------------------------------------------

png(file.path(out_dir, "figure2_posterior_recovery.png"), width = 1600, height = 800, res = 180)

par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

plot(
  density(draws[, "S"]),
  lwd = 3,
  xlab = "Survival",
  ylab = "Posterior density",
  main = "Survival"
)
abline(v = truth$S, lty = 2, lwd = 2)
legend("topright", legend = "Truth", lty = 2, lwd = 2, bty = "n")

plot(
  density(draws[, "median_T"]),
  lwd = 3,
  xlab = "Median travel time (days)",
  ylab = "Posterior density",
  main = "Travel time"
)
abline(v = truth$median_T, lty = 2, lwd = 2)
legend("topright", legend = "Truth", lty = 2, lwd = 2, bty = "n")

dev.off()

# ----------------------------------------------------------------------------
# 9. QUICK OBSERVATION-SPECIFIC SURVIVAL ESTIMATES
# ----------------------------------------------------------------------------
# These are not Bayesian fits. They simply show what each observation system
# implies when p_AT and p_CWT are known.

S_AT_simple  <- n_AT_obs  / (N_AT  * truth$p_AT)
S_CWT_simple <- n_CWT_obs / (N_CWT * truth$p_CWT)

cat("\nSIMPLE SOURCE-SPECIFIC ESTIMATES\n")
cat("AT only :", round(S_AT_simple, 3), "\n")
cat("CWT only:", round(S_CWT_simple, 3), "\n")
cat("Joint posterior mean:", round(mean(draws[, "S"]), 3), "\n")

# ----------------------------------------------------------------------------
# IMPORTANT IDENTIFIABILITY POINT
# ----------------------------------------------------------------------------
# With only ONE endpoint and UNKNOWN p_AT or p_CWT, the data identify:
#
#        S * p_AT
# and    S * p_CWT
#
# not S and p separately.
#
# That is why this first simulation treats p_AT and p_CWT as known.
# A later model can estimate p_AT using a downstream acoustic receiver and
# p_CWT using trawl-efficiency trials or an informative external estimate.
#
# Likewise, premature acoustic-tag failure is omitted here. A next version can
# replace p_AT with p_AT * q(T), where q(T) is the probability the transmitter
# is still operating at travel time T, following the logic of Holbrook et al.
# ============================================================================

cat("\nOutputs written to:\n", normalizePath(out_dir, winslash = "/", mustWork = FALSE), "\n")
