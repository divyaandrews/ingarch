# Pilot check (about 5 minutes) ------------------------------------------------
#
# Run this first: source("R/00_pilot.R")
# It compiles the Stan model, fits one simulated series by CMLE and by NUTS,
# and checks that both are close to the true values. Do not start the long
# scripts until this passes.

source("R/00_functions.R")
set.seed(2026)

truth <- list(alpha0 = 1, alpha1 = 0.4, beta1 = 0.2, phi = 0.35, nu = NA)
y <- sim_ingarch11(500, truth, family = 1L)

cmle <- fit_cmle(y, family = 1L)
message("CMLE estimates (truth 1, 0.4, 0.2, 0.35):")
print(round(cmle$est, 3))
print(round(cmle$se, 3))

fit <- fit_bayes(y, family = 1L, cores = 4)
hmc <- hmc_diagnostics(fit, family = 1L)
message("Posterior summary:")
print(hmc$summary)
message(
  "Divergences: ", hmc$n_divergent,
  " | max R-hat: ", round(hmc$max_rhat, 3),
  " | min bulk ESS: ", round(hmc$min_ess_bulk),
  " | time (s): ", round(hmc$time_sec, 1)
)

post_mean <- hmc$summary$mean
gap <- abs(post_mean - cmle$est)
if (any(gap > 0.25) || hmc$max_rhat > 1.01) {
  warning("Pilot check failed: inspect the output above before continuing.")
} else {
  message("Pilot check passed.")
}
