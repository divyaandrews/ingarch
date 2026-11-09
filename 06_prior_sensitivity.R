# Prior sensitivity analysis (Table S4) ----------------------------------------
#
# Reviewer 2 (majors 4 and 5); Reviewer 1 (minor 2).
#   (a) Prior-only runs: what each prior implies on the original scale.
#   (b) Three prior sets on the polio and hepatitis-B training data.
#   (c) Three prior sets on 50 simulated series from scenarios III and IV
#       (n = 200): average estimate, bias, RMSE and coverage.

source("R/00_functions.R")
library(future.apply)

priors <- list(
  default = default_prior,
  diffuse = modifyList(
    default_prior,
    list(a0_sd = 2, w_conc = c(0.8, 0.8, 0.8), phi_a = 0.5, phi_b = 0.5)
  ),
  concentrated = modifyList(
    default_prior,
    list(a0_sd = 0.5, w_conc = c(2, 2, 2), phi_a = 2, phi_b = 2)
  )
)

# (a) Prior-only draws ---------------------------------------------------------
set.seed(6)
y_dummy <- sim_ingarch11(
  200,
  list(alpha0 = 1, alpha1 = 0.4, beta1 = 0.2, phi = 0.35, nu = NA)
)
prior_tab <- imap(priors, function(prior, name) {
  fit <- fit_bayes(y_dummy, 1L, prior = prior, prior_only = 1L)
  hmc_diagnostics(fit, 1L)$summary |>
    mutate(prior = name)
}) |>
  list_rbind()
write_csv(prior_tab, proj_path("results", "table_prior_only.csv"))

# (b) Real data ----------------------------------------------------------------
polio <- read_csv(proj_path("data", "polio.csv"), show_col_types = FALSE)$x
hepb <- read_csv(proj_path("data", "hepb.csv"), show_col_types = FALSE)$data
train <- list(polio = polio[1:132], hepb = hepb[1:110])

real_sens <- expand_grid(data = names(train), prior = names(priors)) |>
  mutate(summ = map2(data, prior, function(d, p) {
    fit <- fit_bayes(train[[d]], 1L, prior = priors[[p]], cores = 4)
    hmc_diagnostics(fit, 1L)$summary
  })) |>
  unnest(summ)
write_csv(real_sens, proj_path("results", "table_prior_sens_real.csv"))

# (c) Simulated data -----------------------------------------------------------
truths <- list(
  III = list(alpha0 = 1, alpha1 = 0.4, beta1 = 0.2, phi = 0.55, nu = NA),
  IV = list(alpha0 = 1, alpha1 = 0.4, beta1 = 0.2, phi = 0.35, nu = NA)
)

invisible(get_stan_model())
root <- proj_root
plan(multisession, workers = max(1, parallel::detectCores() - 1))
sim_sens <- future_lapply(seq_len(100), function(i) {
  setwd(root)
  source("R/00_functions.R")
  scen <- if (i <= 50) "III" else "IV"
  y <- sim_ingarch11(200, truths[[scen]])

  imap(priors, function(prior, name) {
    fit <- fit_bayes(y, 1L, prior = prior, seed = i)
    hmc_diagnostics(fit, 1L)$summary |>
      mutate(
        prior = name,
        scenario = scen,
        rep = i,
        truth = unlist(truths[[scen]][variable])
      )
  }) |>
    list_rbind()
}, future.seed = TRUE) |>
  list_rbind()
plan(sequential)

table_sim_sens <- sim_sens |>
  group_by(scenario, prior, variable, truth) |>
  summarise(
    avg_est = mean(mean),
    bias = mean(mean - truth),
    rmse = sqrt(mean((mean - truth)^2)),
    coverage = mean(q2.5 <= truth & truth <= q97.5),
    .groups = "drop"
  )
write_csv(table_sim_sens, proj_path("results", "table_prior_sens_sim.csv"))
latex_rows(
  table_sim_sens,
  file = proj_path("results", "table_prior_sens_sim.tex")
)

print(prior_tab)
print(real_sens, n = Inf)
print(table_sim_sens, n = Inf)
