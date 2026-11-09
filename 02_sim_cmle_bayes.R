# Simulation study: CMLE and Bayesian (NUTS) estimation ------------------------
#
# Finite-sample study of NoGe-INGARCH(1,1) estimation with the corrected Stan
# model. Reviewer 2 (majors 4, 6 and 8); Reviewer 1 (majors 6 and 8).
#   * Identical criteria for both methods: bias, mean absolute bias, MSE,
#     RMSE, 95% coverage (Wald for CMLE, equal-tailed credible interval for
#     Bayes) and interval width.
#   * Bayesian point estimate: posterior mean (posterior median also stored).
#   * Per-replicate HMC diagnostics: R-hat, bulk/tail ESS, divergences,
#     tree-depth saturation, E-BFMI and wall time.
#   * Resumable: one .rds file per replicate in results/sim/; finished
#     replicates are skipped when the script is run again.
#
# Run from a terminal:  Rscript R/02_sim_cmle_bayes.R
# Pilot: set n_reps_bayes <- 5 first and look at results/table_time.csv.

source("R/00_functions.R")
library(future.apply)

# Settings ---------------------------------------------------------------------
n_reps_cmle <- 1000
n_reps_bayes <- 500
sizes <- c(50, 200, 500)
n_workers <- max(1, parallel::detectCores() - 1)
out_dir <- proj_path("results", "sim")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scenarios <- tribble(
  ~scenario, ~alpha0, ~alpha1, ~beta1, ~phi,
  "I",       1.0,     0.2,     0.1,    0.05,
  "II",      1.0,     0.3,     0.1,    0.05,
  "III",     1.0,     0.4,     0.2,    0.55,
  "IV",      1.0,     0.4,     0.2,    0.35,
  "V",       0.5,     0.2,     0.7,    0.20
)

# Compile once in the main session; workers then load the cached model.
invisible(get_stan_model())

# One replicate ----------------------------------------------------------------
run_replicate <- function(s, n, rep_id, do_bayes, root) {
  setwd(root)
  source("R/00_functions.R")

  scen <- scenarios[s, ]
  file <- file.path(
    "results", "sim",
    sprintf("%s_n%d_r%04d.rds", scen$scenario, n, rep_id)
  )
  if (file.exists(file)) {
    return(invisible(NULL))
  }

  seed <- 100000 * s + 1000 * n + rep_id
  set.seed(seed)
  truth <- list(alpha0 = scen$alpha0, alpha1 = scen$alpha1,
                beta1 = scen$beta1, phi = scen$phi, nu = NA)
  y <- sim_ingarch11(n, truth, family = 1L)

  start <- Sys.time()
  ml <- fit_cmle(y, family = 1L)
  ml_time <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  res <- tibble(
    param = names(ml$est),
    est = ml$est,
    post_median = NA_real_,
    se = ml$se,
    lower = ml$est - 1.96 * ml$se,
    upper = ml$est + 1.96 * ml$se,
    method = "CMLE",
    time_sec = ml_time,
    rhat = NA_real_,
    ess_bulk = NA_real_,
    ess_tail = NA_real_,
    n_div = NA_real_,
    n_td = NA_real_,
    ebfmi = NA_real_
  )

  if (do_bayes) {
    fit <- fit_bayes(y, family = 1L, seed = seed)
    hmc <- hmc_diagnostics(fit, family = 1L)
    bayes <- hmc$summary |>
      transmute(
        param = variable,
        est = mean,
        post_median = median,
        se = sd,
        lower = q2.5,
        upper = q97.5,
        method = "Bayes",
        time_sec = hmc$time_sec,
        rhat,
        ess_bulk,
        ess_tail,
        n_div = hmc$n_divergent,
        n_td = hmc$n_max_treedepth,
        ebfmi = hmc$min_ebfmi
      )
    res <- bind_rows(res, bayes)
  }

  res <- res |>
    mutate(
      scenario = scen$scenario,
      n = n,
      rep = rep_id,
      truth = unlist(truth[param])
    )
  saveRDS(res, file)
  invisible(NULL)
}

# Run --------------------------------------------------------------------------
grid <- expand_grid(
  s = seq_len(nrow(scenarios)),
  n = sizes,
  rep = seq_len(n_reps_cmle)
) |>
  mutate(do_bayes = rep <= n_reps_bayes) |>
  arrange(desc(do_bayes))

root <- proj_root
plan(multisession, workers = n_workers)
t_start <- Sys.time()
invisible(future_lapply(
  seq_len(nrow(grid)),
  function(i) {
    run_replicate(grid$s[i], grid$n[i], grid$rep[i], grid$do_bayes[i], root)
  },
  future.seed = TRUE,
  future.scheduling = 5
))
plan(sequential)
message("Total time: ", format(Sys.time() - t_start))

# Summaries --------------------------------------------------------------------
all_res <- list.files(out_dir, full.names = TRUE) |>
  map(readRDS) |>
  list_rbind()

table_sim <- all_res |>
  group_by(scenario, n, method, param, truth) |>
  summarise(
    avg_est = mean(est),
    bias = mean(est - truth),
    abs_bias = mean(abs(est - truth)),
    mse = mean((est - truth)^2),
    rmse = sqrt(mse),
    coverage = mean(lower <= truth & truth <= upper, na.rm = TRUE),
    ci_width = mean(upper - lower, na.rm = TRUE),
    pct_se_na = 100 * mean(is.na(se)),
    n_reps = n(),
    .groups = "drop"
  ) |>
  mutate(param = factor(param, c("alpha0", "alpha1", "beta1", "phi"))) |>
  arrange(scenario, n, param, method)

table_conv <- all_res |>
  filter(method == "Bayes") |>
  group_by(scenario, n, rep) |>
  summarise(
    max_rhat = max(rhat),
    min_bulk = min(ess_bulk),
    min_tail = min(ess_tail),
    n_div = first(n_div),
    n_td = first(n_td),
    ebfmi = first(ebfmi),
    time_sec = first(time_sec),
    .groups = "drop"
  ) |>
  group_by(scenario, n) |>
  summarise(
    pct_rhat_gt_101 = 100 * mean(max_rhat > 1.01),
    med_min_ess_bulk = median(min_bulk),
    med_min_ess_tail = median(min_tail),
    pct_reps_div = 100 * mean(n_div > 0),
    pct_reps_td = 100 * mean(n_td > 0),
    pct_low_ebfmi = 100 * mean(ebfmi < 0.3),
    med_time_sec = median(time_sec),
    med_ess_per_sec = median(min_bulk / time_sec),
    n_reps = n(),
    .groups = "drop"
  )

table_time <- all_res |>
  distinct(scenario, n, rep, method, time_sec) |>
  group_by(scenario, n, method) |>
  summarise(median_time_sec = median(time_sec), .groups = "drop")

write_csv(table_sim, proj_path("results", "table_simulation.csv"))
write_csv(table_conv, proj_path("results", "table_convergence.csv"))
write_csv(table_time, proj_path("results", "table_time.csv"))

# LaTeX rows in the layout of Table 2 (one file per scenario):
# n & parameter & CMLE bias, RMSE, CP & Bayes bias, RMSE, CP
table_tex <- table_sim |>
  select(scenario, n, param, method, bias, rmse, coverage) |>
  pivot_wider(
    names_from = method,
    values_from = c(bias, rmse, coverage)
  ) |>
  transmute(
    scenario,
    n,
    param,
    bias_CMLE, rmse_CMLE, coverage_CMLE,
    bias_Bayes, rmse_Bayes, coverage_Bayes
  )

for (s in unique(table_tex$scenario)) {
  table_tex |>
    filter(scenario == s) |>
    select(-scenario) |>
    latex_rows(file = proj_path("results", sprintf("table_sim_%s.tex", s)))
}
latex_rows(
  table_conv,
  digits = 2,
  file = proj_path("results", "table_convergence.tex")
)

print(table_sim, n = Inf)
print(table_conv, n = Inf)
