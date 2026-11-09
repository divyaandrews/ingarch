# NUTS versus tuned random-walk Metropolis (Table S3) --------------------------
#
# Reviewer 2 (major 7).
#   * The random-walk Metropolis (RWM) sampler evaluates rstan::log_prob() on
#     Stan's unconstrained scale, so its target (likelihood, priors and
#     Jacobians) is identical to the one NUTS samples.
#   * RWM tuning: 5000 adaptive burn-in iterations (proposal covariance learnt
#     on the fly, scale 2.38^2 / d; Roberts and Rosenthal, 2001), then a fixed
#     Gaussian proposal.
#   * Efficiency: minimum bulk and tail ESS over (alpha0, alpha1, beta1, phi),
#     ESS per second of sampling and ESS per 1000 draws.

source("R/00_functions.R")
set.seed(11)

# Settings ---------------------------------------------------------------------
settings <- expand_grid(scenario = c("III", "V"), n = c(200, 500))
truths <- list(
  III = list(alpha0 = 1.0, alpha1 = 0.4, beta1 = 0.2, phi = 0.55, nu = NA),
  V = list(alpha0 = 0.5, alpha1 = 0.2, beta1 = 0.7, phi = 0.20, nu = NA)
)
n_rep <- 20
n_burn <- 5000
n_keep <- 20000
n_chains <- 4
par_names <- c("alpha0", "alpha1", "beta1", "phi")

# Helpers ----------------------------------------------------------------------
# Map an unconstrained point to (alpha0, alpha1, beta1, phi) via Stan.
constrained_draw <- function(fit, u) {
  cp <- constrain_pars(fit, u)
  phi <- cp$phi_arr[1]
  beta1 <- cp$w[2]
  alpha0 <- (1 - phi) * (1 - beta1) + cp$alpha0_ex
  c(alpha0, cp$w[1], beta1, phi)
}

rwm_chain <- function(fit, n_burn, n_keep) {
  d <- get_num_upars(fit)
  u <- unconstrain_pars(fit, get_inits(fit)[[1]])
  lp <- log_prob(fit, u)
  sigma <- diag(0.01, d)
  scale <- 2.38^2 / d
  burn_draws <- matrix(NA_real_, n_burn, d)

  # Adaptive tuning phase
  for (i in seq_len(n_burn)) {
    if (i > 500 && i %% 100 == 0) {
      recent <- burn_draws[(i %/% 2):(i - 1), , drop = FALSE]
      sigma <- cov(recent) + diag(1e-8, d)
    }
    proposal <- u + drop(rnorm(d) %*% chol(scale * sigma))
    lp_prop <- tryCatch(log_prob(fit, proposal), error = function(e) -Inf)
    if (is.finite(lp_prop) && log(runif(1)) < lp_prop - lp) {
      u <- proposal
      lp <- lp_prop
    }
    burn_draws[i, ] <- u
  }

  # Sampling phase with a fixed proposal
  chol_prop <- chol(scale * sigma)
  out <- matrix(NA_real_, n_keep, length(par_names))
  n_accept <- 0
  start <- Sys.time()
  for (i in seq_len(n_keep)) {
    proposal <- u + drop(rnorm(d) %*% chol_prop)
    lp_prop <- tryCatch(log_prob(fit, proposal), error = function(e) -Inf)
    if (is.finite(lp_prop) && log(runif(1)) < lp_prop - lp) {
      u <- proposal
      lp <- lp_prop
      n_accept <- n_accept + 1
    }
    out[i, ] <- constrained_draw(fit, u)
  }

  list(
    draws = out,
    time = as.numeric(difftime(Sys.time(), start, units = "secs")),
    acc_rate = n_accept / n_keep
  )
}

efficiency <- function(draws) {
  summ <- summarise_draws(as_draws_array(draws), "rhat", "ess_bulk",
                          "ess_tail")
  c(
    min_bulk = min(summ$ess_bulk),
    min_tail = min(summ$ess_tail),
    max_rhat = max(summ$rhat)
  )
}

compare_once <- function(scenario, n, rep_id) {
  y <- sim_ingarch11(n, truths[[scenario]])

  fit <- fit_bayes(y, 1L, seed = rep_id, chains = n_chains)
  elapsed <- get_elapsed_time(fit)
  nuts <- as.array(fit, pars = par_names)
  eff_nuts <- efficiency(nuts)

  chains <- lapply(seq_len(n_chains), function(k) {
    rwm_chain(fit, n_burn, n_keep)
  })
  rwm <- simplify2array(lapply(chains, function(ch) ch$draws)) |>
    aperm(c(1, 3, 2))
  dimnames(rwm) <- list(NULL, NULL, par_names)
  eff_rwm <- efficiency(rwm)

  tibble(
    scenario = scenario,
    n = n,
    rep = rep_id,
    sampler = c("NUTS", "RWM"),
    time_sampling = c(
      sum(elapsed[, "sample"]),
      sum(vapply(chains, function(ch) ch$time, numeric(1)))
    ),
    draws = c(prod(dim(nuts)[1:2]), n_chains * n_keep),
    min_bulk = c(eff_nuts[["min_bulk"]], eff_rwm[["min_bulk"]]),
    min_tail = c(eff_nuts[["min_tail"]], eff_rwm[["min_tail"]]),
    max_rhat = c(eff_nuts[["max_rhat"]], eff_rwm[["max_rhat"]]),
    acc_rate = c(NA, mean(vapply(chains, function(ch) ch$acc_rate,
                                 numeric(1))))
  ) |>
    mutate(
      bulk_per_sec = min_bulk / time_sampling,
      bulk_per_1000 = 1000 * min_bulk / draws
    )
}

# Run --------------------------------------------------------------------------
res <- pmap(settings, function(scenario, n) {
  map(seq_len(n_rep), function(r) compare_once(scenario, n, r)) |>
    list_rbind()
}) |>
  list_rbind()

table_hmc_mh <- res |>
  group_by(scenario, n, sampler) |>
  summarise(
    across(
      c(time_sampling, min_bulk, min_tail, bulk_per_sec, bulk_per_1000,
        max_rhat, acc_rate),
      median
    ),
    .groups = "drop"
  )

write_csv(table_hmc_mh, proj_path("results", "table_hmc_vs_mh.csv"))
latex_rows(
  table_hmc_mh,
  digits = 2,
  file = proj_path("results", "table_hmc_vs_mh.tex")
)
print(table_hmc_mh)

# Report the result as it comes out. With four parameters a tuned RWM can be
# competitive in ESS per second; the case for NUTS rests on the higher ESS per
# draw, no manual tuning, and diagnostics that flag boundary problems.
