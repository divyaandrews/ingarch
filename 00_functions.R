# Shared functions for the NoGe-INGARCH revision (COST-D-26-00400) -------------
#
# 1. Probability functions of the six conditional laws
# 2. Simulation of INGARCH(1,1) paths
# 3. Conditional log-likelihood and conditional maximum likelihood (CMLE)
# 4. Closed-form NoGe-INGARCH(1,1) moments (Corollary 2.1)
# 5. Predictive distributions, scores and PIT
# 6. Stan helpers and diagnostics
# 7. Output helpers
#
# Usage: set the working directory to the project folder (the folder that
# contains R/, stan/ and data/) and run source("R/00_functions.R").
# Code style: tidyverse style guide (https://style.tidyverse.org).

stopifnot(getRversion() >= "4.1.0")

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstan)
  library(posterior)
  library(loo)
})

rstan_options(auto_write = TRUE)

proj_root <- getOption("noge.root", default = getwd())

proj_path <- function(...) {
  file.path(proj_root, ...)
}

dir.create(proj_path("results"), showWarnings = FALSE)
dir.create(proj_path("figures"), showWarnings = FALSE)

families <- c(NoGe = 1L, Poisson = 2L, NB = 3L, GP = 4L, ZIP = 5L, ZINB = 6L)

has_phi <- function(family) {
  family %in% c(1L, 4L, 5L, 6L)
}

has_nu <- function(family) {
  family %in% c(3L, 6L)
}

n_par <- function(family) {
  3L + has_phi(family) + has_nu(family)
}


# 1. Conditional laws ----------------------------------------------------------

# NoGe(lambda, phi): P(0) = phi and P(x) = (1 - phi) theta (1 - theta)^(x - 1)
# for x >= 1, with theta = (1 - phi) / lambda (mean parameterisation).
dnoge <- function(x, lambda, phi, log = FALSE) {
  theta <- (1 - phi) / lambda
  theta_c <- pmin(theta, 1)
  tail_term <- ifelse(
    x > 1,
    (x - 1) * log1p(-pmin(theta_c, 1 - 1e-15)),
    0
  )
  log_pos <- log1p(-phi) + log(theta_c) + tail_term
  # theta > 1 is outside the admissible space: no mass above one
  log_pos <- ifelse(theta >= 1 & x > 1, -Inf, log_pos)
  out <- ifelse(x == 0, log(phi), log_pos)

  if (log) out else exp(out)
}

# CDF: P(X <= q) = 1 - (1 - phi) (1 - theta)^q for q >= 0
pnoge <- function(q, lambda, phi) {
  theta <- pmin((1 - phi) / lambda, 1)
  ifelse(q < 0, 0, 1 - (1 - phi) * (1 - theta)^floor(q))
}

rnoge <- function(n, lambda, phi) {
  theta <- pmin(rep_len((1 - phi) / lambda, n), 1)
  ifelse(runif(n) < phi, 0, rgeom(n, theta) + 1)
}

dgpois <- function(x, lambda, kappa, log = FALSE) {
  eta <- (1 - kappa) * lambda
  out <- log(eta) + (x - 1) * log(eta + kappa * x) - eta - kappa * x -
    lgamma(x + 1)

  if (log) out else exp(out)
}

# Inversion sampler for one generalized Poisson draw
rgpois_one <- function(lambda, kappa) {
  u <- runif(1)
  cdf <- 0
  x <- 0
  repeat {
    cdf <- cdf + dgpois(x, lambda, kappa)
    if (u <= cdf || x > 1e5) {
      return(x)
    }
    x <- x + 1
  }
}

# Log-pmf of any family; vectorised over x, lambda, phi and nu
dfam <- function(x, lambda, family, phi = NA, nu = NA, log = FALSE) {
  out <- switch(
    as.character(family),
    "1" = dnoge(x, lambda, phi, log = TRUE),
    "2" = dpois(x, lambda, log = TRUE),
    "3" = dnbinom(x, mu = lambda, size = nu, log = TRUE),
    "4" = dgpois(x, lambda, phi, log = TRUE),
    "5" = ifelse(
      x == 0,
      log(phi + (1 - phi) * dpois(0, lambda)),
      log1p(-phi) + dpois(x, lambda, log = TRUE)
    ),
    "6" = ifelse(
      x == 0,
      log(phi + (1 - phi) * dnbinom(0, mu = lambda, size = nu)),
      log1p(-phi) + dnbinom(x, mu = lambda, size = nu, log = TRUE)
    )
  )

  if (log) out else exp(out)
}

# One draw from any family
rfam_one <- function(lambda, family, phi = NA, nu = NA) {
  switch(
    as.character(family),
    "1" = rnoge(1, lambda, phi),
    "2" = rpois(1, lambda),
    "3" = rnbinom(1, mu = lambda, size = nu),
    "4" = rgpois_one(lambda, phi),
    "5" = if (runif(1) < phi) 0 else rpois(1, lambda),
    "6" = if (runif(1) < phi) 0 else rnbinom(1, mu = lambda, size = nu)
  )
}

# Stationary mean of lambda_t; for ZIP/ZINB, E[X_t | F_{t-1}] = (1 - omega)
# lambda_t, which changes the stationary mean of the recursion.
lambda_start <- function(pars, family) {
  f <- if (family %in% c(5L, 6L)) 1 - pars$phi else 1
  pars$alpha0 / (1 - pars$alpha1 * f - pars$beta1)
}


# 2. Simulation ----------------------------------------------------------------

sim_ingarch11 <- function(n, pars, family = 1L, burn = 300) {
  lambda <- lambda_start(pars, family)
  x <- numeric(n + burn)

  for (t in seq_len(n + burn)) {
    x[t] <- rfam_one(lambda, family, pars$phi, pars$nu)
    lambda <- pars$alpha0 + pars$alpha1 * x[t] + pars$beta1 * lambda
  }

  as.integer(x[-seq_len(burn)])
}


# 3. Conditional log-likelihood and CMLE ---------------------------------------

# The unconstrained vector z is mapped to the admissible parameter space,
# exactly as in the Stan program:
#   (alpha1, beta1, 1 - alpha1 - beta1) = softmax(z1, z2, 0)
#   phi = plogis(z3), nu = exp(z4)
#   alpha0 = lb + exp(z_last), lb = (1 - phi)(1 - beta1) for NoGe, else 0.
z_to_pars <- function(z, family) {
  e <- exp(c(z[1], z[2], 0) - max(z[1], z[2], 0))
  w <- e / sum(e)
  i <- 3
  phi <- NA
  nu <- NA

  if (has_phi(family)) {
    phi <- plogis(z[i])
    i <- i + 1
  }
  if (has_nu(family)) {
    nu <- exp(z[i])
    i <- i + 1
  }
  lb <- if (family == 1L) (1 - phi) * (1 - w[2]) else 0

  list(
    alpha0 = lb + exp(z[i]),
    alpha1 = w[1],
    beta1 = w[2],
    phi = phi,
    nu = nu
  )
}

pars_to_z <- function(pars, family) {
  w3 <- 1 - pars$alpha1 - pars$beta1
  z <- c(log(pars$alpha1 / w3), log(pars$beta1 / w3))

  if (has_phi(family)) {
    z <- c(z, qlogis(pars$phi))
  }
  if (has_nu(family)) {
    z <- c(z, log(pars$nu))
  }
  lb <- if (family == 1L) (1 - pars$phi) * (1 - pars$beta1) else 0

  c(z, log(pars$alpha0 - lb))
}

# lambda_1 is the stationary mean; stats::filter() runs the recursion
# lambda_t = alpha0 + alpha1 * y_{t-1} + beta1 * lambda_{t-1} in C.
lambda_path <- function(y, pars, family) {
  lambda_1 <- lambda_start(pars, family)
  u <- pars$alpha0 + pars$alpha1 * y[-length(y)]
  u[1] <- u[1] + pars$beta1 * lambda_1
  rest <- stats::filter(u, pars$beta1, method = "recursive")

  c(lambda_1, as.numeric(rest))
}

loglik_ingarch11 <- function(pars, y, family) {
  lambda <- lambda_path(y, pars, family)
  sum(dfam(y[-1], lambda[-1], family, pars$phi, pars$nu, log = TRUE))
}

pars_vec <- function(pars, family) {
  v <- c(alpha0 = pars$alpha0, alpha1 = pars$alpha1, beta1 = pars$beta1)
  if (has_phi(family)) {
    v <- c(v, phi = pars$phi)
  }
  if (has_nu(family)) {
    v <- c(v, nu = pars$nu)
  }
  v
}

# Delta-method standard errors from the observed information in z
cmle_se <- function(z_hat, nll, family) {
  k <- length(z_hat)
  se <- rep(NA_real_, k)

  hess <- tryCatch(optimHess(z_hat, nll), error = function(e) NULL)
  if (is.null(hess)) {
    return(se)
  }
  vcov_z <- tryCatch(solve(hess), error = function(e) NULL)
  if (is.null(vcov_z)) {
    return(se)
  }

  jac <- matrix(NA_real_, k, k)
  h <- 1e-6
  for (j in seq_len(k)) {
    step <- replace(numeric(k), j, h)
    upper <- pars_vec(z_to_pars(z_hat + step, family), family)
    lower <- pars_vec(z_to_pars(z_hat - step, family), family)
    jac[, j] <- (upper - lower) / (2 * h)
  }
  vcov_theta <- jac %*% vcov_z %*% t(jac)

  if (all(is.finite(vcov_theta)) && all(diag(vcov_theta) > 0)) {
    se <- sqrt(diag(vcov_theta))
  }
  se
}

# CMLE by BFGS on the unconstrained scale from several starting values
fit_cmle <- function(y,
                     family = 1L,
                     se = TRUE,
                     starts = list(c(0.2, 0.3), c(0.1, 0.6), c(0.4, 0.2))) {
  y <- as.numeric(y)

  nll <- function(z) {
    value <- suppressWarnings(
      -loglik_ingarch11(z_to_pars(z, family), y, family)
    )
    if (is.finite(value)) value else 1e10
  }

  phi_0 <- switch(
    as.character(family),
    "1" = min(max(mean(y == 0), 0.02), 0.95),
    "4" = 0.2,
    0.1
  )

  fits <- lapply(starts, function(ab) {
    lb <- if (family == 1L) (1 - phi_0) * (1 - ab[2]) else 0
    alpha0_0 <- max(mean(y) * (1 - sum(ab)), lb + 0.05)
    z_0 <- pars_to_z(
      list(alpha0 = alpha0_0, alpha1 = ab[1], beta1 = ab[2],
           phi = phi_0, nu = 5),
      family
    )
    tryCatch(
      optim(z_0, nll, method = "BFGS",
            control = list(maxit = 1000, reltol = 1e-12)),
      error = function(e) NULL
    )
  })
  fits <- Filter(Negate(is.null), fits)
  best <- fits[[which.min(vapply(fits, function(f) f$value, numeric(1)))]]

  pars <- z_to_pars(best$par, family)
  est <- pars_vec(pars, family)
  se_vec <- if (se) {
    cmle_se(best$par, nll, family)
  } else {
    rep(NA_real_, length(est))
  }
  names(se_vec) <- names(est)

  k <- n_par(family)
  loglik <- -best$value

  list(
    family = family,
    pars = pars,
    est = est,
    se = se_vec,
    loglik = loglik,
    aic = -2 * loglik + 2 * k,
    bic = -2 * loglik + log(length(y) - 1) * k,
    convergence = best$convergence
  )
}


# 4. NoGe-INGARCH(1,1) moments (Corollary 2.1) ---------------------------------

# zeta = (1 + phi) / (1 - phi), Var(X | F) = lambda (zeta lambda - 1)
# D = 1 - 2 alpha1^2 / (1 - phi) - 2 alpha1 beta1 - beta1^2 (D > 0 iff
# second-order stationary)
noge_moments <- function(alpha0, alpha1, beta1, phi) {
  zeta <- (1 + phi) / (1 - phi)
  mu <- alpha0 / (1 - alpha1 - beta1)
  c0 <- zeta * mu^2 - mu
  d <- 1 - 2 * alpha1^2 / (1 - phi) - 2 * alpha1 * beta1 - beta1^2
  var_x <- ifelse(d > 0, c0 * (1 - 2 * alpha1 * beta1 - beta1^2) / d, Inf)
  rho_1 <- alpha1 * (1 - alpha1 * beta1 - beta1^2) /
    (1 - 2 * alpha1 * beta1 - beta1^2)

  tibble(
    mu = mu,
    var_x = var_x,
    disp_index = var_x / mu,
    cond_disp_at_mu = zeta * mu - 1,
    rho_1 = rho_1,
    d = d,
    admissible = alpha0 >= (1 - phi) * (1 - beta1),
    p_zero = phi
  )
}


# 5. Predictive distributions, scores and PIT ----------------------------------

# Posterior-mixture predictive pmf on 0..k_max at one time point;
# lambda, phi and nu are vectors of posterior draws.
pred_pmf <- function(lambda, family, phi, nu, k_max) {
  grid <- 0:k_max
  m <- length(lambda)
  dens <- dfam(
    rep(grid, each = m),
    rep(lambda, k_max + 1),
    family,
    rep(phi, k_max + 1),
    rep(nu, k_max + 1)
  )
  p <- colMeans(matrix(dens, nrow = m))
  p / sum(p)
}

pmf_summary <- function(p, x, level = 0.95) {
  cdf <- cumsum(p)
  grid <- seq_along(p) - 1
  a <- (1 - level) / 2

  tibble(
    mean = sum(grid * p),
    median = grid[which(cdf >= 0.5)[1]],
    lower = grid[which(cdf >= a)[1]],
    upper = grid[which(cdf >= 1 - a)[1]],
    log_score = -log(p[x + 1]),
    crps = sum((cdf - as.numeric(grid >= x))^2),
    p_lo = if (x == 0) 0 else cdf[x],
    p_hi = cdf[x + 1],
    p_zero = p[1]
  )
}

# Non-randomised PIT histogram for counts (Czado, Gneiting and Held, 2009)
pit_nonrandom <- function(p_lo, p_hi, n_bins = 10) {
  u <- seq(0, 1, length.out = n_bins + 1)
  f_bar <- vapply(u, function(uu) {
    mean(ifelse(
      uu <= p_lo,
      0,
      ifelse(uu >= p_hi, 1, (uu - p_lo) / (p_hi - p_lo))
    ))
  }, numeric(1))

  tibble(
    bin_lower = u[-length(u)],
    bin_upper = u[-1],
    density = diff(f_bar) * n_bins
  )
}

pit_random <- function(p_lo, p_hi) {
  p_lo + runif(length(p_lo)) * (p_hi - p_lo)
}


# 6. Stan helpers --------------------------------------------------------------

get_stan_model <- function() {
  if (!exists(".stan_model_cache", envir = globalenv())) {
    model <- stan_model(proj_path("stan", "ingarch11.stan"))
    assign(".stan_model_cache", model, envir = globalenv())
  }
  get(".stan_model_cache", envir = globalenv())
}

# Default weakly informative priors, all on the original scale:
#   alpha0 ~ LN(0, 1) truncated at lb; (alpha1, beta1, 1 - alpha1 - beta1)
#   ~ Dirichlet(1, 1, 1); phi ~ Beta(1, 1); nu ~ Gamma(2, 0.1)
default_prior <- list(
  a0_mu = 0,
  a0_sd = 1,
  w_conc = c(1, 1, 1),
  phi_a = 1,
  phi_b = 1,
  nu_shape = 2,
  nu_rate = 0.1
)

make_stan_data <- function(y,
                           family,
                           y_test = integer(0),
                           h_max = 0,
                           gq_full = 0L,
                           prior = default_prior,
                           prior_only = 0L) {
  data_list <- list(
    N = length(y),
    y = as.integer(y),
    N_test = length(y_test),
    y_test = as.integer(y_test),
    family = as.integer(family),
    H = as.integer(h_max),
    gq_full = as.integer(gq_full),
    prior_only = as.integer(prior_only)
  )
  c(data_list, prior)
}

fit_bayes <- function(y,
                      family = 1L,
                      y_test = integer(0),
                      h_max = 0,
                      gq_full = 0L,
                      prior = default_prior,
                      chains = 4,
                      iter = 2000,
                      warmup = 1000,
                      cores = 1,
                      seed = 2026,
                      adapt_delta = 0.95,
                      prior_only = 0L) {
  sampling(
    get_stan_model(),
    data = make_stan_data(
      y, family, y_test, h_max, gq_full, prior, prior_only
    ),
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    seed = seed,
    refresh = 0,
    control = list(adapt_delta = adapt_delta, max_treedepth = 12)
  )
}

model_pars <- function(family) {
  c(
    "alpha0", "alpha1", "beta1",
    if (has_phi(family)) "phi",
    if (has_nu(family)) "nu"
  )
}

# Quantitative HMC diagnostics: rank-normalised R-hat, bulk and tail ESS,
# divergences, tree-depth saturation, E-BFMI and wall time.
hmc_diagnostics <- function(fit, family, max_treedepth = 12) {
  draws <- as_draws_array(as.array(fit, pars = model_pars(family)))
  summ <- summarise_draws(
    draws,
    "mean", "median", "sd",
    ~ quantile2(.x, probs = c(0.025, 0.975)),
    "rhat", "ess_bulk", "ess_tail"
  )

  sampler <- get_sampler_params(fit, inc_warmup = FALSE)
  ebfmi <- vapply(sampler, function(s) {
    energy <- s[, "energy__"]
    sum(diff(energy)^2) / sum((energy - mean(energy))^2)
  }, numeric(1))
  n_div <- sum(vapply(sampler, function(s) sum(s[, "divergent__"]),
                      numeric(1)))
  n_td <- sum(vapply(sampler, function(s) {
    sum(s[, "treedepth__"] >= max_treedepth)
  }, numeric(1)))

  list(
    summary = summ,
    n_divergent = n_div,
    n_max_treedepth = n_td,
    min_ebfmi = min(ebfmi),
    max_rhat = max(summ$rhat),
    min_ess_bulk = min(summ$ess_bulk),
    min_ess_tail = min(summ$ess_tail),
    time_sec = sum(get_elapsed_time(fit))
  )
}


# 7. Output helpers ------------------------------------------------------------

# Write a data frame as LaTeX tabular rows ("a & b & c \\") for pasting into
# the manuscript tables. Whole numbers are printed without decimals.
latex_rows <- function(df, digits = 3, file = NULL) {
  fmt <- function(v) {
    if (!is.numeric(v)) {
      return(as.character(v))
    }
    if (all(v == round(v), na.rm = TRUE)) {
      return(as.character(v))
    }
    out <- formatC(v, format = "f", digits = digits)
    ifelse(is.na(v), "--", out)
  }

  cells <- as.data.frame(lapply(df, fmt))
  rows <- paste0(apply(cells, 1, paste, collapse = " & "), " \\\\")

  if (!is.null(file)) {
    writeLines(rows, file)
  }
  invisible(rows)
}
