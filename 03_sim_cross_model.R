# Competing-model simulation (Table 3) -----------------------------------------
#
# Reviewer 1 (majors 1, 3 and 6).
#   Data-generating processes: NoGe (zero-inflated), NoGe (zero-deflated
#   relative to a geometric law with the same mean), ZINB, GP and Poisson.
#   Fitted models: all six families by CMLE on the first 80% of each series.
#   Output: AIC/BIC selection frequencies and the mean one-step hold-out log
#   score and CRPS on the last 20%.
# CMLE is used to keep the run time realistic; the Bayesian comparison of the
# six families is carried out on the real data (05_real_data.R).

source("R/00_functions.R")
library(future.apply)

# Settings ---------------------------------------------------------------------
n_reps <- 300
sizes <- c(200, 500)
n_workers <- max(1, parallel::detectCores() - 1)
dir.create(proj_path("results", "cross"), recursive = TRUE,
           showWarnings = FALSE)

base_pars <- list(alpha0 = 1, alpha1 = 0.4, beta1 = 0.2, phi = NA, nu = NA)
dgps <- list(
  "NoGe (zero-inflated)" = list(
    family = 1L,
    pars = modifyList(base_pars, list(phi = 0.35))
  ),
  "NoGe (zero-deflated)" = list(
    family = 1L,
    pars = modifyList(base_pars, list(alpha1 = 0.3, beta1 = 0.3, phi = 0.10))
  ),
  "ZINB" = list(
    family = 6L,
    pars = modifyList(base_pars, list(phi = 0.30, nu = 2))
  ),
  "GP" = list(
    family = 4L,
    pars = modifyList(base_pars, list(phi = 0.30))
  ),
  "Poisson" = list(
    family = 2L,
    pars = base_pars
  )
)

# Helpers ----------------------------------------------------------------------
score_holdout <- function(y, n_train, fit) {
  pars <- fit$pars
  test_idx <- (n_train + 1):length(y)
  lambda <- lambda_path(y, pars, fit$family)[test_idx]
  k_max <- max(y) * 3 + 60

  map2(lambda, y[test_idx], function(lam, x) {
    p <- dfam(0:k_max, lam, fit$family, pars$phi, pars$nu)
    pmf_summary(p / sum(p), x)
  }) |>
    list_rbind() |>
    summarise(log_score = mean(log_score), crps = mean(crps))
}

run_replicate <- function(d, n, rep_id, root) {
  setwd(root)
  source("R/00_functions.R")

  file <- file.path(
    "results", "cross",
    sprintf("%d_n%d_r%04d.rds", d, n, rep_id)
  )
  if (file.exists(file)) {
    return(invisible(NULL))
  }

  set.seed(rep_id + 1e4 * n + 1e7 * d)
  dgp <- dgps[[d]]
  y <- sim_ingarch11(n, dgp$pars, dgp$family)
  n_train <- floor(0.8 * n)

  res <- imap(families, function(fam, fam_name) {
    fit <- tryCatch(
      fit_cmle(y[seq_len(n_train)], fam, se = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(NULL)
    }
    bind_cols(
      tibble(
        dgp = names(dgps)[d],
        n = n,
        rep = rep_id,
        fitted = fam_name,
        aic = fit$aic,
        bic = fit$bic
      ),
      score_holdout(y, n_train, fit)
    )
  }) |>
    list_rbind()

  saveRDS(res, file)
  invisible(NULL)
}

# Run --------------------------------------------------------------------------
grid <- expand_grid(d = seq_along(dgps), n = sizes, rep = seq_len(n_reps))

root <- proj_root
plan(multisession, workers = n_workers)
invisible(future_lapply(
  seq_len(nrow(grid)),
  function(i) run_replicate(grid$d[i], grid$n[i], grid$rep[i], root),
  future.seed = TRUE
))
plan(sequential)

# Summaries --------------------------------------------------------------------
all_res <- list.files(proj_path("results", "cross"), full.names = TRUE) |>
  map(readRDS) |>
  list_rbind()

table_selection <- all_res |>
  group_by(dgp, n, rep) |>
  summarise(
    AIC = fitted[which.min(aic)],
    BIC = fitted[which.min(bic)],
    .groups = "drop"
  ) |>
  pivot_longer(c(AIC, BIC), names_to = "criterion", values_to = "picked") |>
  count(dgp, n, criterion, picked, name = "n_picked") |>
  group_by(dgp, n, criterion) |>
  mutate(pct = round(100 * n_picked / sum(n_picked), 1)) |>
  ungroup() |>
  select(-n_picked) |>
  pivot_wider(names_from = picked, values_from = pct, values_fill = 0)

table_scores <- all_res |>
  group_by(dgp, n, fitted) |>
  summarise(
    log_score = mean(log_score),
    crps = mean(crps),
    .groups = "drop"
  )

write_csv(table_selection, proj_path("results", "table_cross_selection.csv"))
write_csv(table_scores, proj_path("results", "table_cross_scores.csv"))
print(table_selection, n = Inf)
print(table_scores, n = Inf)
