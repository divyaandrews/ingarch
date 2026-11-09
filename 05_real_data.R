# Real-data analysis: six INGARCH(1,1) families x three data sets --------------
#
#   polio        US poliomyelitis cases per month, Jan 1970 - Dec 1983
#                (n = 168); training 1970-1980 (132), test 1981-1983 (36)
#   hepb         weekly hepatitis-B cases (Bremen, Hamburg, Saxony-Anhalt),
#                RKI SurvStat; first 137 weeks, training 110, test 27
#   transactions Ericsson B transactions per minute, 2 July 2002 (n = 460);
#                training 09:35-15:42 (368), test 15:43-17:14 (92)
#
# For every family:
#   CMLE: estimates, standard errors, log-likelihood, AIC, BIC.
#   Bayes: posterior mean, sd, 95% CrI, R-hat, ESS, divergences, P(D > 0).
#   In-sample: PSIS-LOO, -sum log CPO, mean CRPS, non-randomised PIT,
#     Pearson-residual ACF and Ljung-Box test.
#   Zeros: observed proportion, model-implied mean P(X_t = 0 | F_{t-1}) and a
#     posterior-predictive p-value for the number of zeros (Reviewer 2, 11).
#   Hold-out: h = 1 from the exact posterior-mixture pmf (PRMSE, PMAD, log
#     score, CRPS, 95% coverage) and h = 1, 2, 3 from simulated
#     posterior-predictive paths (Reviewer 2, 9).

source("R/00_functions.R")
library(matrixStats)
set.seed(2026)

# Data -------------------------------------------------------------------------
read_polio <- function() {
  if (requireNamespace("gamlss.data", quietly = TRUE)) {
    env <- new.env()
    data("polio", package = "gamlss.data", envir = env)
    return(as.integer(env$polio))
  }
  read_csv(proj_path("data", "polio.csv"), show_col_types = FALSE)$x
}

hepb <- read_csv(proj_path("data", "hepb.csv"), show_col_types = FALSE)
datasets <- list(
  polio = list(y = read_polio(), n_train = 132),
  hepb = list(y = hepb$data[1:137], n_train = 110),
  transactions = list(
    y = scan(proj_path("data", "EricssonB_Jul2.txt"), quiet = TRUE),
    n_train = 368
  )
)
stopifnot(
  length(datasets$polio$y) == 168,
  length(datasets$transactions$y) == 460
)

h_max <- 3
m_score <- 1000   # thinned posterior draws used for predictive pmfs

# Analysis of one family on one data set ---------------------------------------
analyse_one <- function(y_all, n_train, family) {
  y <- y_all[seq_len(n_train)]
  y_test <- y_all[-seq_len(n_train)]
  k_max <- max(y_all) * 3 + 60

  ml <- fit_cmle(y, family)

  fit <- fit_bayes(
    y, family,
    y_test = y_test,
    h_max = h_max,
    gq_full = 1L,
    cores = 4
  )
  hmc <- hmc_diagnostics(fit, family)
  draws <- rstan::extract(fit)
  m <- length(draws$alpha0)
  idx <- round(seq(1, m, length.out = m_score))
  phi_d <- if (has_phi(family)) draws$phi[idx] else rep(NA_real_, m_score)
  nu_d <- if (has_nu(family)) draws$nu[idx] else rep(NA_real_, m_score)

  # PSIS-LOO and CPO from the one-step conditional log-likelihood terms
  ll <- draws$log_lik
  ll_arr <- extract_log_lik(fit, "log_lik", merge_chains = FALSE)
  loo_res <- suppressWarnings(
    loo(ll_arr, r_eff = relative_eff(exp(ll_arr)))
  )
  log_cpo <- -(colLogSumExps(-ll) - log(m))

  # In-sample one-step predictive distributions
  in_sample <- map(2:n_train, function(t) {
    p <- pred_pmf(draws$lambda[idx, t], family, phi_d, nu_d, k_max)
    grid <- seq_along(p) - 1
    pmf_summary(p, y[t]) |>
      mutate(var = sum(grid^2 * p) - sum(grid * p)^2)
  }) |>
    list_rbind()
  pit <- pit_nonrandom(in_sample$p_lo, in_sample$p_hi)
  pearson <- (y[-1] - in_sample$mean) / sqrt(in_sample$var)

  zeros <- tibble(
    zero_obs = mean(y[-1] == 0),
    p0_model = mean(in_sample$p_zero),
    zero_ppc_p = mean(rowMeans(draws$y_rep == 0) >= mean(y[-1] == 0))
  )

  # Hold-out, h = 1, exact posterior-mixture pmf
  out_h1 <- map(seq_along(y_test), function(k) {
    p <- pred_pmf(draws$lambda_test[idx, k], family, phi_d, nu_d, k_max)
    pmf_summary(p, y_test[k])
  }) |>
    list_rbind() |>
    mutate(obs = y_test, t = n_train + seq_along(y_test))

  # Hold-out, h = 1..h_max, simulated posterior-predictive paths
  out_h <- map(seq_len(h_max), function(h) {
    map(seq_len(length(y_test) - h + 1), function(k) {
      sims <- draws$y_path[, k, h]
      target <- y_test[k + h - 1]
      tibble(
        h = h,
        obs = target,
        mean = mean(sims),
        median = median(sims),
        lower = quantile(sims, 0.025, type = 1),
        upper = quantile(sims, 0.975, type = 1),
        crps = scoringRules::crps_sample(target, sims)
      )
    }) |>
      list_rbind()
  }) |>
    list_rbind()

  prob_d_pos <- if (family == 1L) {
    mean(draws$second_order_D > 0)
  } else {
    NA_real_
  }

  list(
    cmle = tibble(param = names(ml$est), est = ml$est, se = ml$se),
    bayes = hmc$summary,
    diag = tibble(
      max_rhat = hmc$max_rhat,
      min_ess_bulk = hmc$min_ess_bulk,
      min_ess_tail = hmc$min_ess_tail,
      n_div = hmc$n_divergent,
      n_td = hmc$n_max_treedepth,
      min_ebfmi = hmc$min_ebfmi,
      time_sec = hmc$time_sec,
      prob_d_pos = prob_d_pos
    ),
    fit = tibble(
      loglik = ml$loglik,
      aic = ml$aic,
      bic = ml$bic,
      elpd_loo = loo_res$estimates["elpd_loo", "Estimate"],
      se_elpd_loo = loo_res$estimates["elpd_loo", "SE"],
      neg_sum_log_cpo = -sum(log_cpo),
      mean_crps_in = mean(in_sample$crps),
      pit_ks_p = ks.test(
        pit_random(in_sample$p_lo, in_sample$p_hi), "punif"
      )$p.value,
      ljung_box_p = Box.test(pearson, lag = 10, type = "Ljung-Box")$p.value,
      pearson_var = var(pearson)
    ),
    zeros = zeros,
    pit = pit,
    resid_acf = tibble(
      lag = 1:10,
      acf = acf(pearson, lag.max = 10, plot = FALSE)$acf[-1]
    ),
    fc1 = out_h1 |>
      summarise(
        prmse = sqrt(mean((obs - mean)^2)),
        pmad = mean(abs(obs - median)),
        log_score = mean(log_score),
        crps = mean(crps),
        cover95 = mean(obs >= lower & obs <= upper),
        width95 = mean(upper - lower)
      ),
    fch = out_h |>
      group_by(h) |>
      summarise(
        prmse = sqrt(mean((obs - mean)^2)),
        pmad = mean(abs(obs - median)),
        crps = mean(crps),
        cover95 = mean(obs >= lower & obs <= upper),
        .groups = "drop"
      ),
    fc1_detail = out_h1,
    loo = loo_res
  )
}

# Run all models ---------------------------------------------------------------
results <- imap(datasets, function(d, data_name) {
  imap(families, function(fam, fam_name) {
    message(data_name, " / ", fam_name)
    analyse_one(d$y, d$n_train, fam)
  })
})
saveRDS(results, proj_path("results", "real_data_results.rds"))

# Tables -----------------------------------------------------------------------
collect <- function(slot) {
  imap(results, function(res_data, data_name) {
    imap(res_data, function(x, fam_name) {
      x[[slot]] |> mutate(data = data_name, family = fam_name)
    }) |>
      list_rbind()
  }) |>
    list_rbind()
}

slots <- c("fit", "zeros", "fc1", "fch", "diag", "bayes", "cmle")
tabs <- set_names(map(slots, collect), slots)
iwalk(tabs, function(tab, name) {
  write_csv(tab, proj_path("results", paste0("real_", name, ".csv")))
})

table_loo <- imap(results, function(res_data, data_name) {
  comparison <- loo_compare(map(res_data, "loo"))
  tibble(
    data = data_name,
    model = rownames(comparison),
    elpd_diff = comparison[, "elpd_diff"],
    se_diff = comparison[, "se_diff"]
  )
}) |>
  list_rbind()
write_csv(table_loo, proj_path("results", "real_loo_compare.csv"))

# LaTeX rows for the manuscript tables
for (data_name in names(datasets)) {
  estimates <- tabs$cmle |>
    filter(data == data_name) |>
    left_join(
      tabs$bayes |>
        filter(data == data_name) |>
        rename(param = variable),
      by = c("data", "family", "param")
    ) |>
    transmute(
      family,
      param,
      cmle = sprintf("%.3f (%.3f)", est, se),
      bayes = sprintf("%.3f [%.3f]", mean, sd),
      cri = sprintf("(%.3f, %.3f)", q2.5, q97.5)
    )
  latex_rows(
    estimates,
    file = proj_path("results", sprintf("tex_%s_estimates.tex", data_name))
  )

  criteria <- tabs$fit |>
    filter(data == data_name) |>
    left_join(filter(tabs$zeros, data == data_name),
              by = c("data", "family")) |>
    left_join(filter(tabs$fc1, data == data_name),
              by = c("data", "family")) |>
    transmute(
      family, aic, bic, elpd_loo, neg_sum_log_cpo, mean_crps_in, p0_model,
      prmse, pmad, log_score, crps, cover95
    )
  latex_rows(
    criteria,
    file = proj_path("results", sprintf("tex_%s_criteria.tex", data_name))
  )
}

# Figures ----------------------------------------------------------------------
for (data_name in names(datasets)) {
  y <- datasets[[data_name]]$y
  n_train <- datasets[[data_name]]$n_train

  plot_ts <- ggplot(tibble(t = seq_along(y), y = y), aes(t, y)) +
    geom_line(colour = "grey30") +
    geom_point(size = 0.6) +
    geom_vline(xintercept = n_train + 0.5, linetype = 2) +
    labs(x = "Time", y = "Count") +
    theme_bw()
  ggsave(proj_path("figures", sprintf("%s_ts.pdf", data_name)),
         plot_ts, width = 7, height = 2.6)

  pit_df <- imap(results[[data_name]], function(x, fam_name) {
    x$pit |> mutate(family = fam_name)
  }) |>
    list_rbind() |>
    mutate(family = factor(family, names(families)))
  plot_pit <- ggplot(
    pit_df,
    aes(xmin = bin_lower, xmax = bin_upper, ymin = 0, ymax = density)
  ) +
    geom_rect(fill = "grey75", colour = "grey30") +
    geom_hline(yintercept = 1, linetype = 2) +
    facet_wrap(~family, nrow = 1) +
    labs(x = "PIT", y = "Density") +
    theme_bw()
  ggsave(proj_path("figures", sprintf("%s_pit.pdf", data_name)),
         plot_pit, width = 9, height = 2.4)

  best_comp <- tabs$fc1 |>
    filter(data == data_name, family != "NoGe") |>
    slice_min(log_score, n = 1, with_ties = FALSE) |>
    pull(family)
  fc_df <- bind_rows(
    results[[data_name]]$NoGe$fc1_detail |> mutate(model = "NoGe"),
    results[[data_name]][[best_comp]]$fc1_detail |> mutate(model = best_comp)
  )
  plot_fc <- ggplot(fc_df, aes(t)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey85") +
    geom_line(aes(y = mean, linetype = "Predictive mean")) +
    geom_step(aes(y = median, linetype = "Predictive median"),
              direction = "mid") +
    geom_point(aes(y = obs), size = 1) +
    facet_wrap(~model) +
    labs(x = "Time", y = "Count", linetype = NULL) +
    theme_bw() +
    theme(legend.position = "bottom")
  ggsave(proj_path("figures", sprintf("%s_forecast.pdf", data_name)),
         plot_fc, width = 8, height = 3)
}

print(tabs$fit, width = Inf)
print(tabs$zeros, width = Inf)
print(tabs$fc1, width = Inf)
print(table_loo)
