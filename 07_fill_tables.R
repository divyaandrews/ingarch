# Paste-ready LaTeX rows and in-text numbers for the TeX files -----------------
#
# Run after 02, 03, 04, 05 and 06 (any script that has not been run yet is
# skipped). Writes, for every table that contains blue [bullet] placeholders,
# the complete table body in exactly the layout of the TeX file:
#
#   results/fill/table2_body.tex       Table 2   (replace everything between
#                                                 \midrule and \bottomrule)
#   results/fill/table3_body.tex       Table 3
#   results/fill/table4_body.tex       Table 4
#   results/fill/table5_body.tex       Table 5
#   results/fill/table6_body.tex       Table 6
#   results/fill/table7_body.tex       Table 7
#   results/fill/tableS2_body.tex      Table S2
#   results/fill/tableS3_body.tex      Table S3
#   results/fill/tableS4_body.tex      Table S4
#   results/fill/numbers_for_text.txt  every in-text placeholder, labelled by
#                                      document, section and line
#
# The CMLE numbers in these files come from R and replace the Python values.

source("R/00_functions.R")

fill_dir <- proj_path("results", "fill")
dir.create(fill_dir, recursive = TRUE, showWarnings = FALSE)

# Formatting helpers -----------------------------------------------------------
num <- function(x, digits = 3) {
  out <- formatC(x, format = "f", digits = digits)
  out <- ifelse(is.na(x), "--", out)
  gsub("^-", "$-$", out)
}

tex_par <- c(
  alpha0 = "$\\alpha_0$", alpha1 = "$\\alpha_1$",
  beta1 = "$\\beta_1$", phi = "$\\phi$"
)
par_order <- c("alpha0", "alpha1", "beta1", "phi")
model_order <- names(families)

write_body <- function(lines, name) {
  writeLines(lines, file.path(fill_dir, name))
  message("wrote ", name)
}

have <- function(...) file.exists(proj_path("results", ...))

text_lines <- character(0)
add_text <- function(where, value) {
  text_lines <<- c(text_lines, sprintf("%-58s %s", where, value))
}

# Simulation results (02) ------------------------------------------------------
if (dir.exists(proj_path("results", "sim"))) {
  sim <- list.files(proj_path("results", "sim"), full.names = TRUE) |>
    map(readRDS) |>
    list_rbind()

  sim_summary <- sim |>
    group_by(scenario, n, method, param, truth) |>
    summarise(
      avg_est = mean(est),
      bias = mean(est - truth),
      mab = mean(abs(est - truth)),
      mse = mean((est - truth)^2),
      rmse = sqrt(mse),
      cp = mean(lower <= truth & truth <= upper, na.rm = TRUE),
      width = mean(upper - lower, na.rm = TRUE),
      pct_na = 100 * mean(is.na(se)),
      avg_median = mean(post_median),
      .groups = "drop"
    )

  # Table 2: one row per scenario x n x method, 4 parameters x 3 criteria
  table2 <- character(0)
  scen_names <- sort(unique(sim_summary$scenario))
  scen_names <- intersect(c("I", "II", "III", "IV", "V"), scen_names)
  for (s in scen_names) {
    for (j in seq_along(c(50, 200, 500))) {
      n_j <- c(50, 200, 500)[j]
      for (m in c("CMLE", "Bayes")) {
        row <- sim_summary |>
          filter(scenario == s, n == n_j, method == m) |>
          arrange(match(param, par_order))
        cells <- paste(num(row$bias), num(row$rmse), num(row$cp, 2),
                       sep = " & ", collapse = " & ")
        lead <- if (j == 1 && m == "CMLE") {
          sprintf("\\multirow{6}{*}{%s}", s)
        } else {
          ""
        }
        n_cell <- if (m == "CMLE") {
          sprintf("\\multirow{2}{*}{%d}", n_j)
        } else {
          ""
        }
        table2 <- c(table2,
                    sprintf("%s & %s & %s & %s \\\\", lead, n_cell, m, cells))
      }
    }
    if (s != "V") table2 <- c(table2, "\\midrule")
  }
  write_body(table2, "table2_body.tex")

  # Table S2: average estimate, MAB, MSE, width (and NA% / posterior median)
  truths_txt <- c(
    I = "1, 0.2, 0.1, 0.05", II = "1, 0.3, 0.1, 0.05",
    III = "1, 0.4, 0.2, 0.55", IV = "1, 0.4, 0.2, 0.35",
    V = "0.5, 0.2, 0.7, 0.2"
  )
  table_s2 <- character(0)
  for (s in scen_names) {
    header <- paste0(
      "\\multicolumn{12}{l}{\\textit{Scenario ", s, "}: ",
      "$(\\alpha_0,\\alpha_1,\\beta_1,\\phi)=(", truths_txt[[s]], ")$}\\\\"
    )
    table_s2 <- c(table_s2, header)
    for (n_j in c(50, 200, 500)) {
      for (p in par_order) {
        cm <- filter(sim_summary, scenario == s, n == n_j,
                     method == "CMLE", param == p)
        by <- filter(sim_summary, scenario == s, n == n_j,
                     method == "Bayes", param == p)
        n_cell <- if (p == "alpha0") as.character(n_j) else ""
        table_s2 <- c(table_s2, paste0(
          n_cell, " & ", tex_par[[p]], " & ",
          num(cm$avg_est), " & ", num(cm$mab), " & ", num(cm$mse, 4), " & ",
          num(cm$width), " & ", num(cm$pct_na, 1), " & ",
          num(by$avg_est), " & ", num(by$mab), " & ", num(by$mse, 4), " & ",
          num(by$width), " & ", num(by$avg_median), " \\\\"
        ))
      }
    }
  }
  write_body(table_s2, "tableS2_body.tex")

  # Table 4: HMC diagnostics pooled over scenarios I-V, by n
  per_fit <- sim |>
    filter(method == "Bayes") |>
    group_by(scenario, n, rep) |>
    summarise(
      max_rhat = max(rhat),
      min_bulk = min(ess_bulk),
      min_tail = min(ess_tail),
      n_div = first(n_div),
      time_sec = first(time_sec),
      .groups = "drop"
    )
  table4 <- per_fit |>
    group_by(n) |>
    summarise(
      pct_rhat = 100 * mean(max_rhat > 1.01),
      bulk = median(min_bulk),
      tail = median(min_tail),
      pct_div = 100 * mean(n_div > 0),
      time = median(time_sec),
      .groups = "drop"
    )
  write_body(
    sprintf("I--V & %d & %s & %s & %s & %s & %s\\\\",
            table4$n, num(table4$pct_rhat, 1), num(table4$bulk, 0),
            num(table4$tail, 0), num(table4$pct_div, 1), num(table4$time, 1)),
    "table4_body.tex"
  )

  # In-text numbers, Sections 4.2 and 4.4 (manuscript) and the reply
  bayes_cp <- filter(sim_summary, method == "Bayes")$cp
  add_text("Manuscript 4.2: Bayes CP range, lower", num(min(bayes_cp), 2))
  add_text("Manuscript 4.2: Bayes CP range, upper", num(max(bayes_cp), 2))
  rmse_ratio <- sim_summary |>
    select(scenario, n, param, method, rmse) |>
    pivot_wider(names_from = method, values_from = rmse) |>
    mutate(ratio = Bayes / CMLE)
  add_text("Manuscript 4.2: median RMSE ratio Bayes/CMLE, n = 50",
           num(median(filter(rmse_ratio, n == 50)$ratio), 2))
  add_text("Manuscript 4.2: median RMSE ratio Bayes/CMLE, n = 500",
           num(median(filter(rmse_ratio, n == 500)$ratio), 2))
  add_text("Manuscript 4.4 / Reply R1-8: % fits max R-hat > 1.01",
           num(100 * mean(per_fit$max_rhat > 1.01), 1))
  add_text("Manuscript 4.4 / Reply R1-8: median min bulk ESS",
           num(median(per_fit$min_bulk), 0))
  add_text("Manuscript 4.4 / Reply R1-8: % fits with divergences",
           num(100 * mean(per_fit$n_div > 0), 1))
  add_text("Manuscript 4.4: median time per Bayes fit, n = 50 (s)",
           num(median(filter(per_fit, n == 50)$time_sec), 1))
  add_text("Manuscript 4.4: median time per Bayes fit, n = 500 (s)",
           num(median(filter(per_fit, n == 500)$time_sec), 1))
  cmle_time <- sim |>
    filter(method == "CMLE", n == 500) |>
    distinct(scenario, rep, time_sec)
  add_text("Reply R1-minor 3: median CMLE time, n = 500 (s)",
           num(median(cmle_time$time_sec), 3))
}

# Competing-model simulation (03) ----------------------------------------------
if (dir.exists(proj_path("results", "cross"))) {
  cross <- list.files(proj_path("results", "cross"), full.names = TRUE) |>
    map(readRDS) |>
    list_rbind()
  dgp_order <- c("NoGe (zero-inflated)", "NoGe (zero-deflated)", "ZINB",
                 "GP", "Poisson")
  table3 <- character(0)
  for (g in dgp_order) {
    for (n_j in c(200, 500)) {
      d <- filter(cross, dgp == g, n == n_j)
      picked <- d |>
        group_by(rep) |>
        slice_min(aic, n = 1, with_ties = FALSE) |>
        ungroup() |>
        dplyr::count(fitted, name = "k")
      pct <- vapply(model_order, function(f) {
        k <- picked$k[picked$fitted == f]
        if (length(k) == 0) 0 else 100 * k / sum(picked$k)
      }, numeric(1))
      ls <- d |>
        group_by(fitted) |>
        summarise(log_score = mean(log_score), .groups = "drop")
      ls <- ls$log_score[match(model_order, ls$fitted)]
      ls_txt <- num(ls)
      best <- which.min(ls)
      ls_txt[best] <- paste0("\\textbf{", ls_txt[best], "}")
      lead <- if (n_j == 200) g else ""
      table3 <- c(table3, paste0(
        lead, " & ", n_j, " & ",
        paste(num(pct, 1), collapse = " & "), " & ",
        paste(ls_txt, collapse = " & "), " \\\\"
      ))
    }
  }
  write_body(table3, "table3_body.tex")
}

# NUTS vs RWM (04) -------------------------------------------------------------
if (have("table_hmc_vs_mh.csv")) {
  hm <- read_csv(proj_path("results", "table_hmc_vs_mh.csv"),
                 show_col_types = FALSE)
  rows <- hm |>
    arrange(match(scenario, c("III", "V")), n, sampler) |>
    mutate(
      lead = if_else(sampler == "NUTS", scenario, ""),
      n_cell = if_else(sampler == "NUTS", as.character(n), "")
    )
  write_body(
    sprintf("%s & %s & %s & %s & %s & %s & %s & %s\\\\",
            rows$lead, rows$n_cell, rows$sampler, num(rows$time_sampling, 1),
            num(rows$min_bulk, 0), num(rows$min_tail, 0),
            num(rows$bulk_per_sec, 1), num(rows$bulk_per_1000, 1)),
    "tableS3_body.tex"
  )
  med <- hm |>
    group_by(sampler) |>
    summarise(per_1000 = median(bulk_per_1000), per_sec = median(bulk_per_sec))
  add_text("Manuscript 4.4 / Reply R2-7: NUTS ESS per 1000 draws",
           num(med$per_1000[med$sampler == "NUTS"], 1))
  add_text("Manuscript 4.4 / Reply R2-7: RWM ESS per 1000 draws",
           num(med$per_1000[med$sampler == "RWM"], 1))
  add_text("Manuscript 4.4 / Reply R2-7: NUTS ESS per second",
           num(med$per_sec[med$sampler == "NUTS"], 1))
  add_text("Manuscript 4.4 / Reply R2-7: RWM ESS per second",
           num(med$per_sec[med$sampler == "RWM"], 1))
}

# Prior sensitivity (06) -------------------------------------------------------
if (have("table_prior_sens_sim.csv")) {
  ps <- read_csv(proj_path("results", "table_prior_sens_sim.csv"),
                 show_col_types = FALSE)
  table_s4 <- character(0)
  for (s in c("III", "IV")) {
    for (p in par_order) {
      cells <- map_chr(c("default", "diffuse", "concentrated"), function(pr) {
        r <- filter(ps, scenario == s, variable == p, prior == pr)
        paste0(num(r$avg_est), " (", num(r$rmse), ") & ", num(r$coverage, 2))
      })
      lead <- if (p == "alpha0") s else ""
      table_s4 <- c(table_s4, paste0(lead, " & ", tex_par[[p]], " & ",
                                     paste(cells, collapse = " & "), "\\\\"))
    }
  }
  write_body(table_s4, "tableS4_body.tex")

  change <- ps |>
    group_by(scenario, variable) |>
    summarise(
      max_change = max(abs(avg_est - avg_est[prior == "default"])),
      .groups = "drop"
    )
  add_text("Manuscript 4.4 / Reply R2-5: max change of posterior means",
           num(max(change$max_change), 3))
}
if (have("table_prior_only.csv")) {
  po <- read_csv(proj_path("results", "table_prior_only.csv"),
                 show_col_types = FALSE) |>
    filter(prior == "default") |>
    arrange(match(variable, par_order))
  add_text("Supplement S2: prior-only means (alpha0, alpha1, beta1, phi)",
           paste(num(po$mean), collapse = ", "))
}

# Real data (05) ---------------------------------------------------------------
if (have("real_data_results.rds")) {
  real <- readRDS(proj_path("results", "real_data_results.rds"))
  fit_tab <- read_csv(proj_path("results", "real_fit.csv"),
                      show_col_types = FALSE)
  fc1 <- read_csv(proj_path("results", "real_fc1.csv"), show_col_types = FALSE)
  fch <- read_csv(proj_path("results", "real_fch.csv"), show_col_types = FALSE)
  diag_tab <- read_csv(proj_path("results", "real_diag.csv"),
                       show_col_types = FALSE)
  bayes <- read_csv(proj_path("results", "real_bayes.csv"),
                    show_col_types = FALSE)

  data_y <- list(
    polio = list(
      y = read_csv(proj_path("data", "polio.csv"), show_col_types = FALSE)$x,
      n_train = 132
    ),
    hepb = list(
      y = read_csv(proj_path("data", "hepb.csv"),
                   show_col_types = FALSE)$data[1:137],
      n_train = 110
    )
  )

  # CMLE-implied mean zero probability over the training period
  cmle_p0 <- function(y, family) {
    fit <- fit_cmle(y, family, se = FALSE)
    lam <- lambda_path(as.numeric(y), fit$pars, family)
    mean(dfam(rep(0, length(lam) - 1), lam[-1], family, fit$pars$phi, fit$pars$nu))
  }

  cell_est <- function(data_name, fam, p) {
    cm <- real[[data_name]][[fam]]$cmle |> filter(param == p)
    by <- bayes |> filter(data == data_name, family == fam, variable == p)
    if (nrow(cm) == 0) {
      return("-- & --")
    }
    paste0(num(cm$est), " (", num(cm$se), ") & ",
           num(by$mean), " [", num(by$sd), "]")
  }

  # Table 5: polio estimates
  table5 <- map_chr(model_order, function(fam) {
    paste0(fam, " & ", paste(
      map_chr(c(par_order, "nu"), function(p) cell_est("polio", fam, p)),
      collapse = " & "
    ), " \\\\")
  })
  write_body(table5, "table5_body.tex")

  crit_row <- function(data_name, fam, full) {
    f <- filter(fit_tab, data == data_name, family == fam)
    h <- filter(fc1, data == data_name, family == fam)
    d <- data_y[[data_name]]
    p0 <- cmle_p0(d$y[seq_len(d$n_train)], families[[fam]])
    best_aic <- min(filter(fit_tab, data == data_name)$aic)
    aic <- num(f$aic, 2)
    if (isTRUE(all.equal(f$aic, best_aic))) {
      aic <- paste0("\\textbf{", aic, "}")
    }
    if (full) {
      paste(fam, num(f$loglik, 2), aic, num(f$bic, 2), num(p0),
            num(f$elpd_loo, 2), num(f$neg_sum_log_cpo, 2),
            num(f$mean_crps_in), num(h$prmse), num(h$pmad), num(h$log_score),
            num(h$crps), num(h$cover95, 2), sep = " & ")
    } else {
      paste(fam, num(f$loglik, 2), aic, num(f$bic, 2), num(p0),
            num(f$elpd_loo, 2), num(h$prmse), num(h$pmad),
            num(h$log_score), num(h$crps), sep = " & ")
    }
  }
  write_body(paste0(map_chr(model_order, crit_row, data_name = "polio",
                            full = TRUE), " \\\\"), "table6_body.tex")
  write_body(paste0(map_chr(model_order, crit_row, data_name = "hepb",
                            full = FALSE), " \\\\"), "table7_body.tex")

  # In-text numbers, Sections 6.1 and 6.2 and the reply
  for (dn in c("polio", "hepb", "transactions")) {
    add_text(sprintf("Manuscript 6 / Reply R2-3: P(D > 0 | data), %s", dn),
             num(filter(diag_tab, data == dn, family == "NoGe")$prob_d_pos))
  }
  phi_polio <- filter(bayes, data == "polio", family == "NoGe",
                      variable == "phi")
  add_text("Reply R2-11: 95% CrI of phi (polio, NoGe)",
           sprintf("(%s, %s)", num(phi_polio$q2.5), num(phi_polio$q97.5)))

  for (dn in c("polio", "hepb")) {
    ranks <- filter(fc1, data == dn) |>
      mutate(
        rank_logs = rank(log_score),
        rank_crps = rank(crps),
        rank_pmad = rank(pmad),
        rank_prmse = rank(prmse)
      ) |>
      filter(family == "NoGe")
    add_text(
      sprintf("Manuscript 6.%d: NoGe hold-out ranks (of 6), %s",
              if (dn == "polio") 1 else 2, dn),
      sprintf("LogS %s, CRPS %s, PMAD %s, PRMSE %s", ranks$rank_logs,
              ranks$rank_crps, ranks$rank_pmad, ranks$rank_prmse)
    )
    add_text(
      sprintf("Manuscript 6: NoGe hold-out values, %s", dn),
      sprintf("LogS %s, CRPS %s, PMAD %s, PRMSE %s", num(ranks$log_score),
              num(ranks$crps), num(ranks$pmad), num(ranks$prmse))
    )
  }

  best_comp <- fc1 |>
    filter(data == "polio", family != "NoGe") |>
    slice_min(log_score, n = 1, with_ties = FALSE) |>
    pull(family)
  for (h_j in 2:3) {
    add_text(sprintf("Supplement S5: polio CRPS h = %d, NoGe", h_j),
             num(filter(fch, data == "polio", family == "NoGe", h == h_j)$crps))
    add_text(sprintf("Supplement S5: polio CRPS h = %d, %s", h_j, best_comp),
             num(filter(fch, data == "polio", family == best_comp,
                        h == h_j)$crps))
  }
}

add_text("Manuscript Declarations / Reply minor 8: code URL",
         "your GitHub or Zenodo link")
writeLines(text_lines, file.path(fill_dir, "numbers_for_text.txt"))
message("wrote numbers_for_text.txt")
cat(text_lines, sep = "\n")
