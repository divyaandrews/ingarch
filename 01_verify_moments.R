# Monte Carlo check of the NoGe-INGARCH(1,1) moments (Corollary 2.1) -----------
#
# Characterises the simulation scenarios (Table 1) and checks the closed-form
# moments against long simulated paths (Table S1). About 1 minute.
# Reviewer 1 (majors 2 and 6); Reviewer 2 (major 3).
# Output: results/table_scenarios.csv, results/table_scenarios.tex

source("R/00_functions.R")
set.seed(20261109)

scenarios <- tribble(
  ~scenario, ~alpha0, ~alpha1, ~beta1, ~phi,
  "I",       1.0,     0.2,     0.1,    0.05,
  "II",      1.0,     0.3,     0.1,    0.05,
  "III",     1.0,     0.4,     0.2,    0.55,
  "IV",      1.0,     0.4,     0.2,    0.35,
  "V",       0.5,     0.2,     0.7,    0.20
)

# Theoretical values -----------------------------------------------------------
theory <- scenarios |>
  mutate(moments = pmap(list(alpha0, alpha1, beta1, phi), noge_moments)) |>
  unnest(moments)

# Monte Carlo values -----------------------------------------------------------
simulate_summary <- function(alpha0, alpha1, beta1, phi) {
  pars <- list(alpha0 = alpha0, alpha1 = alpha1, beta1 = beta1,
               phi = phi, nu = NA)
  x <- sim_ingarch11(3e5, pars)

  tibble(
    mean_mc = mean(x),
    var_mc = var(x),
    rho_1_mc = cor(x[-1], x[-length(x)]),
    zero_mc = mean(x == 0)
  )
}

empirical <- scenarios |>
  mutate(mc = pmap(list(alpha0, alpha1, beta1, phi), simulate_summary)) |>
  unnest(mc) |>
  select(scenario, mean_mc:zero_mc)

# Output -----------------------------------------------------------------------
table_scenarios <- theory |>
  left_join(empirical, by = "scenario")

print(table_scenarios, width = Inf)
write_csv(table_scenarios, proj_path("results", "table_scenarios.csv"))

table_scenarios |>
  select(scenario, mu, var_x, disp_index, rho_1, p_zero, d) |>
  latex_rows(digits = 3, file = proj_path("results", "table_scenarios.tex"))

# In scenarios with d close to 0 (III: d = 0.089) the sample variance
# converges slowly because of the heavy right tail, so var_mc understates
# var_x; the mean, lag-one autocorrelation and zero proportion agree closely.
