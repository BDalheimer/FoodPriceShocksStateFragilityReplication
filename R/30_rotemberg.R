# =============================================================================
# 30_rotemberg.R
# Goldsmith-Pinkham, Sorkin & Swift (2020) Rotemberg weight decomposition
#
# The aggregate consumption shock B_{it} = Σ_g B_{igt} is a sum of commodity-
# group components. Following GPS (2020), we decompose the LP coefficient as:
#
#   β̂_OLS = Σ_g α̂_g · β̂_g^IV
#
# where:
#   α̂_g = Cov(B̃_g, B̃) / Var(B̃)   [Rotemberg weight; sums to 1]
#   β̂_g^IV = just-identified IV using only B_{igt} as instrument for B_{it}
#
# Residualization (tildes): year FE + two shock lags, matching the LP spec.
#
# Focus: consumption exposure (our main finding), horizons h=1 and h=4.
# Output: output/tables/tab_rotemberg.tex
# =============================================================================

source(here::here("R", "00_setup.R"))

panel <- readRDS(file.path(DATA_PRO, "panel.rds"))

message("Panel: ", nrow(panel), " obs | ",
        n_distinct(panel$iso3c), " countries | ",
        min(panel$year), "–", max(panel$year))

# ---- Constants --------------------------------------------------------------
GROUPS <- c("cereals", "oils", "animal", "sugar_fruit", "cash_crops")

GROUP_LABELS <- c(
  cereals     = "Cereals",
  oils        = "Oilseeds \\& oils",
  animal      = "Animal products",
  sugar_fruit = "Sugar \\& fruit",
  cash_crops  = "Cash crops"
)

HORIZONS_REPORT <- c(1, 4)

# Group shock columns are already in panel.rds (built by 07_panel.R)
panel_r <- panel

grp_cols <- paste0("shock_cons_", GROUPS)
message("Group shock columns available: ", paste(grp_cols, collapse = ", "))

# ---- Helpers ----------------------------------------------------------------
stars_fn <- function(p) {
  case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.10  ~ "$\\cdot$",
    TRUE      ~ ""
  )
}

write_tex <- function(lines, filename) {
  path <- file.path(OUT_TABS, filename)
  writeLines(lines, path)
  message("  Saved: ", filename)
}

# =============================================================================
# Section 1: Rotemberg weights and just-identified IV coefficients
# =============================================================================
message("\n=== Section 1: Rotemberg decomposition ===")

# Compute Rotemberg weight and just-identified IV coefficient for one group + horizon
run_rotemberg <- function(data, group, h) {
  gc    <- paste0("shock_cons_", group)
  y_col <- paste0("dfsi_h", h)

  d <- data |>
    filter(!is.na(.data[[y_col]]),
           !is.na(shock_cons),
           !is.na(.data[[gc]]),
           !is.na(shock_cons_l1),
           !is.na(shock_cons_l2))

  if (nrow(d) < 20) return(NULL)

  # Residualize aggregate shock on year FE + 2 lags (same controls as LP)
  fit_agg <- feols(shock_cons ~ shock_cons_l1 + shock_cons_l2 | year,
                   data = d, notes = FALSE, warn = FALSE)
  r_agg   <- residuals(fit_agg)

  # Residualize group shock on same controls
  fml_g  <- as.formula(paste0(gc, " ~ shock_cons_l1 + shock_cons_l2 | year"))
  fit_g  <- feols(fml_g, data = d, notes = FALSE, warn = FALSE)
  r_g    <- residuals(fit_g)

  # Rotemberg weight: Cov(r_g, r_agg) / Var(r_agg)
  alpha_g <- sum(r_g * r_agg, na.rm = TRUE) / sum(r_agg^2, na.rm = TRUE)

  # Just-identified IV: instrument shock_cons with shock_cons_{group}
  fml_iv <- as.formula(paste0(
    y_col, " ~ shock_cons_l1 + shock_cons_l2 | year | shock_cons ~ ", gc
  ))
  fit_iv <- tryCatch(
    feols(fml_iv, data = d, cluster = ~iso3c, notes = FALSE, warn = FALSE),
    error = function(e) NULL
  )

  if (is.null(fit_iv)) {
    return(tibble(group = group, h = h, alpha = alpha_g,
                  est = NA_real_, se = NA_real_, pval = NA_real_, n_obs = nrow(d)))
  }

  ct        <- summary(fit_iv)$coeftable
  coef_name <- grep("^fit_shock_cons", rownames(ct), value = TRUE)[1]

  tibble(
    group = group,
    h     = h,
    alpha = alpha_g,
    est   = ct[coef_name, "Estimate"],
    se    = ct[coef_name, "Std. Error"],
    pval  = ct[coef_name, "Pr(>|t|)"],
    n_obs = nobs(fit_iv)
  )
}

# Run over all groups × horizons
results <- map_dfr(HORIZONS_REPORT, function(h) {
  map_dfr(GROUPS, \(g) run_rotemberg(panel_r, g, h))
})

# ---- OLS baseline for verification ------------------------------------------
baseline <- map_dfr(HORIZONS_REPORT, function(h) {
  y_col <- paste0("dfsi_h", h)
  d <- panel_r |>
    filter(!is.na(.data[[y_col]]), !is.na(shock_cons),
           !is.na(shock_cons_l1), !is.na(shock_cons_l2))
  fit <- feols(as.formula(paste0(
    y_col, " ~ shock_cons + shock_cons_l1 + shock_cons_l2 | year"
  )), data = d, cluster = ~iso3c, notes = FALSE, warn = FALSE)
  ct <- summary(fit)$coeftable
  tibble(h = h, beta_ols = ct["shock_cons", "Estimate"], n_ols = nobs(fit))
})

# ---- Verification: Σ α̂_g · β̂_g^IV ≈ β̂_OLS --------------------------------
verification <- results |>
  group_by(h) |>
  summarise(
    implied_beta = sum(alpha * est, na.rm = TRUE),
    sum_alpha    = sum(alpha, na.rm = TRUE)
  ) |>
  left_join(baseline, by = "h")

message("\nVerification (Σ α̂_g · β̂_g^IV vs OLS):")
print(verification)

# ---- Console output ---------------------------------------------------------
cat("\nRotemberg weights and just-identified IV coefficients (consumption, h=1 and h=4):\n")
cat(sprintf("%-20s  %7s  %8s  %6s  %5s\n", "Group", "alpha", "beta_IV", "se", "p"))
cat(strrep("-", 55), "\n")

for (h in HORIZONS_REPORT) {
  b_ols <- filter(baseline, h == !!h)$beta_ols
  cat(sprintf("h = %d  (OLS baseline: %.3f)\n", h, b_ols))
  res_h <- filter(results, h == !!h) |> arrange(desc(alpha))
  for (i in seq_len(nrow(res_h))) {
    r <- res_h[i, ]
    stars <- ifelse(is.na(r$pval), "", ifelse(r$pval < 0.001, "***",
                    ifelse(r$pval < 0.01, "**", ifelse(r$pval < 0.05, "*",
                    ifelse(r$pval < 0.10, ".", "")))))
    cat(sprintf("  %-18s  %+7.3f  %+8.3f  %6.3f  %.3f %s\n",
                GROUP_LABELS[r$group], r$alpha, r$est, r$se, r$pval, stars))
  }
  cat("\n")
}

# =============================================================================
# Section 2: LaTeX table
# tab:rotemberg
# =============================================================================
message("\n=== Section 2: LaTeX table ===")

tab <- c(
  "% Rotemberg weight decomposition of consumption LP — GPS (2020)",
  "",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Rotemberg weight decomposition of the consumption shock (GPS 2020)}",
  "\\label{tab:rotemberg}",
  "\\begin{adjustbox}{max width=\\textwidth}",
  "\\begin{threeparttable}\\small",
  "\\begin{tabular}{l rr@{\\quad}rr}",
  "\\toprule",
  " & \\multicolumn{2}{c}{$h=1$} & \\multicolumn{2}{c}{$h=4$} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  paste0("Commodity group & $\\hat{\\alpha}_g$ & $\\hat{\\beta}_g^{\\text{IV}}$ ",
         "& $\\hat{\\alpha}_g$ & $\\hat{\\beta}_g^{\\text{IV}}$ \\\\"),
  "\\midrule"
)

for (g in GROUPS) {
  res_g <- filter(results, group == g)
  r1    <- filter(res_g, h == 1)
  r4    <- filter(res_g, h == 4)

  if (nrow(r1) == 0 || nrow(r4) == 0) next

  est_row <- sprintf(
    "%s & $%+.3f$ & $%+.3f$%s & $%+.3f$ & $%+.3f$%s \\\\",
    GROUP_LABELS[g],
    r1$alpha, r1$est, stars_fn(r1$pval),
    r4$alpha, r4$est, stars_fn(r4$pval)
  )
  se_row <- sprintf(
    "   & & $(%.3f)$ & & $(%.3f)$ \\\\[4pt]",
    r1$se, r4$se
  )
  tab <- c(tab, est_row, se_row)
}

# OLS baseline row
b1 <- filter(baseline, h == 1)$beta_ols
b4 <- filter(baseline, h == 4)$beta_ols
n1 <- filter(baseline, h == 1)$n_ols
n4 <- filter(baseline, h == 4)$n_ols

tab <- c(tab,
  "\\midrule",
  sprintf(
    "\\textit{OLS baseline (weighted sum)} & $1.000$ & $%+.3f$ & $1.000$ & $%+.3f$ \\\\",
    b1, b4
  ),
  "\\midrule",
  sprintf("$N$ & \\multicolumn{2}{c}{%s} & \\multicolumn{2}{c}{%s} \\\\",
          format(n1, big.mark = ","), format(n4, big.mark = ",")),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0(
    "\\item \\textit{Notes:} Rotemberg weight decomposition following ",
    "\\citet{GoldsmithPinkham2020}. ",
    "The aggregate consumption shock $B_{it} = \\sum_g B_{igt}$ is decomposed ",
    "into five commodity-group components. ",
    "$\\hat{\\alpha}_g = \\widetilde{\\text{Cov}}(B_g, B) / \\widetilde{\\text{Var}}(B)$ ",
    "is the Rotemberg weight, measuring group $g$'s contribution to aggregate ",
    "shock variation after partialing out year fixed effects and two shock lags. ",
    "Weights sum to one; a negative weight indicates the group's shock variation ",
    "is negatively correlated with the aggregate after conditioning. ",
    "$\\hat{\\beta}_g^{\\text{IV}}$ is the just-identified IV coefficient from using ",
    "only $B_{igt}$ as the instrument for $B_{it}$ in the LP; ",
    "standard errors clustered by country in parentheses. ",
    "The OLS baseline equals $\\sum_g \\hat{\\alpha}_g \\hat{\\beta}_g^{\\text{IV}}$ by construction. ",
    "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."
  ),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{adjustbox}",
  "\\end{table}"
)

write_tex(tab, "tab_rotemberg.tex")

# =============================================================================
# Section 3: Commodity-level Rotemberg weights
# Build individual commodity shocks from shares × price changes, then
# run the same decomposition for each of the 22 commodities.
# Output: output/tables/tab_rotemberg_commodity.tex
# =============================================================================
message("\n=== Section 3: Commodity-level Rotemberg weights ===")

shares <- readRDS(file.path(DATA_PRO, "fao_shares.rds"))
prices <- readRDS(file.path(DATA_PRO, "pinksheet_annual.rds"))

COMMODITY_LABELS <- c(
  wheat        = "Wheat",
  rice         = "Rice",
  maize        = "Maize",
  barley       = "Barley",
  sorghum      = "Sorghum",
  soybeans     = "Soybeans",
  soybean_oil  = "Soybean oil",
  palm_oil     = "Palm oil",
  groundnut_oil = "Groundnut oil",
  sunflower_oil = "Sunflower oil",
  coconut_oil  = "Coconut oil",
  rapeseed_oil = "Rapeseed oil",
  sugar        = "Sugar",
  bananas      = "Bananas",
  oranges      = "Oranges",
  beef         = "Beef",
  chicken      = "Chicken",
  lamb         = "Lamb",
  shrimp       = "Shrimp",
  coffee       = "Coffee",
  tea          = "Tea",
  cocoa        = "Cocoa"
)

# Group membership for table organisation
COMMODITY_GROUPS <- c(
  wheat = "Cereals", rice = "Cereals", maize = "Cereals",
  barley = "Cereals", sorghum = "Cereals",
  soybeans = "Oils", soybean_oil = "Oils", palm_oil = "Oils",
  groundnut_oil = "Oils", sunflower_oil = "Oils",
  coconut_oil = "Oils", rapeseed_oil = "Oils",
  sugar = "Sugar \\& fruit", bananas = "Sugar \\& fruit", oranges = "Sugar \\& fruit",
  beef = "Animal", chicken = "Animal", lamb = "Animal", shrimp = "Animal",
  coffee = "Cash crops", tea = "Cash crops", cocoa = "Cash crops"
)

# ---- Build individual commodity shocks --------------------------------------
# s_cons is the caloric consumption share (time-varying)
# commodity shock c: s_{ic,t} × dlog_p_{c,t}
commodity_shocks_long <- shares |>
  select(iso3c, year, pink_var, s_cons) |>
  mutate(s_cons = replace_na(s_cons, 0)) |>
  inner_join(prices |> select(year, pink_var, dlog_p), by = c("pink_var", "year")) |>
  filter(!is.na(dlog_p), year >= SAMPLE_START, year <= SAMPLE_END) |>
  mutate(shock_c = s_cons * dlog_p)

COMMODITIES <- unique(commodity_shocks_long$pink_var)
message("Commodities: ", paste(sort(COMMODITIES), collapse = ", "))

# Pivot to wide: one column per commodity
commodity_shocks_wide <- commodity_shocks_long |>
  select(iso3c, year, pink_var, shock_c) |>
  pivot_wider(names_from = pink_var, values_from = shock_c,
              names_prefix = "sc_")

# Merge into panel
panel_c <- panel |>
  left_join(commodity_shocks_wide, by = c("iso3c", "year"))

sc_cols <- paste0("sc_", COMMODITIES)
message("Commodity shock columns: ", length(sc_cols))

# ---- Rotemberg decomposition at commodity level ------------------------------
run_rotemberg_c <- function(data, commodity, h) {
  sc  <- paste0("sc_", commodity)
  y_col <- paste0("dfsi_h", h)

  d <- data |>
    filter(!is.na(.data[[y_col]]),
           !is.na(shock_cons),
           !is.na(.data[[sc]]),
           !is.na(shock_cons_l1),
           !is.na(shock_cons_l2))

  if (nrow(d) < 20) return(NULL)

  fit_agg <- feols(shock_cons ~ shock_cons_l1 + shock_cons_l2 | year,
                   data = d, notes = FALSE, warn = FALSE)
  r_agg   <- residuals(fit_agg)

  fml_c  <- as.formula(paste0(sc, " ~ shock_cons_l1 + shock_cons_l2 | year"))
  fit_c  <- feols(fml_c, data = d, notes = FALSE, warn = FALSE)
  r_c    <- residuals(fit_c)

  alpha_c <- sum(r_c * r_agg, na.rm = TRUE) / sum(r_agg^2, na.rm = TRUE)

  fml_iv <- as.formula(paste0(
    y_col, " ~ shock_cons_l1 + shock_cons_l2 | year | shock_cons ~ ", sc
  ))
  fit_iv <- tryCatch(
    feols(fml_iv, data = d, cluster = ~iso3c, notes = FALSE, warn = FALSE),
    error = function(e) NULL
  )

  if (is.null(fit_iv)) {
    return(tibble(commodity = commodity, h = h, alpha = alpha_c,
                  est = NA_real_, se = NA_real_, pval = NA_real_, n_obs = nrow(d)))
  }

  ct        <- summary(fit_iv)$coeftable
  coef_name <- grep("^fit_shock_cons", rownames(ct), value = TRUE)[1]

  tibble(
    commodity = commodity,
    h         = h,
    alpha     = alpha_c,
    est       = ct[coef_name, "Estimate"],
    se        = ct[coef_name, "Std. Error"],
    pval      = ct[coef_name, "Pr(>|t|)"],
    n_obs     = nobs(fit_iv)
  )
}

results_c <- map_dfr(HORIZONS_REPORT, function(h) {
  message("  Running h=", h, " ...")
  map_dfr(COMMODITIES, \(c) run_rotemberg_c(panel_c, c, h))
})

# ---- Verification -----------------------------------------------------------
verif_c <- results_c |>
  group_by(h) |>
  summarise(
    implied_beta = sum(alpha * est, na.rm = TRUE),
    sum_alpha    = sum(alpha, na.rm = TRUE)
  ) |>
  left_join(baseline, by = "h")

message("\nVerification — commodity level (Σ α̂_c · β̂_c^IV vs OLS):")
print(verif_c)

# ---- Console output ---------------------------------------------------------
cat("\nCommodity-level Rotemberg weights (h=4, sorted by |alpha|):\n")
cat(sprintf("%-16s  %-14s  %7s  %8s  %6s\n",
            "Commodity", "Group", "alpha", "beta_IV", "p"))
cat(strrep("-", 58), "\n")

res_c4 <- filter(results_c, h == 4) |>
  arrange(desc(abs(alpha))) |>
  mutate(grp = COMMODITY_GROUPS[commodity])

for (i in seq_len(nrow(res_c4))) {
  r <- res_c4[i, ]
  stars <- ifelse(is.na(r$pval), "",
           ifelse(r$pval < 0.001, "***", ifelse(r$pval < 0.01, "**",
           ifelse(r$pval < 0.05, "*", ifelse(r$pval < 0.10, ".", "")))))
  cat(sprintf("%-16s  %-14s  %+7.4f  %+8.2f  %.3f %s\n",
              COMMODITY_LABELS[r$commodity], r$grp,
              r$alpha, r$est, r$pval, stars))
}

# ---- LaTeX table: commodity level, h=4 (sorted by alpha desc) ---------------
message("\n=== Section 4: LaTeX table — commodity level ===")

# Order: within each group, by descending alpha
GROUP_ORDER <- c("Cereals", "Oils", "Sugar \\& fruit", "Animal", "Cash crops")

res_c_tab <- results_c |>
  filter(h == 4) |>
  mutate(
    grp       = COMMODITY_GROUPS[commodity],
    grp_f     = factor(grp, levels = GROUP_ORDER),
    label     = COMMODITY_LABELS[commodity]
  ) |>
  arrange(grp_f, desc(alpha))

# Also pull h=1 for the same commodity order
res_c1 <- results_c |>
  filter(h == 1) |>
  select(commodity, alpha_h1 = alpha, est_h1 = est, se_h1 = se, pval_h1 = pval)

res_c_tab <- res_c_tab |>
  left_join(res_c1, by = "commodity")

b4c <- filter(verif_c, h == 4)$beta_ols
b1c <- filter(verif_c, h == 1)$beta_ols
n4c <- nrow(filter(panel_c, !is.na(dfsi_h4), !is.na(shock_cons)))
n1c <- nrow(filter(panel_c, !is.na(dfsi_h1), !is.na(shock_cons)))

tab_c <- c(
  "% Rotemberg weight decomposition — individual commodity level",
  "",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Rotemberg weight decomposition by individual commodity (GPS 2020)}",
  "\\label{tab:rotemberg-commodity}",
  "\\begin{adjustbox}{max width=\\textwidth, max totalheight=0.92\\textheight}",
  "\\begin{threeparttable}\\small",
  "\\begin{tabular}{l rr@{\\quad}rr}",
  "\\toprule",
  " & \\multicolumn{2}{c}{$h=1$} & \\multicolumn{2}{c}{$h=4$} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  paste0("Commodity & $\\hat{\\alpha}_c$ & $\\hat{\\beta}_c^{\\text{IV}}$ ",
         "& $\\hat{\\alpha}_c$ & $\\hat{\\beta}_c^{\\text{IV}}$ \\\\"),
  "\\midrule"
)

current_grp <- ""
for (i in seq_len(nrow(res_c_tab))) {
  r <- res_c_tab[i, ]

  # Group header
  if (r$grp != current_grp) {
    if (current_grp != "") tab_c <- c(tab_c, "\\addlinespace")
    tab_c <- c(tab_c,
      sprintf("\\multicolumn{5}{l}{\\textit{%s}} \\\\", r$grp),
      "\\midrule"
    )
    current_grp <- r$grp
  }

  est_row <- sprintf(
    "\\quad %s & $%+.4f$ & $%+.3f$%s & $%+.4f$ & $%+.3f$%s \\\\",
    r$label,
    r$alpha_h1, r$est_h1, stars_fn(r$pval_h1),
    r$alpha,    r$est,    stars_fn(r$pval)
  )
  se_row <- sprintf(
    "   & & $(%.3f)$ & & $(%.3f)$ \\\\[2pt]",
    r$se_h1, r$se
  )
  tab_c <- c(tab_c, est_row, se_row)
}

tab_c <- c(tab_c,
  "\\midrule",
  sprintf(
    "\\textit{OLS baseline (weighted sum)} & $1.000$ & $%+.3f$ & $1.000$ & $%+.3f$ \\\\",
    b1c, b4c
  ),
  "\\midrule",
  sprintf("$N$ & \\multicolumn{2}{c}{%s} & \\multicolumn{2}{c}{%s} \\\\",
          format(n1c, big.mark = ","), format(n4c, big.mark = ",")),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0(
    "\\item \\textit{Notes:} Rotemberg weight decomposition following ",
    "\\citet{GoldsmithPinkham2020} at the individual commodity level. ",
    "The aggregate consumption shock $B_{it} = \\sum_c s_{ic,t}\\,\\Delta\\log p_{c,t}$ ",
    "is decomposed into 22 individual commodity components. ",
    "$\\hat{\\alpha}_c$ = Rotemberg weight; ",
    "$\\hat{\\beta}_c^{\\text{IV}}$ = just-identified IV coefficient using only ",
    "commodity $c$'s shock as the instrument for the aggregate. ",
    "Year FE; SE clustered by country in parentheses. ",
    "Commodities sorted by descending $\\hat{\\alpha}_c$ within each group at $h=4$. ",
    "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."
  ),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{adjustbox}",
  "\\end{table}"
)

write_tex(tab_c, "tab_rotemberg_commodity.tex")

# =============================================================================
# Section 5: Within-group Rotemberg decomposition
#
# For each GROUP-SPECIFIC LP model (regressor = shock_cons_{g}), decompose the
# group shock into individual commodity contributions.  This answers: within the
# cereals model, does identification come from rice + maize only, or are all
# cereals contributing?
#
#   α_c^g = Cov(B̃_c, B̃_g) / Var(B̃_g)
#
# where tildes denote residualisation on year FE + 2 lags of the GROUP shock.
# Weights sum to 1 within each group; HHI_g = Σ_c (α_c^g)^2, eff_N_g = 1/HHI_g.
# =============================================================================
message("\n=== Section 5: Within-group Rotemberg decomposition ===")

# Map commodity → group
COMM_TO_GRP <- c(
  wheat = "cereals", rice = "cereals", maize = "cereals",
  barley = "cereals", sorghum = "cereals",
  soybeans = "oils", soybean_oil = "oils", palm_oil = "oils",
  groundnut_oil = "oils", sunflower_oil = "oils",
  coconut_oil = "oils", rapeseed_oil = "oils",
  sugar = "sugar_fruit", bananas = "sugar_fruit", oranges = "sugar_fruit",
  beef = "animal", chicken = "animal", lamb = "animal", shrimp = "animal",
  coffee = "cash_crops", tea = "cash_crops", cocoa = "cash_crops"
)

run_within_group <- function(data, grp, h) {
  grp_col <- paste0("shock_cons_", grp)
  y_col   <- paste0("dfsi_h", h)

  # Commodities in this group
  comms_in_grp <- names(COMM_TO_GRP)[COMM_TO_GRP == grp]
  sc_cols_grp  <- paste0("sc_", comms_in_grp)
  sc_cols_grp  <- sc_cols_grp[sc_cols_grp %in% names(data)]
  if (length(sc_cols_grp) == 0) return(NULL)

  d <- data |>
    filter(!is.na(.data[[y_col]]),
           !is.na(.data[[grp_col]]),
           if_all(all_of(sc_cols_grp), ~ !is.na(.x)))

  if (nrow(d) < 20) return(NULL)

  # Compute 2 lags of the group shock within d
  d <- d |>
    arrange(iso3c, year) |>
    group_by(iso3c) |>
    mutate(
      grp_l1 = lag(.data[[grp_col]], 1),
      grp_l2 = lag(.data[[grp_col]], 2)
    ) |>
    ungroup() |>
    filter(!is.na(grp_l1), !is.na(grp_l2))

  if (nrow(d) < 20) return(NULL)

  # Residualise group shock on year FE + group lags
  fit_grp <- feols(as.formula(paste0(grp_col, " ~ grp_l1 + grp_l2 | year")),
                   data = d, notes = FALSE, warn = FALSE)
  r_grp <- residuals(fit_grp)

  results_within <- map_dfr(sc_cols_grp, function(sc) {
    comm <- sub("^sc_", "", sc)

    fit_c <- feols(as.formula(paste0(sc, " ~ grp_l1 + grp_l2 | year")),
                   data = d, notes = FALSE, warn = FALSE)
    r_c <- residuals(fit_c)

    alpha_c <- sum(r_c * r_grp, na.rm = TRUE) / sum(r_grp^2, na.rm = TRUE)

    fml_iv <- as.formula(paste0(
      y_col, " ~ grp_l1 + grp_l2 | year | ", grp_col, " ~ ", sc
    ))
    fit_iv <- tryCatch(
      feols(fml_iv, data = d, cluster = ~iso3c, notes = FALSE, warn = FALSE),
      error = function(e) NULL
    )

    if (is.null(fit_iv)) {
      return(tibble(grp = grp, commodity = comm, h = h,
                    alpha = alpha_c, est = NA_real_, se = NA_real_,
                    pval = NA_real_, n_obs = nrow(d)))
    }
    ct        <- summary(fit_iv)$coeftable
    coef_name <- grep("^fit_", rownames(ct), value = TRUE)[1]
    tibble(
      grp = grp, commodity = comm, h = h,
      alpha = alpha_c,
      est   = ct[coef_name, "Estimate"],
      se    = ct[coef_name, "Std. Error"],
      pval  = ct[coef_name, "Pr(>|t|)"],
      n_obs = nobs(fit_iv)
    )
  })
  results_within
}

wg_results <- map_dfr(HORIZONS_REPORT, function(h) {
  map_dfr(GROUPS, \(g) run_within_group(panel_c, g, h))
})

# ---- HHI and effective N per group ------------------------------------------
wg_hhi <- wg_results |>
  group_by(grp, h) |>
  summarise(
    hhi   = sum(alpha^2, na.rm = TRUE),
    eff_n = 1 / hhi,
    sum_alpha = sum(alpha, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nWithin-group HHI and effective number of shocks:\n")
cat(sprintf("%-14s  %4s  %6s  %6s  %8s\n", "Group", "h", "HHI", "eff_N", "Σ alpha"))
cat(strrep("-", 44), "\n")
for (i in seq_len(nrow(wg_hhi))) {
  r <- wg_hhi[i, ]
  cat(sprintf("%-14s  %4d  %6.3f  %6.2f  %8.3f\n",
              r$grp, r$h, r$hhi, r$eff_n, r$sum_alpha))
}

cat("\nWithin-group commodity weights (h=4):\n")
wg4 <- wg_results |>
  filter(h == 4) |>
  arrange(grp, desc(alpha))
for (g in GROUPS) {
  cat(sprintf("\n  %s:\n", g))
  wg4g <- filter(wg4, grp == g)
  for (i in seq_len(nrow(wg4g))) {
    r <- wg4g[i, ]
    cat(sprintf("    %-16s  alpha=%+.4f\n",
                COMMODITY_LABELS[r$commodity], r$alpha))
  }
}

# ---- LaTeX table: within-group HHI, eff N, and top commodity weights ---------
GRP_LABELS_CLEAN <- c(
  cereals     = "Cereals",
  oils        = "Oils",
  animal      = "Animal products",
  sugar_fruit = "Sugar \\& fruit",
  cash_crops  = "Beverage crops"
)

# Build one row per group with h=1 and h=4 HHI/eff_N + top commodity at h=4
wg_tab <- map_dfr(GROUPS, function(g) {
  h1 <- filter(wg_hhi, grp == g, h == 1)
  h4 <- filter(wg_hhi, grp == g, h == 4)
  top_c <- wg_results |>
    filter(grp == g, h == 4) |>
    slice_max(alpha, n = 1)
  tibble(
    grp        = g,
    hhi_h1     = if (nrow(h1) > 0) h1$hhi   else NA_real_,
    effn_h1    = if (nrow(h1) > 0) h1$eff_n  else NA_real_,
    hhi_h4     = if (nrow(h4) > 0) h4$hhi   else NA_real_,
    effn_h4    = if (nrow(h4) > 0) h4$eff_n  else NA_real_,
    top_comm   = if (nrow(top_c) > 0) COMMODITY_LABELS[top_c$commodity[1]] else NA_character_,
    top_alpha  = if (nrow(top_c) > 0) top_c$alpha[1] else NA_real_
  )
})

tab_wg <- c(
  "% Within-group Rotemberg weight decomposition",
  "",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Within-group Rotemberg weight concentration by commodity group}",
  "\\label{tab:rotemberg-within}",
  "\\begin{adjustbox}{max width=\\textwidth}",
  "\\begin{threeparttable}\\small",
  "\\begin{tabular}{l cc cc lc}",
  "\\toprule",
  " & \\multicolumn{2}{c}{$h = 1$} & \\multicolumn{2}{c}{$h = 4$} & \\multicolumn{2}{c}{Dominant commodity ($h=4$)} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
  "Group & HHI & Eff.\\ $N$ & HHI & Eff.\\ $N$ & Commodity & $\\hat{\\alpha}$ \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(wg_tab))) {
  r <- wg_tab[i, ]
  h1_str <- if (!is.na(r$hhi_h1))
    sprintf("%.3f & %.2f", r$hhi_h1, r$effn_h1) else "--- & ---"
  h4_str <- if (!is.na(r$hhi_h4))
    sprintf("%.3f & %.2f", r$hhi_h4, r$effn_h4) else "--- & ---"
  dom_str <- if (!is.na(r$top_comm))
    sprintf("%s & $%+.3f$", r$top_comm, r$top_alpha) else "--- & ---"

  tab_wg <- c(tab_wg,
    sprintf("%s & %s & %s & %s \\\\",
            GRP_LABELS_CLEAN[r$grp], h1_str, h4_str, dom_str))
}

# Add aggregate reference row for comparison
agg_hhi_h4 <- filter(results, h == 4) |>
  summarise(hhi = sum(alpha^2, na.rm = TRUE), eff_n = 1 / sum(alpha^2, na.rm = TRUE))

tab_wg <- c(tab_wg,
  "\\midrule",
  sprintf(
    "\\textit{Aggregate model} & \\multicolumn{2}{c}{---} & %.3f & %.2f & Cereals & ---  \\\\",
    agg_hhi_h4$hhi, agg_hhi_h4$eff_n
  ),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0(
    "\\item \\textit{Notes:} Within-group Rotemberg weight decomposition following \\citet{goldsmith2020bartik}. ",
    "For each group-specific local projection (outcome: $\\Delta$FSI; regressor: group consumption shock $B^g_{it}$), ",
    "the within-group weight $\\hat{\\alpha}_c^g = \\text{Cov}(\\tilde{B}_c, \\tilde{B}_g) / \\text{Var}(\\tilde{B}_g)$ ",
    "measures the contribution of commodity $c$ to identification within group $g$. ",
    "Residualisation ($\\tilde{\\cdot}$) removes year fixed effects and two lags of the group shock. ",
    "HHI $= \\sum_c (\\hat{\\alpha}_c^g)^2$; effective $N = 1/$HHI. ",
    "The bottom row reproduces the aggregate-model HHI for reference. ",
    "Effective $N > 1$ indicates identification is not fully concentrated in a single commodity."
  ),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{adjustbox}",
  "\\end{table}"
)

write_tex(tab_wg, "tab_rotemberg_within.tex")
message("\nDone.")

