# =============================================================================
# 29_gov_revenue_lp.R
# Fiscal channel robustness — LP on WDI fiscal indicators.
#
# Reviewer concern: food price shocks may affect FSI through a fiscal channel
# (shock → revenues/expenditure → public services → fragility).
# Tests the reduced-form relationship between each Bartik shock and five fiscal
# outcomes. A null result rules out a large aggregate fiscal response.
#
# Fiscal outcomes (all from WDI, sourced in 06_other_variables.R):
#   gov_rev   — Revenue excl. grants, % GDP              (GC.REV.XGRT.GD.ZS)
#   trade_tax — Taxes on international trade, % revenue  (GC.TAX.INTT.RV.ZS)
#   gs_tax    — Goods & services tax revenue, % revenue  (GC.TAX.GSRV.RV.ZS)
#   gov_exp   — General govt expenditure, % GDP          (GC.XPN.TOTL.GD.ZS)
#   (GC.BAL.CASH.GD.ZS and DT.DOD.DECT.GD.ZS not available via WDI API)
#
# NOTE: Run 06_other_variables.R (with download block un-commented) first if
# the new fiscal indicators are not yet in wdi_covariates.csv.
#
# Output:
#   output/tables/tab_fiscal_lp.tex   — all outcomes × cons + prod, h=1 and h=4
# =============================================================================

source(here::here("R", "00_setup.R"))

panel <- readRDS(file.path(DATA_PRO, "panel.rds"))
covariates <- readRDS(file.path(DATA_PRO, "covariates.rds"))

message("Panel: ", nrow(panel), " obs | ",
        n_distinct(panel$iso3c), " countries")

# Fiscal indicators to test
FISCAL_VARS <- c("gov_rev", "gov_exp")

FISCAL_LABELS <- c(
  gov_rev = "Revenue excl. grants (\\% GDP)",
  gov_exp = "Govt expenditure (\\% GDP)"
)

# Keep only fiscal vars that are present in covariates
FISCAL_VARS <- intersect(FISCAL_VARS, names(covariates))
if (length(FISCAL_VARS) == 0)
  stop("No fiscal variables found in covariates.rds — re-run 06_other_variables.R with new indicators.")

cat("\nFiscal variable coverage:\n")
covariates |>
  select(any_of(FISCAL_VARS)) |>
  summarise(across(everything(), \(x) sum(!is.na(x)))) |>
  print()

# ---- Merge fiscal vars into panel -------------------------------------------
# Drop any FISCAL_VARS already present in panel to avoid join suffixes
new_fiscal <- setdiff(FISCAL_VARS, names(panel))
panel_fis <- panel |>
  left_join(covariates |> select(iso3c, year, any_of(new_fiscal)),
            by = c("iso3c", "year"))

# ---- Build LP horizon outcomes for each fiscal variable ---------------------
panel_fis <- panel_fis |>
  arrange(iso3c, year) |>
  group_by(iso3c) |>
  mutate(across(
    all_of(FISCAL_VARS),
    list(
      lag1 = \(x) lag(x, 1),
      h0   = \(x) x               - lag(x, 1),
      h1   = \(x) lead(x, 1)      - lag(x, 1),
      h2   = \(x) lead(x, 2)      - lag(x, 1),
      h3   = \(x) lead(x, 3)      - lag(x, 1),
      h4   = \(x) lead(x, 4)      - lag(x, 1)
    ),
    .names = "{.col}_{.fn}"
  )) |>
  ungroup()

cat("\nHorizon outcome coverage (h=0 / h=4):\n")
for (v in FISCAL_VARS) {
  n0 <- sum(!is.na(panel_fis[[paste0(v, "_h0")]]))
  n4 <- sum(!is.na(panel_fis[[paste0(v, "_h4")]]))
  cat(sprintf("  %-12s  h=0: %d  |  h=4: %d\n", v, n0, n4))
}

# ---- LP runner --------------------------------------------------------------
CI_Z <- 1.645

stars_fn <- function(p)
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", p < 0.10 ~ "$\\cdot$", TRUE ~ "")

run_lp <- function(data, shock_var, outcome_stem, horizons = LP_HORIZONS) {
  lag1 <- paste0(shock_var, "_l1")
  lag2 <- paste0(shock_var, "_l2")
  map_dfr(horizons, function(h) {
    y_col <- paste0(outcome_stem, "_h", h)
    if (!y_col %in% names(data) || !shock_var %in% names(data)) return(NULL)
    d <- data |>
      rename(.y = all_of(y_col), .s = all_of(shock_var)) |>
      filter(!is.na(.y), !is.na(.s))
    if (nrow(d) < 30 || n_distinct(d$iso3c) < 5) return(NULL)
    if (lag1 %in% names(d)) d <- rename(d, .s_l1 = all_of(lag1))
    if (lag2 %in% names(d)) d <- rename(d, .s_l2 = all_of(lag2))
    rhs_lags <- intersect(c(".s_l1", ".s_l2"), names(d))
    fml <- as.formula(paste0(
      ".y ~ .s",
      if (length(rhs_lags)) paste0(" + ", paste(rhs_lags, collapse = " + ")),
      " | year"
    ))
    fit <- tryCatch(
      feols(fml, data = d, cluster = ~iso3c, warn = FALSE, notes = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    ct <- summary(fit)$coeftable
    if (!".s" %in% rownames(ct)) return(NULL)
    tibble(h     = h,
           est   = ct[".s", "Estimate"],
           se    = ct[".s", "Std. Error"],
           pval  = ct[".s", "Pr(>|t|)"],
           n_obs = nobs(fit))
  })
}

# ---- Estimate: cons + prod × all fiscal outcomes ----------------------------
TREAT_TYPES  <- c("cons", "prod")
TREAT_LABELS <- c(cons = "Consumption", prod = "Production")

message("\n=== Fiscal LP ===")
irf_fis <- map_dfr(FISCAL_VARS, function(fv) {
  map_dfr(TREAT_TYPES, function(tt) {
    sv <- paste0("shock_", tt)
    run_lp(panel_fis, sv, fv) |>
      mutate(fiscal_var = fv, treat = tt)
  })
})

# Console summary
cat("\n=== Results: h=1 and h=4 ===\n")
cat(sprintf("%-40s  %-12s  %8s  %6s  %8s  %6s\n",
            "Outcome", "Treat", "h=1 est", "p", "h=4 est", "p"))
cat(strrep("-", 84), "\n")
for (fv in FISCAL_VARS) {
  for (tt in TREAT_TYPES) {
    r1 <- filter(irf_fis, fiscal_var == fv, treat == tt, h == 1)
    r4 <- filter(irf_fis, fiscal_var == fv, treat == tt, h == 4)
    cat(sprintf("%-40s  %-12s  %+8.3f  %6.3f  %+8.3f  %6.3f%s\n",
                gsub("\\\\", "", FISCAL_LABELS[fv]), TREAT_LABELS[tt],
                if (nrow(r1)) r1$est else NA, if (nrow(r1)) r1$pval else NA,
                if (nrow(r4)) r4$est else NA, if (nrow(r4)) r4$pval else NA,
                if (nrow(r4)) stars_fn(r4$pval) else ""))
  }
  cat("\n")
}

# =============================================================================
# LaTeX table: rows = fiscal outcomes, columns = cons (h=1, h=4) + prod (h=1, h=4)
# =============================================================================
fmt_cell <- function(fv, tt, h_) {
  r <- filter(irf_fis, fiscal_var == fv, treat == tt, h == h_)
  if (nrow(r) == 0 || is.na(r$est))
    return(list(e = "---", s = ""))
  list(e = sprintf("$%+.3f%s$", r$est, stars_fn(r$pval)),
       s = sprintf("$(%.3f)$", r$se))
}

tab <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Fiscal channel robustness: LP on WDI fiscal indicators}",
  "\\label{tab:fiscal-lp}",
  "\\begin{adjustbox}{max width=\\textwidth}",
  "\\begin{threeparttable}\\small",
  "\\begin{tabular}{lrrrr}",
  "\\toprule",
  " & \\multicolumn{2}{c}{Consumption} & \\multicolumn{2}{c}{Production} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "Fiscal outcome & $h=1$ & $h=4$ & $h=1$ & $h=4$ \\\\",
  "\\midrule"
)

for (fv in FISCAL_VARS) {
  n_obs_cons <- filter(irf_fis, fiscal_var == fv, treat == "cons", h == 0)$n_obs
  n_obs_prod <- filter(irf_fis, fiscal_var == fv, treat == "prod", h == 0)$n_obs
  n_str <- if (length(n_obs_cons)) format(n_obs_cons[1], big.mark = ",") else "---"

  tab <- c(tab,
    paste0(FISCAL_LABELS[fv], " & ",
           fmt_cell(fv, "cons", 1)$e, " & ", fmt_cell(fv, "cons", 4)$e, " & ",
           fmt_cell(fv, "prod", 1)$e, " & ", fmt_cell(fv, "prod", 4)$e, " \\\\"),
    paste0("& ",
           fmt_cell(fv, "cons", 1)$s, " & ", fmt_cell(fv, "cons", 4)$s, " & ",
           fmt_cell(fv, "prod", 1)$s, " & ", fmt_cell(fv, "prod", 4)$s, " \\\\"),
    paste0("& \\multicolumn{4}{l}{\\footnotesize $N = ", n_str, "$} \\\\[4pt]")
  )
}

tab <- c(tab,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0("\\item \\textit{Notes:} Each row runs the baseline LP ",
         "($\\Delta_h Y_{it} = Y_{i,t+h} - Y_{i,t-1}$) with a WDI fiscal ",
         "indicator as the outcome in place of FSI. Specification: year FE, ",
         "two lags of the shock, SE clustered by country. ",
         "A null result across outcomes indicates that food price shocks do not ",
         "operate through an aggregate fiscal channel. ",
         "Data: World Bank WDI. ",
         "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{adjustbox}",
  "\\end{table}"
)

writeLines(tab, file.path(OUT_TABS, "tab_fiscal_lp.tex"))
message("\nSaved: tab_fiscal_lp.tex")

# =============================================================================
# By income group: cons + prod × all fiscal outcomes, h=1 and h=4
# One table per fiscal variable, rows = income groups, cols = cons/prod × h
# =============================================================================
message("\n=== Fiscal LP by income group ===")

INCOME_LEVELS <- c("Low", "Lower-middle", "Upper-middle", "High")
INCOME_LABELS_SHORT <- c(
  Low            = "Low income",
  `Lower-middle` = "Lower-middle",
  `Upper-middle` = "Upper-middle",
  High           = "High income"
)

panel_by_inc <- setNames(
  map(INCOME_LEVELS, ~ filter(panel_fis, income_group == .x)),
  INCOME_LEVELS
)

irf_fis_inc <- map_dfr(INCOME_LEVELS, function(ig) {
  d <- panel_by_inc[[ig]]
  message("  ", ig, " (", n_distinct(d$iso3c), " countries)")
  map_dfr(FISCAL_VARS, function(fv) {
    map_dfr(TREAT_TYPES, function(tt) {
      sv <- paste0("shock_", tt)
      run_lp(d, sv, fv) |>
        mutate(fiscal_var = fv, treat = tt, income_group = ig)
    })
  })
})

# Console summary — gov_rev only for brevity
cat("\n=== By income group: gov_rev, cons (h=1 and h=4) ===\n")
cat(sprintf("%-16s  %+8s  %6s  %+8s  %6s\n", "Income", "h=1", "p", "h=4", "p"))
cat(strrep("-", 50), "\n")
for (ig in INCOME_LEVELS) {
  for (tt in TREAT_TYPES) {
    r1 <- filter(irf_fis_inc, income_group == ig, fiscal_var == "gov_rev",
                 treat == tt, h == 1)
    r4 <- filter(irf_fis_inc, income_group == ig, fiscal_var == "gov_rev",
                 treat == tt, h == 4)
    cat(sprintf("%-16s  %-12s  %+8.3f  %6.3f  %+8.3f  %6.3f%s\n",
                ig, TREAT_LABELS[tt],
                if (nrow(r1)) r1$est else NA, if (nrow(r1)) r1$pval else NA,
                if (nrow(r4)) r4$est else NA, if (nrow(r4)) r4$pval else NA,
                if (nrow(r4)) stars_fn(r4$pval) else ""))
  }
  cat("\n")
}

# LaTeX: one table per fiscal variable
# Rows: income groups; cols: cons (h=1, h=4) + prod (h=1, h=4)
fmt_inc <- function(fv, ig, tt, h_) {
  r <- filter(irf_fis_inc, fiscal_var == fv, income_group == ig,
              treat == tt, h == h_)
  if (nrow(r) == 0 || is.na(r$est)) return(list(e = "---", s = ""))
  list(e = sprintf("$%+.3f%s$", r$est, stars_fn(r$pval)),
       s = sprintf("$(%.3f)$", r$se))
}

for (fv in FISCAL_VARS) {
  tab_inc <- c(
    sprintf("\\begin{table}[htbp]"),
    "\\centering",
    sprintf("\\caption{Fiscal LP by income group: %s}", gsub("\\\\", "", FISCAL_LABELS[fv])),
    sprintf("\\label{tab:fiscal-inc-%s}", fv),
    "\\begin{adjustbox}{max width=\\textwidth}",
    "\\begin{threeparttable}\\small",
    "\\begin{tabular}{lrrrr}",
    "\\toprule",
    " & \\multicolumn{2}{c}{Consumption} & \\multicolumn{2}{c}{Production} \\\\",
    "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
    "Income group & $h=1$ & $h=4$ & $h=1$ & $h=4$ \\\\",
    "\\midrule"
  )
  for (ig in INCOME_LEVELS) {
    n_ <- filter(irf_fis_inc, fiscal_var == fv, income_group == ig,
                 treat == "cons", h == 0)$n_obs
    n_str <- if (length(n_)) format(n_[1], big.mark = ",") else "---"
    tab_inc <- c(tab_inc,
      paste0(INCOME_LABELS_SHORT[ig], " & ",
             fmt_inc(fv, ig, "cons", 1)$e, " & ", fmt_inc(fv, ig, "cons", 4)$e, " & ",
             fmt_inc(fv, ig, "prod", 1)$e, " & ", fmt_inc(fv, ig, "prod", 4)$e, " \\\\"),
      paste0("& ",
             fmt_inc(fv, ig, "cons", 1)$s, " & ", fmt_inc(fv, ig, "cons", 4)$s, " & ",
             fmt_inc(fv, ig, "prod", 1)$s, " & ", fmt_inc(fv, ig, "prod", 4)$s, " \\\\"),
      paste0("& \\multicolumn{4}{l}{\\footnotesize $N = ", n_str, "$} \\\\[4pt]")
    )
  }
  tab_inc <- c(tab_inc,
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{tablenotes}\\small",
    paste0("\\item \\textit{Notes:} LP with ", gsub("\\\\", "", FISCAL_LABELS[fv]),
           " as the outcome ($\\Delta_h Y_{it} = Y_{i,t+h} - Y_{i,t-1}$). ",
           "Year FE; SE clustered by country. Data: World Bank WDI. ",
           "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."),
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\end{adjustbox}",
    "\\end{table}"
  )
  fname <- sprintf("tab_fiscal_inc_%s.tex", fv)
  writeLines(tab_inc, file.path(OUT_TABS, fname))
  message("  Saved: ", fname)
}

# =============================================================================
# Combined table: all fiscal outcomes × all income groups in one table
# Panels: one per income group; rows = fiscal outcomes; cols = cons/prod × h=1,h=4
# =============================================================================
tab_all <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Fiscal channel LP by income group: all outcomes}",
  "\\label{tab:fiscal-inc-all}",
  "\\begin{adjustbox}{max width=\\textwidth}",
  "\\begin{threeparttable}\\small",
  "\\begin{tabular}{lrrrr}",
  "\\toprule",
  " & \\multicolumn{2}{c}{Consumption} & \\multicolumn{2}{c}{Production} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "Fiscal outcome & $h=1$ & $h=4$ & $h=1$ & $h=4$ \\\\",
  "\\midrule"
)

for (ig in INCOME_LEVELS) {
  tab_all <- c(tab_all,
    sprintf("\\multicolumn{5}{l}{\\textit{%s}} \\\\", INCOME_LABELS_SHORT[ig]),
    "\\midrule"
  )
  for (fv in FISCAL_VARS) {
    n_ <- filter(irf_fis_inc, fiscal_var == fv, income_group == ig,
                 treat == "cons", h == 0)$n_obs
    n_str <- if (length(n_)) format(n_[1], big.mark = ",") else "---"
    tab_all <- c(tab_all,
      paste0("\\quad ", FISCAL_LABELS[fv], " & ",
             fmt_inc(fv, ig, "cons", 1)$e, " & ", fmt_inc(fv, ig, "cons", 4)$e, " & ",
             fmt_inc(fv, ig, "prod", 1)$e, " & ", fmt_inc(fv, ig, "prod", 4)$e, " \\\\"),
      paste0("& ",
             fmt_inc(fv, ig, "cons", 1)$s, " & ", fmt_inc(fv, ig, "cons", 4)$s, " & ",
             fmt_inc(fv, ig, "prod", 1)$s, " & ", fmt_inc(fv, ig, "prod", 4)$s,
             " \\\\[2pt]")
    )
  }
  tab_all <- c(tab_all, "\\addlinespace")
}

tab_all <- c(tab_all,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0("\\item \\textit{Notes:} Each cell reports the LP coefficient ",
         "($\\Delta_h Y_{it} = Y_{i,t+h} - Y_{i,t-1}$) for the indicated ",
         "fiscal outcome and income group. Specification: year FE, two lags ",
         "of the shock, SE clustered by country. Data: World Bank WDI. ",
         "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{adjustbox}",
  "\\end{table}"
)

writeLines(tab_all, file.path(OUT_TABS, "tab_fiscal_inc_all.tex"))
message("  Saved: tab_fiscal_inc_all.tex")
message("Done.")
