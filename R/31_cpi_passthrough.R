# =============================================================================
# 31_cpi_passthrough.R
# Food price pass-through: does the Bartik shock reach domestic food CPI,
# and does pass-through differ systematically across income groups?
#
# Outcome: annual log-change in FAO food CPI (country-level)
# Treatment: shift-share consumption shock (aggregate + 5 commodity groups)
#
# Specification:
#   Δlog food_CPI_{it} = α_t + β shock_{it} + ε_{it}
#   Year FE; SE clustered by country.
#
# A null result on income-group interactions validates the Bartik shock as
# a relevant first stage across all country subsamples.
#
# Output:
#   output/tables/tab_cpi_passthrough.tex   — aggregate + group, by income group
# =============================================================================

source(here::here("R", "00_setup.R"))

# ---- Load FAO CPI -----------------------------------------------------------
CPI_RAW <- file.path(DATA_RAW, "fao", "ConsumerPriceIndices_E_All_Data_(Normalized).csv")
message("Reading: ", CPI_RAW)
cpi_raw <- read_csv(CPI_RAW, show_col_types = FALSE)

# Keep food CPI and general CPI; annual average from monthly values
cpi_annual <- cpi_raw |>
  filter(Item %in% c("Consumer Prices, Food Indices (2015 = 100)",
                     "Consumer Prices, General Indices (2015 = 100)"),
         !is.na(Value)) |>
  mutate(
    cpi_type = if_else(grepl("Food", Item), "food_cpi", "gen_cpi"),
    year     = as.integer(Year),
    Area     = as.character(Area)
  ) |>
  filter(year >= SAMPLE_START - 1, year <= SAMPLE_END) |>
  group_by(Area, year, cpi_type) |>
  summarise(cpi = mean(Value, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = cpi_type, values_from = cpi)

message("CPI annual: ", nrow(cpi_annual), " obs | ",
        n_distinct(cpi_annual$Area), " areas | ",
        min(cpi_annual$year), "–", max(cpi_annual$year))

# ---- Country name → iso3c ---------------------------------------------------
cpi_annual <- cpi_annual |>
  mutate(
    iso3c = countrycode(Area, "country.name", "iso3c",
                        warn = FALSE, nomatch = NA_character_)
  ) |>
  filter(!is.na(iso3c), nchar(iso3c) == 3)

message("After iso3c match: ", n_distinct(cpi_annual$iso3c), " countries")

# ---- Compute log changes ----------------------------------------------------
cpi_lc <- cpi_annual |>
  arrange(iso3c, year) |>
  group_by(iso3c) |>
  mutate(
    dlog_food_cpi = log(food_cpi) - log(lag(food_cpi, 1)),
    dlog_gen_cpi  = log(gen_cpi)  - log(lag(gen_cpi,  1))
  ) |>
  ungroup() |>
  filter(year >= SAMPLE_START, year <= SAMPLE_END)

message("Log-change CPI: ", nrow(cpi_lc), " obs | ",
        sum(!is.na(cpi_lc$dlog_food_cpi)), " non-NA food CPI changes")

# ---- Merge with panel -------------------------------------------------------
panel <- readRDS(file.path(DATA_PRO, "panel.rds"))

panel_cpi <- panel |>
  left_join(cpi_lc |> select(iso3c, year, dlog_food_cpi, dlog_gen_cpi),
            by = c("iso3c", "year"))

message("Panel + CPI: ", nrow(panel_cpi), " obs | ",
        sum(!is.na(panel_cpi$dlog_food_cpi)), " non-NA food CPI obs")

# ---- Constants --------------------------------------------------------------
INCOME_LEVELS <- c("Low", "Lower-middle", "Upper-middle", "High")
INCOME_LABELS <- c("Low", "Lower-middle", "Upper-middle", "High income")

GROUPS      <- c("cereals", "oils", "animal", "sugar_fruit", "cash_crops")
GRP_LABELS  <- c(
  cereals     = "Cereals",
  oils        = "Oils",
  animal      = "Animal products",
  sugar_fruit = "Sugar \\& fruit",
  cash_crops  = "Beverage crops"
)

stars_math <- function(p)
  case_when(is.na(p) ~ "", p < 0.001 ~ "^{***}", p < 0.01 ~ "^{**}",
            p < 0.05 ~ "^{*}", p < 0.10 ~ "^{\\cdot}", TRUE ~ "")

write_tex <- function(lines, filename) {
  path <- file.path(OUT_TABS, filename)
  writeLines(lines, path)
  message("  Saved: ", filename)
}

# ---- Helper: run pass-through regression ------------------------------------
# Regress dlog_food_cpi on shock_var, year FE, clustered by country
run_pt <- function(data, shock_var) {
  d <- data |>
    filter(!is.na(dlog_food_cpi), !is.na(.data[[shock_var]]))
  if (nrow(d) < 20 || n_distinct(d$iso3c) < 5) return(NULL)

  fml <- as.formula(paste0("dlog_food_cpi ~ ", shock_var, " | year"))
  fit <- tryCatch(
    feols(fml, data = d, cluster = ~iso3c, notes = FALSE, warn = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  ct <- summary(fit)$coeftable
  if (!shock_var %in% rownames(ct)) return(NULL)

  tibble(
    est   = ct[shock_var, "Estimate"],
    se    = ct[shock_var, "Std. Error"],
    pval  = ct[shock_var, "Pr(>|t|)"],
    n_obs = nobs(fit),
    n_cty = n_distinct(d$iso3c)
  )
}

# =============================================================================
# Section 1: Pass-through by income group — aggregate and commodity groups
# =============================================================================
message("\n=== Section 1: Pass-through regressions ===")

SHOCKS <- c("shock_cons", paste0("shock_cons_", GROUPS))
SHOCK_LABELS <- c(
  shock_cons             = "\\textit{Aggregate}",
  shock_cons_cereals     = "\\quad Cereals",
  shock_cons_oils        = "\\quad Oils",
  shock_cons_animal      = "\\quad Animal products",
  shock_cons_sugar_fruit = "\\quad Sugar \\& fruit",
  shock_cons_cash_crops  = "\\quad Beverage crops"
)

results_pt <- map_dfr(SHOCKS, function(sv) {
  # All countries
  r_all <- run_pt(panel_cpi, sv)
  # By income group
  r_inc <- map_dfr(INCOME_LEVELS, function(ig) {
    d_ig <- filter(panel_cpi, income_group == ig)
    run_pt(d_ig, sv) |> mutate(income_group = ig)
  })
  bind_rows(
    if (!is.null(r_all)) mutate(r_all, income_group = "All"),
    r_inc
  ) |> mutate(shock = sv)
})

# Console output
cat("\nPass-through coefficients (dlog food CPI ~ shock_cons | year FE):\n")
cat(sprintf("%-28s  %5s  %8s  %6s\n", "Shock", "Group", "beta", "p"))
cat(strrep("-", 52), "\n")

for (sv in SHOCKS) {
  cat(SHOCK_LABELS[sv], "\n")
  for (ig in c("All", INCOME_LEVELS)) {
    r <- filter(results_pt, shock == sv, income_group == ig)
    if (nrow(r) == 0) next
    stars <- ifelse(is.na(r$pval), "",
             ifelse(r$pval < 0.001, "***", ifelse(r$pval < 0.01, "**",
             ifelse(r$pval < 0.05, "*", ifelse(r$pval < 0.10, ".", "")))))
    cat(sprintf("  %-24s  %+8.4f  %.3f %s  (N=%d)\n",
                ig, r$est, r$pval, stars, r$n_obs))
  }
  cat("\n")
}

# ---- Test: F-test for equality of pass-through across income groups ---------
message("\n=== Section 2: Equality test across income groups ===")

pt_tests <- map_dfr(SHOCKS, function(sv) {
  d <- panel_cpi |>
    filter(!is.na(dlog_food_cpi), !is.na(.data[[sv]]),
           income_group %in% INCOME_LEVELS) |>
    mutate(income_f = factor(income_group, levels = INCOME_LEVELS))

  # Interact shock with income group, Low income as reference
  # H0: deviations from Low income group are jointly zero = equal pass-through
  fml_int <- as.formula(paste0(
    "dlog_food_cpi ~ i(income_f, ", sv, ", ref = 'Low') | year"
  ))
  fit_int <- tryCatch(
    feols(fml_int, data = d, cluster = ~iso3c, notes = FALSE, warn = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit_int)) return(tibble(shock = sv, fstat = NA, pval_f = NA))

  # Wald test on interaction terms only (deviations from Low income)
  # If H0 not rejected: pass-through is equal across groups
  wt <- tryCatch(
    wald(fit_int, "income_f::"),
    error = function(e) NULL
  )
  if (is.null(wt))
    return(tibble(shock = sv, fstat = NA, pval_f = NA))

  tibble(shock = sv, fstat = wt$stat, pval_f = wt$p)
})

cat("\nWald test — H0: pass-through equal across income groups:\n")
print(pt_tests)

# =============================================================================
# Section 3: LaTeX table
# Rows: aggregate + 5 commodity groups
# Columns: All | Low | Lower-middle | Upper-middle | High
# =============================================================================
message("\n=== Section 3: LaTeX table ===")

GROUPS_DISPLAY <- c("All", INCOME_LEVELS)
GROUPS_HEADER  <- c("All", "Low", "Lower-middle", "Upper-middle", "High income")

tab <- c(
  "% CPI pass-through by income group",
  "",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Food price pass-through to domestic food CPI by income group}",
  "\\label{tab:cpi-passthrough}",
  "\\begin{adjustbox}{max width=\\textwidth}",
  "\\begin{threeparttable}\\small",
  "\\begin{tabular}{l rrrrr}",
  "\\toprule",
  paste0(" & ", paste(sprintf("\\multicolumn{1}{c}{%s}", GROUPS_HEADER),
                      collapse = " & "), " \\\\"),
  "\\midrule"
)

prev_block <- ""
for (sv in SHOCKS) {
  lbl <- SHOCK_LABELS[sv]

  # Add midrule before group entries
  is_group <- grepl("^\\\\quad", lbl)
  if (!is_group && sv != "shock_cons") tab <- c(tab, "\\addlinespace")

  est_row <- paste(
    map_chr(GROUPS_DISPLAY, function(ig) {
      r <- filter(results_pt, shock == sv, income_group == ig)
      if (nrow(r) == 0 || is.na(r$est)) return("---")
      sprintf("$%+.4f%s$", r$est, stars_math(r$pval))
    }),
    collapse = " & "
  )
  se_row <- paste(
    map_chr(GROUPS_DISPLAY, function(ig) {
      r <- filter(results_pt, shock == sv, income_group == ig)
      if (nrow(r) == 0 || is.na(r$se)) return("")
      sprintf("$(%.4f)$", r$se)
    }),
    collapse = " & "
  )
  n_row <- paste(
    map_chr(GROUPS_DISPLAY, function(ig) {
      r <- filter(results_pt, shock == sv, income_group == ig)
      if (nrow(r) == 0) return("")
      format(r$n_obs, big.mark = ",")
    }),
    collapse = " & "
  )

  tab <- c(tab,
    sprintf("%s & %s \\\\", lbl, est_row),
    sprintf("   & %s \\\\", se_row),
    sprintf("   & %s \\\\[3pt]", n_row)
  )
}

# Wald test row for aggregate
wt_agg <- filter(pt_tests, shock == "shock_cons")
tab <- c(tab,
  "\\midrule",
  sprintf(
    "\\multicolumn{6}{l}{\\textit{Wald test — equal pass-through across income groups: $F = %.2f$, $p = %.3f$}} \\\\",
    wt_agg$fstat, wt_agg$pval_f
  ),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0(
    "\\item \\textit{Notes:} OLS regression of annual log-change in country food CPI ",
    "(FAO, base 2015$=$100) on the shift-share consumption shock, estimated with year ",
    "fixed effects and standard errors clustered by country. ",
    "Commodity-group rows decompose the aggregate shock into five components; ",
    "coefficients measure pass-through from each group's global price change to the ",
    "aggregate domestic food price index. ",
    "The Wald test (bottom row) tests whether the aggregate pass-through coefficient ",
    "is equal across income groups; a $p$-value above 0.10 supports the null of ",
    "homogeneous pass-through. ",
    "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."
  ),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{adjustbox}",
  "\\end{table}"
)

write_tex(tab, "tab_cpi_passthrough.tex")
message("\nDone.")
