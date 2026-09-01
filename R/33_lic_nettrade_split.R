# =============================================================================
# 33_lic_nettrade_split.R
# Robustness check: LP split by income group × net food trade position.
#
# Motivation: A reviewer asks whether the stabilising LIC result could merely
# reflect that LICs are net food exporters gaining export revenue when prices
# rise.  This table allows three separate discriminating tests:
#
#  (1) Net-importer LICs stabilise → rules out export-revenue channel for LICs
#  (2) Net-importer LICs stabilise (not destabilise) → rules out relative
#      deprivation as the LIC mechanism: relative deprivation would predict
#      net importers become MORE fragile when prices rise, since import costs
#      squeeze consumer budgets.  Finding stabilisation among net importers
#      is consistent only with the household farm-income / opportunity-cost
#      channel, which operates through producer income regardless of aggregate
#      trade balance.
#  (3) Pattern across income groups: if opportunity cost dominates in LICs
#      (net importers stabilise) but relative deprivation dominates in HICs
#      (net importers destabilise), that directly maps to the income-group
#      channel heterogeneity already documented in the main results.
#
# Classification:
#   Net food exporter = median aggregate net-import share < 0 over the first
#   five sample years (predetermined).  s_net_ma = (imports - exports)/supply,
#   summed over all commodities; negative = exports exceed imports in caloric
#   terms.
#
# Specification:
#   Consumption shock, year FE, 2 lags, clustered SE by country.
#   Horizons h=1 and h=4.
#
# Output:
#   output/tables/tab_nettrade_split.tex   — income group × net trade position
# =============================================================================

source(here::here("R", "00_setup.R"))

panel      <- readRDS(file.path(DATA_PRO, "panel.rds"))
shares_raw <- readRDS(file.path(DATA_PRO, "fao_shares.rds"))

message("Panel: ", nrow(panel), " obs | ",
        n_distinct(panel$iso3c), " countries | ",
        min(panel$year), "–", max(panel$year))

CI_Z <- 1.645
HORIZONS_REPORT <- c(1, 4)
INCOME_LEVELS   <- c("Low", "Lower-middle", "Upper-middle", "High")
INCOME_LABELS   <- c(
  "Low"           = "Low income",
  "Lower-middle"  = "Lower-middle income",
  "Upper-middle"  = "Upper-middle income",
  "High"          = "High income"
)

write_tex <- function(lines, filename) {
  path <- file.path(OUT_TABS, filename)
  writeLines(lines, path)
  message("  Saved: ", filename)
}

# =============================================================================
# Section 1: Classify ALL countries by net food trade position
# =============================================================================
message("\n=== Section 1: Net-trade classification (all countries) ===")

net_pos <- shares_raw |>
  filter(!is.na(s_net_ma),
         year >= SAMPLE_START, year <= SAMPLE_START + 4) |>
  group_by(iso3c, year) |>
  summarise(net_share = sum(s_net_ma, na.rm = TRUE), .groups = "drop") |>
  group_by(iso3c) |>
  summarise(net_share_baseline = median(net_share, na.rm = TRUE),
            .groups = "drop") |>
  mutate(net_exporter = net_share_baseline < 0)

panel <- panel |>
  left_join(net_pos |> select(iso3c, net_share_baseline, net_exporter),
            by = "iso3c")

cat("\nNet-trade classification by income group:\n")
panel |>
  filter(!duplicated(iso3c), !is.na(income_group)) |>
  count(income_group, net_exporter) |>
  pivot_wider(names_from = net_exporter, values_from = n,
              names_prefix = "net_exp_") |>
  rename(net_importers = net_exp_FALSE, net_exporters = net_exp_TRUE) |>
  print()

# =============================================================================
# Section 2: LP helper
# =============================================================================
run_lp <- function(data, shock_var = "shock_cons",
                   fd_stem = "dfsi", horizons = LP_HORIZONS) {
  lag1 <- paste0(shock_var, "_l1")
  lag2 <- paste0(shock_var, "_l2")
  map_dfr(horizons, function(h) {
    y_col <- paste0(fd_stem, "_h", h)
    if (!y_col %in% names(data)) return(NULL)
    d <- data |>
      rename(.y = all_of(y_col), .s = all_of(shock_var)) |>
      filter(!is.na(.y), !is.na(.s))
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
    tibble(
      h     = h,
      est   = ct[".s", "Estimate"],
      se    = ct[".s", "Std. Error"],
      pval  = ct[".s", "Pr(>|t|)"],
      ci_lo = est - CI_Z * se,
      ci_hi = est + CI_Z * se,
      n_obs = nobs(fit),
      n_cty = n_distinct(d$iso3c)
    )
  })
}

# =============================================================================
# Section 3: Run LP for each income group × net-trade position
# =============================================================================
message("\n=== Section 3: LP by income group × net-trade position ===")

results <- map_dfr(INCOME_LEVELS, function(ig) {
  d_all <- filter(panel, income_group == ig)
  d_exp <- filter(d_all, net_exporter == TRUE)
  d_imp <- filter(d_all, net_exporter == FALSE)

  groups <- list(
    all      = d_all,
    exporter = d_exp,
    importer = d_imp
  )
  map_dfr(names(groups), function(grp) {
    d <- groups[[grp]]
    if (n_distinct(d$iso3c) < 5) {
      message("  Skipping ", ig, " / ", grp,
              " (only ", n_distinct(d$iso3c), " countries)")
      return(tibble(h = HORIZONS_REPORT, est = NA_real_, se = NA_real_,
                    pval = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                    n_obs = NA_integer_, n_cty = n_distinct(d$iso3c),
                    income_group = ig, trade_pos = grp))
    }
    run_lp(d) |>
      filter(h %in% HORIZONS_REPORT) |>
      mutate(income_group = ig, trade_pos = grp)
  })
})

# Console output
cat("\nConsumption LP by income group × net trade position:\n")
cat(sprintf("%-18s  %-10s  %2s  %7s  %6s  %5s  %4s\n",
            "Income group", "Trade pos", "h", "est", "se", "p", "N_c"))
cat(strrep("-", 60), "\n")
for (ig in INCOME_LEVELS) {
  for (tp in c("all", "exporter", "importer")) {
    for (h in HORIZONS_REPORT) {
      r <- filter(results, income_group == ig, trade_pos == tp, h == !!h)
      if (nrow(r) == 0 || is.na(r$est)) {
        cat(sprintf("%-18s  %-10s  %2d  %7s\n", ig, tp, h, "(insuff N)"))
        next
      }
      stars <- ifelse(r$pval < 0.001, "***", ifelse(r$pval < 0.01, "**",
               ifelse(r$pval < 0.05, "*",   ifelse(r$pval < 0.10, ".", " "))))
      cat(sprintf("%-18s  %-10s  %2d  %+7.2f  %6.2f  %.3f%s  %3d\n",
                  ig, tp, h, r$est, r$se, r$pval, stars, r$n_cty))
    }
  }
  cat("\n")
}

# =============================================================================
# Section 4: LaTeX table
#
# Layout: one row-block per income group; within each block three rows:
#   (a) All countries, (b) Net exporters, (c) Net importers
# Columns: h=1 (est, SE) | h=4 (est, SE) | N_countries
# =============================================================================
message("\n=== Section 4: LaTeX table ===")

stars_math <- function(p) case_when(
  is.na(p)  ~ "",
  p < 0.001 ~ "^{***}",
  p < 0.01  ~ "^{**}",
  p < 0.05  ~ "^{*}",
  p < 0.10  ~ "^{\\cdot}",
  TRUE      ~ ""
)

fmt_cell <- function(r, h_val) {
  row <- filter(r, h == h_val)
  if (nrow(row) == 0 || is.na(row$est))
    return(list(est = "---", se = ""))
  list(
    est = sprintf("$%+.2f%s$", row$est, stars_math(row$pval)),
    se  = sprintf("$(%.2f)$", row$se)
  )
}

TRADE_LABELS <- c(
  all      = "\\quad All",
  exporter = "\\quad Net exporters",
  importer = "\\quad Net importers"
)

tab <- c(
  "% LP by income group x net food trade position: consumption shock",
  "",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Consumption shock LP by income group and net food trade position}",
  "\\label{tab:nettrade-split}",
  "\\begin{adjustbox}{max width=\\textwidth}",
  "\\begin{threeparttable}\\small",
  "\\begin{tabular}{l rr rr r}",
  "\\toprule",
  " & \\multicolumn{2}{c}{$h = 1$} & \\multicolumn{2}{c}{$h = 4$} & \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "Sample & Est. & SE & Est. & SE & Countries \\\\",
  "\\midrule"
)

for (ig in INCOME_LEVELS) {
  tab <- c(tab, sprintf("\\textit{%s} \\\\", INCOME_LABELS[ig]))

  for (tp in c("all", "exporter", "importer")) {
    r_tp <- filter(results, income_group == ig, trade_pos == tp)
    c1 <- fmt_cell(r_tp, 1)
    c4 <- fmt_cell(r_tp, 4)

    n_cty_val <- r_tp$n_cty[1]
    n_str <- if (!is.na(n_cty_val)) as.character(n_cty_val) else "---"

    tab <- c(tab,
      sprintf("%s & %s & %s & %s & %s & %s \\\\",
              TRADE_LABELS[tp], c1$est, c1$se, c4$est, c4$se, n_str),
      sprintf(" & & %s & & %s & \\\\", c1$se, c4$se)
    )
  }
  tab <- c(tab, "\\addlinespace")
}

# Compute baseline share stats for table note
n_exp_by_group <- panel |>
  filter(!duplicated(iso3c), !is.na(income_group), !is.na(net_exporter)) |>
  group_by(income_group) |>
  summarise(n_exp = sum(net_exporter, na.rm = TRUE),
            n_tot = n(), .groups = "drop")

note_str <- paste(
  map_chr(INCOME_LEVELS, function(ig) {
    r <- filter(n_exp_by_group, income_group == ig)
    if (nrow(r) == 0) return("")
    sprintf("%s: %d/%d net exporters", INCOME_LABELS[ig], r$n_exp, r$n_tot)
  }),
  collapse = "; "
)

tab <- c(tab,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0(
    "\\item \\textit{Notes:} LP: ",
    "$\\text{FSI}_{i,t+h} - \\text{FSI}_{i,t-1} = \\beta\\,B^{\\text{cons}}_{it} + ",
    "\\sum_{\\ell=1}^{2}\\gamma_\\ell B^{\\text{cons}}_{i,t-\\ell} + \\alpha_t + \\varepsilon_{iht}$. ",
    "Year FE; SE clustered by country. ",
    "\\textit{Net exporter} = median aggregate net-import share ",
    "($s^{\\text{net}} = (\\text{imports} - \\text{exports})/\\text{supply}$, summed over all commodities) ",
    "negative over the first five sample years (predetermined). ",
    "Country counts by group: ", note_str, ". ",
    "Key discriminating predictions: the \\textit{relative deprivation} channel predicts net importers ",
    "become \\textit{more} fragile when prices rise (import costs compress consumer budgets); ",
    "the \\textit{opportunity cost} channel predicts net importers still \\textit{stabilise} ",
    "(smallholder producer income rises regardless of aggregate trade balance). ",
    "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."
  ),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{adjustbox}",
  "\\end{table}"
)

write_tex(tab, "tab_nettrade_split.tex")
message("\nDone.")
