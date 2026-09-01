# =============================================================================
# 10_income_lp.R
# LP models by World Bank income group
#
# Specification: same as 08_lp.R (year FE, 2 lags, cluster by country)
# Sample splits: Low / Lower-middle / Upper-middle / High income
#                plus combined LMIC (Low + Lower-middle)
#
# Outputs:
#   output/tables/tab_income_lp.tex   — all exposures × income groups, h=1 & h=4
#   output/figures/fig3_income_irf.pdf/.png
#     — cons + prod exposures, facet_grid(treat ~ income_group), 180×100 mm
# =============================================================================

source(here::here("R", "00_setup.R"))

panel <- readRDS(file.path(DATA_PRO, "panel.rds"))

message("Panel: ", nrow(panel), " obs | ",
        n_distinct(panel$iso3c), " countries | ",
        min(panel$year), "–", max(panel$year))

# ---- Constants --------------------------------------------------------------
TREAT_TYPES  <- c("cons", "prod", "imp", "exp", "net")
TREAT_LABELS <- c(cons = "Consumption",
                  prod = "Production",
                  imp  = "Imports",
                  exp  = "Exports",
                  net  = "Net imports")

INCOME_LEVELS <- c("Low", "Lower-middle", "Upper-middle", "High")
INCOME_LABELS <- c("Low\nincome", "Lower-middle\nincome",
                   "Upper-middle\nincome", "High\nincome")

CI_Z <- 1.645

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

stars_console <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
  ifelse(p < 0.05,  "*",   ifelse(p < 0.10, ".",  ""))))
}

write_tex <- function(lines, filename) {
  path <- file.path(OUT_TABS, filename)
  writeLines(lines, path)
  message("  Saved: ", filename)
}

run_lp <- function(data, shock_var, fd_stem = "dfsi", horizons = LP_HORIZONS) {
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
      n_obs = nobs(fit)
    )
  })
}

# =============================================================================
# Section 1: Run LP for each income group × exposure type
# =============================================================================
message("\n=== Section 1: Income group LP ===")

# Check income group column
if (!"income_group" %in% names(panel)) {
  stop("income_group not in panel — run 07_panel.R with income_groups.rds present.")
}

cat("\nIncome group distribution:\n")
print(table(panel$income_group[!duplicated(panel$iso3c)], useNA = "ifany"))

# Build LMIC subset (Low + Lower-middle combined)
panel_lmic <- panel |>
  filter(income_group %in% c("Low", "Lower-middle")) |>
  mutate(income_group = "LMIC")

# All subsets: named list of data frames
panel_by_group <- c(
  setNames(
    map(INCOME_LEVELS, ~ filter(panel, income_group == .x)),
    INCOME_LEVELS
  ),
  list(LMIC = panel_lmic)
)

# Report country counts per group
walk2(names(panel_by_group), panel_by_group, function(grp, d) {
  n_c <- n_distinct(d$iso3c)
  message("  ", grp, ": ", n_c, " countries")
})

# Run LP for all groups × all exposures
irf_inc <- map_dfr(names(panel_by_group), function(grp) {
  d_grp <- panel_by_group[[grp]]
  n_c   <- n_distinct(d_grp$iso3c)
  if (n_c < 5) {
    message("  Skipping ", grp, " (only ", n_c, " countries)")
    return(NULL)
  }
  map_dfr(TREAT_TYPES, function(tt) {
    run_lp(d_grp, paste0("shock_", tt)) |>
      mutate(treat = tt, income_group = grp)
  })
})

# Console summary: h=1 and h=4
SHOW_H <- c(1, 4)
cat(sprintf("\n%-14s  %-14s  %5s  %8s  %6s  %5s\n",
            "Exposure", "Income group", "h", "est", "se", "p"))
cat(strrep("-", 58), "\n")
for (tt in TREAT_TYPES) {
  for (grp in INCOME_LEVELS) {
    for (h_ in SHOW_H) {
      r <- filter(irf_inc, treat == tt, income_group == grp, h == h_)
      if (nrow(r) == 0) next
      cat(sprintf("%-14s  %-14s  h=%d  %+8.3f  %6.3f  %4.3f%s\n",
                  TREAT_LABELS[tt], grp, h_,
                  r$est, r$se, r$pval, stars_console(r$pval)))
    }
  }
  cat("\n")
}

# =============================================================================
# Section 2: LaTeX table
# tab:income-lp
# Rows: exposure types (coeff + SE)
# Column blocks: Low | Lower-middle | Upper-middle | High (× h=1, h=4)
# =============================================================================
message("\n=== Section 2: LaTeX table (tab_income_lp.tex) ===")

# 4 income groups × 2 horizons = 8 data columns
# col spec: l rr rr rr rr
col_spec <- "l rr rr rr rr"

tab_inc <- c(
  "% Income group LP: all exposure types, h=1 and h=4",
  "",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Local projections by income group (FSI outcome)}",
  "\\label{tab:app-income-lp}",
  "\\begin{threeparttable}",
  "\\small",
  sprintf("\\begin{tabular}{%s}", col_spec),
  "\\toprule",
  # Income group spanning headers
  paste0(
    " & ",
    paste(
      sprintf("\\multicolumn{2}{c}{%s}", INCOME_LEVELS),
      collapse = " & "
    ),
    " \\\\"
  ),
  # Horizontal rules under each span
  paste0(
    paste(sprintf("\\cmidrule(lr){%d-%d}", seq(2, 9, 2), seq(3, 10, 2)),
          collapse = " "),
    ""
  ),
  # h sub-headers
  paste0(
    " & ",
    paste(rep("$h=1$ & $h=4$", 4), collapse = " & "),
    " \\\\"
  ),
  "\\midrule"
)

for (tt in TREAT_TYPES) {
  d_tt <- filter(irf_inc, treat == tt)

  est_cells <- map_chr(INCOME_LEVELS, function(grp) {
    r1 <- filter(d_tt, income_group == grp, h == 1)
    r4 <- filter(d_tt, income_group == grp, h == 4)
    c1 <- if (nrow(r1) > 0 && !is.na(r1$est))
            sprintf("%.3f%s", r1$est, stars_fn(r1$pval)) else "---"
    c4 <- if (nrow(r4) > 0 && !is.na(r4$est))
            sprintf("%.3f%s", r4$est, stars_fn(r4$pval)) else "---"
    paste(c1, c4, sep = " & ")
  })

  se_cells <- map_chr(INCOME_LEVELS, function(grp) {
    r1 <- filter(d_tt, income_group == grp, h == 1)
    r4 <- filter(d_tt, income_group == grp, h == 4)
    s1 <- if (nrow(r1) > 0 && !is.na(r1$se)) sprintf("(%.3f)", r1$se) else ""
    s4 <- if (nrow(r4) > 0 && !is.na(r4$se)) sprintf("(%.3f)", r4$se) else ""
    paste(s1, s4, sep = " & ")
  })

  tab_inc <- c(tab_inc,
    sprintf("%s & %s \\\\", TREAT_LABELS[tt], paste(est_cells, collapse = " & ")),
    sprintf("   & %s \\\\[4pt]", paste(se_cells, collapse = " & "))
  )
}

# N row: one N per (income group × horizon), using consumption exposure
n_cells <- map_chr(INCOME_LEVELS, function(grp) {
  r1 <- filter(irf_inc, treat == "cons", income_group == grp, h == 1)
  r4 <- filter(irf_inc, treat == "cons", income_group == grp, h == 4)
  n1 <- if (nrow(r1) > 0) format(r1$n_obs, big.mark = ",") else "---"
  n4 <- if (nrow(r4) > 0) format(r4$n_obs, big.mark = ",") else "---"
  paste(n1, n4, sep = " & ")
})

tab_inc <- c(tab_inc, "\\midrule",
             sprintf("$N$ & %s \\\\", paste(n_cells, collapse = " & ")))

tab_inc <- c(
  tab_inc,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0(
    "\\item \\textit{Notes:} LP: $\\text{FSI}_{i,t+h} - \\text{FSI}_{i,t-1} = ",
    "\\beta\\, x_{it} + \\sum_{\\ell=1}^{2}\\gamma_\\ell x_{i,t-\\ell} + \\alpha_t + \\varepsilon_{iht}$. ",
    "Samples split by World Bank income classification (2010, fixed). ",
    "Year fixed effects; standard errors clustered by country. ",
    "$N$ from consumption exposure; other exposures may differ slightly. ",
    "Significance: $***$ $p<0.001$, $**$ $p<0.01$, $*$ $p<0.05$, $\\cdot$ $p<0.10$."
  ),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)

write_tex(tab_inc, "tab_app_income_lp.tex")

# =============================================================================
# Section 3: Income group IRF figure
# Consumption + production exposures × 4 income groups
# facet_grid(treat ~ income_group), 180×100 mm
# =============================================================================
message("\n=== Section 3: Income group IRF figure ===")

COL_NEG <- "#2166ac"
COL_POS <- "#d6604d"
COL_RIB <- "grey70"

theme_nature <- function(base_size = 8) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid.minor    = element_blank(),
      panel.grid.major.x  = element_blank(),
      panel.grid.major.y  = element_line(colour = "grey92", linewidth = 0.3),
      axis.line.x         = element_line(colour = "grey40", linewidth = 0.4),
      axis.ticks.x        = element_line(colour = "grey40", linewidth = 0.4),
      axis.ticks.length   = unit(2, "pt"),
      strip.text          = element_text(face = "bold", size = base_size),
      axis.title          = element_text(size = base_size),
      axis.text           = element_text(size = base_size - 1),
      legend.position     = "bottom",
      plot.tag            = element_text(face = "bold", size = base_size + 1)
    )
}

fig3_dat <- irf_inc |>
  filter(income_group %in% INCOME_LEVELS) |>
  mutate(
    income_f = factor(income_group,
                      levels = INCOME_LEVELS,
                      labels = INCOME_LABELS),
    treat_f  = factor(treat,
                      levels = TREAT_TYPES,
                      labels = unname(TREAT_LABELS))
  )

# Color each facet panel by dominant sign of its mean IRF
fig3_sign <- fig3_dat |>
  group_by(treat_f, income_f) |>
  summarise(col = if_else(mean(est, na.rm = TRUE) < 0, COL_NEG, COL_POS),
            .groups = "drop")

fig3_dat <- fig3_dat |>
  left_join(fig3_sign, by = c("treat_f", "income_f"))

fig3 <- ggplot(fig3_dat, aes(x = h, y = est)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.35) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              fill = COL_RIB, alpha = 0.25) +
  geom_line(aes(colour = col), linewidth = 0.7) +
  geom_point(aes(colour = col), size = 1.6) +
  scale_colour_identity() +
  scale_x_continuous(breaks = 0:4, name = "Horizon (years)") +
  scale_y_continuous(name = "FSI change (index points)") +
  facet_grid(treat_f ~ income_f, scales = "free_y") +
  theme_nature() +
  theme(strip.text.y = element_text(angle = 0, hjust = 0))

ggsave(file.path(OUT_FIGS, "fig_app_income_irf.pdf"),
       fig3, width = 180, height = 150, units = "mm")
ggsave(file.path(OUT_FIGS, "fig_app_income_irf.png"),
       fig3, width = 180, height = 150, units = "mm", dpi = 300)

message("  Saved: fig_app_income_irf.pdf / .png")

# =============================================================================
# Section 4: Per-treatment appendix tables — income groups × all horizons
# One table per treatment type; rows = income group, cols = h=0..4
# Labels: tab:app-income-cons, tab:app-income-prod, etc.
# =============================================================================
message("\n=== Section 4: Per-treatment appendix tables ===")

INCOME_LABELS_TEX <- c(
  Low            = "Low income",
  `Lower-middle` = "Lower-middle",
  `Upper-middle` = "Upper-middle",
  High           = "High income"
)

TREAT_TITLES <- c(
  cons = "Consumption",
  prod = "Production",
  imp  = "Import",
  exp  = "Export",
  net  = "Net import"
)

for (tt in TREAT_TYPES) {
  d_tt <- filter(irf_inc, treat == tt, income_group %in% INCOME_LEVELS)

  tab <- c(
    sprintf("%% Income-group LP: %s exposure, all horizons", TREAT_TITLES[tt]),
    "",
    "\\begin{table}[htbp]",
    "\\centering",
    sprintf("\\caption{Income-group local projections: %s exposure (FSI outcome)}",
            TREAT_TITLES[tt]),
    sprintf("\\label{tab:app-income-%s}", tt),
    "\\begin{threeparttable}",
    "\\small",
    "\\begin{tabular}{l rrrrr}",
    "\\toprule",
    paste0("Income group & $h=0$ & $h=1$ & $h=2$ & $h=3$ & $h=4$ \\\\"),
    "\\midrule"
  )

  for (ig in INCOME_LEVELS) {
    est_row <- map_chr(LP_HORIZONS, function(h_) {
      r <- filter(d_tt, income_group == ig, h == h_)
      if (nrow(r) == 0 || is.na(r$est)) return("---")
      sprintf("$%+.3f%s$", r$est, stars_fn(r$pval))
    })
    se_row <- map_chr(LP_HORIZONS, function(h_) {
      r <- filter(d_tt, income_group == ig, h == h_)
      if (nrow(r) == 0 || is.na(r$se)) return("")
      sprintf("$(%.3f)$", r$se)
    })
    tab <- c(tab,
      paste0(INCOME_LABELS_TEX[ig], " & ",
             paste(est_row, collapse = " & "), " \\\\"),
      paste0("& ", paste(se_row, collapse = " & "), " \\\\[3pt]")
    )
  }

  # N row: one N per horizon, from this treatment type
  n_row <- map_chr(LP_HORIZONS, function(h_) {
    r <- filter(d_tt, income_group == "Low", h == h_)
    if (nrow(r) == 0) return("---")
    format(r$n_obs, big.mark = ",")
  })

  tab <- c(tab,
    "\\midrule",
    paste0("$N$ (Low) & ", paste(n_row, collapse = " & "), " \\\\"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{tablenotes}\\small",
    paste0("\\item \\textit{Notes:} LP estimates of ",
           "$\\Delta_h\\text{FSI}_{it} = \\beta_h x_{it} + ",
           "\\sum_{\\ell=1}^{2}\\gamma_\\ell x_{i,t-\\ell} + \\alpha_i + \\alpha_t + \\varepsilon_{iht}$, ",
           "where $\\Delta_h\\text{FSI}_{it} = \\text{FSI}_{i,t+h} - \\text{FSI}_{i,t-1}$. ",
           "Country and year fixed effects; SE clustered by country. ",
           "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."),
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\end{table}"
  )

  fname <- sprintf("tab_app_income_%s.tex", tt)
  write_tex(tab, fname)
}

message("\nDone.")
