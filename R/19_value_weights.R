# =============================================================================
# 19_value_weights.R
# Robustness: value-weighted exposure shocks vs caloric-share baseline.
#
# The main specification weights commodity price changes by caloric shares
# (consumption) or quantity shares (production, trade). This script runs the
# same LP using value-weighted shares instead:
#   shock_val_cons: FAO FBS monetary value shares × Δlog price
#   shock_val_prod: FAO production value shares × Δlog price
#
# Both alternative shocks are already in panel.rds (built by 04_shocks.R).
#
# Outputs:
#   output/tables/tab_value_weights.tex
#     Aggregate LP — caloric vs value, cons + prod, h = 0..4
#   output/tables/tab_value_weights_income.tex
#     By income group — h = 1 and h = 4 (main results horizons)
# =============================================================================

source(here::here("R", "00_setup.R"))

panel <- readRDS(file.path(DATA_PRO, "panel.rds"))

message("Panel: ", nrow(panel), " obs | ", min(panel$year), "–", max(panel$year))

INCOME_LEVELS <- c("Low", "Lower-middle", "Upper-middle", "High")
INCOME_LABELS <- c("Low income", "Lower-middle", "Upper-middle", "High income")

# ---- LP runner (reuses style from 08_lp.R) ----------------------------------
run_lp <- function(data, shock_var, horizons = LP_HORIZONS) {
  lag1 <- paste0(shock_var, "_l1")
  lag2 <- paste0(shock_var, "_l2")
  map_dfr(horizons, function(h) {
    y_col <- paste0("dfsi_h", h)
    if (!y_col %in% names(data)) return(NULL)
    d <- data |>
      rename(.y = all_of(y_col), .s = all_of(shock_var)) |>
      filter(!is.na(.y), !is.na(.s))
    if (lag1 %in% names(d)) d <- rename(d, .sl1 = all_of(lag1))
    if (lag2 %in% names(d)) d <- rename(d, .sl2 = all_of(lag2))
    rhs <- intersect(c(".sl1", ".sl2"), names(d))
    fml <- as.formula(paste0(".y ~ .s",
      if (length(rhs)) paste0(" + ", paste(rhs, collapse = " + ")), " | year"))
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

stars_fn <- function(p)
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", p < 0.10 ~ "$\\cdot$", TRUE ~ "")

# =============================================================================
# Section 1: Aggregate LP — caloric vs value, cons + prod
# =============================================================================
message("\n=== Section 1: Aggregate comparison (caloric vs value) ===")

specs <- list(
  list(sv = "shock_cons",     label = "Consumption (caloric)"),
  list(sv = "shock_val_cons", label = "Consumption (value)"),
  list(sv = "shock_prod",     label = "Production (caloric)"),
  list(sv = "shock_val_prod", label = "Production (value)")
)

irf_agg <- map_dfr(specs, function(s) {
  run_lp(panel, s$sv) |> mutate(label = s$label, shock_var = s$sv)
})

cat("\nAggregate LP: caloric vs value (h=0..4)\n")
cat(sprintf("%-28s  %5s  %8s  %6s  %6s\n", "Spec", "h", "est", "se", "p"))
cat(strrep("-", 58), "\n")
for (s in specs) {
  for (h_ in LP_HORIZONS) {
    r <- filter(irf_agg, shock_var == s$sv, h == h_)
    if (nrow(r) == 0) next
    cat(sprintf("%-28s  h=%d  %+8.4f  %6.4f  %6.3f%s\n",
                s$label, h_, r$est, r$se, r$pval, stars_fn(r$pval)))
  }
  cat("\n")
}

# LaTeX: aggregate table
fmt_cell <- function(est, se, pval) {
  if (is.na(est)) return("---")
  sprintf("$%+.3f%s$ \\\\ $(%.3f)$", est, stars_fn(pval), se)
}

agg_n <- map_int(specs, ~ {
  r <- filter(irf_agg, shock_var == .x$sv, h == 0)
  if (nrow(r) == 0) NA_integer_ else as.integer(r$n_obs)
})

make_two_row <- function(sv, lbl) {
  r <- filter(irf_agg, shock_var == sv)
  est_row <- paste(sapply(LP_HORIZONS, function(h_) {
    rx <- filter(r, h == h_)
    if (nrow(rx) == 0) return("---")
    sprintf("$%+.3f%s$", rx$est, stars_fn(rx$pval))
  }), collapse = " & ")
  se_row <- paste(sapply(LP_HORIZONS, function(h_) {
    rx <- filter(r, h == h_)
    if (nrow(rx) == 0) return("")
    sprintf("$(%.3f)$", rx$se)
  }), collapse = " & ")
  n_ <- filter(r, h == 0)$n_obs
  c(paste0(lbl, " & ", est_row, " \\\\"),
    paste0("& ", se_row, " \\\\"),
    paste0("& \\multicolumn{5}{l}{$N = ", ifelse(length(n_) > 0, n_[1], "---"), "$} \\\\[2pt]"))
}

tex_agg <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Value-weighted vs caloric-weighted shocks: aggregate LP}",
  "\\label{tab:value-weights-agg}",
  "\\begin{threeparttable}\\small",
  "\\begin{tabular}{lrrrrr}",
  "\\toprule",
  "Specification & $h=0$ & $h=1$ & $h=2$ & $h=3$ & $h=4$ \\\\",
  "\\midrule",
  "\\textit{Consumption exposure} \\\\",
  make_two_row("shock_cons",     "\\quad Caloric shares"),
  make_two_row("shock_val_cons", "\\quad Value shares"),
  "\\textit{Production exposure} \\\\",
  make_two_row("shock_prod",     "\\quad Caloric shares"),
  make_two_row("shock_val_prod", "\\quad Value shares"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0("\\item \\textit{Notes:} LP specification: ",
         "$\\Delta_{h}\\text{FSI}_{it} = \\beta x_{it} + \\sum_{k=1}^{2}\\gamma_k x_{i,t-k} + \\alpha_t + \\varepsilon_{iht}$. ",
         "Year FE; SE clustered by country. ",
         "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."),
  "\\end{tablenotes}\\end{threeparttable}\\end{table}"
)
writeLines(tex_agg, file.path(OUT_TABS, "tab_value_weights.tex"))
message("  Saved: tab_value_weights.tex")

# =============================================================================
# Section 2: By income group — h=1 and h=4
# =============================================================================
message("\n=== Section 2: Value weights by income group (h=1, h=4) ===")

panel_by_inc <- c(
  list(All = panel),
  setNames(map(INCOME_LEVELS, ~ filter(panel, income_group == .x)), INCOME_LEVELS)
)

irf_inc <- map_dfr(names(panel_by_inc), function(ig) {
  d <- panel_by_inc[[ig]]
  map_dfr(specs, function(s) {
    run_lp(d, s$sv) |> mutate(income_group = ig, label = s$label, shock_var = s$sv)
  })
})

cat("\nBy income group: h=1 and h=4\n")
irf_inc |>
  filter(h %in% c(1, 4)) |>
  mutate(cell = sprintf("%+.3f (%.3f)%s", est, se,
                        ifelse(pval < 0.05, "*", ifelse(pval < 0.10, ".", "")))) |>
  select(income_group, label, h, cell) |>
  pivot_wider(names_from = h, values_from = cell) |>
  arrange(income_group, label) |>
  print(n = 50)

# LaTeX: income group table (h=1 and h=4)
inc_groups <- c("All", INCOME_LEVELS)
fmt_inc <- function(ig, sv, h_) {
  r <- filter(irf_inc, income_group == ig, shock_var == sv, h == h_)
  if (nrow(r) == 0) return(list(e = "---", s = ""))
  list(e = sprintf("$%+.3f%s$", r$est, stars_fn(r$pval)),
       s = sprintf("$(%.3f)$", r$se))
}

tex_inc_lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Value-weighted vs caloric-weighted shocks by income group ($h=1$ and $h=4$)}",
  "\\label{tab:value-weights-income}",
  "\\begin{threeparttable}\\small",
  sprintf("\\begin{tabular}{ll%s}", paste(rep("rr", length(inc_groups)), collapse="")),
  "\\toprule",
  paste0(" & & ",
    paste(sapply(inc_groups, function(ig)
      paste0("\\multicolumn{2}{c}{", if (ig == "All") "All" else ig, "}")),
      collapse = " & "), " \\\\"),
  paste0("\\cmidrule(lr){3-4}",
    paste(sprintf("\\cmidrule(lr){%d-%d}", seq(5, 5 + 2*(length(inc_groups)-1), 2),
                  seq(6, 6 + 2*(length(inc_groups)-1), 2)), collapse = "")),
  paste0("Exposure & Weight & ",
    paste(rep("$h=1$ & $h=4$", length(inc_groups)), collapse = " & "), " \\\\"),
  "\\midrule"
)
for (base_sv in c("shock_cons", "shock_prod")) {
  lbl <- if (grepl("cons", base_sv)) "Consumption" else "Production"
  val_sv <- paste0("shock_val_", sub("shock_", "", base_sv))
  tex_inc_lines <- c(tex_inc_lines,
    paste0("\\multirow{4}{*}{", lbl, "} & Caloric &",
      paste(sapply(inc_groups, function(ig)
        paste0(fmt_inc(ig, base_sv, 1)$e, " & ", fmt_inc(ig, base_sv, 4)$e)),
        collapse = " & "), " \\\\"),
    paste0("& &",
      paste(sapply(inc_groups, function(ig)
        paste0(fmt_inc(ig, base_sv, 1)$s, " & ", fmt_inc(ig, base_sv, 4)$s)),
        collapse = " & "), " \\\\"),
    paste0("& Value &",
      paste(sapply(inc_groups, function(ig)
        paste0(fmt_inc(ig, val_sv, 1)$e, " & ", fmt_inc(ig, val_sv, 4)$e)),
        collapse = " & "), " \\\\"),
    paste0("& &",
      paste(sapply(inc_groups, function(ig)
        paste0(fmt_inc(ig, val_sv, 1)$s, " & ", fmt_inc(ig, val_sv, 4)$s)),
        collapse = " & "), " \\\\[4pt]")
  )
}
tex_inc_lines <- c(tex_inc_lines,
  "\\bottomrule", "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0("\\item \\textit{Notes:} Each cell shows the LP coefficient at the indicated horizon. ",
         "Year FE; SE clustered by country in parentheses. ",
         "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."),
  "\\end{tablenotes}\\end{threeparttable}\\end{table}"
)
writeLines(tex_inc_lines, file.path(OUT_TABS, "tab_value_weights_income.tex"))
message("  Saved: tab_value_weights_income.tex")
message("\nDone.")
