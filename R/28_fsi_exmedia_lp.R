# =============================================================================
# 28_fsi_exmedia_lp.R
# Baseline LP with media-independent FSI sub-index as outcome.
#
# Concern: the cohesion indicators (C1–C3) and the political indicators P1 and
# P3 are scored primarily via CAST media content analysis, so their scores may
# correlate with food-price shocks through coverage intensity rather than true
# political effects.
#
# This script drops C1, C2, C3, P1, and P3 from the FSI total:
#   fsi_exmedia = fsi_total - C1 - C2 - C3 - P1 - P3   (7-indicator sub-index)
#   Outcome: fsi_exmedia_{t+h} - fsi_exmedia_{t-1}
#
# Outputs:
#   output/tables/tab_exmedia_lp.tex    — side-by-side vs. baseline, cons + prod
# =============================================================================

source(here::here("R", "00_setup.R"))

# ---- Load data --------------------------------------------------------------
panel   <- readRDS(file.path(DATA_PRO, "panel.rds"))
fsi_raw <- readRDS(file.path(DATA_PRO, "fsi.rds"))

message("Panel: ", nrow(panel), " obs | ",
        n_distinct(panel$iso3c), " countries | ",
        min(panel$year), "–", max(panel$year))

# ---- Build fsi_exmedia ------------------------------------------------------
fsi_em <- fsi_raw |>
  select(iso3c, year, fsi_total,
         c1 = c1_security_apparatus,
         c2 = c2_factionalized_elites,
         c3 = c3_group_grievance,
         p1 = p1_state_legitimacy,
         p3 = p3_human_rights) |>
  mutate(fsi_exmedia = fsi_total - c1 - c2 - c3 - p1 - p3) |>
  arrange(iso3c, year) |>
  group_by(iso3c) |>
  mutate(
    fsi_em_lag1  = lag(fsi_exmedia, 1),
    dfsi_em_h0   = fsi_exmedia                - fsi_em_lag1,
    dfsi_em_h1   = lead(fsi_exmedia, 1)       - fsi_em_lag1,
    dfsi_em_h2   = lead(fsi_exmedia, 2)       - fsi_em_lag1,
    dfsi_em_h3   = lead(fsi_exmedia, 3)       - fsi_em_lag1,
    dfsi_em_h4   = lead(fsi_exmedia, 4)       - fsi_em_lag1
  ) |>
  ungroup() |>
  select(iso3c, year, fsi_exmedia, fsi_em_lag1,
         dfsi_em_h0, dfsi_em_h1, dfsi_em_h2, dfsi_em_h3, dfsi_em_h4)

cat("fsi_exmedia range:", round(range(fsi_em$fsi_exmedia, na.rm=TRUE), 2),
    "| non-NA h=0:", sum(!is.na(fsi_em$dfsi_em_h0)),
    "| non-NA h=4:", sum(!is.na(fsi_em$dfsi_em_h4)), "\n\n")

# Merge into panel
panel_em <- panel |>
  left_join(fsi_em |> select(iso3c, year, starts_with("dfsi_em")),
            by = c("iso3c", "year"))

message("Panel after merge: ", nrow(panel_em),
        " rows | dfsi_em_h0 non-NA: ", sum(!is.na(panel_em$dfsi_em_h0)))

# ---- LP runner (generic outcome stem) ----------------------------------------
CI_Z <- 1.645

stars_fn <- function(p)
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", p < 0.10 ~ "$\\cdot$", TRUE ~ "")

run_lp <- function(data, shock_var, fd_stem, horizons = LP_HORIZONS) {
  lag1 <- paste0(shock_var, "_l1")
  lag2 <- paste0(shock_var, "_l2")
  map_dfr(horizons, function(h) {
    y_col <- paste0(fd_stem, "_h", h)
    if (!y_col %in% names(data)) return(NULL)
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
           ci_lo = est - CI_Z * se,
           ci_hi = est + CI_Z * se,
           n_obs = nobs(fit))
  })
}

# ---- Estimate: baseline (dfsi) vs. ex-media (dfsi_em) -----------------------
TREAT_TYPES  <- c("cons", "prod", "imp", "exp", "net")
TREAT_LABELS <- c(cons = "Consumption", prod = "Production",
                  imp  = "Imports",     exp  = "Exports",
                  net  = "Net imports")

message("\n=== Running baseline LP ===")
irf_base <- map_dfr(TREAT_TYPES, function(tt) {
  sv <- paste0("shock_", tt)
  run_lp(panel_em, sv, "dfsi") |> mutate(treat = tt, spec = "baseline")
})

message("=== Running ex-media LP ===")
irf_em <- map_dfr(TREAT_TYPES, function(tt) {
  sv <- paste0("shock_", tt)
  run_lp(panel_em, sv, "dfsi_em") |> mutate(treat = tt, spec = "exmedia")
})

irf_all <- bind_rows(irf_base, irf_em)

# Console summary
cat("\n=== Results: baseline vs. ex-media (h=1 and h=4) ===\n")
cat(sprintf("%-14s  %-10s  %8s  %8s  %8s  %8s\n",
            "Treat", "Spec", "h=1 est", "h=1 p", "h=4 est", "h=4 p"))
cat(strrep("-", 70), "\n")
for (tt in TREAT_TYPES) {
  for (sp in c("baseline", "exmedia")) {
    r1 <- filter(irf_all, treat == tt, spec == sp, h == 1)
    r4 <- filter(irf_all, treat == tt, spec == sp, h == 4)
    cat(sprintf("%-14s  %-10s  %+8.3f  %8.3f  %+8.3f  %8.3f\n",
                TREAT_LABELS[tt], sp,
                if (nrow(r1)) r1$est else NA,
                if (nrow(r1)) r1$pval else NA,
                if (nrow(r4)) r4$est else NA,
                if (nrow(r4)) r4$pval else NA))
  }
  cat("\n")
}

# =============================================================================
# LaTeX table: cons + prod, all horizons, baseline vs. ex-media side by side
# =============================================================================
make_rows <- function(tt, spec_tag, spec_label) {
  r <- filter(irf_all, treat == tt, spec == spec_tag)
  est_row <- paste(map_chr(LP_HORIZONS, function(h_) {
    rx <- filter(r, h == h_)
    if (nrow(rx) == 0) return("---")
    sprintf("$%+.3f%s$", rx$est, stars_fn(rx$pval))
  }), collapse = " & ")
  se_row <- paste(map_chr(LP_HORIZONS, function(h_) {
    rx <- filter(r, h == h_)
    if (nrow(rx) == 0) return("")
    sprintf("$(%.3f)$", rx$se)
  }), collapse = " & ")
  n_ <- filter(r, h == 0)$n_obs
  c(paste0("\\quad ", spec_label, " & ", est_row, " \\\\"),
    paste0("& ", se_row, " \\\\[2pt]"))
}

tab <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Baseline LP vs.\\ media-independent FSI sub-index}",
  "\\label{tab:exmedia-lp}",
  "\\begin{threeparttable}\\small",
  "\\begin{tabular}{lrrrrr}",
  "\\toprule",
  "Exposure & $h=0$ & $h=1$ & $h=2$ & $h=3$ & $h=4$ \\\\",
  "\\midrule",
  "\\multicolumn{6}{l}{\\textit{Consumption exposure}} \\\\",
  "\\midrule",
  make_rows("cons", "baseline", "Full FSI"),
  make_rows("cons", "exmedia",  "Ex media (FSI$-$C1$-$C2$-$C3$-$P1$-$P3)"),
  sprintf("\\multicolumn{6}{l}{$N = %s$ (both specs)}  \\\\[6pt]",
          format(filter(irf_all, treat=="cons", spec=="baseline", h==0)$n_obs,
                 big.mark=",")),
  "\\midrule",
  "\\multicolumn{6}{l}{\\textit{Production exposure}} \\\\",
  "\\midrule",
  make_rows("prod", "baseline", "Full FSI"),
  make_rows("prod", "exmedia",  "Ex media (FSI$-$C1$-$C2$-$C3$-$P1$-$P3)"),
  sprintf("\\multicolumn{6}{l}{$N = %s$ (both specs)}  \\\\[6pt]",
          format(filter(irf_all, treat=="prod", spec=="baseline", h==0)$n_obs,
                 big.mark=",")),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0("\\item \\textit{Notes:} Each panel shows the baseline LP (outcome: ",
         "$\\Delta_h\\text{FSI}_{it}$) alongside the LP on the media-independent ",
         "sub-index ($\\text{FSI} - \\text{C1} - \\text{C2} - \\text{C3} - \\text{P1} - \\text{P3}$), ",
         "which excludes all cohesion indicators (C1--C3) and the political indicators ",
         "State Legitimacy (P1) and Human Rights (P3), all scored primarily via the ",
         "Fund for Peace CAST news content analysis tool. The remaining seven indicators ",
         "draw predominantly on administrative and survey data. ",
         "Specification: year FE, two lags of the shock; SE clustered by country. ",
         "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."),
  "\\end{tablenotes}\\end{threeparttable}\\end{table}"
)

out_path <- file.path(OUT_TABS, "tab_exmedia_lp.tex")
writeLines(tab, out_path)
message("\nSaved: tab_exmedia_lp.tex")

# =============================================================================
# By income group: baseline vs. ex-media, cons + prod, h=1 and h=4
# =============================================================================
message("\n=== Running by income group ===")

INCOME_LEVELS <- c("Low", "Lower-middle", "Upper-middle", "High")
INCOME_LABELS <- c("Low income", "Lower-middle", "Upper-middle", "High income")

panel_by_inc <- setNames(
  map(INCOME_LEVELS, ~ filter(panel_em, income_group == .x)),
  INCOME_LEVELS
)

irf_inc_base <- map_dfr(names(panel_by_inc), function(ig) {
  d <- panel_by_inc[[ig]]
  map_dfr(c("cons", "prod"), function(tt) {
    sv <- paste0("shock_", tt)
    run_lp(d, sv, "dfsi") |> mutate(treat = tt, spec = "baseline", income_group = ig)
  })
})

irf_inc_em <- map_dfr(names(panel_by_inc), function(ig) {
  d <- panel_by_inc[[ig]]
  map_dfr(c("cons", "prod"), function(tt) {
    sv <- paste0("shock_", tt)
    run_lp(d, sv, "dfsi_em") |> mutate(treat = tt, spec = "exmedia", income_group = ig)
  })
})

irf_inc_all <- bind_rows(irf_inc_base, irf_inc_em)

# Console summary at h=1 and h=4
cat("\n=== By income group: baseline vs. ex-media (cons, h=1 and h=4) ===\n")
cat(sprintf("%-16s  %-10s  %8s  %8s  %8s  %8s\n",
            "Income", "Spec", "h=1 est", "h=1 p", "h=4 est", "h=4 p"))
cat(strrep("-", 72), "\n")
for (ig in INCOME_LEVELS) {
  for (sp in c("baseline", "exmedia")) {
    r1 <- filter(irf_inc_all, income_group == ig, treat == "cons", spec == sp, h == 1)
    r4 <- filter(irf_inc_all, income_group == ig, treat == "cons", spec == sp, h == 4)
    cat(sprintf("%-16s  %-10s  %+8.3f  %8.3f  %+8.3f  %8.3f\n",
                ig, sp,
                if (nrow(r1)) r1$est else NA, if (nrow(r1)) r1$pval else NA,
                if (nrow(r4)) r4$est else NA, if (nrow(r4)) r4$pval else NA))
  }
  cat("\n")
}

cat("\n=== By income group: baseline vs. ex-media (prod, h=1 and h=4) ===\n")
cat(sprintf("%-16s  %-10s  %8s  %8s  %8s  %8s\n",
            "Income", "Spec", "h=1 est", "h=1 p", "h=4 est", "h=4 p"))
cat(strrep("-", 72), "\n")
for (ig in INCOME_LEVELS) {
  for (sp in c("baseline", "exmedia")) {
    r1 <- filter(irf_inc_all, income_group == ig, treat == "prod", spec == sp, h == 1)
    r4 <- filter(irf_inc_all, income_group == ig, treat == "prod", spec == sp, h == 4)
    cat(sprintf("%-16s  %-10s  %+8.3f  %8.3f  %+8.3f  %8.3f\n",
                ig, sp,
                if (nrow(r1)) r1$est else NA, if (nrow(r1)) r1$pval else NA,
                if (nrow(r4)) r4$est else NA, if (nrow(r4)) r4$pval else NA))
  }
  cat("\n")
}

# =============================================================================
# LaTeX: income group table, cons only, h=1 and h=4, baseline vs. ex-media
# =============================================================================
fmt_cell <- function(ig, tt, sp, h_) {
  r <- filter(irf_inc_all, income_group == ig, treat == tt, spec == sp, h == h_)
  if (nrow(r) == 0 || is.na(r$est)) return(list(e = "---", s = ""))
  list(e = sprintf("$%+.3f%s$", r$est, stars_fn(r$pval)),
       s = sprintf("$(%.3f)$", r$se))
}

make_inc_panel <- function(tt, treat_lbl) {
  rows <- c(sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", 1 + 2*length(INCOME_LEVELS), treat_lbl),
            "\\midrule")
  for (sp in c("baseline", "exmedia")) {
    sp_lbl <- if (sp == "baseline") "Full FSI" else "Ex media"
    est_cells <- paste(map_chr(INCOME_LEVELS, function(ig)
      paste0(fmt_cell(ig, tt, sp, 1)$e, " & ", fmt_cell(ig, tt, sp, 4)$e)),
      collapse = " & ")
    se_cells <- paste(map_chr(INCOME_LEVELS, function(ig)
      paste0(fmt_cell(ig, tt, sp, 1)$s, " & ", fmt_cell(ig, tt, sp, 4)$s)),
      collapse = " & ")
    rows <- c(rows,
      paste0("\\quad ", sp_lbl, " & ", est_cells, " \\\\"),
      paste0("& ", se_cells, " \\\\[3pt]"))
  }
  c(rows, "\\addlinespace")
}

n_cols <- 1 + 2 * length(INCOME_LEVELS)
col_spec <- paste0("l", paste(rep("rr", length(INCOME_LEVELS)), collapse = ""))

tab_inc <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Baseline vs.\\ ex-media LP by income group ($h=1$ and $h=4$)}",
  "\\label{tab:exmedia-lp-income}",
  "\\begin{adjustbox}{max width=\\textwidth}",
  "\\begin{threeparttable}\\small",
  sprintf("\\begin{tabular}{%s}", col_spec),
  "\\toprule",
  paste0(" & ",
    paste(map_chr(INCOME_LABELS, ~ sprintf("\\multicolumn{2}{c}{%s}", .x)),
          collapse = " & "), " \\\\"),
  paste(sprintf("\\cmidrule(lr){%d-%d}",
                2 + 2*(0:(length(INCOME_LEVELS)-1)),
                3 + 2*(0:(length(INCOME_LEVELS)-1))),
        collapse = ""),
  paste0("Spec & ",
    paste(rep("$h=1$ & $h=4$", length(INCOME_LEVELS)), collapse = " & "), " \\\\"),
  "\\midrule",
  make_inc_panel("cons", "Consumption exposure"),
  make_inc_panel("prod", "Production exposure"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0("\\item \\textit{Notes:} Each cell pair shows $h=1$ and $h=4$ LP coefficients. ",
         "Full FSI outcome: $\\Delta_h\\text{FSI}_{it}$. ",
         "Ex media outcome: FSI minus C1, C2, C3, P1, and P3 (CAST-scored indicators). ",
         "Year FE; SE clustered by country. ",
         "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."),
  "\\end{tablenotes}\\end{threeparttable}",
  "\\end{adjustbox}\\end{table}"
)

writeLines(tab_inc, file.path(OUT_TABS, "tab_exmedia_lp_income.tex"))
message("  Saved: tab_exmedia_lp_income.tex")
message("Done.")
