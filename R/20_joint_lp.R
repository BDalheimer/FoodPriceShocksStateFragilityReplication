# =============================================================================
# 20_joint_lp.R
# Joint consumer + producer LP: both exposure types in one equation.
#
# In the main specifications (08_lp.R), consumption and production shocks are
# estimated in separate regressions. This script includes both simultaneously:
#
#   ΔFSI_{i,t+h} = β_c shock_cons_{it} + β_p shock_prod_{it}
#                + Σ_{k=1}^{2} [γ_c^k shock_cons_{i,t-k} + γ_p^k shock_prod_{i,t-k}]
#                + α_t + ε_{iht}
#
# Motivation: tests whether consumption and production effects are independent
# or whether one absorbs the other. In agricultural economies the two shocks
# are correlated (countries that produce cereals also consume them), so the
# joint specification separates the channels.
#
# Output:
#   output/tables/tab_joint_lp.tex   — aggregate, all horizons
# =============================================================================

source(here::here("R", "00_setup.R"))

panel <- readRDS(file.path(DATA_PRO, "panel.rds"))

message("Panel: ", nrow(panel), " obs | ",
        n_distinct(panel$iso3c), " countries | ",
        min(panel$year), "–", max(panel$year))

stars_fn <- function(p)
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", p < 0.10 ~ "$\\cdot$", TRUE ~ "")

# ---- Joint LP runner --------------------------------------------------------
run_joint_lp <- function(data, horizons = LP_HORIZONS) {
  map_dfr(horizons, function(h) {
    y_col <- paste0("dfsi_h", h)
    if (!y_col %in% names(data)) return(NULL)
    d <- data |>
      rename(.y = all_of(y_col)) |>
      filter(!is.na(.y), !is.na(shock_cons), !is.na(shock_prod))
    if (nrow(d) < 30 || n_distinct(d$iso3c) < 5) return(NULL)
    fit <- tryCatch(
      feols(.y ~ shock_cons + shock_cons_l1 + shock_cons_l2 +
                 shock_prod + shock_prod_l1 + shock_prod_l2 | year,
            data = d, cluster = ~iso3c, warn = FALSE, notes = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    ct <- summary(fit)$coeftable
    bind_rows(
      if ("shock_cons" %in% rownames(ct))
        tibble(h = h, term = "cons",
               est   = ct["shock_cons", "Estimate"],
               se    = ct["shock_cons", "Std. Error"],
               pval  = ct["shock_cons", "Pr(>|t|)"],
               n_obs = nobs(fit)),
      if ("shock_prod" %in% rownames(ct))
        tibble(h = h, term = "prod",
               est   = ct["shock_prod", "Estimate"],
               se    = ct["shock_prod", "Std. Error"],
               pval  = ct["shock_prod", "Pr(>|t|)"],
               n_obs = nobs(fit))
    )
  })
}

# =============================================================================
# Aggregate joint LP
# =============================================================================
message("\n=== Aggregate joint LP ===")
irf_joint <- run_joint_lp(panel)

cat(sprintf("%-14s  %5s  %8s  %6s  %6s\n", "Treatment", "h", "est", "se", "p"))
cat(strrep("-", 48), "\n")
for (tt in c("cons", "prod")) {
  for (h_ in LP_HORIZONS) {
    r <- filter(irf_joint, term == tt, h == h_)
    if (nrow(r) == 0) next
    lab <- if (tt == "cons") "Consumption" else "Production"
    cat(sprintf("%-14s  h=%d  %+8.4f  %6.4f  %6.3f%s\n",
                lab, h_, r$est, r$se, r$pval, stars_fn(r$pval)))
  }
  cat("\n")
}

# =============================================================================
# LaTeX table — aggregate, all horizons
# =============================================================================
make_two_row_joint <- function(term_tag, term_lbl) {
  r <- filter(irf_joint, term == term_tag)
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
  c(paste0(term_lbl, " & ", est_row, " \\\\"),
    paste0("& ", se_row, " \\\\"),
    paste0("& \\multicolumn{5}{l}{$N = ",
           ifelse(length(n_) > 0, n_[1], "---"), "$} \\\\[2pt]"))
}

tex_agg <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Joint consumption and production LP}",
  "\\label{tab:joint-lp}",
  "\\begin{threeparttable}\\small",
  "\\begin{tabular}{lrrrrr}",
  "\\toprule",
  "Exposure & $h=0$ & $h=1$ & $h=2$ & $h=3$ & $h=4$ \\\\",
  "\\midrule",
  make_two_row_joint("cons", "Consumption"),
  make_two_row_joint("prod", "Production"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0("\\item \\textit{Notes:} Consumption and production Bartik shocks estimated jointly ",
         "in the same LP equation with two lags of each. Outcome: ",
         "$\\text{FSI}_{i,t+h} - \\text{FSI}_{i,t-1}$. Year FE; SE clustered by country. ",
         "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."),
  "\\end{tablenotes}\\end{threeparttable}\\end{table}"
)
writeLines(tex_agg, file.path(OUT_TABS, "tab_joint_lp.tex"))
message("  Saved: tab_joint_lp.tex")
message("\nDone.")
