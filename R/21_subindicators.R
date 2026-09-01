# =============================================================================
# 21_subindicators.R
# LP on FSI sub-indicators to identify which fragility dimension food price
# shocks affect most.
#
# The FSI comprises 12 sub-indicators across four dimensions:
#   Cohesion (C): C1 Security apparatus, C2 Factionalized elites,
#                 C3 Group grievance
#   Economic (E): E1 Economy, E2 Economic inequality, E3 Human flight
#   Political (P): P1 State legitimacy, P2 Public services, P3 Human rights
#   Social (S):   S1 Demographic pressures, S2 Refugees/IDPs,
#                 X1 External intervention
#
# We focus on the consumption exposure shock (main identification argument)
# and run the same LP specification as in 08_lp.R for each sub-indicator
# outcome, across the full sample and by income group.
#
# Interpretation guide:
#   E1 (Economy) / E2 (Inequality)  → income / distributional channels
#   P1 (State legitimacy) / P2 (Public services) → state capacity channel
#   C1 (Security apparatus) / C3 (Grievance)     → conflict/rapacity channel
#
# Outputs:
#   output/tables/tab_subindicators.tex    (aggregate, h=0..4)
#   output/tables/tab_subind_income.tex    (by income group, h=1 and h=4)
# =============================================================================

source(here::here("R", "00_setup.R"))

panel <- readRDS(file.path(DATA_PRO, "panel.rds"))
fsi   <- readRDS(file.path(DATA_PRO, "fsi.rds"))

message("Panel: ", nrow(panel), " obs | ", min(panel$year), "–", max(panel$year))
message("FSI sub-indicators: ", paste(grep("^[cepsx][0-9]", names(fsi), value = TRUE), collapse = ", "))

INCOME_LEVELS <- c("Low", "Lower-middle", "Upper-middle", "High")
INCOME_LABELS <- c("Low income", "Lower-middle", "Upper-middle", "High income")

SUB_VARS <- c(
  c1 = "c1_security_apparatus",
  c2 = "c2_factionalized_elites",
  c3 = "c3_group_grievance",
  e1 = "e1_economy",
  e2 = "e2_economic_inequality",
  e3 = "e3_human_flight_and_brain_drain",
  p1 = "p1_state_legitimacy",
  p2 = "p2_public_services",
  p3 = "p3_human_rights",
  s1 = "s1_demographic_pressures",
  s2 = "s2_refugees_and_id_ps",
  x1 = "x1_external_intervention"
)
SUB_VARS <- SUB_VARS[SUB_VARS %in% names(fsi)]

SUB_LABELS <- c(
  c1 = "C1 Security apparatus",
  c2 = "C2 Factionalized elites",
  c3 = "C3 Group grievance",
  e1 = "E1 Economy",
  e2 = "E2 Economic inequality",
  e3 = "E3 Human flight \\& brain drain",
  p1 = "P1 State legitimacy",
  p2 = "P2 Public services",
  p3 = "P3 Human rights",
  s1 = "S1 Demographic pressures",
  s2 = "S2 Refugees \\& IDPs",
  x1 = "X1 External intervention"
)
SUB_LABELS <- SUB_LABELS[names(SUB_VARS)]

# ---- Merge sub-indicators and build forward-difference outcomes --------------
fsi_sub <- fsi |>
  select(iso3c, year, all_of(unname(SUB_VARS))) |>
  filter(!is.na(iso3c))

panel_sub <- panel |>
  left_join(fsi_sub, by = c("iso3c", "year"))

# Build FSI_{t+h} - FSI_{t-1} for each sub-indicator and each horizon
message("Building forward-difference outcomes for ", length(SUB_VARS), " sub-indicators...")
for (sv in unname(SUB_VARS)) {
  panel_sub <- panel_sub |>
    group_by(iso3c) |>
    arrange(year) |>
    mutate(
      !!paste0("lag1_", sv) := lag(.data[[sv]], 1),
      !!paste0("d_", sv, "_h0") := .data[[sv]] - lag(.data[[sv]], 1),
      !!paste0("d_", sv, "_h1") := lead(.data[[sv]], 1) - lag(.data[[sv]], 1),
      !!paste0("d_", sv, "_h2") := lead(.data[[sv]], 2) - lag(.data[[sv]], 1),
      !!paste0("d_", sv, "_h3") := lead(.data[[sv]], 3) - lag(.data[[sv]], 1),
      !!paste0("d_", sv, "_h4") := lead(.data[[sv]], 4) - lag(.data[[sv]], 1)
    ) |>
    ungroup()
}
message("Done building outcomes.")

# ---- LP runner --------------------------------------------------------------
run_lp_sub <- function(data, y_col, shock_var = "shock_cons") {
  lag1 <- paste0(shock_var, "_l1")
  lag2 <- paste0(shock_var, "_l2")
  if (!y_col %in% names(data) || !shock_var %in% names(data)) return(NULL)
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
  tibble(est = ct[".s","Estimate"], se = ct[".s","Std. Error"],
         pval = ct[".s","Pr(>|t|)"], n_obs = nobs(fit))
}

stars_fn <- function(p)
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", p < 0.10 ~ "$\\cdot$", TRUE ~ "")

# =============================================================================
# Section 1: Aggregate LP on all sub-indicators
# =============================================================================
message("\n=== Section 1: Aggregate sub-indicator LP ===")
irf_sub_agg <- map_dfr(names(SUB_VARS), function(sk) {
  sv <- SUB_VARS[sk]
  map_dfr(LP_HORIZONS, function(h) {
    r <- run_lp_sub(panel_sub, paste0("d_", sv, "_h", h))
    if (is.null(r)) return(NULL)
    r |> mutate(h = h, sub_key = sk, sub_label = SUB_LABELS[sk])
  })
})

cat("\nConsumption exposure LP on FSI sub-indicators (h=0 and h=4):\n")
cat(sprintf("%-38s  %5s  %8s  %6s  %6s\n", "Sub-indicator", "h", "est", "se", "p"))
cat(strrep("-", 66), "\n")
for (sk in names(SUB_VARS)) {
  for (h_ in c(0, 2, 4)) {
    r <- filter(irf_sub_agg, sub_key == sk, h == h_)
    if (nrow(r) == 0) next
    cat(sprintf("%-38s  h=%d  %+8.4f  %6.4f  %6.3f%s\n",
                SUB_LABELS[sk], h_, r$est, r$se, r$pval, stars_fn(r$pval)))
  }
  cat("\n")
}

# Memo: FSI total for comparison
cat("FSI total (consumption, h=0 and h=4):\n")
for (h_ in c(0, 4)) {
  r <- run_lp_sub(panel_sub, paste0("dfsi_h", h_))
  if (!is.null(r))
    cat(sprintf("  FSI total  h=%d  %+8.4f  (p=%.3f%s)\n",
                h_, r$est, r$pval, stars_fn(r$pval)))
}

# =============================================================================
# Section 2: LP by income group (consumption, h=1 and h=4)
# =============================================================================
message("\n=== Section 2: Sub-indicator LP by income group ===")
panel_by_inc <- c(
  list(All = panel_sub),
  setNames(map(INCOME_LEVELS, ~ filter(panel_sub, income_group == .x)), INCOME_LEVELS)
)

irf_sub_inc <- map_dfr(names(panel_by_inc), function(ig) {
  d <- panel_by_inc[[ig]]
  message("  ", ig, " (", n_distinct(d$iso3c), " countries)")
  map_dfr(names(SUB_VARS), function(sk) {
    sv <- SUB_VARS[sk]
    map_dfr(c(1, 4), function(h) {
      r <- run_lp_sub(d, paste0("d_", sv, "_h", h))
      if (is.null(r)) return(NULL)
      r |> mutate(h = h, income_group = ig, sub_key = sk, sub_label = SUB_LABELS[sk])
    })
  })
})

# =============================================================================
# LaTeX: aggregate table (h = 0, 2, 4)
# =============================================================================
dim_labels <- c(
  "\\textit{C — Cohesion}", "c1", "c2", "c3",
  "\\textit{E — Economic}", "e1", "e2", "e3",
  "\\textit{P — Political}", "p1", "p2", "p3",
  "\\textit{S/X — Social}", "s1", "s2", "x1"
)

fmt_h <- function(sk, h_) {
  r <- filter(irf_sub_agg, sub_key == sk, h == h_)
  if (nrow(r) == 0) return("---")
  sprintf("$%+.3f%s$", r$est, stars_fn(r$pval))
}
fmt_h_se <- function(sk, h_) {
  r <- filter(irf_sub_agg, sub_key == sk, h == h_)
  if (nrow(r) == 0) return("")
  sprintf("$(%.3f)$", r$se)
}

make_sub_rows <- function(sk) {
  c(
    paste0("\\quad ", SUB_LABELS[sk], " & ",
      paste(sapply(c(0, 2, 4), function(h_) fmt_h(sk, h_)), collapse = " & "),
      " \\\\"),
    paste0("& ",
      paste(sapply(c(0, 2, 4), function(h_) fmt_h_se(sk, h_)), collapse = " & "),
      " \\\\")
  )
}

dim_order <- list(
  C = c("c1","c2","c3"),
  E = c("e1","e2","e3"),
  P = c("p1","p2","p3"),
  S = intersect(c("s1","s2","x1"), names(SUB_VARS))
)
dim_labels_map <- c(C = "C --- Cohesion", E = "E --- Economic",
                    P = "P --- Political", S = "S/X --- Social")

tex_body <- character(0)
for (dim in names(dim_order)) {
  tex_body <- c(tex_body,
    sprintf("\\multicolumn{4}{l}{\\textit{%s}} \\\\", dim_labels_map[dim]))
  for (sk in dim_order[[dim]]) {
    if (sk %in% names(SUB_VARS))
      tex_body <- c(tex_body, make_sub_rows(sk))
  }
  tex_body <- c(tex_body, "\\addlinespace")
}

# N from e1 at h=0 as representative
n_rep <- filter(irf_sub_agg, sub_key == "e1", h == 0)$n_obs
n_rep <- if (length(n_rep) > 0) n_rep[1] else "---"

tex_sub <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{LP on FSI sub-indicators: consumption exposure (full sample)}",
  "\\label{tab:subindicators}",
  "\\begin{threeparttable}\\small",
  "\\begin{tabular}{p{5.5cm} rrr}",
  "\\toprule",
  "Sub-indicator & $h=0$ & $h=2$ & $h=4$ \\\\",
  "\\midrule",
  tex_body,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0("\\item \\textit{Notes:} Each cell shows the LP coefficient ",
         "from regressing the cumulative change in the sub-indicator ",
         "($\\text{Sub}_{i,t+h} - \\text{Sub}_{i,t-1}$) on the consumption-weighted ",
         "shift-share shock, with two shock lags and year FE. ",
         "SE clustered by country. Representative $N \\approx ", n_rep, "$. ",
         "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."),
  "\\end{tablenotes}\\end{threeparttable}\\end{table}"
)
writeLines(tex_sub, file.path(OUT_TABS, "tab_subindicators.tex"))
message("  Saved: tab_subindicators.tex")

# =============================================================================
# LaTeX: income group table — key sub-indicators (e1, p2, c1, c3), h=1 and h=4
# =============================================================================
KEY_SUBS <- intersect(c("e1","p2","c1","c3","p1"), names(SUB_VARS))
inc_groups_all <- c("All", INCOME_LEVELS)

fmt_inc_sub <- function(ig, sk, h_) {
  r <- filter(irf_sub_inc, income_group == ig, sub_key == sk, h == h_)
  if (nrow(r) == 0) return(list(e = "---", s = ""))
  list(e = sprintf("$%+.3f%s$", r$est, stars_fn(r$pval)),
       s = sprintf("$(%.3f)$", r$se))
}

tex_inc <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{LP on key FSI sub-indicators by income group ($h=1$ and $h=4$)}",
  "\\label{tab:subind-income}",
  "\\begin{threeparttable}\\small",
  sprintf("\\begin{tabular}{l%s}", paste(rep("rr", length(inc_groups_all)), collapse="")),
  "\\toprule",
  paste0(" & ",
    paste(sapply(inc_groups_all, function(ig)
      paste0("\\multicolumn{2}{c}{", if (ig == "All") "All" else ig, "}")),
      collapse = " & "), " \\\\"),
  paste0("\\cmidrule(lr){2-3}",
    paste(sprintf("\\cmidrule(lr){%d-%d}",
                  seq(4, 4 + 2*(length(inc_groups_all)-1), 2),
                  seq(5, 5 + 2*(length(inc_groups_all)-1), 2)),
          collapse = "")),
  paste0("Sub-indicator & ",
    paste(rep("$h=1$ & $h=4$", length(inc_groups_all)), collapse = " & "), " \\\\"),
  "\\midrule"
)
for (sk in KEY_SUBS) {
  tex_inc <- c(tex_inc,
    paste0("\\quad ", SUB_LABELS[sk], " &",
      paste(sapply(inc_groups_all, function(ig)
        paste0(fmt_inc_sub(ig, sk, 1)$e, " & ", fmt_inc_sub(ig, sk, 4)$e)),
        collapse = " & "), " \\\\"),
    paste0("&",
      paste(sapply(inc_groups_all, function(ig)
        paste0(fmt_inc_sub(ig, sk, 1)$s, " & ", fmt_inc_sub(ig, sk, 4)$s)),
        collapse = " & "), " \\\\[4pt]")
  )
}
tex_inc <- c(tex_inc,
  "\\bottomrule", "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0("\\item \\textit{Notes:} Consumption exposure LP on sub-indicator outcomes. ",
         "Year FE; SE clustered by country. ",
         "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."),
  "\\end{tablenotes}\\end{threeparttable}\\end{table}"
)
writeLines(tex_inc, file.path(OUT_TABS, "tab_subind_income.tex"))
message("  Saved: tab_subind_income.tex")
message("\nDone.")
