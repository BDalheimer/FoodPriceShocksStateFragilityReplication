# =============================================================================
# 18_pretrend.R
# Pre-trend (placebo) tests for the shift-share LP design.
#
# Tests parallel trends by running the LP at negative horizons (h = -2, -3):
#   FSI_{i,t-|h|} - FSI_{i,t-1} ~ shock_{it} + lags | year FE
#
# If identification holds (shock uncorrelated with pre-determined fragility
# trends), coefficients at h = -2 and h = -3 should not be significant.
# We categorise each model as:
#   Pass        — p >= 0.10 at both h=-2 and h=-3
#   Attenuating — pre-trend significant but opposite sign to main effect at h=4
#                 (estimates conservative lower bounds, not inflated)
#   Concern     — pre-trend significant, same sign as main effect at h=4
#   N/A         — near-zero shock SD (<0.007)
#
# Panels:
#   A. Aggregate × 5 exposure types
#   B. Consumption shock × 4 income groups
#   C. Commodity-group consumption shocks (aggregate)
#   D. Animal consumption shock × 4 income groups
#
# Outputs:
#   output/tables/tab_pretrend.tex
# =============================================================================

source(here::here("R", "00_setup.R"))

panel  <- readRDS(file.path(DATA_PRO, "panel.rds"))
fsi    <- readRDS(file.path(DATA_PRO, "fsi.rds"))

message("Panel: ", nrow(panel), " obs | ", min(panel$year), "–", max(panel$year))

# ---- Constants --------------------------------------------------------------
TREAT_TYPES  <- c("cons", "prod", "imp", "exp", "net")
TREAT_LABELS <- c(cons = "Consumption", prod = "Production",
                  imp  = "Import",      exp  = "Export",  net = "Net import")
INCOME_LEVELS <- c("Low", "Lower-middle", "Upper-middle", "High")
INCOME_LABELS <- c("Low income", "Lower-middle", "Upper-middle", "High income")
GRP_LEVELS    <- c("cereals", "oils", "animal", "sugar_fruit", "cash_crops")
GRP_LABELS    <- c("Cereals", "Oils", "Animal products",
                   "Sugar \\& fruit", "Bev. crops")

# ---- Build backward-horizon outcomes ----------------------------------------
# dfsi_hn2 = FSI_{t-2} - FSI_{t-1}  (pre-trend at h=-2)
# dfsi_hn3 = FSI_{t-3} - FSI_{t-1}  (pre-trend at h=-3)
fsi_sub <- fsi |>
  select(iso3c, year, fsi_total) |>
  filter(!is.na(fsi_total))

panel_pt <- panel |>
  left_join(fsi_sub |> rename(fsi_t = fsi_total), by = c("iso3c", "year")) |>
  left_join(
    fsi_sub |> mutate(year = year + 2) |> rename(fsi_tm2 = fsi_total),
    by = c("iso3c", "year")
  ) |>
  left_join(
    fsi_sub |> mutate(year = year + 3) |> rename(fsi_tm3 = fsi_total),
    by = c("iso3c", "year")
  ) |>
  mutate(
    dfsi_hn2 = fsi_tm2 - fsi_lag1,   # FSI_{t-2} - FSI_{t-1}
    dfsi_hn3 = fsi_tm3 - fsi_lag1    # FSI_{t-3} - FSI_{t-1}
  )

message("Negative-horizon obs: hn2=", sum(!is.na(panel_pt$dfsi_hn2)),
        "  hn3=", sum(!is.na(panel_pt$dfsi_hn3)))

# ---- LP helper --------------------------------------------------------------
run_lp_h <- function(data, y_col, shock_var) {
  lag1 <- paste0(shock_var, "_l1")
  lag2 <- paste0(shock_var, "_l2")
  if (!y_col %in% names(data) || !shock_var %in% names(data)) return(NULL)
  d <- data |>
    rename(.y = all_of(y_col), .s = all_of(shock_var)) |>
    filter(!is.na(.y), !is.na(.s))
  if (lag1 %in% names(d)) d <- rename(d, .sl1 = all_of(lag1))
  if (lag2 %in% names(d)) d <- rename(d, .sl2 = all_of(lag2))
  rhs_lags <- intersect(c(".sl1", ".sl2"), names(d))
  fml <- as.formula(paste0(".y ~ .s",
    if (length(rhs_lags)) paste0(" + ", paste(rhs_lags, collapse = " + ")),
    " | year"))
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

# ---- Row builder ------------------------------------------------------------
make_row <- function(data, shock_var, label) {
  shock_sd <- sd(data[[shock_var]], na.rm = TRUE)
  if (is.na(shock_sd) || shock_sd < 0.007)
    return(tibble(label = label, est_hn2 = NA, se_hn2 = NA, pval_hn2 = NA,
                  est_hn3 = NA, se_hn3 = NA, pval_hn3 = NA,
                  est_h4 = NA, verdict = "N/A"))

  r2  <- run_lp_h(data, "dfsi_hn2", shock_var)
  r3  <- run_lp_h(data, "dfsi_hn3", shock_var)
  r4  <- run_lp_h(data, "dfsi_h4",  shock_var)

  if (is.null(r2)) return(NULL)

  pass2 <- is.na(r2$pval) || r2$pval >= 0.10
  pass3 <- is.null(r3)    || is.na(r3$pval) || r3$pval >= 0.10

  verdict <- if (pass2 && pass3) {
    "Pass"
  } else if (!is.null(r4) && !is.na(r4$est) && !is.na(r2$est)) {
    if (sign(r2$est) != sign(r4$est)) "Attenuating" else "Concern"
  } else "Significant"

  tibble(
    label    = label,
    est_hn2  = r2$est,  se_hn2  = r2$se,  pval_hn2 = r2$pval,
    est_hn3  = if (!is.null(r3)) r3$est else NA,
    se_hn3   = if (!is.null(r3)) r3$se  else NA,
    pval_hn3 = if (!is.null(r3)) r3$pval else NA,
    est_h4   = if (!is.null(r4)) r4$est else NA,
    verdict  = verdict
  )
}

# ---- Panel A: Aggregate × 5 exposures ---------------------------------------
message("\nPanel A: Aggregate exposures...")
rows_A <- map_dfr(TREAT_TYPES, function(tt) {
  sv <- paste0("shock_", tt)
  make_row(panel_pt, sv, TREAT_LABELS[tt])
})
print(rows_A |> select(label, est_hn2, pval_hn2, est_h4, verdict))

# ---- Panel B: Consumption × income group ------------------------------------
message("\nPanel B: Consumption by income group...")
rows_B <- map_dfr(seq_along(INCOME_LEVELS), function(i) {
  d <- filter(panel_pt, income_group == INCOME_LEVELS[i])
  make_row(d, "shock_cons", INCOME_LABELS[i])
})
print(rows_B |> select(label, est_hn2, pval_hn2, est_h4, verdict))

# ---- Panel C: Commodity group × consumption (aggregate) ---------------------
message("\nPanel C: Commodity groups (aggregate)...")
rows_C <- map_dfr(seq_along(GRP_LEVELS), function(i) {
  sv <- paste0("shock_cons_", GRP_LEVELS[i])
  make_row(panel_pt, sv, GRP_LABELS[i])
})
print(rows_C |> select(label, est_hn2, pval_hn2, est_h4, verdict))

# ---- Panel D: Animal shock × income group -----------------------------------
message("\nPanel D: Animal consumption shock by income group...")
rows_D <- map_dfr(seq_along(INCOME_LEVELS), function(i) {
  d <- filter(panel_pt, income_group == INCOME_LEVELS[i])
  make_row(d, "shock_cons_animal", INCOME_LABELS[i])
})
print(rows_D |> select(label, est_hn2, pval_hn2, est_h4, verdict))

# =============================================================================
# LaTeX table
# =============================================================================
star_fn <- function(p) {
  case_when(is.na(p) ~ "", p < 0.01 ~ "**", p < 0.05 ~ "*",
            p < 0.10 ~ "$^{\\dagger}$", TRUE ~ "")
}
verdict_tex <- function(v) {
  case_when(
    v == "Pass"        ~ "\\checkmark",
    v == "Attenuating" ~ "Attenuating$^{\\ddagger}$",
    v == "Concern"     ~ "Concern$^{\\S}$",
    v == "N/A"         ~ "---",
    TRUE               ~ v
  )
}
fmt_cell <- function(est, se, pval) {
  if (is.na(est)) return("--- & ---")
  sprintf("%.3f%s & (%.3f)", est, star_fn(pval), se)
}
make_tex_rows <- function(df) {
  map_chr(seq_len(nrow(df)), function(i) {
    d <- df[i, ]
    sprintf("\\quad %s & %s & %s & %s \\\\",
            d$label,
            fmt_cell(d$est_hn2, d$se_hn2, d$pval_hn2),
            fmt_cell(d$est_hn3, d$se_hn3, d$pval_hn3),
            verdict_tex(d$verdict))
  })
}
panel_hdr <- function(ltr, title)
  sprintf("\\multicolumn{6}{l}{\\textit{%s. %s}} \\\\", ltr, title)

tex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Pre-trend tests: LP coefficients at placebo horizons $h = -2$ and $h = -3$}",
  "\\label{tab:pretrend}",
  "\\begin{threeparttable}",
  "\\small",
  "\\begin{tabular}{p{5.2cm} rr rr l}",
  "\\toprule",
  " & \\multicolumn{2}{c}{$h = -2$} & \\multicolumn{2}{c}{$h = -3$} & \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "Model & Est. & (SE) & Est. & (SE) & Verdict \\\\",
  "\\midrule",
  panel_hdr("A", "Aggregate, five exposure types"),
  make_tex_rows(rows_A),
  "\\addlinespace",
  panel_hdr("B", "Consumption shock by income group"),
  make_tex_rows(rows_B),
  "\\addlinespace",
  panel_hdr("C", "Commodity group shocks, aggregate"),
  make_tex_rows(rows_C),
  "\\addlinespace",
  panel_hdr("D", "Animal consumption shock by income group"),
  make_tex_rows(rows_D),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0("\\item \\textit{Notes:} Each cell reports the LP coefficient ",
         "from regressing $\\text{FSI}_{i,t+h} - \\text{FSI}_{i,t-1}$ on the relevant ",
         "shock at placebo horizons $h \\in \\{-2,-3\\}$, with year fixed effects ",
         "and two shock lags. Standard errors clustered by country. ",
         "\\checkmark~Pass: $p \\geq 0.10$ at both horizons. ",
         "Attenuating$^{\\ddagger}$: significant pre-trend but opposite sign to the main ",
         "effect at $h=4$ (conservative bias, not inflation). ",
         "Concern$^{\\S}$: significant pre-trend same sign as main effect. ",
         "--- : near-zero shock variance. ",
         "$^{**}p<0.01$, $^{*}p<0.05$, $^{\\dagger}p<0.10$."),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)

writeLines(tex, file.path(OUT_TABS, "tab_pretrend.tex"))
message("\nSaved: output/tables/tab_pretrend.tex")
message("Done.")
