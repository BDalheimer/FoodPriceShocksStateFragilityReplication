# =============================================================================
# 22_agemployment.R
# Agricultural employment interaction: tests whether the opportunity cost
# channel scales with the share of the population working in agriculture.
#
# The opportunity cost mechanism predicts that food price rises stabilize
# most strongly where more people earn income from agriculture. We test this
# by interacting the consumption and production shocks with a country's
# agricultural employment share.
#
# Data: rural_pct (rural population share, %) from WDI is in the panel as a
# proxy. We use it because agricultural employment data (WDI SL.AGR.EMPL.ZS)
# is available for fewer countries in the sample period.
#
# Two complementary approaches:
#   1. Sample split: bottom tercile (low rural/ag) vs top tercile (high rural/ag)
#      separately estimated LP, reported as IRF overlays
#   2. Continuous interaction: shock × rural_pct (demeaned) in a single LP
#      — identifies whether the slope of the IRF varies linearly with ag exposure
#
# Both approaches use consumption and production shocks.
#
# Outputs:
#   output/tables/tab_agemployment_split.tex   (split-sample LP, h=1 and h=4)
#   output/tables/tab_agemployment_interact.tex (interaction LP, h=1 and h=4)
#   output/figures/fig_agemployment_irf.pdf     (IRF overlay by tercile)
# =============================================================================

source(here::here("R", "00_setup.R"))

panel <- readRDS(file.path(DATA_PRO, "panel.rds"))

message("Panel: ", nrow(panel), " obs | ", min(panel$year), "–", max(panel$year))
message("rural_pct coverage: ", sum(!is.na(panel$rural_pct)), " obs")

INCOME_LEVELS <- c("Low", "Lower-middle", "Upper-middle", "High")

COL_LO  <- "#2166ac"   # low rural — blue
COL_MED <- "#4dac26"   # medium rural — green
COL_HI  <- "#d6604d"   # high rural — red
CI_Z    <- 1.645       # 90% CI

stars_fn <- function(p)
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", p < 0.10 ~ "$\\cdot$", TRUE ~ "")

# ---- Assign tercile based on country-level median rural_pct -----------------
# Use median over sample period per country so split is time-stable
rural_median <- panel |>
  group_by(iso3c) |>
  summarise(rural_med = median(rural_pct, na.rm = TRUE), .groups = "drop") |>
  filter(!is.na(rural_med))

breaks <- quantile(rural_median$rural_med, c(1/3, 2/3), na.rm = TRUE)
rural_median <- rural_median |>
  mutate(rural_tercile = cut(rural_med, breaks = c(-Inf, breaks, Inf),
                             labels = c("Low rural", "Mid rural", "High rural")))

panel_ag <- panel |>
  left_join(rural_median, by = "iso3c") |>
  filter(!is.na(rural_med)) |>
  group_by(iso3c) |>
  mutate(rural_dm = rural_pct - mean(rural_pct, na.rm = TRUE)) |>
  ungroup()

cat("Tercile distribution:\n")
print(table(rural_median$rural_tercile))
cat("Median rural_pct by tercile:\n")
rural_median |>
  group_by(rural_tercile) |>
  summarise(median_rural = median(rural_med), .groups = "drop") |>
  print()

# ---- Standard LP runner -----------------------------------------------------
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
           ci_lo = est - CI_Z * se,
           ci_hi = est + CI_Z * se,
           n_obs = nobs(fit),
           n_cty = n_distinct(d$iso3c))
  })
}

# =============================================================================
# Section 1: Split-sample LP by rural tercile
# =============================================================================
message("\n=== Section 1: Split-sample LP by rural tercile ===")

TERCILE_LEVELS <- c("Low rural", "Mid rural", "High rural")
panel_by_tercile <- setNames(
  map(TERCILE_LEVELS, ~ filter(panel_ag, rural_tercile == .x)),
  TERCILE_LEVELS
)
walk2(names(panel_by_tercile), panel_by_tercile,
      ~ message("  ", .x, ": ", n_distinct(.y$iso3c), " countries"))

irf_split <- map_dfr(TERCILE_LEVELS, function(terc) {
  d <- panel_by_tercile[[terc]]
  map_dfr(c("cons", "prod"), function(tt) {
    run_lp(d, paste0("shock_", tt)) |>
      mutate(tercile = terc, treat = tt)
  })
})

cat("\nSplit-sample LP (h=1 and h=4):\n")
irf_split |>
  filter(h %in% c(1, 4)) |>
  mutate(cell = sprintf("%+.3f (%.3f)%s", est, se, stars_fn(pval))) |>
  select(tercile, treat, h, cell) |>
  pivot_wider(names_from = h, values_from = cell) |>
  arrange(tercile, treat) |>
  print()

# =============================================================================
# Section 2: Continuous interaction LP
# shock × rural_dm (demeaned within country)
# =============================================================================
message("\n=== Section 2: Interaction LP (shock × rural_pct demeaned) ===")

run_interact_lp <- function(data, shock_var, horizons = LP_HORIZONS) {
  lag1 <- paste0(shock_var, "_l1")
  lag2 <- paste0(shock_var, "_l2")
  map_dfr(horizons, function(h) {
    y_col <- paste0("dfsi_h", h)
    if (!y_col %in% names(data)) return(NULL)
    d <- data |>
      rename(.y = all_of(y_col), .s = all_of(shock_var)) |>
      filter(!is.na(.y), !is.na(.s), !is.na(rural_dm)) |>
      mutate(.sx = .s * rural_dm)
    if (lag1 %in% names(d)) d <- rename(d, .sl1 = all_of(lag1))
    if (lag2 %in% names(d)) d <- rename(d, .sl2 = all_of(lag2))
    rhs_lags <- intersect(c(".sl1", ".sl2"), names(d))
    fml <- as.formula(paste0(".y ~ .s + .sx + rural_dm",
      if (length(rhs_lags)) paste0(" + ", paste(rhs_lags, collapse = " + ")),
      " | year"))
    fit <- tryCatch(
      feols(fml, data = d, cluster = ~iso3c, warn = FALSE, notes = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    ct <- summary(fit)$coeftable
    get_coef <- function(v) {
      if (!v %in% rownames(ct)) return(tibble(est = NA, se = NA, pval = NA))
      tibble(est = ct[v,"Estimate"], se = ct[v,"Std. Error"], pval = ct[v,"Pr(>|t|)"])
    }
    main  <- get_coef(".s")
    inter <- get_coef(".sx")
    tibble(h = h,
           est_main  = main$est,  se_main  = main$se,  pval_main  = main$pval,
           est_inter = inter$est, se_inter = inter$se, pval_inter = inter$pval,
           n_obs = nobs(fit), n_cty = n_distinct(d$iso3c))
  })
}

irf_interact <- map_dfr(c("cons", "prod"), function(tt) {
  run_interact_lp(panel_ag, paste0("shock_", tt)) |> mutate(treat = tt)
})

cat("\nInteraction LP — main effect and interaction (h=1 and h=4):\n")
irf_interact |>
  filter(h %in% c(1, 4)) |>
  mutate(
    main  = sprintf("%+.4f (%.4f)%s", est_main,  se_main,  stars_fn(pval_main)),
    inter = sprintf("%+.4f (%.4f)%s", est_inter, se_inter, stars_fn(pval_inter))
  ) |>
  select(treat, h, main, inter, n_obs, n_cty) |>
  print()

# =============================================================================
# Figure: IRF overlay by tercile (cons and prod)
# =============================================================================
message("\n=== Section 3: IRF figure ===")

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
      legend.position     = "bottom",
      legend.key.width    = unit(20, "pt"),
      legend.text         = element_text(size = base_size - 1)
    )
}

TERC_COLS   <- c(`Low rural` = COL_LO, `Mid rural` = COL_MED, `High rural` = COL_HI)
TERC_LTYPE  <- c(`Low rural` = "dashed", `Mid rural` = "dotted", `High rural` = "solid")

fig_dat <- irf_split |>
  filter(treat %in% c("cons", "prod")) |>
  mutate(
    treat_f  = factor(treat, levels = c("cons","prod"),
                      labels = c("Consumption exposure","Production exposure")),
    tercile_f = factor(tercile, levels = TERCILE_LEVELS)
  )

fig_ag <- ggplot(fig_dat, aes(x = h, y = est, colour = tercile_f, fill = tercile_f,
                               linetype = tercile_f)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.35) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  scale_colour_manual(values = TERC_COLS, name = "Rural population tercile") +
  scale_fill_manual(values   = TERC_COLS, name = "Rural population tercile") +
  scale_linetype_manual(values = TERC_LTYPE, name = "Rural population tercile") +
  scale_x_continuous(breaks = 0:4, name = "Horizon (years)") +
  scale_y_continuous(name = "FSI change (index points)") +
  facet_wrap(~ treat_f, nrow = 1, scales = "free_y") +
  theme_nature() +
  theme(strip.text.y = element_text(angle = 0, hjust = 0))

ggsave(file.path(OUT_FIGS, "fig_agemployment_irf.pdf"),
       fig_ag, width = 180, height = 80, units = "mm")
ggsave(file.path(OUT_FIGS, "fig_agemployment_irf.png"),
       fig_ag, width = 180, height = 80, units = "mm", dpi = 300)
message("  Saved: fig_agemployment_irf.pdf / .png")

# =============================================================================
# LaTeX: split-sample table
# =============================================================================
fmt_cell <- function(terc, tt, h_) {
  r <- filter(irf_split, tercile == terc, treat == tt, h == h_)
  if (nrow(r) == 0) return(list(e = "---", s = ""))
  list(e = sprintf("$%+.3f%s$", r$est, stars_fn(r$pval)),
       s = sprintf("$(%.3f)$", r$se))
}

TERC_LABELS_TEX <- c(
  `Low rural`  = "Low rural",
  `Mid rural`  = "Mid rural",
  `High rural` = "High rural"
)

tex_split <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Food price LP by rural population tercile ($h = 1$ and $h = 4$)}",
  "\\label{tab:ag-split}",
  "\\begin{threeparttable}\\small",
  "\\begin{tabular}{llrrrrrr}",
  "\\toprule",
  paste0(" & & \\multicolumn{2}{c}{Low rural} & \\multicolumn{2}{c}{Mid rural}",
         " & \\multicolumn{2}{c}{High rural} \\\\"),
  "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}",
  "Exposure & & $h=1$ & $h=4$ & $h=1$ & $h=4$ & $h=1$ & $h=4$ \\\\",
  "\\midrule"
)
for (tt in c("cons", "prod")) {
  lbl <- if (tt == "cons") "Consumption" else "Production"
  tex_split <- c(tex_split,
    paste0(lbl, " &",
      paste(sapply(TERCILE_LEVELS, function(terc)
        paste0(fmt_cell(terc, tt, 1)$e, " & ", fmt_cell(terc, tt, 4)$e)),
        collapse = " & "), " \\\\"),
    paste0("&",
      paste(sapply(TERCILE_LEVELS, function(terc)
        paste0(fmt_cell(terc, tt, 1)$s, " & ", fmt_cell(terc, tt, 4)$s)),
        collapse = " & "), " \\\\[4pt]")
  )
}
n_cty <- irf_split |>
  filter(h == 0) |>
  group_by(tercile) |>
  summarise(n = first(n_cty), .groups = "drop")
n_row <- paste(sapply(TERCILE_LEVELS, function(t) {
  n <- filter(n_cty, tercile == t)$n
  paste0("\\multicolumn{2}{c}{$N_{c}=", ifelse(length(n) > 0, n[1], "---"), "$}")
}), collapse = " & ")

tex_split <- c(tex_split,
  "\\midrule",
  paste0("Countries & & ", n_row, " \\\\"),
  "\\bottomrule", "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0("\\item \\textit{Notes:} Countries split by terciles of median rural ",
         "population share (World Bank WDI, SP.RUR.TOTL.ZS). Tercile thresholds: ",
         sprintf("%.1f\\%% and %.1f\\%%", breaks[1], breaks[2]), ". ",
         "LP with year FE and two shock lags; SE clustered by country. ",
         "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."),
  "\\end{tablenotes}\\end{threeparttable}\\end{table}"
)
writeLines(tex_split, file.path(OUT_TABS, "tab_agemployment_split.tex"))
message("  Saved: tab_agemployment_split.tex")

# =============================================================================
# LaTeX: interaction table
# =============================================================================
fmt_int <- function(tt, term, h_) {
  r <- filter(irf_interact, treat == tt, h == h_)
  if (nrow(r) == 0) return(list(e = "---", s = ""))
  est_ <- if (term == "main") r$est_main else r$est_inter
  se_  <- if (term == "main") r$se_main  else r$se_inter
  pval_<- if (term == "main") r$pval_main else r$pval_inter
  list(e = sprintf("$%+.4f%s$", est_, stars_fn(pval_)),
       s = sprintf("$(%.4f)$", se_))
}

tex_int <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Food price LP with rural exposure interaction ($h = 1$ and $h = 4$)}",
  "\\label{tab:ag-interact}",
  "\\begin{threeparttable}\\small",
  "\\begin{tabular}{llrrrr}",
  "\\toprule",
  paste0(" & & \\multicolumn{2}{c}{Consumption} & \\multicolumn{2}{c}{Production} \\\\"),
  "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
  "Term & & $h=1$ & $h=4$ & $h=1$ & $h=4$ \\\\",
  "\\midrule"
)
for (term in c("main","inter")) {
  lbl <- if (term == "main") "Shock ($\\beta$)" else "Shock $\\times$ Rural (\\%)"
  tex_int <- c(tex_int,
    paste0(lbl, " &",
      paste(sapply(c("cons","prod"), function(tt)
        paste0(fmt_int(tt, term, 1)$e, " & ", fmt_int(tt, term, 4)$e)),
        collapse = " & "), " \\\\"),
    paste0("&",
      paste(sapply(c("cons","prod"), function(tt)
        paste0(fmt_int(tt, term, 1)$s, " & ", fmt_int(tt, term, 4)$s)),
        collapse = " & "), " \\\\[4pt]")
  )
}
n_int <- filter(irf_interact, treat == "cons", h == 0)$n_obs
tex_int <- c(tex_int,
  "\\midrule",
  paste0("$N$ & & \\multicolumn{2}{c}{", ifelse(length(n_int)>0,n_int[1],"---"),
         "} & \\multicolumn{2}{c}{",
         ifelse(length(filter(irf_interact,treat=="prod",h==0)$n_obs)>0,
                filter(irf_interact,treat=="prod",h==0)$n_obs[1],"---"), "} \\\\"),
  "\\bottomrule", "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0("\\item \\textit{Notes:} The model adds shock $\\times$ rural population share ",
         "(country-demeaned) as an interaction. A negative interaction coefficient for ",
         "consumption indicates that the stabilising effect is stronger in more rural ",
         "(agricultural) economies, consistent with the opportunity cost channel. ",
         "Year FE; SE clustered by country. ",
         "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."),
  "\\end{tablenotes}\\end{threeparttable}\\end{table}"
)
writeLines(tex_int, file.path(OUT_TABS, "tab_agemployment_interact.tex"))
message("  Saved: tab_agemployment_interact.tex")
message("\nDone.")
