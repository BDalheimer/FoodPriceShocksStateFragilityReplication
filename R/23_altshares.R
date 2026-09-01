# =============================================================================
# 23_altshares.R
# Robustness: alternative consumption exposure using 1999–2009 base-period shares.
#
# The main specification (08_lp.R, 12_group_lp.R) uses caloric shares averaged
# over 2010–2023, which overlap with the outcome sample period (2007–2023).
# Here we use shares averaged over 1999–2009 — fully predetermined relative to
# the estimation sample — strengthening the shares-based identification argument.
#
# Data: data/processed/fao_shares_cons_9909.rds (iso3c × pink_var × grp × s_cons)
#       Same structure as the main shares; built by 02_fao_fbs.R.
#
# Procedure:
#   1. Cross the fixed 9909 shares with all sample years, join Δlog p, sum.
#   2. Compute group-level shocks by summing within commodity group.
#   3. Add lags (1–4).
#   4. Merge with panel outcome variables.
#   5. Run aggregate LP, income-group LP, and commodity-group heatmap LP.
#   6. Produce comparison table: main (2010–2023 shares) vs alt (1999–2009).
#
# Outputs:
#   output/tables/tab_altshares.tex         (aggregate + income, h=1 and h=4)
#   output/tables/tab_altshares_heatmap.tex (group × income, h=4, main vs alt)
#   output/figures/fig_altshares_irf.pdf    (aggregate IRF overlay, cons+prod)
# =============================================================================

source(here::here("R", "00_setup.R"))

panel  <- readRDS(file.path(DATA_PRO, "panel.rds"))
prices <- readRDS(file.path(DATA_PRO, "pinksheet_annual.rds"))
shares_9909 <- readRDS(file.path(DATA_PRO, "fao_shares_cons_9909.rds"))

message("Panel: ",     nrow(panel),        " obs | ", min(panel$year), "–", max(panel$year))
message("Alt shares: ", nrow(shares_9909), " rows | ",
        n_distinct(shares_9909$iso3c), " countries | ",
        n_distinct(shares_9909$pink_var), " commodities")

INCOME_LEVELS <- c("Low", "Lower-middle", "Upper-middle", "High")
INCOME_LABELS <- c("Low income", "Lower-middle", "Upper-middle", "High income")
GRP_LEVELS    <- c("cereals", "oils", "animal", "sugar_fruit", "cash_crops")
GRP_LABELS    <- c("Cereals", "Oils", "Animal prods", "Sugar \\& fruit", "Beverage crops")

CI_Z <- 1.645   # 90% CI for figures

stars_fn <- function(p)
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", p < 0.10 ~ "$\\cdot$", TRUE ~ "")

# =============================================================================
# 1. Build alternative consumption shocks
# =============================================================================
message("\n=== Building alternative shocks (1999–2009 shares) ===")

sample_years <- sort(unique(panel$year))

prices_filt <- prices |>
  filter(!is.na(dlog_p), year %in% sample_years) |>
  select(year, pink_var, dlog_p)

# Aggregate shock: sum over all commodities
shock_alt_agg <- shares_9909 |>
  select(iso3c, pink_var, s_cons) |>
  crossing(year = sample_years) |>
  left_join(prices_filt, by = c("pink_var", "year")) |>
  filter(!is.na(dlog_p)) |>
  group_by(iso3c, year) |>
  summarise(shock_cons_alt = sum(s_cons * dlog_p, na.rm = TRUE),
            .groups = "drop")

message("Aggregate alt shock: ", nrow(shock_alt_agg), " obs | ",
        n_distinct(shock_alt_agg$iso3c), " countries")

# Group-level shocks: one per commodity group
shock_alt_grp <- map_dfr(GRP_LEVELS, function(grp) {
  shares_9909 |>
    filter(grp == .env$grp) |>
    select(iso3c, pink_var, s_cons) |>
    crossing(year = sample_years) |>
    left_join(prices_filt, by = c("pink_var", "year")) |>
    filter(!is.na(dlog_p)) |>
    group_by(iso3c, year) |>
    summarise(shock = sum(s_cons * dlog_p, na.rm = TRUE), .groups = "drop") |>
    mutate(grp = grp)
}) |>
  pivot_wider(names_from = grp, values_from = shock,
              names_prefix = "shock_cons_alt_")

message("Group alt shocks built.")

# =============================================================================
# 2. Merge with panel and add lags
# =============================================================================
panel_alt <- panel |>
  left_join(shock_alt_agg, by = c("iso3c", "year")) |>
  left_join(shock_alt_grp, by = c("iso3c", "year")) |>
  group_by(iso3c) |>
  arrange(year) |>
  mutate(
    shock_cons_alt_l1 = lag(shock_cons_alt, 1),
    shock_cons_alt_l2 = lag(shock_cons_alt, 2),
    across(paste0("shock_cons_alt_", GRP_LEVELS),
           list(l1 = ~ lag(.x, 1), l2 = ~ lag(.x, 2)),
           .names = "{.col}_{.fn}")
  ) |>
  ungroup()

message("Alt shock coverage: ", sum(!is.na(panel_alt$shock_cons_alt)),
        " obs (of ", nrow(panel_alt), " total)")

# =============================================================================
# 3. LP runner
# =============================================================================
run_lp <- function(data, shock_var, horizons = LP_HORIZONS) {
  lag1 <- paste0(shock_var, "_l1")
  lag2 <- paste0(shock_var, "_l2")
  map_dfr(horizons, function(h) {
    y_col <- paste0("dfsi_h", h)
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
    tibble(h     = h,
           est   = ct[".s", "Estimate"],
           se    = ct[".s", "Std. Error"],
           pval  = ct[".s", "Pr(>|t|)"],
           ci_lo = est - CI_Z * se,
           ci_hi = est + CI_Z * se,
           n_obs = nobs(fit))
  })
}

# =============================================================================
# 4. Aggregate LP: main vs alt
# =============================================================================
message("\n=== Section 4: Aggregate LP — main vs alt shares ===")

irf_main_agg <- run_lp(panel_alt, "shock_cons")     |> mutate(spec = "Main (2010–23)")
irf_alt_agg  <- run_lp(panel_alt, "shock_cons_alt") |> mutate(spec = "Alt (1999–09)")
irf_agg_both <- bind_rows(irf_main_agg, irf_alt_agg)

cat("\nAggregate LP comparison (h=0..4):\n")
irf_agg_both |>
  mutate(cell = sprintf("%+.4f (%.4f)%s", est, se, stars_fn(pval))) |>
  select(spec, h, cell) |>
  pivot_wider(names_from = h, values_from = cell) |>
  print()

# =============================================================================
# 5. Income-group LP: main vs alt
# =============================================================================
message("\n=== Section 5: Income-group LP — main vs alt shares ===")

panel_by_inc <- c(
  list(All = panel_alt),
  setNames(map(INCOME_LEVELS, ~ filter(panel_alt, income_group == .x)), INCOME_LEVELS)
)

irf_inc_both <- map_dfr(names(panel_by_inc), function(ig) {
  d <- panel_by_inc[[ig]]
  message("  ", ig, " (", n_distinct(d$iso3c), " countries)")
  bind_rows(
    run_lp(d, "shock_cons")     |> mutate(income_group = ig, spec = "main"),
    run_lp(d, "shock_cons_alt") |> mutate(income_group = ig, spec = "alt")
  )
})

cat("\nIncome LP comparison (h=1 and h=4):\n")
irf_inc_both |>
  filter(h %in% c(1, 4)) |>
  mutate(cell = sprintf("%+.3f (%.3f)%s", est, se, stars_fn(pval))) |>
  select(income_group, spec, h, cell) |>
  pivot_wider(names_from = c(spec, h), values_from = cell) |>
  arrange(factor(income_group, levels = c("All", INCOME_LEVELS))) |>
  print()

# =============================================================================
# 6. Commodity-group heatmap LP (h=4): main vs alt, by income group
# =============================================================================
message("\n=== Section 6: Group LP heatmap — main vs alt shares ===")

GRP_SD_MAIN <- map_dbl(setNames(GRP_LEVELS, GRP_LEVELS),
                        ~ sd(panel_alt[[paste0("shock_cons_", .x)]], na.rm = TRUE))
GRP_SD_ALT  <- map_dbl(setNames(GRP_LEVELS, GRP_LEVELS),
                        ~ sd(panel_alt[[paste0("shock_cons_alt_", .x)]], na.rm = TRUE))

cat("\nGroup shock SDs — main vs alt:\n")
print(data.frame(grp = GRP_LEVELS,
                 sd_main = round(GRP_SD_MAIN, 4),
                 sd_alt  = round(GRP_SD_ALT,  4)))

irf_grp_both <- map_dfr(c("All", INCOME_LEVELS), function(ig) {
  d <- if (ig == "All") panel_alt else filter(panel_alt, income_group == ig)
  if (n_distinct(d$iso3c) < 5) return(NULL)
  map_dfr(GRP_LEVELS, function(grp) {
    bind_rows(
      run_lp(d, paste0("shock_cons_",     grp)) |>
        mutate(income_group = ig, grp = grp, spec = "main"),
      run_lp(d, paste0("shock_cons_alt_", grp)) |>
        mutate(income_group = ig, grp = grp, spec = "alt")
    )
  })
})

cat("\nGroup LP at h=4, All countries (scaled by alt SD):\n")
irf_grp_both |>
  filter(income_group == "All", h == 4) |>
  mutate(est_sc = est * if_else(spec == "main", GRP_SD_MAIN[grp], GRP_SD_ALT[grp]),
         sig    = if_else(pval < 0.05, "*", ""),
         cell   = sprintf("%+.2f%s", est_sc, sig)) |>
  select(spec, grp, cell) |>
  pivot_wider(names_from = grp, values_from = cell) |>
  print()

# =============================================================================
# 7. Figure: aggregate IRF overlay
# =============================================================================
message("\n=== Section 7: Figure ===")

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
      legend.key.width    = unit(20, "pt"),
      legend.text         = element_text(size = base_size - 1)
    )
}

COL_MAIN <- "#2166ac"
COL_ALT  <- "#d6604d"

fig_dat <- irf_agg_both |>
  mutate(
    spec_f = factor(spec, levels = c("Main (2010–23)", "Alt (1999–09)"))
  )

fig_alt <- ggplot(fig_dat, aes(x = h, y = est, colour = spec_f, fill = spec_f,
                                linetype = spec_f)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.35) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c(`Main (2010–23)` = COL_MAIN,
                                  `Alt (1999–09)`  = COL_ALT),
                      name = "Base-period shares") +
  scale_fill_manual(values   = c(`Main (2010–23)` = COL_MAIN,
                                  `Alt (1999–09)`  = COL_ALT),
                    name = "Base-period shares") +
  scale_linetype_manual(values = c(`Main (2010–23)` = "solid",
                                    `Alt (1999–09)`  = "dashed"),
                        name = "Base-period shares") +
  scale_x_continuous(breaks = 0:4, name = "Horizon (years)") +
  scale_y_continuous(name = "FSI change (index points)") +
  ggtitle("Consumption exposure: main vs alternative base-period shares") +
  theme_nature()

ggsave(file.path(OUT_FIGS, "fig_altshares_irf.pdf"),
       fig_alt, width = 120, height = 75, units = "mm")
ggsave(file.path(OUT_FIGS, "fig_altshares_irf.png"),
       fig_alt, width = 120, height = 75, units = "mm", dpi = 300)
message("  Saved: fig_altshares_irf.pdf / .png")

# =============================================================================
# 8. LaTeX: comparison table (aggregate + income group, h=1 and h=4)
# =============================================================================
inc_groups_all <- c("All", INCOME_LEVELS)

fmt_comp <- function(ig, spec_, h_) {
  r <- filter(irf_inc_both, income_group == ig, spec == spec_, h == h_)
  if (nrow(r) == 0) return(list(e = "---", s = ""))
  list(e = sprintf("$%+.3f%s$", r$est, stars_fn(r$pval)),
       s = sprintf("$(%.3f)$", r$se))
}

tex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Consumption LP: main (2010--23 shares) vs alternative (1999--09 shares)}",
  "\\label{tab:altshares}",
  "\\begin{threeparttable}\\small",
  sprintf("\\begin{tabular}{ll%s}", paste(rep("rr", length(inc_groups_all)), collapse="")),
  "\\toprule",
  paste0(" & & ",
    paste(sapply(inc_groups_all, function(ig)
      sprintf("\\multicolumn{2}{c}{%s}", if (ig == "All") "All" else ig)),
      collapse = " & "), " \\\\"),
  paste0("\\cmidrule(lr){3-4}",
    paste(sprintf("\\cmidrule(lr){%d-%d}",
                  seq(5, 5 + 2*(length(inc_groups_all)-1), 2),
                  seq(6, 6 + 2*(length(inc_groups_all)-1), 2)),
          collapse = "")),
  paste0("Specification & & ",
    paste(rep("$h=1$ & $h=4$", length(inc_groups_all)), collapse = " & "), " \\\\"),
  "\\midrule"
)

for (spec_ in c("main", "alt")) {
  lbl <- if (spec_ == "main") "Main (2010--23)" else "Alt (1999--09)"
  tex <- c(tex,
    paste0(lbl, " &",
      paste(sapply(inc_groups_all, function(ig)
        paste0(fmt_comp(ig, spec_, 1)$e, " & ", fmt_comp(ig, spec_, 4)$e)),
        collapse = " & "), " \\\\"),
    paste0("&",
      paste(sapply(inc_groups_all, function(ig)
        paste0(fmt_comp(ig, spec_, 1)$s, " & ", fmt_comp(ig, spec_, 4)$s)),
        collapse = " & "), " \\\\[4pt]")
  )
}

n_main <- filter(irf_inc_both, income_group == "All", spec == "main", h == 0)$n_obs
n_alt  <- filter(irf_inc_both, income_group == "All", spec == "alt",  h == 0)$n_obs

tex <- c(tex,
  "\\midrule",
  paste0("$N$ (main) & & \\multicolumn{", 2*length(inc_groups_all), "}{l}{",
         ifelse(length(n_main) > 0, n_main[1], "---"), "} \\\\"),
  paste0("$N$ (alt)  & & \\multicolumn{", 2*length(inc_groups_all), "}{l}{",
         ifelse(length(n_alt)  > 0, n_alt[1],  "---"), "} \\\\"),
  "\\bottomrule", "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0("\\item \\textit{Notes:} Both rows estimate the same LP specification ",
         "$\\Delta_{h}\\text{FSI}_{it} = \\beta x_{it} + \\sum_{k=1}^{2}\\gamma_k x_{i,t-k} + \\alpha_t + \\varepsilon_{iht}$ ",
         "but differ in how the consumption-weighted Bartik shock is constructed. ",
         "\\textit{Main}: caloric shares averaged over 2010--2023 (overlap with outcome sample). ",
         "\\textit{Alt}: caloric shares averaged over 1999--2009 (fully predetermined). ",
         "Year FE; SE clustered by country. ",
         "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."),
  "\\end{tablenotes}\\end{threeparttable}\\end{table}"
)
writeLines(tex, file.path(OUT_TABS, "tab_altshares.tex"))
message("  Saved: tab_altshares.tex")

# =============================================================================
# 9. LaTeX: group heatmap comparison (h=4, INCOME_LEVELS, main vs alt)
# =============================================================================
heat_dat <- irf_grp_both |>
  filter(income_group %in% INCOME_LEVELS, h == 4) |>
  mutate(
    grp_sd    = if_else(spec == "main", GRP_SD_MAIN[grp], GRP_SD_ALT[grp]),
    est_sc    = est * grp_sd,
    sig       = if_else(pval < 0.05, "*", ""),
    cell      = sprintf("%+.2f%s", est_sc, sig)
  )

fmt_heat <- function(ig, grp_, spec_) {
  r <- filter(heat_dat, income_group == ig, grp == grp_, spec == spec_)
  if (nrow(r) == 0) return("---")
  r$cell[1]
}

tex_heat <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Group LP heatmap ($h=4$): main vs alternative shares, scaled by group shock SD}",
  "\\label{tab:altshares-heatmap}",
  "\\begin{threeparttable}\\small",
  sprintf("\\begin{tabular}{ll%s}", paste(rep("r", length(GRP_LEVELS)), collapse="")),
  "\\toprule",
  paste0(" & & ", paste(GRP_LABELS, collapse = " & "), " \\\\"),
  "\\midrule"
)
for (ig in INCOME_LEVELS) {
  tex_heat <- c(tex_heat,
    paste0("\\multirow{2}{*}{", ig, "} & Main &",
      paste(sapply(GRP_LEVELS, function(g) fmt_heat(ig, g, "main")),
            collapse = " & "), " \\\\"),
    paste0("& Alt &",
      paste(sapply(GRP_LEVELS, function(g) fmt_heat(ig, g, "alt")),
            collapse = " & "), " \\\\[3pt]")
  )
}
tex_heat <- c(tex_heat,
  "\\bottomrule", "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0("\\item \\textit{Notes:} Each cell reports LP coefficient $\\times$ SD ",
         "of the group shock (FSI points per 1-SD shock). ",
         "Main: 2010--23 shares. Alt: 1999--09 shares. $^{*}p<0.05$."),
  "\\end{tablenotes}\\end{threeparttable}\\end{table}"
)
writeLines(tex_heat, file.path(OUT_TABS, "tab_altshares_heatmap.tex"))
message("  Saved: tab_altshares_heatmap.tex")
message("\nDone.")
