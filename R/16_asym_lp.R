# =============================================================================
# 16_asym_lp.R
# Asymmetric LP by commodity group × income group
#
# Splits each group shock into positive (price surges) and negative (price
# crashes) components and estimates their effects jointly:
#
#   ΔFSI_{i,t+h} = α_t + β⁺ shock_pos_{i,t} + β⁻ shock_neg_{i,t}
#                       + γ₁⁺ shock_pos_{i,t-1} + γ₁⁻ shock_neg_{i,t-1}
#                       + γ₂⁺ shock_pos_{i,t-2} + γ₂⁻ shock_neg_{i,t-2}
#                       + ε_{iht}
#
# Sign convention in heatmap:
#   Surge  panel  → β⁺ × SD(group shock)      positive = destabilising
#   Crash  panel  → −β⁻ × SD(group shock)     positive = destabilising
#   (shock_neg ≤ 0, so multiplying by −1 aligns the colour direction)
#
# Significance markers: p < 0.05 at h = 4.
#
# Outputs:
#   output/figures/fig9_asym_heatmap.pdf/.png
#     — 2×2 facet_grid(direction ~ exposure), cons + prod, h=4
#   data/processed/asym_lp_results.rds
# =============================================================================

source(here::here("R", "00_setup.R"))

panel <- readRDS(file.path(DATA_PRO, "panel.rds"))

message("Panel: ", nrow(panel), " obs | ",
        n_distinct(panel$iso3c), " countries | ",
        min(panel$year), "–", max(panel$year))

# ---- Constants --------------------------------------------------------------
TREAT_TYPES  <- c("cons", "prod", "imp", "exp", "net")
TREAT_LABELS <- c(cons = "Consumption", prod = "Production",
                  imp  = "Imports",     exp  = "Exports",
                  net  = "Net imports")

GRP_LEVELS <- c("cereals", "oils", "animal", "sugar_fruit", "cash_crops")
GRP_LABELS <- c("Cereals", "Oils", "Animal\nprods", "Sugar &\nfruit", "Beverage\ncrops")

INCOME_LEVELS <- c("Low", "Lower-middle", "Upper-middle", "High")
INCOME_LABELS <- c("Low\nincome", "Lower-middle\nincome",
                   "Upper-middle\nincome", "High\nincome")

CI_Z     <- 1.96    # 95% CI for asymmetric models (matches original paper)
HEAT_LIM <- 2

# ---- Group shock SDs --------------------------------------------------------
GRP_SD <- map_dbl(setNames(GRP_LEVELS, GRP_LEVELS),
                  ~ sd(panel[[paste0("shock_cons_", .x)]], na.rm = TRUE))
message("Group shock SDs: ")
print(round(GRP_SD, 4))

# ---- Verify columns ---------------------------------------------------------
grp_shock_cols <- as.vector(outer(paste0("shock_", TREAT_TYPES, "_"), GRP_LEVELS, paste0))
missing <- setdiff(grp_shock_cols, names(panel))
if (length(missing) > 0)
  stop("Missing group shock columns — run 05_group_shocks.R:\n  ",
       paste(missing, collapse = ", "))

# =============================================================================
# Asymmetric LP runner
# Takes data subset and a shock column name; computes pos/neg splits and
# their lags internally, then runs the joint LP for each horizon.
# Returns one row per (h, component) where component ∈ {pos, neg}.
# =============================================================================
run_lp_asym <- function(data, shock_var, fd_stem = "dfsi", horizons = LP_HORIZONS) {
  if (n_distinct(data$iso3c) < 5) return(NULL)

  d_prep <- data |>
    arrange(iso3c, year) |>
    group_by(iso3c) |>
    mutate(
      .sp   = pmax(.data[[shock_var]], 0),
      .sn   = pmin(.data[[shock_var]], 0),
      .sp_l1 = lag(.sp, 1), .sp_l2 = lag(.sp, 2),
      .sn_l1 = lag(.sn, 1), .sn_l2 = lag(.sn, 2)
    ) |>
    ungroup()

  map_dfr(horizons, function(h) {
    y_col <- paste0(fd_stem, "_h", h)
    if (!y_col %in% names(d_prep)) return(NULL)
    d <- d_prep |>
      rename(.y = all_of(y_col)) |>
      filter(!is.na(.y), !is.na(.sp), !is.na(.sn))
    if (nrow(d) < 20) return(NULL)

    fit <- tryCatch(
      feols(.y ~ .sp + .sn + .sp_l1 + .sn_l1 + .sp_l2 + .sn_l2 | year,
            data = d, cluster = ~iso3c, warn = FALSE, notes = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    ct <- summary(fit)$coeftable

    get_coef <- function(v, comp) {
      if (!v %in% rownames(ct)) return(NULL)
      tibble(h, component = comp,
             est   = ct[v, "Estimate"],
             se    = ct[v, "Std. Error"],
             pval  = ct[v, "Pr(>|t|)"],
             ci_lo = est - CI_Z * se,
             ci_hi = est + CI_Z * se,
             n_obs = nobs(fit))
    }
    bind_rows(get_coef(".sp", "pos"), get_coef(".sn", "neg"))
  })
}

# =============================================================================
# Section 1: Run for all (income group × exposure × commodity group)
# =============================================================================
message("\n=== Section 1: Asymmetric LP ===")

panel_by_group <- c(
  list(All = panel),
  setNames(map(INCOME_LEVELS, ~ filter(panel, income_group == .x)), INCOME_LEVELS)
)

irf_asym <- map_dfr(names(panel_by_group), function(inc_grp) {
  d_inc <- panel_by_group[[inc_grp]]
  message("  ", inc_grp, " (", n_distinct(d_inc$iso3c), " countries)")
  map_dfr(TREAT_TYPES, function(tt) {
    map_dfr(GRP_LEVELS, function(grp) {
      sv <- paste0("shock_", tt, "_", grp)
      run_lp_asym(d_inc, sv) |>
        mutate(treat = tt, grp = grp, income_group = inc_grp)
    })
  })
})

message("Asymmetric IRFs: ", nrow(irf_asym), " rows")

# Console: cons exposure, h=4, key cells
cat("\n=== h=4, consumption: β⁺ and −β⁻ × SD ===\n")
irf_asym |>
  filter(treat == "cons", h == 4, income_group %in% INCOME_LEVELS) |>
  mutate(
    est_disp = if_else(component == "neg", -est, est),
    est_sc   = est_disp * GRP_SD[grp],
    sig      = if_else(pval < 0.05, "*", ""),
    cell     = sprintf("%+.2f%s", est_sc, sig)
  ) |>
  select(component, income_group, grp, cell) |>
  pivot_wider(names_from = grp, values_from = cell) |>
  arrange(component, factor(income_group, levels = INCOME_LEVELS)) |>
  print(n = 40)

# =============================================================================
# Section 2: Heatmap — income × commodity, cons + prod, h=4
# Rows: surge (β⁺) | crash (−β⁻); Columns: cons | prod
# =============================================================================
message("\n=== Section 2: Asymmetric heatmap ===")

COL_NEG <- "#2166ac"
COL_POS <- "#d6604d"

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
      legend.key.width    = unit(14, "pt"),
      legend.text         = element_text(size = base_size - 1),
      plot.tag            = element_text(face = "bold", size = base_size + 1)
    )
}

# Shared data prep: significance markers + scaled estimates
# sig_label shows comma-separated horizons where p < 0.05 (matching scripts 12/13)
sig_labels <- irf_asym |>
  filter(income_group %in% INCOME_LEVELS,
         treat %in% c("cons", "prod"),
         !is.na(pval), pval < 0.05) |>
  group_by(income_group, grp, treat, component) |>
  summarise(sig_label = paste(sort(h), collapse = ","), .groups = "drop")

heat_dat <- irf_asym |>
  filter(income_group %in% INCOME_LEVELS,
         treat %in% c("cons", "prod"),
         h == 4) |>
  left_join(sig_labels, by = c("income_group", "grp", "treat", "component")) |>
  mutate(
    sig_label   = replace_na(sig_label, ""),
    grp_sd      = GRP_SD[grp],
    est_disp    = if_else(component == "neg", -est, est),  # flip crash sign
    est_scaled  = est_disp * grp_sd,
    income_f    = factor(income_group,
                         levels = rev(INCOME_LEVELS),
                         labels = rev(INCOME_LABELS)),
    grp_f       = factor(grp, levels = GRP_LEVELS, labels = GRP_LABELS),
    direction_f = factor(component,
                         levels = c("pos", "neg"),
                         labels = c("Price rises (shock > 0)",
                                    "Price falls (shock < 0)"))
  )

# Shared heatmap builder
make_heatmap <- function(dat) {
  ggplot(dat, aes(x = grp_f, y = income_f, fill = est_scaled)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sig_label), size = 2.5,
              vjust = 0.75, colour = "grey20") +
    scale_fill_gradient2(
      low      = COL_NEG,
      mid      = "white",
      high     = COL_POS,
      midpoint = 0,
      limits   = c(-HEAT_LIM, HEAT_LIM),
      oob      = scales::squish,
      name     = "FSI change per 1 SD shock (index points)",
      guide    = guide_colourbar(barwidth = 12, barheight = 0.4,
                                 title.position = "top", title.hjust = 0.5)
    ) +
    facet_wrap(~ direction_f, nrow = 1) +
    scale_x_discrete(name = NULL) +
    scale_y_discrete(name = NULL) +
    theme_nature() +
    theme(panel.grid = element_blank(), axis.ticks = element_blank())
}

# Figure 9a: consumption exposure
fig9a <- make_heatmap(filter(heat_dat, treat == "cons"))
ggsave(file.path(OUT_FIGS, "fig9a_asym_cons.pdf"),
       fig9a, width = 180, height = 100, units = "mm")
ggsave(file.path(OUT_FIGS, "fig9a_asym_cons.png"),
       fig9a, width = 180, height = 100, units = "mm", dpi = 300)
message("  Saved: fig9a_asym_cons.pdf / .png")

# Figure 9b: production exposure
fig9b <- make_heatmap(filter(heat_dat, treat == "prod"))
ggsave(file.path(OUT_FIGS, "fig9b_asym_prod.pdf"),
       fig9b, width = 180, height = 100, units = "mm")
ggsave(file.path(OUT_FIGS, "fig9b_asym_prod.png"),
       fig9b, width = 180, height = 100, units = "mm", dpi = 300)
message("  Saved: fig9b_asym_prod.pdf / .png")

# =============================================================================
# Section 3: LaTeX tables — one per treatment type
# Format: Panel A (β+, price rises) | Panel B (β−, price falls)
# Rows: commodity groups; Columns: income groups × horizons (h = 0..4)
# =============================================================================
message("\n=== Section 3: Asymmetric LP tables ===")

stars_fn <- function(p)
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", p < 0.10 ~ "$\\cdot$", TRUE ~ "")

# Stars as math superscripts — for use inside $...$
stars_math <- function(p)
  case_when(is.na(p) ~ "", p < 0.001 ~ "^{***}", p < 0.01 ~ "^{**}",
            p < 0.05 ~ "^{*}", p < 0.10 ~ "^{\\cdot}", TRUE ~ "")

GRP_LABELS_TEX <- c(
  cereals     = "Cereals",
  oils        = "Oils",
  animal      = "Animal prods",
  sugar_fruit = "Sugar \\& fruit",
  cash_crops  = "Beverage crops"
)

INCOME_LABELS_SHORT <- c(
  Low            = "Low income",
  `Lower-middle` = "Lower-middle",
  `Upper-middle` = "Upper-middle",
  High           = "High income"
)

# Column spec: label + 4 × (5 horizons separated by \quad)
col_spec <- paste0("l",
  paste(rep("rrrrr", length(INCOME_LEVELS)), collapse = "@{\\quad}"))

make_panel_rows <- function(tt, comp_tag) {
  rows <- character(0)
  for (grp in GRP_LEVELS) {
    # Estimate row
    est_cells <- map_chr(INCOME_LEVELS, function(ig) {
      paste(map_chr(LP_HORIZONS, function(h_) {
        r <- filter(irf_asym, treat == tt, grp == .env$grp,
                    income_group == ig, component == comp_tag, h == h_)
        if (nrow(r) == 0 || is.na(r$est)) return("---")
        sprintf("$%+.2f%s$", r$est, stars_math(r$pval))
      }), collapse = " & ")
    })
    # SE row
    se_cells <- map_chr(INCOME_LEVELS, function(ig) {
      paste(map_chr(LP_HORIZONS, function(h_) {
        r <- filter(irf_asym, treat == tt, grp == .env$grp,
                    income_group == ig, component == comp_tag, h == h_)
        if (nrow(r) == 0 || is.na(r$se)) return("")
        sprintf("$(%.2f)$", r$se)
      }), collapse = " & ")
    })
    rows <- c(rows,
      paste0(GRP_LABELS_TEX[grp], " & ",
             paste(est_cells, collapse = " & "), " \\\\"),
      paste0("  & ", paste(se_cells, collapse = " & "), " \\\\[3pt]")
    )
  }
  rows
}

make_asym_table <- function(tt) {
  treat_lbl <- TREAT_LABELS[tt]

  # Header rows: income group labels spanning 5 cols each
  n_inc  <- length(INCOME_LEVELS)
  n_h    <- length(LP_HORIZONS)
  cmidr  <- paste(sprintf("\\cmidrule(lr){%d-%d}",
                           2 + n_h * (0:(n_inc-1)),
                           1 + n_h * (1:n_inc)),
                  collapse = " ")
  inc_header <- paste0(" & ",
    paste(sprintf("\\multicolumn{%d}{c}{%s}", n_h, INCOME_LABELS_SHORT[INCOME_LEVELS]),
          collapse = " & "), " \\\\")
  h_header <- paste0("  & ",
    paste(rep(paste(paste0("$", LP_HORIZONS, "$"), collapse = "& "), n_inc),
          collapse = " & "), " \\\\")

  panel_header <- function(panel_lbl) {
    c(sprintf("\\multicolumn{%d}{l}{\\textbf{%s}} \\\\",
              1 + n_inc * n_h, panel_lbl),
      "\\midrule",
      inc_header,
      cmidr,
      h_header,
      "\\midrule")
  }

  tab <- c(
    sprintf("%% Asymmetric LP table: %s exposure", treat_lbl),
    "",
    "\\begin{table}[p]",
    "\\centering",
    sprintf("\\caption{Asymmetric local projection: \\textit{%s} exposure}", treat_lbl),
    sprintf("\\label{tab:asym:%s}", tt),
    "\\begin{adjustbox}{max width=\\textwidth, max totalheight=\\textheight}",
    "\\begin{threeparttable}",
    "\\small",
    sprintf("\\begin{tabular}{%s}", col_spec),
    "\\toprule",
    panel_header(sprintf("Panel A: Price rises ($\\hat{\\beta}^+$, shock $> 0$)")),
    make_panel_rows(tt, "pos"),
    "\\midrule",
    panel_header(sprintf("Panel B: Price falls ($\\hat{\\beta}^-$, shock $< 0$)")),
    make_panel_rows(tt, "neg"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{tablenotes}\\small",
    paste0("\\item \\textit{Notes:} ",
           "LP: $\\text{FSI}_{i,t+h} - \\text{FSI}_{i,t-1} = ",
           "\\hat{\\beta}^+ x^+_{it} + \\hat{\\beta}^- x^-_{it} + ",
           "\\sum_{\\ell=1}^{2}(\\gamma_\\ell x^+_{i,t-\\ell} + ",
           "\\delta_\\ell x^-_{i,t-\\ell}) + \\alpha_t + \\varepsilon_{iht}$, ",
           "where $x^+_{it}=\\max(x_{it},0)$ and $x^-_{it}=\\min(x_{it},0)$. ",
           sprintf("Shock: %s exposure by commodity group. ", treat_lbl),
           "Raw coefficients; multiply by commodity-group SD for standardised effects. ",
           "Year FE; SE clustered by country. ",
           "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."),
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\end{adjustbox}",
    "\\end{table}"
  )

  fname <- sprintf("tab_asym_%s.tex", tt)
  writeLines(tab, file.path(OUT_TABS, fname))
  message("  Saved: ", fname)
}

for (tt in TREAT_TYPES) make_asym_table(tt)

# =============================================================================
# Section 4: Save results
# =============================================================================
saveRDS(list(irf_asym = irf_asym, grp_sd = GRP_SD),
        file.path(DATA_PRO, "asym_lp_results.rds"))
message("\nSaved: data/processed/asym_lp_results.rds")
