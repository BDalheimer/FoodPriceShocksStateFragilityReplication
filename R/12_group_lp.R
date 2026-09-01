# =============================================================================
# 12_group_lp.R
# LP models by commodity group × income group
#
# For each combination of income group and commodity group, runs the LP using
# the additive group-level shocks from shocks_group.rds:
#   shock_{treat}_{grp} = Σ_{c ∈ grp} share_{i,c,t} × Δlog p_{c,t}
#
# Cells with fewer than 5 countries are skipped.
#
# Outputs:
#   output/figures/fig5_heatmap.pdf/.png
#     — income groups (y) × commodity groups (x), faceted by exposure (cons + net)
#     — coefficient at h=4, scaled by SD of group shock; significance text in cells
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

CI_Z <- 1.645

# ---- Compute group shock SDs (for heatmap normalisation) --------------------
# SD of each group-level consumption shock across all country-years
GRP_SD <- map_dbl(
  setNames(GRP_LEVELS, GRP_LEVELS),
  ~ sd(panel[[paste0("shock_cons_", .x)]], na.rm = TRUE)
)
message("\nGroup shock SDs (consumption):")
print(round(GRP_SD, 4))

# ---- Helpers ----------------------------------------------------------------
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

write_tex <- function(lines, filename) {
  writeLines(lines, file.path(OUT_TABS, filename))
  message("  Saved: ", filename)
}

# =============================================================================
# Section 1: Run LP for all (income group × treatment × commodity group)
# =============================================================================
message("\n=== Section 1: Group LP by income group ===")

if (!"income_group" %in% names(panel)) {
  stop("income_group not in panel — run 07_panel.R with income_groups.rds present.")
}

panel_lmic <- panel |>
  filter(income_group %in% c("Low", "Lower-middle")) |>
  mutate(income_group = "LMIC")

panel_by_group <- c(
  list(All = panel),
  setNames(map(INCOME_LEVELS, ~ filter(panel, income_group == .x)), INCOME_LEVELS),
  list(LMIC = panel_lmic)
)

# Verify group shock columns exist
grp_shock_cols <- outer(paste0("shock_", TREAT_TYPES, "_"), GRP_LEVELS, paste0) |> c()
missing <- setdiff(grp_shock_cols, names(panel))
if (length(missing) > 0) {
  stop("Missing group shock columns — run 05_group_shocks.R first:\n  ",
       paste(missing, collapse = ", "))
}

# Run LP for each combination
irf_grp <- map_dfr(names(panel_by_group), function(inc_grp) {
  d_inc <- panel_by_group[[inc_grp]]
  n_c   <- n_distinct(d_inc$iso3c)
  if (n_c < 5) {
    message("  Skipping ", inc_grp, " (only ", n_c, " countries)")
    return(NULL)
  }
  message("  ", inc_grp, " (", n_c, " countries)")
  map_dfr(TREAT_TYPES, function(tt) {
    map_dfr(GRP_LEVELS, function(grp) {
      sv <- paste0("shock_", tt, "_", grp)
      run_lp(d_inc, sv) |>
        mutate(treat = tt, grp = grp, income_group = inc_grp)
    })
  })
})

message("\nIRFs computed: ", nrow(irf_grp), " rows")

# Console snapshot: h=4, consumption, all income groups
cat("\n=== h=4 snapshot: consumption exposure ===\n")
irf_grp |>
  filter(treat == "cons", h == 4, income_group %in% INCOME_LEVELS) |>
  mutate(cell = sprintf("%+.4f (p=%.3f)", est, pval)) |>
  select(income_group, grp, cell) |>
  pivot_wider(names_from = grp, values_from = cell) |>
  print()

# =============================================================================
# Section 2: Heatmap — income × commodity, cons + net, h=4
# Normalised: est × SD(group shock) → FSI points per 1-SD group shock
# Significance text: comma-separated horizons where p < 0.05
# =============================================================================
message("\n=== Section 2: Heatmap figure ===")

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

# Significance text: horizons with p < 0.05 for each (income × commodity × treat) cell
sig_labels <- irf_grp |>
  filter(income_group %in% INCOME_LEVELS,
         treat %in% c("cons", "prod", "net"),
         !is.na(pval), pval < 0.05) |>
  group_by(income_group, grp, treat) |>
  summarise(sig_label = paste(sort(h), collapse = ","), .groups = "drop")

# Heatmap data: h=4, scale by group shock SD
fig5_dat <- irf_grp |>
  filter(income_group %in% INCOME_LEVELS,
         treat %in% c("cons", "prod", "net"),
         h == 4) |>
  left_join(sig_labels, by = c("income_group", "grp", "treat")) |>
  mutate(
    sig_label  = replace_na(sig_label, ""),
    grp_sd     = GRP_SD[grp],
    est_scaled = est * grp_sd,
    income_f   = factor(income_group,
                        levels = rev(INCOME_LEVELS),
                        labels = rev(INCOME_LABELS)),
    grp_f      = factor(grp, levels = GRP_LEVELS, labels = GRP_LABELS),
    treat_f    = factor(treat,
                        levels = c("cons", "prod", "net"),
                        labels = c("Consumption exposure",
                                   "Production exposure",
                                   "Net-import exposure"))
  )

# Colour scale truncated at ±2: values beyond are squished to the extreme colour.
# Prevents one large insignificant outlier from washing out the rest of the scale.
HEAT_LIM <- 2

fig5 <- ggplot(fig5_dat, aes(x = grp_f, y = income_f, fill = est_scaled)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sig_label), size = 2.0,
            vjust = 0.75, colour = "grey20") +
  scale_fill_gradient2(
    low      = COL_NEG,
    mid      = "white",
    high     = COL_POS,
    midpoint = 0,
    limits   = c(-HEAT_LIM, HEAT_LIM),
    oob      = scales::squish,
    name     = "FSI change per 1 SD shock (index points)",
    guide    = guide_colourbar(barwidth = 10, barheight = 0.4,
                               title.position = "top", title.hjust = 0.5)
  ) +
  facet_wrap(~ treat_f, ncol = 3) +
  scale_x_discrete(name = NULL) +
  scale_y_discrete(name = NULL) +
  theme_nature() +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(file.path(OUT_FIGS, "fig5_heatmap.pdf"),
       fig5, width = 180, height = 80, units = "mm")
ggsave(file.path(OUT_FIGS, "fig5_heatmap.png"),
       fig5, width = 180, height = 80, units = "mm", dpi = 300)

message("  Saved: fig5_heatmap.pdf / .png")

# =============================================================================
# Section 3: Appendix tables — one per commodity group
# Rows: income groups; columns: h=0..4
# Two panels per table: consumption and production exposure.
# Labels: tab:app-grp-cereals, tab:app-grp-oils, etc.
# =============================================================================
message("\n=== Section 3: Appendix tables (one per commodity group) ===")

stars_fn <- function(p)
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", p < 0.10 ~ "$\\cdot$", TRUE ~ "")

GRP_LABELS_TEX <- c(
  cereals     = "Cereals",
  oils        = "Oils",
  animal      = "Animal products",
  sugar_fruit = "Sugar \\& fruit",
  cash_crops  = "Beverage crops"
)

INCOME_LABELS_TEX <- c(
  Low            = "Low income",
  `Lower-middle` = "Lower-middle",
  `Upper-middle` = "Upper-middle",
  High           = "High income"
)

for (grp in GRP_LEVELS) {
  grp_lbl <- GRP_LABELS_TEX[grp]

  make_panel <- function(treat_tag, treat_lbl) {
    d_treat <- filter(irf_grp, grp == .env$grp, treat == treat_tag,
                      income_group %in% INCOME_LEVELS)
    rows <- character(0)
    for (ig in INCOME_LEVELS) {
      est_row <- map_chr(LP_HORIZONS, function(h_) {
        r <- filter(d_treat, income_group == ig, h == h_)
        if (nrow(r) == 0 || is.na(r$est)) return("---")
        sprintf("$%+.3f%s$", r$est, stars_fn(r$pval))
      })
      se_row <- map_chr(LP_HORIZONS, function(h_) {
        r <- filter(d_treat, income_group == ig, h == h_)
        if (nrow(r) == 0 || is.na(r$se)) return("")
        sprintf("$(%.3f)$", r$se)
      })
      rows <- c(rows,
        paste0(INCOME_LABELS_TEX[ig], " & ",
               paste(est_row, collapse = " & "), " \\\\"),
        paste0("& ", paste(se_row, collapse = " & "), " \\\\[3pt]")
      )
    }
    n_low <- filter(d_treat, income_group == "Low", h == 0)$n_obs
    n_row <- paste0("$N$ (Low) & ",
                    paste(map_chr(LP_HORIZONS, function(h_) {
                      r <- filter(d_treat, income_group == "Low", h == h_)
                      if (nrow(r) == 0) "---" else format(r$n_obs, big.mark = ",")
                    }), collapse = " & "), " \\\\")
    c(sprintf("\\multicolumn{6}{l}{\\textit{%s}} \\\\", treat_lbl),
      "\\midrule",
      rows,
      "\\midrule",
      n_row,
      "\\addlinespace")
  }

  tab <- c(
    sprintf("%% Group LP appendix table: %s", grp_lbl),
    "",
    "\\begin{table}[htbp]",
    "\\centering",
    sprintf("\\caption{Commodity-group LP: %s (FSI outcome, by income group)}", grp_lbl),
    sprintf("\\label{tab:app-grp-%s}", grp),
    "\\begin{adjustbox}{max width=\\textwidth}",
    "\\begin{threeparttable}",
    "\\small",
    "\\begin{tabular}{l rrrrr}",
    "\\toprule",
    "Income group & $h=0$ & $h=1$ & $h=2$ & $h=3$ & $h=4$ \\\\",
    "\\toprule",
    make_panel("cons", "Consumption exposure"),
    make_panel("prod", "Production exposure"),
    make_panel("imp",  "Import exposure"),
    make_panel("exp",  "Export exposure"),
    make_panel("net",  "Net-import exposure"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{tablenotes}\\small",
    paste0("\\item \\textit{Notes:} LP estimates of ",
           "$\\Delta_h\\text{FSI}_{it} = \\beta_h x^g_{it} + ",
           "\\sum_{\\ell=1}^{2}\\gamma_\\ell x^g_{i,t-\\ell} + \\alpha_t + \\varepsilon_{iht}$, ",
           "where $x^g_{it}$ is the ", grp_lbl, " Bartik shock. ",
           "Year FE; SE clustered by country. ",
           "Coefficients are raw LP estimates in FSI index points per unit shock; ",
           sprintf("multiply by the group shock SD ($%.4f$) to convert to ", GRP_SD[grp]),
           "FSI points per 1-SD shock. ",
           "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."),
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\end{adjustbox}",
    "\\end{table}"
  )

  fname <- sprintf("tab_app_grp_%s.tex", grp)
  write_tex(tab, fname)
}

# =============================================================================
# Section 4: Save IRF results for use in other scripts
# =============================================================================
saveRDS(list(irf_grp = irf_grp, grp_sd = GRP_SD),
        file.path(DATA_PRO, "group_lp_results.rds"))
message("\nSaved: data/processed/group_lp_results.rds")
