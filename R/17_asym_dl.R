# =============================================================================
# 17_asym_dl.R
# Asymmetric distributed lag (DL) models by commodity group × income group
#
# Mirrors 16_asym_lp.R but uses the DL specification from 09_dl.R.
# Shock is split into positive (surges) and negative (crashes) components;
# all lags of both are included jointly in a single regression on dfsi_h0:
#
#   ΔFSI_{i,t} = α_t + Σ_{k=0}^{H} [β⁺_k shock_pos_{i,t-k}
#                                   + β⁻_k shock_neg_{i,t-k}] + ε_{it}
#
# IRF at horizon h: cumulative multiplier B⁺_h = Σ_{k=0}^{h} β⁺_k (and B⁻_h).
# Delta-method SEs from the clustered vcov sub-matrix.
#
# Sign convention in heatmap:
#   Surge  panel  → B⁺_h × SD(group shock)      positive = destabilising
#   Crash  panel  → −B⁻_h × SD(group shock)     positive = destabilising
#
# Significance markers: comma-separated horizons where p < 0.05.
#
# Outputs:
#   output/figures/fig10a_asym_dl_cons.pdf/.png
#   output/figures/fig10b_asym_dl_prod.pdf/.png
#   data/processed/asym_dl_results.rds
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

CI_Z     <- 1.96    # 95% CI for asymmetric models
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
# Asymmetric DL runner
# Splits shock into pos/neg, adds 4 lags of each, runs one regression on
# dfsi_h0, then returns cumulative multipliers for each component at each h.
# =============================================================================
run_dl_asym <- function(data, shock_var, fd_stem = "dfsi", horizons = LP_HORIZONS) {
  if (n_distinct(data$iso3c) < 5) return(NULL)
  max_h <- max(horizons)
  y_col <- paste0(fd_stem, "_h0")
  if (!y_col %in% names(data)) return(NULL)

  d_prep <- data |>
    arrange(iso3c, year) |>
    group_by(iso3c) |>
    mutate(
      .sp    = pmax(.data[[shock_var]], 0),
      .sn    = pmin(.data[[shock_var]], 0),
      .sp_l1 = lag(.sp, 1), .sp_l2 = lag(.sp, 2),
      .sp_l3 = lag(.sp, 3), .sp_l4 = lag(.sp, 4),
      .sn_l1 = lag(.sn, 1), .sn_l2 = lag(.sn, 2),
      .sn_l3 = lag(.sn, 3), .sn_l4 = lag(.sn, 4)
    ) |>
    ungroup()

  sp_vars <- c(".sp",  paste0(".sp_l", seq_len(max_h)))
  sn_vars <- c(".sn",  paste0(".sn_l", seq_len(max_h)))
  all_vars <- c(sp_vars, sn_vars)
  all_vars <- all_vars[all_vars %in% names(d_prep)]

  d <- d_prep |>
    filter(!is.na(.data[[y_col]]), !is.na(.sp), !is.na(.sn))
  if (nrow(d) < 20) return(NULL)

  fml <- as.formula(paste0(y_col, " ~ ", paste(all_vars, collapse = " + "), " | year"))
  fit <- tryCatch(
    feols(fml, data = d, cluster = ~iso3c, warn = FALSE, notes = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  ct         <- summary(fit)$coeftable
  vc         <- vcov(fit)
  n          <- nobs(fit)
  sp_coefs   <- sp_vars[sp_vars %in% rownames(ct)]
  sn_coefs   <- sn_vars[sn_vars %in% rownames(ct)]

  get_cum <- function(coef_seq, comp, h) {
    idx <- coef_seq[seq_len(min(h + 1L, length(coef_seq)))]
    if (length(idx) == 0) return(NULL)
    B    <- sum(ct[idx, "Estimate"])
    se   <- sqrt(sum(vc[idx, idx]))
    pval <- 2 * pnorm(-abs(B / se))
    tibble(h = h, component = comp,
           est   = B, se = se, pval = pval,
           ci_lo = B - CI_Z * se, ci_hi = B + CI_Z * se,
           n_obs = n)
  }

  map_dfr(horizons, function(h) {
    bind_rows(get_cum(sp_coefs, "pos", h),
              get_cum(sn_coefs, "neg", h))
  })
}

# =============================================================================
# Section 1: Run for all (income group × exposure × commodity group)
# =============================================================================
message("\n=== Section 1: Asymmetric DL ===")

panel_by_group <- c(
  list(All = panel),
  setNames(map(INCOME_LEVELS, ~ filter(panel, income_group == .x)), INCOME_LEVELS)
)

irf_asym_dl <- map_dfr(names(panel_by_group), function(inc_grp) {
  d_inc <- panel_by_group[[inc_grp]]
  message("  ", inc_grp, " (", n_distinct(d_inc$iso3c), " countries)")
  map_dfr(TREAT_TYPES, function(tt) {
    map_dfr(GRP_LEVELS, function(grp) {
      sv <- paste0("shock_", tt, "_", grp)
      run_dl_asym(d_inc, sv) |>
        mutate(treat = tt, grp = grp, income_group = inc_grp)
    })
  })
})

message("Asymmetric DL IRFs: ", nrow(irf_asym_dl), " rows")

# Console: cons exposure, h=4, key cells
cat("\n=== h=4, consumption: B⁺ and −B⁻ × SD ===\n")
irf_asym_dl |>
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
# Section 2: Heatmaps — income × commodity, cons + prod, h=4
# =============================================================================
message("\n=== Section 2: Asymmetric DL heatmaps ===")

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

sig_labels <- irf_asym_dl |>
  filter(income_group %in% INCOME_LEVELS,
         treat %in% c("cons", "prod"),
         !is.na(pval), pval < 0.05) |>
  group_by(income_group, grp, treat, component) |>
  summarise(sig_label = paste(sort(h), collapse = ","), .groups = "drop")

heat_dat <- irf_asym_dl |>
  filter(income_group %in% INCOME_LEVELS,
         treat %in% c("cons", "prod"),
         h == 4) |>
  left_join(sig_labels, by = c("income_group", "grp", "treat", "component")) |>
  mutate(
    sig_label   = replace_na(sig_label, ""),
    grp_sd      = GRP_SD[grp],
    est_disp    = if_else(component == "neg", -est, est),
    est_scaled  = est_disp * grp_sd,
    income_f    = factor(income_group,
                         levels = rev(INCOME_LEVELS),
                         labels = rev(INCOME_LABELS)),
    grp_f       = factor(grp, levels = GRP_LEVELS, labels = GRP_LABELS),
    direction_f = factor(component,
                         levels = c("pos", "neg"),
                         labels = c("Price surges (B\u207a)",
                                    "Price crashes (\u2212B\u207b)"))
  )

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
    facet_wrap(~ direction_f, ncol = 1) +
    scale_x_discrete(name = NULL) +
    scale_y_discrete(name = NULL) +
    theme_nature() +
    theme(panel.grid = element_blank(), axis.ticks = element_blank())
}

# Figure 10a: consumption exposure
fig10a <- make_heatmap(filter(heat_dat, treat == "cons"))
ggsave(file.path(OUT_FIGS, "fig10a_asym_dl_cons.pdf"),
       fig10a, width = 180, height = 100, units = "mm")
ggsave(file.path(OUT_FIGS, "fig10a_asym_dl_cons.png"),
       fig10a, width = 180, height = 100, units = "mm", dpi = 300)
message("  Saved: fig10a_asym_dl_cons.pdf / .png")

# Figure 10b: production exposure
fig10b <- make_heatmap(filter(heat_dat, treat == "prod"))
ggsave(file.path(OUT_FIGS, "fig10b_asym_dl_prod.pdf"),
       fig10b, width = 180, height = 100, units = "mm")
ggsave(file.path(OUT_FIGS, "fig10b_asym_dl_prod.png"),
       fig10b, width = 180, height = 100, units = "mm", dpi = 300)
message("  Saved: fig10b_asym_dl_prod.pdf / .png")

# =============================================================================
# Section 3: Save results
# =============================================================================
saveRDS(list(irf_asym_dl = irf_asym_dl, grp_sd = GRP_SD),
        file.path(DATA_PRO, "asym_dl_results.rds"))
message("\nSaved: data/processed/asym_dl_results.rds")
