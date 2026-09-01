# =============================================================================
# 13_group_dl.R
# Distributed lag models by commodity group × income group
#
# Mirrors 12_group_lp.R but uses the DL specification from 09_dl.R.
# IRF at horizon h: cumulative multiplier B_h = Σ_{k=0}^{h} β_k,
# delta-method SEs from the clustered vcov.
#
# Outputs:
#   output/figures/fig6_heatmap_dl.pdf/.png
#   data/processed/group_dl_results.rds
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

CI_Z     <- 1.645
HEAT_LIM <- 2   # colour scale truncated at ±2 (squished, not clipped)

# ---- Group shock SDs (same as 12_group_lp.R) --------------------------------
GRP_SD <- map_dbl(
  setNames(GRP_LEVELS, GRP_LEVELS),
  ~ sd(panel[[paste0("shock_cons_", .x)]], na.rm = TRUE)
)
message("\nGroup shock SDs (consumption):")
print(round(GRP_SD, 4))

# ---- Helpers ----------------------------------------------------------------
run_dl <- function(data, shock_var, fd_stem = "dfsi", horizons = LP_HORIZONS) {
  max_h    <- max(horizons)
  y_col    <- paste0(fd_stem, "_h0")
  lag_vars <- c(shock_var, paste0(shock_var, "_l", seq_len(max_h)))
  lag_vars <- lag_vars[lag_vars %in% names(data)]

  d <- data |> filter(!is.na(.data[[y_col]]), !is.na(.data[[shock_var]]))
  if (nrow(d) < 10) return(NULL)

  fml <- as.formula(paste0(
    y_col, " ~ ", paste(lag_vars, collapse = " + "), " | year"
  ))
  fit <- tryCatch(
    feols(fml, data = d, cluster = ~iso3c, warn = FALSE, notes = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  ct         <- summary(fit)$coeftable
  vc         <- vcov(fit)
  coef_names <- lag_vars[lag_vars %in% rownames(ct)]
  n          <- nobs(fit)

  map_dfr(horizons, function(h) {
    idx  <- coef_names[seq_len(min(h + 1L, length(coef_names)))]
    if (length(idx) == 0) return(NULL)
    B_h  <- sum(ct[idx, "Estimate"])
    se_h <- sqrt(sum(vc[idx, idx]))
    pval <- 2 * pnorm(-abs(B_h / se_h))
    tibble(h = h, est = B_h, se = se_h, pval = pval,
           ci_lo = B_h - CI_Z * se_h,
           ci_hi = B_h + CI_Z * se_h,
           n_obs = n)
  })
}

# =============================================================================
# Section 1: Run DL for all (income group × treatment × commodity group)
# =============================================================================
message("\n=== Section 1: Group DL by income group ===")

if (!"income_group" %in% names(panel)) {
  stop("income_group not in panel — run 07_panel.R with income_groups.rds present.")
}

grp_shock_cols <- outer(paste0("shock_", TREAT_TYPES, "_"), GRP_LEVELS, paste0) |> c()
missing <- setdiff(grp_shock_cols, names(panel))
if (length(missing) > 0) {
  stop("Missing group shock columns — run 05_group_shocks.R first:\n  ",
       paste(missing, collapse = ", "))
}

panel_lmic <- panel |>
  filter(income_group %in% c("Low", "Lower-middle")) |>
  mutate(income_group = "LMIC")

panel_by_group <- c(
  list(All = panel),
  setNames(map(INCOME_LEVELS, ~ filter(panel, income_group == .x)), INCOME_LEVELS),
  list(LMIC = panel_lmic)
)

irf_grp_dl <- map_dfr(names(panel_by_group), function(inc_grp) {
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
      run_dl(d_inc, sv) |>
        mutate(treat = tt, grp = grp, income_group = inc_grp)
    })
  })
})

message("\nIRFs computed: ", nrow(irf_grp_dl), " rows")

# Console snapshot: h=4, consumption, all income groups
cat("\n=== h=4 snapshot: consumption exposure ===\n")
irf_grp_dl |>
  filter(treat == "cons", h == 4, income_group %in% INCOME_LEVELS) |>
  mutate(cell = sprintf("%+.4f (p=%.3f)", est, pval)) |>
  select(income_group, grp, cell) |>
  pivot_wider(names_from = grp, values_from = cell) |>
  print()

# =============================================================================
# Section 2: Heatmap — income × commodity, cons + net, h=4
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

sig_labels <- irf_grp_dl |>
  filter(income_group %in% INCOME_LEVELS,
         treat %in% c("cons", "net"),
         !is.na(pval), pval < 0.05) |>
  group_by(income_group, grp, treat) |>
  summarise(sig_label = paste(sort(h), collapse = ","), .groups = "drop")

fig6_dat <- irf_grp_dl |>
  filter(income_group %in% INCOME_LEVELS,
         treat %in% c("cons", "net"),
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
                        levels = c("cons", "net"),
                        labels = c("Consumption exposure",
                                   "Net-import exposure"))
  )

fig6 <- ggplot(fig6_dat, aes(x = grp_f, y = income_f, fill = est_scaled)) +
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
  facet_wrap(~ treat_f, ncol = 2) +
  scale_x_discrete(name = NULL) +
  scale_y_discrete(name = NULL) +
  theme_nature() +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(file.path(OUT_FIGS, "fig6_heatmap_dl.pdf"),
       fig6, width = 180, height = 80, units = "mm")
ggsave(file.path(OUT_FIGS, "fig6_heatmap_dl.png"),
       fig6, width = 180, height = 80, units = "mm", dpi = 300)

message("  Saved: fig6_heatmap_dl.pdf / .png")

# =============================================================================
# Section 3: Save results
# =============================================================================
saveRDS(list(irf_grp_dl = irf_grp_dl, grp_sd = GRP_SD),
        file.path(DATA_PRO, "group_dl_results.rds"))
message("\nSaved: data/processed/group_dl_results.rds")
