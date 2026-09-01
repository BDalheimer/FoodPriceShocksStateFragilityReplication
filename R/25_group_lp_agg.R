# =============================================================================
# 25_group_lp_agg.R
# LP by commodity group — aggregate (full sample, no income split).
#
# Mirrors 10_income_lp.R but replaces the income-group dimension with
# commodity group. Runs separate LP for each group-level Bartik shock
# (already in panel.rds) across all five exposure types.
#
# Specification: country + year FEs, 2 lags, SE clustered by country.
# This is consistent with 10_income_lp.R.
#
# Outputs:
#   output/tables/tab_group_lp_agg.tex        — groups × treatments, h=1 & h=4
#   output/figures/fig_group_lp_agg_irf.pdf/.png
#     — cons + prod × 5 commodity groups, facet_grid(grp ~ treat), 120×170 mm
# =============================================================================

source(here::here("R", "00_setup.R"))

panel <- readRDS(file.path(DATA_PRO, "panel.rds"))

message("Panel: ", nrow(panel), " obs | ",
        n_distinct(panel$iso3c), " countries | ",
        min(panel$year), "–", max(panel$year))

# ---- Constants --------------------------------------------------------------
TREAT_TYPES  <- c("cons", "prod", "imp", "exp", "net")
TREAT_LABELS <- c(cons = "Consumption", prod = "Production",
                  imp  = "Imports",     exp  = "Exports", net = "Net imports")

GRP_LEVELS <- c("cereals", "oils", "animal", "sugar_fruit", "cash_crops")
GRP_LABELS <- c(cereals     = "Cereals",
                oils        = "Oils",
                animal      = "Animal products",
                sugar_fruit = "Sugar \\& fruit",
                cash_crops  = "Beverage crops")
GRP_LABELS_FIG <- c(cereals     = "Cereals",
                    oils        = "Oils",
                    animal      = "Animal products",
                    sugar_fruit = "Sugar & fruit",
                    cash_crops  = "Beverage crops")

CI_Z <- 1.645

stars_fn <- function(p)
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", p < 0.10 ~ "$\\cdot$", TRUE ~ "")

stars_console <- function(p)
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
  ifelse(p < 0.05,  "*",   ifelse(p < 0.10, ".", ""))))

write_tex <- function(lines, filename) {
  writeLines(lines, file.path(OUT_TABS, filename))
  message("  Saved: ", filename)
}

# ---- LP runner (country + year FEs, matching 10_income_lp.R) ----------------
run_lp <- function(data, shock_var, fd_stem = "dfsi", horizons = LP_HORIZONS) {
  lag1 <- paste0(shock_var, "_l1")
  lag2 <- paste0(shock_var, "_l2")
  map_dfr(horizons, function(h) {
    y_col <- paste0(fd_stem, "_h", h)
    if (!y_col %in% names(data) || !shock_var %in% names(data)) return(NULL)
    d <- data |>
      rename(.y = all_of(y_col), .s = all_of(shock_var)) |>
      filter(!is.na(.y), !is.na(.s))
    if (nrow(d) < 30) return(NULL)
    if (lag1 %in% names(d)) d <- rename(d, .s_l1 = all_of(lag1))
    if (lag2 %in% names(d)) d <- rename(d, .s_l2 = all_of(lag2))
    rhs_lags <- intersect(c(".s_l1", ".s_l2"), names(d))
    fml <- as.formula(paste0(
      ".y ~ .s",
      if (length(rhs_lags)) paste0(" + ", paste(rhs_lags, collapse = " + ")),
      " | iso3c + year"
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

# =============================================================================
# Section 1: Run LP for each commodity group × exposure type
# =============================================================================
message("\n=== Section 1: Group LP (aggregate) ===")

irf_grp <- map_dfr(GRP_LEVELS, function(grp) {
  map_dfr(TREAT_TYPES, function(tt) {
    sv <- paste0("shock_", tt, "_", grp)
    if (!sv %in% names(panel)) {
      message("  Missing: ", sv)
      return(NULL)
    }
    run_lp(panel, sv) |> mutate(grp = grp, treat = tt)
  })
})

# Console summary
cat(sprintf("\n%-16s  %-14s  %5s  %8s  %6s  %5s\n",
            "Group", "Exposure", "h", "est", "se", "p"))
cat(strrep("-", 60), "\n")
for (grp in GRP_LEVELS) {
  for (tt in TREAT_TYPES) {
    for (h_ in c(1, 4)) {
      r <- filter(irf_grp, grp == .env$grp, treat == tt, h == h_)
      if (nrow(r) == 0) next
      cat(sprintf("%-16s  %-14s  h=%d  %+8.3f  %6.3f  %4.3f%s\n",
                  GRP_LABELS[grp], TREAT_LABELS[tt], h_,
                  r$est, r$se, r$pval, stars_console(r$pval)))
    }
  }
  cat("\n")
}

# =============================================================================
# Section 2: LaTeX table — groups × treatments, h=1 and h=4
# Rows: commodity groups; column blocks: treatment types
# =============================================================================
message("\n=== Section 2: LaTeX table ===")

# 5 treatments × 2 horizons = 10 data columns → landscape recommended
col_spec <- paste0("l ", paste(rep("rr", length(TREAT_TYPES)), collapse = " "))

tab_lines <- c(
  "% Commodity group LP: all exposure types, h=1 and h=4",
  "",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Local projections by commodity group (aggregate sample)}",
  "\\label{tab:group-lp-agg}",
  "\\begin{threeparttable}",
  "\\small",
  sprintf("\\begin{tabular}{%s}", col_spec),
  "\\toprule",
  paste0(" & ",
    paste(sapply(TREAT_TYPES, function(tt)
      sprintf("\\multicolumn{2}{c}{%s}", TREAT_LABELS[tt])),
      collapse = " & "), " \\\\"),
  paste0(paste(sprintf("\\cmidrule(lr){%d-%d}",
                       seq(2, 2 + 2*(length(TREAT_TYPES)-1), 2),
                       seq(3, 3 + 2*(length(TREAT_TYPES)-1), 2)),
               collapse = " ")),
  paste0("Commodity group & ",
    paste(rep("$h=1$ & $h=4$", length(TREAT_TYPES)), collapse = " & "), " \\\\"),
  "\\midrule"
)

for (grp in GRP_LEVELS) {
  est_cells <- map_chr(TREAT_TYPES, function(tt) {
    r1 <- filter(irf_grp, grp == .env$grp, treat == tt, h == 1)
    r4 <- filter(irf_grp, grp == .env$grp, treat == tt, h == 4)
    c1 <- if (nrow(r1) > 0 && !is.na(r1$est))
            sprintf("%.3f%s", r1$est, stars_fn(r1$pval)) else "---"
    c4 <- if (nrow(r4) > 0 && !is.na(r4$est))
            sprintf("%.3f%s", r4$est, stars_fn(r4$pval)) else "---"
    paste(c1, c4, sep = " & ")
  })
  se_cells <- map_chr(TREAT_TYPES, function(tt) {
    r1 <- filter(irf_grp, grp == .env$grp, treat == tt, h == 1)
    r4 <- filter(irf_grp, grp == .env$grp, treat == tt, h == 4)
    s1 <- if (nrow(r1) > 0 && !is.na(r1$se)) sprintf("(%.3f)", r1$se) else ""
    s4 <- if (nrow(r4) > 0 && !is.na(r4$se)) sprintf("(%.3f)", r4$se) else ""
    paste(s1, s4, sep = " & ")
  })
  n_ <- filter(irf_grp, grp == .env$grp, treat == "cons", h == 1)$n_obs
  tab_lines <- c(tab_lines,
    sprintf("%s & %s \\\\", GRP_LABELS[grp], paste(est_cells, collapse = " & ")),
    sprintf("   & %s \\\\[4pt]", paste(se_cells, collapse = " & "))
  )
}

# N row from consumption h=1
n_cells <- map_chr(TREAT_TYPES, function(tt) {
  r1 <- filter(irf_grp, grp == "cereals", treat == tt, h == 1)
  r4 <- filter(irf_grp, grp == "cereals", treat == tt, h == 4)
  n1 <- if (nrow(r1) > 0) format(r1$n_obs, big.mark = ",") else "---"
  n4 <- if (nrow(r4) > 0) format(r4$n_obs, big.mark = ",") else "---"
  paste(n1, n4, sep = " & ")
})

tab_lines <- c(tab_lines, "\\midrule",
               sprintf("$N$ & %s \\\\", paste(n_cells, collapse = " & ")),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0("\\item \\textit{Notes:} Each cell reports the LP coefficient from ",
         "$\\text{FSI}_{i,t+h} - \\text{FSI}_{i,t-1} = \\beta x^g_{it} + ",
         "\\sum_{\\ell=1}^{2}\\gamma_\\ell x^g_{i,t-\\ell} + \\alpha_i + \\alpha_t + \\varepsilon_{iht}$, ",
         "where $x^g_{it}$ is the commodity-group Bartik shock for group $g$. ",
         "Country and year fixed effects; SE clustered by country. ",
         "$N$ from cereals group (other groups similar). ",
         "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)

write_tex(tab_lines, "tab_group_lp_agg.tex")

# =============================================================================
# Section 3: Figure — cons + prod × 5 commodity groups
# =============================================================================
message("\n=== Section 3: Figure ===")

COL_NEG <- "#2166ac"
COL_POS <- "#d6604d"
COL_RIB <- "grey70"

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
      legend.position     = "bottom"
    )
}

fig_dat <- irf_grp |>
  filter(treat %in% c("cons", "prod")) |>
  mutate(
    grp_f   = factor(grp,   levels = GRP_LEVELS, labels = unname(GRP_LABELS_FIG)),
    treat_f = factor(treat, levels = c("cons", "prod"),
                     labels = c("Consumption exposure", "Production exposure"))
  )

fig_sign <- fig_dat |>
  group_by(treat_f, grp_f) |>
  summarise(col = if_else(mean(est, na.rm = TRUE) < 0, COL_NEG, COL_POS),
            .groups = "drop")

fig_dat <- left_join(fig_dat, fig_sign, by = c("treat_f", "grp_f"))

fig_grp <- ggplot(fig_dat, aes(x = h, y = est)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.35) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              fill = COL_RIB, alpha = 0.25) +
  geom_line(aes(colour = col), linewidth = 0.7) +
  geom_point(aes(colour = col), size = 1.6) +
  scale_colour_identity() +
  scale_x_continuous(breaks = 0:4, name = "Horizon (years)") +
  scale_y_continuous(name = "FSI change (index points)") +
  facet_grid(grp_f ~ treat_f, scales = "free_y") +
  theme_nature() +
  theme(strip.text.y = element_text(angle = 0, hjust = 0))

ggsave(file.path(OUT_FIGS, "fig_group_lp_agg_irf.pdf"),
       fig_grp, width = 120, height = 170, units = "mm")
ggsave(file.path(OUT_FIGS, "fig_group_lp_agg_irf.png"),
       fig_grp, width = 120, height = 170, units = "mm", dpi = 300)
message("  Saved: fig_group_lp_agg_irf.pdf / .png")
message("\nDone.")
