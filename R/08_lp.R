# =============================================================================
# 08_lp.R
# Local projection estimates and paper tables
#
# LP specification:
#   FSI_{i,t+h} - FSI_{i,t-1} = β x_{it} + γ₁ x_{i,t-1} + γ₂ x_{i,t-2}
#                                + α_t + ε_{iht}
#   Year FE; SE clustered by country.
#
# Built step by step. Each section saves its table to output/tables/.
# =============================================================================

source(here::here("R", "00_setup.R"))

panel <- readRDS(file.path(DATA_PRO, "panel.rds"))

message("Panel: ", nrow(panel), " obs | ",
        n_distinct(panel$iso3c), " countries | ",
        min(panel$year), "–", max(panel$year))

# ---- Constants --------------------------------------------------------------
TREAT_TYPES  <- c("cons", "prod", "imp", "exp", "net", "val_cons", "val_prod")
TREAT_LABELS <- c(cons      = "Consumption",
                  prod      = "Production",
                  imp       = "Imports",
                  exp       = "Exports",
                  net       = "Net imports",
                  val_cons  = "Consumption (value)",
                  val_prod  = "Production (value)")

CI_Z <- 1.645   # 90% CI for IRF figures

# ---- Helpers ----------------------------------------------------------------
stars_fn <- function(p) {
  case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.10  ~ "$\\cdot$",
    TRUE      ~ ""
  )
}

stars_console <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
  ifelse(p < 0.05,  "*",   ifelse(p < 0.10, ".",  ""))))
}

# run_lp: runs LP for one shock variable across all horizons
# outcome stem: "dfsi" → looks for dfsi_h0 ... dfsi_h4
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

# write_tex: saves a character vector as a .tex file
write_tex <- function(lines, filename) {
  path <- file.path(OUT_TABS, filename)
  writeLines(lines, path)
  message("  Saved: ", filename)
}

# =============================================================================
# Section 1: Aggregate LP — all five exposure measures, FSI outcome
# tab:agg-lp
# Rows: exposure type; Columns: horizon h = 0..4
# =============================================================================
message("\n=== Section 1: Aggregate LP (tab:agg-lp) ===")

irf_agg <- map_dfr(TREAT_TYPES, function(tt) {
  run_lp(panel, paste0("shock_", tt)) |>
    mutate(treat = tt)
})

# Console summary
cat(sprintf("\n%-14s  %5s  %8s  %6s  %6s\n", "Exposure", "h", "est", "se", "p"))
cat(strrep("-", 44), "\n")
for (tt in TREAT_TYPES) {
  for (h in LP_HORIZONS) {
    r <- filter(irf_agg, treat == tt, h == !!h)
    if (nrow(r) == 0) next
    cat(sprintf("%-14s  h=%d  %+8.3f  %6.3f  %5.3f%s\n",
                TREAT_LABELS[tt], h, r$est, r$se, r$pval,
                stars_console(r$pval)))
  }
  cat("\n")
}

# LaTeX table
# Rows = exposure types; Columns = horizons
n_h      <- length(LP_HORIZONS)
col_spec <- sprintf("l%s", paste(rep("r", n_h), collapse = ""))
h_header <- paste(
  " ",
  paste(sprintf("& $h=%d$", LP_HORIZONS), collapse = " "),
  "\\\\"
)

tab_agg <- c(
  "% Aggregate LP: FSI outcome, five exposure types",
  "",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Aggregate local projections: all five exposure measures (FSI outcome)}",
  "\\label{tab:agg-lp}",
  "\\begin{threeparttable}",
  "\\small",
  sprintf("\\begin{tabular}{%s}", col_spec),
  "\\toprule",
  h_header,
  "\\midrule"
)

for (tt in TREAT_TYPES) {
  d_tt <- filter(irf_agg, treat == tt)
  est_row <- paste(
    map_chr(LP_HORIZONS, function(h) {
      r <- filter(d_tt, h == !!h)
      if (nrow(r) == 0 || is.na(r$est)) return("---")
      sprintf("%.3f%s", r$est, stars_fn(r$pval))
    }),
    collapse = " & "
  )
  se_row <- paste(
    map_chr(LP_HORIZONS, function(h) {
      r <- filter(d_tt, h == !!h)
      if (nrow(r) == 0 || is.na(r$se)) return("")
      sprintf("(%.3f)", r$se)
    }),
    collapse = " & "
  )
  tab_agg <- c(tab_agg,
    sprintf("%s & %s \\\\", TREAT_LABELS[tt], est_row),
    sprintf("   & %s \\\\[4pt]", se_row)
  )
}

# N row: use consumption as representative (same country-year coverage per horizon)
d_cons <- filter(irf_agg, treat == "cons")
n_row <- paste(
  map_chr(LP_HORIZONS, function(h) {
    r <- filter(d_cons, h == !!h)
    if (nrow(r) == 0) "---" else format(r$n_obs, big.mark = ",")
  }),
  collapse = " & "
)
tab_agg <- c(tab_agg, "\\midrule",
             sprintf("$N$ (consumption) & %s \\\\", n_row))

tab_agg <- c(
  tab_agg,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0(
    "\\item \\textit{Notes:} LP: $\\text{FSI}_{i,t+h} - \\text{FSI}_{i,t-1} = ",
    "\\beta\\, x_{it} + \\sum_{\\ell=1}^{2}\\gamma_\\ell x_{i,t-\\ell} + \\alpha_t + \\varepsilon_{iht}$. ",
    "Year fixed effects; standard errors clustered by country. ",
    "Significance: $***$ $p<0.001$, $**$ $p<0.01$, $*$ $p<0.05$, $\\cdot$ $p<0.10$."
  ),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)

write_tex(tab_agg, "tab_agg_lp.tex")

# =============================================================================
# Figure 1: Aggregate IRF — one facet per exposure type
# fig1_aggregate_irf.pdf
# =============================================================================
message("\n=== Figure 1: Aggregate IRF ===")

COL_NEG <- "#2166ac"    # blue — negative mean impulse response
COL_POS <- "#d6604d"    # red  — positive mean impulse response
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
      legend.position     = "bottom",
      plot.tag            = element_text(face = "bold", size = base_size + 1)
    )
}

# Assign color per exposure based on sign of mean estimate across horizons
fig1_dat <- irf_agg |>
  group_by(treat) |>
  mutate(mean_est = mean(est, na.rm = TRUE)) |>
  ungroup() |>
  mutate(
    treat_label = factor(TREAT_LABELS[treat], levels = TREAT_LABELS),
    col         = if_else(mean_est < 0, COL_NEG, COL_POS)
  )

fig1 <- ggplot(fig1_dat, aes(x = h, y = est)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.35) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              fill = COL_RIB, alpha = 0.25) +
  geom_line(aes(colour = col), linewidth = 0.7) +
  geom_point(aes(colour = col), size = 1.6) +
  scale_colour_identity() +
  scale_x_continuous(breaks = 0:4, name = "Horizon (years)") +
  scale_y_continuous(name = "FSI change (index points)") +
  facet_wrap(~ treat_label, nrow = 2, scales = "free_y") +
  theme_nature()

ggsave(file.path(OUT_FIGS, "fig1_aggregate_irf.pdf"),
       fig1, width = 180, height = 90, units = "mm")
ggsave(file.path(OUT_FIGS, "fig1_aggregate_irf.png"),
       fig1, width = 180, height = 90, units = "mm", dpi = 300)

message("  Saved: fig1_aggregate_irf.pdf / .png")
