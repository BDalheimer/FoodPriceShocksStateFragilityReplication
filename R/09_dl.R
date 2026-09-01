# =============================================================================
# 09_dp.R
# Distributed lag (DP) models: single-equation cumulative multipliers
#
# Specification (per exposure type):
#   ΔFSI_{i,t} = α_t + β₀ shock_{i,t} + β₁ shock_{i,t-1} + ⋯ + β_H shock_{i,t-H}
#                + ε_{i,t}
#   Year FE; SE clustered by country.
#
# IRF at horizon h: cumulative multiplier  B_h = Σ_{k=0}^{h} β_k
# SE(B_h): delta method — sqrt(ι'Σ_hι), ι = ones vector, Σ_h = clustered vcov
#          for the first (h+1) lag coefficients
#
# Comparable to LP in 08_lp.R: same horizons, same sample, same exposure types.
# Key distinction: one regression rather than H+1; cumulative rather than direct
# projection; all-horizon precision comes from off-diagonal covariance terms.
# =============================================================================

source(here::here("R", "00_setup.R"))

panel <- readRDS(file.path(DATA_PRO, "panel.rds"))

message("Panel: ", nrow(panel), " obs | ",
        n_distinct(panel$iso3c), " countries | ",
        min(panel$year), "–", max(panel$year))

# ---- Constants --------------------------------------------------------------
TREAT_TYPES  <- c("cons", "prod", "imp", "exp", "net")
TREAT_LABELS <- c(cons = "Consumption",
                  prod = "Production",
                  imp  = "Imports",
                  exp  = "Exports",
                  net  = "Net imports")

CI_Z <- 1.645   # 90% CI

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

write_tex <- function(lines, filename) {
  path <- file.path(OUT_TABS, filename)
  writeLines(lines, path)
  message("  Saved: ", filename)
}

# run_dl: distributed lag model; returns cumulative multipliers at each horizon
# Outcome: ΔFSI_t = dfsi_h0 (FSI_t - FSI_{t-1}).
# RHS: shock_{i,t}, shock_{i,t-1}, ..., shock_{i,t-H} plus year FE.
run_dl <- function(data, shock_var, fd_stem = "dfsi", horizons = LP_HORIZONS) {
  max_h    <- max(horizons)
  y_col    <- paste0(fd_stem, "_h0")
  lag_vars <- c(shock_var, paste0(shock_var, "_l", seq_len(max_h)))
  lag_vars <- lag_vars[lag_vars %in% names(data)]   # keep only available lags

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
  vc         <- vcov(fit)                              # clustered by iso3c (from feols call)
  coef_names <- lag_vars[lag_vars %in% rownames(ct)]  # ordered lag names in ct
  n          <- nobs(fit)

  map_dfr(horizons, function(h) {
    idx <- coef_names[seq_len(min(h + 1L, length(coef_names)))]
    if (length(idx) == 0) return(NULL)
    B_h  <- sum(ct[idx, "Estimate"])
    se_h <- sqrt(sum(vc[idx, idx]))    # Var(Σ β_k) = ι'Σι
    pval <- 2 * pnorm(-abs(B_h / se_h))
    tibble(
      h     = h,
      est   = B_h,
      se    = se_h,
      pval  = pval,
      ci_lo = B_h - CI_Z * se_h,
      ci_hi = B_h + CI_Z * se_h,
      n_obs = n
    )
  })
}

# =============================================================================
# Section 1: Distributed lag IRFs — all exposure types
# tab:agg-dl
# =============================================================================
message("\n=== Section 1: Distributed lag IRFs (tab:agg-dl) ===")

irf_dl <- map_dfr(TREAT_TYPES, function(tt) {
  run_dl(panel, paste0("shock_", tt)) |>
    mutate(treat = tt)
})

# Console summary
cat(sprintf("\n%-14s  %5s  %8s  %6s  %6s\n", "Exposure", "h", "est", "se", "p"))
cat(strrep("-", 44), "\n")
for (tt in TREAT_TYPES) {
  for (h in LP_HORIZONS) {
    r <- filter(irf_dl, treat == tt, h == !!h)
    if (nrow(r) == 0) next
    cat(sprintf("%-14s  h=%d  %+8.3f  %6.3f  %5.3f%s\n",
                TREAT_LABELS[tt], h, r$est, r$se, r$pval,
                stars_console(r$pval)))
  }
  cat("\n")
}

# LaTeX table
n_h      <- length(LP_HORIZONS)
col_spec <- sprintf("l%s", paste(rep("r", n_h), collapse = ""))
h_header <- paste(
  " ",
  paste(sprintf("& $h=%d$", LP_HORIZONS), collapse = " "),
  "\\\\"
)

tab_dl <- c(
  "% Distributed lag IRFs: FSI outcome, five exposure types",
  "",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Distributed lag models: cumulative impulse responses (FSI outcome)}",
  "\\label{tab:agg-dl}",
  "\\begin{threeparttable}",
  "\\small",
  sprintf("\\begin{tabular}{%s}", col_spec),
  "\\toprule",
  h_header,
  "\\midrule"
)

for (tt in TREAT_TYPES) {
  d_tt <- filter(irf_dl, treat == tt)
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
  tab_dl <- c(tab_dl,
    sprintf("%s & %s \\\\", TREAT_LABELS[tt], est_row),
    sprintf("   & %s \\\\[4pt]", se_row)
  )
}

# N: same across all horizons (single regression per exposure type)
d_cons <- filter(irf_dl, treat == "cons", h == 0)
n_val  <- if (nrow(d_cons) > 0) format(d_cons$n_obs, big.mark = ",") else "---"
n_row  <- paste(rep(n_val, n_h), collapse = " & ")
tab_dl <- c(tab_dl, "\\midrule",
            sprintf("$N$ & %s \\\\", n_row))

tab_dl <- c(
  tab_dl,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0(
    "\\item \\textit{Notes:} Distributed lag specification: ",
    "$\\Delta\\text{FSI}_{i,t} = \\alpha_t + \\sum_{k=0}^{H}\\beta_k\\,x_{i,t-k} + \\varepsilon_{it}$. ",
    "Entries are cumulative multipliers $B_h = \\sum_{k=0}^{h}\\beta_k$ with ",
    "delta-method standard errors in parentheses ",
    "(using the full clustered variance-covariance matrix). ",
    "Year fixed effects; standard errors clustered by country. ",
    "$N$ is the same across all horizons (single regression per row). ",
    "Significance: $***$ $p<0.001$, $**$ $p<0.01$, $*$ $p<0.05$, $\\cdot$ $p<0.10$."
  ),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)

write_tex(tab_dl, "tab_agg_dl.tex")

# =============================================================================
# Figure: Distributed lag IRFs — one facet per exposure type
# fig2_dl_irf.pdf
# =============================================================================
message("\n=== Figure: Distributed lag IRF ===")

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
      legend.position     = "bottom",
      plot.tag            = element_text(face = "bold", size = base_size + 1)
    )
}

fig_dl_dat <- irf_dl |>
  group_by(treat) |>
  mutate(mean_est = mean(est, na.rm = TRUE)) |>
  ungroup() |>
  mutate(
    treat_label = factor(TREAT_LABELS[treat], levels = TREAT_LABELS),
    col         = if_else(mean_est < 0, COL_NEG, COL_POS)
  )

fig_dl <- ggplot(fig_dl_dat, aes(x = h, y = est)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.35) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              fill = COL_RIB, alpha = 0.25) +
  geom_line(aes(colour = col), linewidth = 0.7) +
  geom_point(aes(colour = col), size = 1.6) +
  scale_colour_identity() +
  scale_x_continuous(breaks = 0:4, name = "Horizon (years)") +
  scale_y_continuous(name = "FSI change (index points)") +
  facet_wrap(~ treat_label, nrow = 1, scales = "free_y") +
  theme_nature()

ggsave(file.path(OUT_FIGS, "fig2_dl_irf.pdf"),
       fig_dl, width = 180, height = 55, units = "mm")
ggsave(file.path(OUT_FIGS, "fig2_dl_irf.png"),
       fig_dl, width = 180, height = 55, units = "mm", dpi = 300)

message("  Saved: fig2_dl_irf.pdf / .png")
