# =============================================================================
# 11_income_dl.R
# Distributed lag models by World Bank income group
#
# Mirrors 10_income_lp.R but uses the DL specification from 09_dl.R:
#   ΔFSI_{i,t} = α_t + β₀ shock_{i,t} + ⋯ + β_H shock_{i,t-H} + ε_{i,t}
# IRF at horizon h: cumulative multiplier B_h = Σ_{k=0}^{h} β_k,
# delta-method SEs from the clustered vcov.
#
# Outputs:
#   output/tables/tab_income_dl.tex
#   output/figures/fig4_income_dl_irf.pdf/.png
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

INCOME_LEVELS <- c("Low", "Lower-middle", "Upper-middle", "High")
INCOME_LABELS <- c("Low\nincome", "Lower-middle\nincome",
                   "Upper-middle\nincome", "High\nincome")

CI_Z <- 1.645

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
# Section 1: Run DL for each income group × exposure type
# =============================================================================
message("\n=== Section 1: Income group DL ===")

if (!"income_group" %in% names(panel)) {
  stop("income_group not in panel — run 07_panel.R with income_groups.rds present.")
}

panel_lmic <- panel |>
  filter(income_group %in% c("Low", "Lower-middle")) |>
  mutate(income_group = "LMIC")

panel_by_group <- c(
  setNames(
    map(INCOME_LEVELS, ~ filter(panel, income_group == .x)),
    INCOME_LEVELS
  ),
  list(LMIC = panel_lmic)
)

walk2(names(panel_by_group), panel_by_group, function(grp, d) {
  message("  ", grp, ": ", n_distinct(d$iso3c), " countries")
})

irf_inc_dl <- map_dfr(names(panel_by_group), function(grp) {
  d_grp <- panel_by_group[[grp]]
  if (n_distinct(d_grp$iso3c) < 5) {
    message("  Skipping ", grp, " (too few countries)")
    return(NULL)
  }
  map_dfr(TREAT_TYPES, function(tt) {
    run_dl(d_grp, paste0("shock_", tt)) |>
      mutate(treat = tt, income_group = grp)
  })
})

# Console summary: h=1 and h=4
SHOW_H <- c(1, 4)
cat(sprintf("\n%-14s  %-14s  %5s  %8s  %6s  %5s\n",
            "Exposure", "Income group", "h", "est", "se", "p"))
cat(strrep("-", 58), "\n")
for (tt in TREAT_TYPES) {
  for (grp in INCOME_LEVELS) {
    for (h_ in SHOW_H) {
      r <- filter(irf_inc_dl, treat == tt, income_group == grp, h == h_)
      if (nrow(r) == 0) next
      cat(sprintf("%-14s  %-14s  h=%d  %+8.3f  %6.3f  %4.3f%s\n",
                  TREAT_LABELS[tt], grp, h_,
                  r$est, r$se, r$pval, stars_console(r$pval)))
    }
  }
  cat("\n")
}

# =============================================================================
# Section 2: LaTeX table  (tab:income-dl)
# =============================================================================
message("\n=== Section 2: LaTeX table (tab_income_dl.tex) ===")

col_spec <- "l rr rr rr rr"

tab_inc_dl <- c(
  "% Distributed lag IRFs by income group: all exposure types, h=1 and h=4",
  "",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Distributed lag models by income group (FSI outcome)}",
  "\\label{tab:income-dl}",
  "\\begin{threeparttable}",
  "\\small",
  sprintf("\\begin{tabular}{%s}", col_spec),
  "\\toprule",
  paste0(
    " & ",
    paste(sprintf("\\multicolumn{2}{c}{%s}", INCOME_LEVELS), collapse = " & "),
    " \\\\"
  ),
  paste0(
    paste(sprintf("\\cmidrule(lr){%d-%d}", seq(2, 9, 2), seq(3, 10, 2)),
          collapse = " ")
  ),
  paste0(
    " & ",
    paste(rep("$h=1$ & $h=4$", 4), collapse = " & "),
    " \\\\"
  ),
  "\\midrule"
)

for (tt in TREAT_TYPES) {
  d_tt <- filter(irf_inc_dl, treat == tt)

  est_cells <- map_chr(INCOME_LEVELS, function(grp) {
    r1 <- filter(d_tt, income_group == grp, h == 1)
    r4 <- filter(d_tt, income_group == grp, h == 4)
    c1 <- if (nrow(r1) > 0 && !is.na(r1$est))
            sprintf("%.3f%s", r1$est, stars_fn(r1$pval)) else "---"
    c4 <- if (nrow(r4) > 0 && !is.na(r4$est))
            sprintf("%.3f%s", r4$est, stars_fn(r4$pval)) else "---"
    paste(c1, c4, sep = " & ")
  })

  se_cells <- map_chr(INCOME_LEVELS, function(grp) {
    r1 <- filter(d_tt, income_group == grp, h == 1)
    r4 <- filter(d_tt, income_group == grp, h == 4)
    s1 <- if (nrow(r1) > 0 && !is.na(r1$se)) sprintf("(%.3f)", r1$se) else ""
    s4 <- if (nrow(r4) > 0 && !is.na(r4$se)) sprintf("(%.3f)", r4$se) else ""
    paste(s1, s4, sep = " & ")
  })

  tab_inc_dl <- c(tab_inc_dl,
    sprintf("%s & %s \\\\", TREAT_LABELS[tt], paste(est_cells, collapse = " & ")),
    sprintf("   & %s \\\\[4pt]", paste(se_cells, collapse = " & "))
  )
}

# N: same across h=1 and h=4 within each group (single regression per cell)
n_cells <- map_chr(INCOME_LEVELS, function(grp) {
  r <- filter(irf_inc_dl, treat == "cons", income_group == grp, h == 0)
  n_val <- if (nrow(r) > 0) format(r$n_obs, big.mark = ",") else "---"
  paste(n_val, n_val, sep = " & ")
})

tab_inc_dl <- c(tab_inc_dl, "\\midrule",
                sprintf("$N$ & %s \\\\", paste(n_cells, collapse = " & ")))

tab_inc_dl <- c(
  tab_inc_dl,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0(
    "\\item \\textit{Notes:} Distributed lag specification: ",
    "$\\Delta\\text{FSI}_{i,t} = \\alpha_t + \\sum_{k=0}^{H}\\beta_k\\,x_{i,t-k} + \\varepsilon_{it}$. ",
    "Entries are cumulative multipliers $B_h = \\sum_{k=0}^{h}\\beta_k$ with ",
    "delta-method standard errors in parentheses. ",
    "Samples split by World Bank income classification (2010, fixed). ",
    "Year fixed effects; standard errors clustered by country. ",
    "$N$ is the same across both horizon columns (single regression per cell). ",
    "Significance: $***$ $p<0.001$, $**$ $p<0.01$, $*$ $p<0.05$, $\\cdot$ $p<0.10$."
  ),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)

write_tex(tab_inc_dl, "tab_income_dl.tex")

# =============================================================================
# Section 3: Income group IRF figure (DL)
# =============================================================================
message("\n=== Section 3: Income group DL IRF figure ===")

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

fig4_dat <- irf_inc_dl |>
  filter(treat %in% c("cons", "prod"),
         income_group %in% INCOME_LEVELS) |>
  mutate(
    income_f = factor(income_group,
                      levels = INCOME_LEVELS,
                      labels = INCOME_LABELS),
    treat_f  = factor(treat,
                      levels = c("cons", "prod"),
                      labels = c("Consumption exposure",
                                 "Production exposure"))
  )

fig4_sign <- fig4_dat |>
  group_by(treat_f, income_f) |>
  summarise(col = if_else(mean(est, na.rm = TRUE) < 0, COL_NEG, COL_POS),
            .groups = "drop")

fig4_dat <- fig4_dat |>
  left_join(fig4_sign, by = c("treat_f", "income_f"))

fig4 <- ggplot(fig4_dat, aes(x = h, y = est)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.35) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              fill = COL_RIB, alpha = 0.25) +
  geom_line(aes(colour = col), linewidth = 0.7) +
  geom_point(aes(colour = col), size = 1.6) +
  scale_colour_identity() +
  scale_x_continuous(breaks = 0:4, name = "Horizon (years)") +
  scale_y_continuous(name = "FSI change (index points)") +
  facet_grid(treat_f ~ income_f, scales = "free_y") +
  theme_nature() +
  theme(strip.text.y = element_text(angle = 0, hjust = 0))

ggsave(file.path(OUT_FIGS, "fig4_income_dl_irf.pdf"),
       fig4, width = 180, height = 100, units = "mm")
ggsave(file.path(OUT_FIGS, "fig4_income_dl_irf.png"),
       fig4, width = 180, height = 100, units = "mm", dpi = 300)

message("  Saved: fig4_income_dl_irf.pdf / .png")
