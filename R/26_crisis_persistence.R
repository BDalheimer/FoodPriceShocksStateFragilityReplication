# =============================================================================
# 26_crisis_persistence.R
# Robustness: crisis-episode heterogeneity and shock persistence decomposition.
#
# Addresses reviewer comment on temporal heterogeneity across food price crises
# and Dube & Vargas (2013) observation that persistent vs transitory price
# movements have qualitatively different political economy effects.
#
# Two exercises on the aggregate (full sample) panel:
#
# A. Crisis-period interaction
#    Crisis years: 2007, 2008, 2010, 2011, 2022 (the three major spike episodes)
#    shock_crisis = shock × I(crisis year)
#    shock_normal = shock × I(non-crisis year)
#    LP includes both as separate regressors plus their lags.
#
# B. Persistent vs transitory decomposition
#    Persistent = 3-year backward moving average of the shock
#    Transitory  = shock minus the 3-year MA (deviation from trend)
#    Both lags already available from existing l1–l4 columns.
#
# Both exercises run for consumption and production exposure.
#
# Outputs:
#   output/tables/tab_crisis_interaction.tex   — crisis vs normal, cons+prod, h=0..4
#   output/tables/tab_persistence.tex          — perm vs trans, cons+prod, h=0..4
#   output/figures/fig_crisis_irf.pdf/.png     — IRF overlay: crisis vs normal
#   output/figures/fig_persistence_irf.pdf/.png — IRF overlay: perm vs trans
# =============================================================================

source(here::here("R", "00_setup.R"))

panel <- readRDS(file.path(DATA_PRO, "panel.rds"))

message("Panel: ", nrow(panel), " obs | ",
        n_distinct(panel$iso3c), " countries | ",
        min(panel$year), "–", max(panel$year))

CI_Z <- 1.645

CRISIS_YRS <- c(2007, 2008, 2010, 2011, 2022)
message("Crisis years: ", paste(CRISIS_YRS, collapse = ", "),
        " (", sum(panel$year %in% CRISIS_YRS), " obs)")

stars_fn <- function(p)
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", p < 0.10 ~ "$\\cdot$", TRUE ~ "")

# =============================================================================
# 1. Build new variables
# =============================================================================
panel <- panel |>
  group_by(iso3c) |>
  arrange(year) |>
  mutate(
    crisis_t = as.integer(year %in% CRISIS_YRS),

    # A: Crisis interaction
    shock_cons_crisis = shock_cons * crisis_t,
    shock_cons_normal = shock_cons * (1L - crisis_t),
    shock_prod_crisis = shock_prod * crisis_t,
    shock_prod_normal = shock_prod * (1L - crisis_t),

    shock_cons_crisis_l1 = shock_cons_l1 * lag(crisis_t, 1),
    shock_cons_crisis_l2 = shock_cons_l2 * lag(crisis_t, 2),
    shock_cons_normal_l1 = shock_cons_l1 * (1L - lag(crisis_t, 1)),
    shock_cons_normal_l2 = shock_cons_l2 * (1L - lag(crisis_t, 2)),
    shock_prod_crisis_l1 = shock_prod_l1 * lag(crisis_t, 1),
    shock_prod_crisis_l2 = shock_prod_l2 * lag(crisis_t, 2),
    shock_prod_normal_l1 = shock_prod_l1 * (1L - lag(crisis_t, 1)),
    shock_prod_normal_l2 = shock_prod_l2 * (1L - lag(crisis_t, 2)),

    # B: Persistent (3-year backward MA) and transitory (deviation from MA)
    shock_cons_perm  = (shock_cons + shock_cons_l1 + shock_cons_l2) / 3,
    shock_cons_trans = shock_cons - shock_cons_perm,
    shock_prod_perm  = (shock_prod + shock_prod_l1 + shock_prod_l2) / 3,
    shock_prod_trans = shock_prod - shock_prod_perm,

    shock_cons_perm_l1  = (shock_cons_l1 + shock_cons_l2 + shock_cons_l3) / 3,
    shock_cons_perm_l2  = (shock_cons_l2 + shock_cons_l3 + shock_cons_l4) / 3,
    shock_cons_trans_l1 = shock_cons_l1 - shock_cons_perm_l1,
    shock_cons_trans_l2 = shock_cons_l2 - shock_cons_perm_l2,
    shock_prod_perm_l1  = (shock_prod_l1 + shock_prod_l2 + shock_prod_l3) / 3,
    shock_prod_perm_l2  = (shock_prod_l2 + shock_prod_l3 + shock_prod_l4) / 3,
    shock_prod_trans_l1 = shock_prod_l1 - shock_prod_perm_l1,
    shock_prod_trans_l2 = shock_prod_l2 - shock_prod_perm_l2
  ) |>
  ungroup()

# =============================================================================
# 2. LP runner (joint: two treatments in one equation + their lags)
# =============================================================================
run_lp_joint <- function(data, sv1, sv2, h) {
  y_col <- paste0("dfsi_h", h)
  svs   <- c(sv1, sv2,
             paste0(sv1, "_l1"), paste0(sv2, "_l1"),
             paste0(sv1, "_l2"), paste0(sv2, "_l2"))
  svs   <- svs[svs %in% names(data)]
  if (!y_col %in% names(data)) return(NULL)
  d <- data[!is.na(data[[y_col]]), ]
  fml <- as.formula(paste(y_col, "~",
                          paste(svs, collapse = " + "), "| year"))
  fit <- tryCatch(feols(fml, data = d, cluster = ~iso3c,
                        warn = FALSE, notes = FALSE),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  ct <- summary(fit)$coeftable
  map_dfr(c(sv1, sv2), function(sv) {
    if (!sv %in% rownames(ct)) return(NULL)
    tibble(
      h      = h,
      term   = sv,
      est    = ct[sv, "Estimate"],
      se     = ct[sv, "Std. Error"],
      pval   = ct[sv, "Pr(>|t|)"],
      ci_lo  = est - CI_Z * se,
      ci_hi  = est + CI_Z * se,
      n_obs  = nobs(fit)
    )
  })
}

# =============================================================================
# 3. Run both exercises
# =============================================================================
message("\n=== Exercise A: Crisis interaction ===")

irf_crisis <- map_dfr(LP_HORIZONS, function(h) {
  bind_rows(
    run_lp_joint(panel, "shock_cons_crisis", "shock_cons_normal", h) |>
      mutate(treat = "cons"),
    run_lp_joint(panel, "shock_prod_crisis", "shock_prod_normal", h) |>
      mutate(treat = "prod")
  )
})

cat("\nCrisis interaction LP:\n")
cat(sprintf("  %-26s  %8s  %8s  %8s  %8s  %8s\n",
            "term", "h=0", "h=1", "h=2", "h=3", "h=4"))
for (tm in unique(irf_crisis$term)) {
  cells <- sapply(LP_HORIZONS, function(h_) {
    r <- filter(irf_crisis, term == tm, h == h_)
    if (nrow(r) == 0) return("   ---  ")
    sprintf("%+7.3f%s", r$est, ifelse(r$pval < 0.10, "*", " "))
  })
  cat(sprintf("  %-26s  %s\n", tm, paste(cells, collapse = "  ")))
}

message("\n=== Exercise B: Persistent vs transitory ===")

irf_persist <- map_dfr(LP_HORIZONS, function(h) {
  bind_rows(
    run_lp_joint(panel, "shock_cons_perm", "shock_cons_trans", h) |>
      mutate(treat = "cons"),
    run_lp_joint(panel, "shock_prod_perm", "shock_prod_trans", h) |>
      mutate(treat = "prod")
  )
})

cat("\nPersistence decomposition LP:\n")
cat(sprintf("  %-26s  %8s  %8s  %8s  %8s  %8s\n",
            "term", "h=0", "h=1", "h=2", "h=3", "h=4"))
for (tm in unique(irf_persist$term)) {
  cells <- sapply(LP_HORIZONS, function(h_) {
    r <- filter(irf_persist, term == tm, h == h_)
    if (nrow(r) == 0) return("   ---  ")
    sprintf("%+7.3f%s", r$est, ifelse(r$pval < 0.10, "*", " "))
  })
  cat(sprintf("  %-26s  %s\n", tm, paste(cells, collapse = "  ")))
}

# =============================================================================
# 4. LaTeX table helpers
# =============================================================================
make_two_row <- function(irf_dat, term_tag, term_lbl) {
  r <- filter(irf_dat, term == term_tag)
  est_cells <- sapply(LP_HORIZONS, function(h_) {
    rx <- filter(r, h == h_)
    if (nrow(rx) == 0) return("---")
    sprintf("$%+.3f%s$", rx$est, stars_fn(rx$pval))
  })
  se_cells <- sapply(LP_HORIZONS, function(h_) {
    rx <- filter(r, h == h_)
    if (nrow(rx) == 0) return("")
    sprintf("$(%.3f)$", rx$se)
  })
  n_ <- filter(r, h == 0)$n_obs
  c(paste0(term_lbl, " & ", paste(est_cells, collapse = " & "), " \\\\"),
    paste0("& ", paste(se_cells, collapse = " & "), " \\\\"),
    paste0("& \\multicolumn{5}{l}{$N = ",
           ifelse(length(n_) > 0, n_[1], "---"), "$} \\\\[2pt]"))
}

notes_str <- paste0(
  "Year FE; SE clustered by country. ",
  "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."
)

# =============================================================================
# 5. Table A: Crisis interaction
# =============================================================================
tex_crisis <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Crisis-episode heterogeneity: LP with crisis--period interaction}",
  "\\label{tab:crisis-interaction}",
  "\\begin{threeparttable}\\small",
  "\\begin{tabular}{lrrrrr}",
  "\\toprule",
  "Exposure & $h=0$ & $h=1$ & $h=2$ & $h=3$ & $h=4$ \\\\",
  "\\midrule",
  "\\textit{Consumption exposure} \\\\",
  make_two_row(irf_crisis, "shock_cons_crisis", "\\quad Crisis years"),
  make_two_row(irf_crisis, "shock_cons_normal", "\\quad Non-crisis years"),
  "\\textit{Production exposure} \\\\",
  make_two_row(irf_crisis, "shock_prod_crisis", "\\quad Crisis years"),
  make_two_row(irf_crisis, "shock_prod_normal", "\\quad Non-crisis years"),
  "\\bottomrule", "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0("\\item \\textit{Notes:} Crisis years are ",
         paste(CRISIS_YRS, collapse = ", "),
         " (the 2007--08 spike, the 2010--11 surge, and the 2022 Ukraine episode). ",
         "Both crisis and non-crisis shocks enter jointly with two lags each. ",
         notes_str),
  "\\end{tablenotes}\\end{threeparttable}\\end{table}"
)
writeLines(tex_crisis, file.path(OUT_TABS, "tab_crisis_interaction.tex"))
message("  Saved: tab_crisis_interaction.tex")

# =============================================================================
# 6. Table B: Persistence decomposition
# =============================================================================
tex_persist <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Shock persistence decomposition: persistent vs transitory LP}",
  "\\label{tab:persistence}",
  "\\begin{threeparttable}\\small",
  "\\begin{tabular}{lrrrrr}",
  "\\toprule",
  "Component & $h=0$ & $h=1$ & $h=2$ & $h=3$ & $h=4$ \\\\",
  "\\midrule",
  "\\textit{Consumption exposure} \\\\",
  make_two_row(irf_persist, "shock_cons_perm",  "\\quad Persistent (3-yr MA)"),
  make_two_row(irf_persist, "shock_cons_trans", "\\quad Transitory (deviation)"),
  "\\textit{Production exposure} \\\\",
  make_two_row(irf_persist, "shock_prod_perm",  "\\quad Persistent (3-yr MA)"),
  make_two_row(irf_persist, "shock_prod_trans", "\\quad Transitory (deviation)"),
  "\\bottomrule", "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0("\\item \\textit{Notes:} The persistent component is the 3-year backward ",
         "moving average of the Bartik shock, $\\bar{x}_{it} = (x_{it} + x_{i,t-1} + x_{i,t-2})/3$. ",
         "The transitory component is the deviation $x_{it} - \\bar{x}_{it}$. ",
         "Both enter jointly with two lags each. ",
         notes_str),
  "\\end{tablenotes}\\end{threeparttable}\\end{table}"
)
writeLines(tex_persist, file.path(OUT_TABS, "tab_persistence.tex"))
message("  Saved: tab_persistence.tex")

# =============================================================================
# 7. Figures
# =============================================================================
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

COL_CRISIS <- "#d6604d"
COL_NORMAL <- "#2166ac"
COL_PERM   <- "#1a9850"
COL_TRANS  <- "#984ea3"

TREAT_LABS <- c(cons = "Consumption exposure", prod = "Production exposure")

# Figure A: Crisis interaction
fig_crisis_dat <- irf_crisis |>
  filter(treat %in% c("cons", "prod")) |>
  mutate(
    component = case_when(
      grepl("crisis", term) ~ "Crisis years",
      grepl("normal", term) ~ "Non-crisis years"
    ),
    treat_f = factor(treat, levels = c("cons", "prod"), labels = TREAT_LABS),
    comp_f  = factor(component, levels = c("Crisis years", "Non-crisis years"))
  )

fig_crisis <- ggplot(fig_crisis_dat, aes(x = h, y = est,
                                          colour = comp_f, fill = comp_f,
                                          linetype = comp_f)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.35) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c("Crisis years"     = COL_CRISIS,
                                  "Non-crisis years" = COL_NORMAL),
                      name = NULL) +
  scale_fill_manual(values   = c("Crisis years"     = COL_CRISIS,
                                  "Non-crisis years" = COL_NORMAL),
                    name = NULL) +
  scale_linetype_manual(values = c("Crisis years"     = "dashed",
                                    "Non-crisis years" = "solid"),
                        name = NULL) +
  scale_x_continuous(breaks = 0:4, name = "Horizon (years)") +
  scale_y_continuous(name = "FSI change (index points)") +
  facet_wrap(~ treat_f, ncol = 2, scales = "free_y") +
  theme_nature()

ggsave(file.path(OUT_FIGS, "fig_crisis_irf.pdf"),
       fig_crisis, width = 150, height = 75, units = "mm")
ggsave(file.path(OUT_FIGS, "fig_crisis_irf.png"),
       fig_crisis, width = 150, height = 75, units = "mm", dpi = 300)
message("  Saved: fig_crisis_irf.pdf / .png")

# Figure B: Persistence decomposition
fig_persist_dat <- irf_persist |>
  filter(treat %in% c("cons", "prod")) |>
  mutate(
    component = case_when(
      grepl("perm",  term) ~ "Persistent (3-yr MA)",
      grepl("trans", term) ~ "Transitory (deviation)"
    ),
    treat_f = factor(treat, levels = c("cons", "prod"), labels = TREAT_LABS),
    comp_f  = factor(component,
                     levels = c("Persistent (3-yr MA)", "Transitory (deviation)"))
  )

fig_persist <- ggplot(fig_persist_dat, aes(x = h, y = est,
                                             colour = comp_f, fill = comp_f,
                                             linetype = comp_f)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.35) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c("Persistent (3-yr MA)"   = COL_PERM,
                                  "Transitory (deviation)" = COL_TRANS),
                      name = NULL) +
  scale_fill_manual(values   = c("Persistent (3-yr MA)"   = COL_PERM,
                                  "Transitory (deviation)" = COL_TRANS),
                    name = NULL) +
  scale_linetype_manual(values = c("Persistent (3-yr MA)"   = "solid",
                                    "Transitory (deviation)" = "dashed"),
                        name = NULL) +
  scale_x_continuous(breaks = 0:4, name = "Horizon (years)") +
  scale_y_continuous(name = "FSI change (index points)") +
  facet_wrap(~ treat_f, ncol = 2, scales = "free_y") +
  theme_nature()

ggsave(file.path(OUT_FIGS, "fig_persistence_irf.pdf"),
       fig_persist, width = 150, height = 75, units = "mm")
ggsave(file.path(OUT_FIGS, "fig_persistence_irf.png"),
       fig_persist, width = 150, height = 75, units = "mm", dpi = 300)
message("  Saved: fig_persistence_irf.pdf / .png")
message("\nDone.")
