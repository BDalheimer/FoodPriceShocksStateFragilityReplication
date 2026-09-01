# =============================================================================
# 27_leave_one_out.R
# Leave-one-country-out robustness for the main LP estimates.
#
# Concern: large or politically unusual countries (China, India, USA, …)
# may dominate the aggregate estimate. This script re-runs the baseline LP
# for each treatment type dropping one country at a time and collects the
# resulting distribution of coefficients.
#
# Specification: same as 08_lp.R (year FE, 2 lags, cluster by country).
#
# Outputs:
#   output/figures/fig_loo_cons.pdf/.png   — LOO distribution, cons, h=1 & h=4
#   output/figures/fig_loo_prod.pdf/.png   — LOO distribution, prod, h=1 & h=4
#   output/tables/tab_loo_influence.tex    — top influential countries (by |Δest|)
# =============================================================================

source(here::here("R", "00_setup.R"))

panel <- readRDS(file.path(DATA_PRO, "panel.rds"))

message("Panel: ", nrow(panel), " obs | ",
        n_distinct(panel$iso3c), " countries | ",
        min(panel$year), "–", max(panel$year))

TREAT_TYPES  <- c("cons", "prod")
TREAT_LABELS <- c(cons = "Consumption exposure", prod = "Production exposure")
SHOW_H       <- c(1, 4)      # horizons to report in table
CI_Z         <- 1.645

stars_fn <- function(p)
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", p < 0.10 ~ "$\\cdot$", TRUE ~ "")

# ---- LP runner --------------------------------------------------------------
run_lp_h <- function(data, shock_var, h) {
  y_col <- paste0("dfsi_h", h)
  lag1  <- paste0(shock_var, "_l1")
  lag2  <- paste0(shock_var, "_l2")
  if (!y_col %in% names(data) || !shock_var %in% names(data)) return(NULL)
  d <- data |>
    rename(.y = all_of(y_col), .s = all_of(shock_var)) |>
    filter(!is.na(.y), !is.na(.s))
  if (nrow(d) < 20 || n_distinct(d$iso3c) < 5) return(NULL)
  if (lag1 %in% names(d)) d <- rename(d, .sl1 = all_of(lag1))
  if (lag2 %in% names(d)) d <- rename(d, .sl2 = all_of(lag2))
  rhs <- intersect(c(".sl1", ".sl2"), names(d))
  fml <- as.formula(paste0(".y ~ .s",
    if (length(rhs)) paste0(" + ", paste(rhs, collapse = " + ")),
    " | year"))
  fit <- tryCatch(feols(fml, data = d, cluster = ~iso3c, warn = FALSE, notes = FALSE),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  ct <- summary(fit)$coeftable
  if (!".s" %in% rownames(ct)) return(NULL)
  tibble(est  = ct[".s", "Estimate"],
         se   = ct[".s", "Std. Error"],
         pval = ct[".s", "Pr(>|t|)"],
         n_obs = nobs(fit))
}

# =============================================================================
# 1. Baseline estimates
# =============================================================================
message("\n=== Baseline estimates ===")
baseline <- map_dfr(TREAT_TYPES, function(tt) {
  sv <- paste0("shock_", tt)
  map_dfr(LP_HORIZONS, function(h) {
    r <- run_lp_h(panel, sv, h)
    if (is.null(r)) return(NULL)
    r |> mutate(treat = tt, h = h, iso3c = "FULL")
  })
})

cat(sprintf("%-6s  %-4s  %+9s  %7s  %6s\n", "treat", "h", "est", "se", "p"))
cat(strrep("-", 40), "\n")
baseline |> filter(h %in% SHOW_H) |>
  arrange(treat, h) |>
  with(mapply(function(tt, h_, e, s, p)
    cat(sprintf("%-6s  h=%d  %+9.4f  %7.4f  %5.3f\n", tt, h_, e, s, p)),
    treat, h, est, se, pval))

# =============================================================================
# 2. Leave-one-out loop
# =============================================================================
countries <- sort(unique(panel$iso3c))
message("\n=== LOO loop: ", length(countries), " countries × ",
        length(TREAT_TYPES), " treatments × ", length(LP_HORIZONS), " horizons ===")

pb_step <- max(1L, floor(length(countries) / 20))

loo_results <- map_dfr(seq_along(countries), function(i) {
  iso <- countries[i]
  if (i %% pb_step == 0)
    message("  ", i, "/", length(countries), " — dropping ", iso)
  d_loo <- filter(panel, iso3c != iso)
  map_dfr(TREAT_TYPES, function(tt) {
    sv <- paste0("shock_", tt)
    map_dfr(LP_HORIZONS, function(h) {
      r <- run_lp_h(d_loo, sv, h)
      if (is.null(r)) return(NULL)
      r |> mutate(treat = tt, h = h, iso3c = iso)
    })
  })
})

message("LOO complete: ", nrow(loo_results), " rows")

# =============================================================================
# 3. Influence measure: |est_loo - est_baseline| at each h
# =============================================================================
influence <- loo_results |>
  left_join(baseline |> select(treat, h, est_base = est, pval_base = pval),
            by = c("treat", "h")) |>
  mutate(
    delta    = est - est_base,
    sign_flip = sign(est) != sign(est_base),
    sig_base  = pval_base < 0.10,
    sig_loo   = pval    < 0.10
  )

# Top influential countries at h=4, by |delta|, consumption
cat("\n=== Top 15 influential countries (cons, h=4, |Δest|) ===\n")
influence |>
  filter(treat == "cons", h == 4) |>
  arrange(desc(abs(delta))) |>
  slice_head(n = 15) |>
  mutate(flip = if_else(sign_flip, "FLIP", "")) |>
  with(mapply(function(iso, e, b, d, fl)
    cat(sprintf("  %-6s  loo=%+.4f  base=%+.4f  delta=%+.4f  %s\n",
                iso, e, b, d, fl)),
    iso3c, est, est_base, delta, flip))

cat("\n=== Top 15 influential countries (prod, h=4, |Δest|) ===\n")
influence |>
  filter(treat == "prod", h == 4) |>
  arrange(desc(abs(delta))) |>
  slice_head(n = 15) |>
  mutate(flip = if_else(sign_flip, "FLIP", "")) |>
  with(mapply(function(iso, e, b, d, fl)
    cat(sprintf("  %-6s  loo=%+.4f  base=%+.4f  delta=%+.4f  %s\n",
                iso, e, b, d, fl)),
    iso3c, est, est_base, delta, flip))

# =============================================================================
# 4. Figure: LOO distribution at h=1 and h=4, cons and prod
# =============================================================================
message("\n=== Figures ===")

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
      legend.position     = "none"
    )
}

make_loo_fig <- function(tt) {
  dat <- influence |>
    filter(treat == tt, h %in% SHOW_H) |>
    mutate(
      h_f     = factor(h, levels = SHOW_H,
                       labels = paste0("h = ", SHOW_H)),
      # sort countries by loo estimate at h=4
      iso3c   = factor(iso3c,
                       levels = influence |>
                         filter(treat == tt, h == 4) |>
                         arrange(est) |>
                         pull(iso3c))
    )
  base_dat <- baseline |>
    filter(treat == tt, h %in% SHOW_H) |>
    mutate(h_f = factor(h, levels = SHOW_H, labels = paste0("h = ", SHOW_H)))

  # Identify countries whose exclusion flips sign at h=4
  flip_isos <- influence |>
    filter(treat == tt, h == 4, sign_flip) |>
    pull(iso3c)

  dat <- dat |>
    mutate(col = case_when(
      iso3c %in% flip_isos ~ "#d6604d",
      TRUE                 ~ "grey60"
    ))

  ggplot(dat, aes(x = est, y = after_stat(density))) +
    geom_histogram(bins = 40, fill = "grey70", colour = "white",
                   linewidth = 0.2) +
    geom_vline(data = base_dat, aes(xintercept = est),
               colour = "#2166ac", linewidth = 0.8, linetype = "solid") +
    geom_vline(xintercept = 0, colour = "grey40",
               linewidth = 0.4, linetype = "dashed") +
    scale_x_continuous(name = "LP coefficient (FSI pts)") +
    scale_y_continuous(name = "Density") +
    facet_wrap(~ h_f, scales = "free_x", ncol = 2) +
    ggtitle(TREAT_LABELS[tt]) +
    theme_nature()
}

fig_loo_cons <- make_loo_fig("cons")
fig_loo_prod <- make_loo_fig("prod")

ggsave(file.path(OUT_FIGS, "fig_loo_cons.pdf"),
       fig_loo_cons, width = 120, height = 60, units = "mm")
ggsave(file.path(OUT_FIGS, "fig_loo_cons.png"),
       fig_loo_cons, width = 120, height = 60, units = "mm", dpi = 300)
ggsave(file.path(OUT_FIGS, "fig_loo_prod.pdf"),
       fig_loo_prod, width = 120, height = 60, units = "mm")
ggsave(file.path(OUT_FIGS, "fig_loo_prod.png"),
       fig_loo_prod, width = 120, height = 60, units = "mm", dpi = 300)
message("  Saved: fig_loo_cons/prod .pdf / .png")

# =============================================================================
# 5. Table: top influential countries at h=4, cons + prod
# =============================================================================
N_TOP <- 10

make_top <- function(tt, n = N_TOP) {
  influence |>
    filter(treat == tt, h == 4) |>
    arrange(desc(abs(delta))) |>
    slice_head(n = n) |>
    mutate(
      flip_str = if_else(sign_flip, "yes", "no"),
      sig_str  = case_when(
        sig_base & sig_loo  ~ "both",
        sig_base & !sig_loo ~ "lost",
        !sig_base & sig_loo ~ "gained",
        TRUE                ~ "neither"
      ),
      est_base_str = sprintf("$%+.3f%s$", est_base, stars_fn(pval_base)),
      est_loo_str  = sprintf("$%+.3f%s$", est, stars_fn(pval)),
      delta_str    = sprintf("$%+.3f$", delta)
    ) |>
    select(iso3c, est_base_str, est_loo_str, delta_str, flip_str, sig_str)
}

top_cons <- make_top("cons")
top_prod <- make_top("prod")

tab <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Leave-one-out influence: top countries by $|\\Delta\\hat\\beta|$ at $h=4$}",
  "\\label{tab:loo-influence}",
  "\\begin{threeparttable}\\small",
  "\\begin{tabular}{l rrrll}",
  "\\toprule",
  "Country & Baseline & LOO & $\\Delta$ & Sign flip & Significance \\\\",
  "\\midrule",
  "\\multicolumn{6}{l}{\\textit{Consumption exposure}} \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(top_cons))) {
  r <- top_cons[i, ]
  tab <- c(tab, sprintf("%s & %s & %s & %s & %s & %s \\\\",
    r$iso3c, r$est_base_str, r$est_loo_str, r$delta_str,
    r$flip_str, r$sig_str))
}
tab <- c(tab,
  "\\midrule",
  "\\multicolumn{6}{l}{\\textit{Production exposure}} \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(top_prod))) {
  r <- top_prod[i, ]
  tab <- c(tab, sprintf("%s & %s & %s & %s & %s & %s \\\\",
    r$iso3c, r$est_base_str, r$est_loo_str, r$delta_str,
    r$flip_str, r$sig_str))
}
tab <- c(tab,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  paste0("\\item \\textit{Notes:} Each row reports the baseline $h=4$ ",
         "LP coefficient, the coefficient when the indicated country is ",
         "excluded, and the difference $\\Delta = \\hat\\beta_{\\text{LOO}} - ",
         "\\hat\\beta_{\\text{full}}$. ``Sign flip'': exclusion reverses the ",
         "sign of the estimate. ``Significance'' reports whether the 10\\% ",
         "significance status changes (both / lost / gained / neither). ",
         "Year FE; SE clustered by country. ",
         "$^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\cdot}p<0.10$."),
  "\\end{tablenotes}\\end{threeparttable}\\end{table}"
)
writeLines(tab, file.path(OUT_TABS, "tab_loo_influence.tex"))
message("  Saved: tab_loo_influence.tex")

# Save full LOO results for further inspection
saveRDS(list(baseline = baseline, loo = loo_results, influence = influence),
        file.path(DATA_PRO, "loo_results.rds"))
message("  Saved: data/processed/loo_results.rds")
message("\nDone.")
