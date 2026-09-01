# =============================================================================
# 32_bhj_variance.R
# BHJ (2022) shock-level variance decomposition of the shift-share instrument.
#
# Question: is variation in the Bartik exposure B_it driven by heterogeneity
# in the macro shocks g_jt or by baseline country-level share structures z_ij?
#
# Following Borusyak, Hull & Jaravel (2022), for each year t:
#
#   s_j  = (1/N) sum_i z_ij              -- average share ("aggregate weight")
#   g_bar_t = sum_j s_j * g_jt           -- share-weighted mean shock
#   Var_s(g_t) = sum_j s_j*(g_jt - g_bar_t)^2   -- weighted shock variance
#
# Var_s(g_t) captures how much the macro shocks themselves vary across
# commodities.  If it is low, the Bartik variation is almost entirely share-
# driven (exposure heterogeneity across countries); if it is high, the shocks
# themselves are doing the work.
#
# We also report:
#   - Var(B_it) : observed variance of the Bartik exposure across countries
#   - Cov decomposition: Var(B_it) ~ Var_s(g_t) * Var(z_i) component
#   - Effective number of shocks: 1/HHI  where HHI = sum_j s_j^2
#   - Year-by-year panel of the key statistics
#
# Output:
#   output/tables/tab_bhj_variance.tex
# =============================================================================

source(here::here("R", "00_setup.R"))

# ---- Load data --------------------------------------------------------------
shares_raw  <- readRDS(file.path(DATA_PRO, "fao_shares.rds"))
pinksheet   <- readRDS(file.path(DATA_PRO, "pinksheet_annual.rds"))
panel       <- readRDS(file.path(DATA_PRO, "panel.rds"))

SAMPLE_YRS <- SAMPLE_START:SAMPLE_END

# ---- Shares: time-invariant baseline consumption shares ---------------------
# s_cons is already fixed (identical across years per country-commodity pair).
# Take one row per country-commodity; drop missing shares.
shares <- shares_raw |>
  filter(!is.na(s_cons)) |>
  distinct(iso3c, pink_var, grp, s_cons)

# Restrict to commodities present in the pinksheet
commodities <- intersect(unique(shares$pink_var), unique(pinksheet$pink_var))
shares <- filter(shares, pink_var %in% commodities)

message("Commodities: ", length(commodities))
message("Countries with shares: ", n_distinct(shares$iso3c))

# ---- Shocks: annual log price changes g_jt ----------------------------------
shocks <- pinksheet |>
  filter(year %in% SAMPLE_YRS, pink_var %in% commodities) |>
  left_join(distinct(shares, pink_var, grp), by = "pink_var") |>
  select(year, pink_var, grp, g = dlog_p)

# =============================================================================
# Section 1: Aggregate weights s_j and HHI
# =============================================================================
message("\n=== Section 1: Aggregate weights ===")

N_countries <- n_distinct(shares$iso3c)

agg_weights <- shares |>
  group_by(pink_var, grp) |>
  summarise(s_j = sum(s_cons, na.rm = TRUE) / N_countries, .groups = "drop") |>
  arrange(desc(s_j))

# Normalise so weights sum to 1 (they should already be close, but may not sum
# to exactly 1 if some countries have missing shares for some commodities)
agg_weights <- agg_weights |>
  mutate(s_j_norm = s_j / sum(s_j))

hhi      <- sum(agg_weights$s_j_norm^2)
eff_n    <- 1 / hhi

cat(sprintf("\nAggregate weights (top 10 commodities, normalised):\n"))
cat(sprintf("%-20s  %-12s  %8s\n", "Commodity", "Group", "s_j (norm)"))
cat(strrep("-", 44), "\n")
agg_weights |>
  slice_head(n = 10) |>
  with(invisible(mapply(function(v, g, s)
    cat(sprintf("%-20s  %-12s  %8.4f\n", v, g, s)),
    pink_var, grp, s_j_norm)))

cat(sprintf("\nHHI (share concentration) = %.4f\n", hhi))
cat(sprintf("Effective number of shocks (1/HHI) = %.2f\n", eff_n))

# =============================================================================
# Section 2: Year-by-year BHJ weighted shock variance
# =============================================================================
message("\n=== Section 2: Year-by-year Var_s(g_t) ===")

# For each year: merge weights onto shocks, compute g_bar and Var_s
bhj_annual <- shocks |>
  left_join(agg_weights |> select(pink_var, s_j_norm), by = "pink_var") |>
  filter(!is.na(g), !is.na(s_j_norm)) |>
  group_by(year) |>
  mutate(g_bar = sum(s_j_norm * g, na.rm = TRUE)) |>
  summarise(
    g_bar      = first(g_bar),
    var_s_g    = sum(s_j_norm * (g - g_bar)^2, na.rm = TRUE),   # weighted variance of shocks
    sd_s_g     = sqrt(var_s_g),
    n_comm     = n(),
    .groups    = "drop"
  )

# Observed Bartik variance (across countries, within year)
bartik_obs <- shares |>
  left_join(shocks, by = "pink_var", relationship = "many-to-many") |>
  filter(year %in% SAMPLE_YRS, !is.na(g)) |>
  group_by(iso3c, year) |>
  summarise(B_it = sum(s_cons * g, na.rm = TRUE), .groups = "drop") |>
  group_by(year) |>
  summarise(var_B = var(B_it, na.rm = TRUE),
            sd_B  = sd(B_it, na.rm = TRUE),
            mean_B = mean(B_it, na.rm = TRUE),
            .groups = "drop")

bhj_annual <- bhj_annual |>
  left_join(bartik_obs, by = "year") |>
  mutate(
    # Ratio SD_s(g) / SD(B_it): how much the shocks themselves vary relative to
    # the observed Bartik exposure cross-sectional variation. A high ratio means
    # the share averaging compresses the cross-country variation substantially,
    # i.e. share heterogeneity, not shock variation, is the binding constraint.
    shock_to_bartik = sd_s_g / (sd_B + 1e-12)
  )

cat("\nYear-by-year BHJ shock variance:\n")
cat(sprintf("%-6s  %10s  %10s  %10s  %10s  %8s\n",
            "Year", "Var_s(g)", "SD_s(g)", "Var(B_it)", "SD(B_it)", "SD ratio"))
cat(strrep("-", 60), "\n")
for (i in seq_len(nrow(bhj_annual))) {
  r <- bhj_annual[i, ]
  cat(sprintf("%-6d  %10.6f  %10.6f  %10.6f  %10.6f  %8.2f\n",
              r$year, r$var_s_g, r$sd_s_g, r$var_B, r$sd_B, r$shock_to_bartik))
}

cat(sprintf("\nMean Var_s(g) across years  : %.6f  (SD: %.6f)\n",
            mean(bhj_annual$var_s_g, na.rm = TRUE), sd(bhj_annual$var_s_g, na.rm = TRUE)))
cat(sprintf("Mean SD_s(g) across years   : %.6f\n", mean(bhj_annual$sd_s_g, na.rm = TRUE)))
cat(sprintf("Mean Var(B_it) across years : %.6f  (SD: %.6f)\n",
            mean(bhj_annual$var_B, na.rm = TRUE), sd(bhj_annual$var_B, na.rm = TRUE)))
cat(sprintf("Mean SD_s(g)/SD(B) ratio    : %.2f\n",
            mean(bhj_annual$shock_to_bartik, na.rm = TRUE)))
cat(sprintf("\nInterpretation: SD_s(g)/SD(B) >> 1 means shares compress the cross-country\n"))
cat(sprintf("variation substantially; share heterogeneity is the binding variation source.\n"))

# =============================================================================
# Section 3: Pooled decomposition
# =============================================================================
message("\n=== Section 3: Pooled decomposition ===")

# Pooled: stack all year-commodity observations; pool s_j (time-invariant)
pooled_shocks <- shocks |>
  filter(year %in% SAMPLE_YRS) |>
  left_join(agg_weights |> select(pink_var, s_j_norm), by = "pink_var") |>
  filter(!is.na(g), !is.na(s_j_norm))

# Pooled g_bar (averaging s_j * g over years)
g_bar_pooled <- sum(pooled_shocks$s_j_norm * pooled_shocks$g, na.rm = TRUE) /
                length(SAMPLE_YRS)

var_s_pooled <- pooled_shocks |>
  group_by(pink_var) |>
  summarise(s_j_norm = first(s_j_norm),
            mean_g   = mean(g, na.rm = TRUE), .groups = "drop") |>
  summarise(var_s = sum(s_j_norm * (mean_g - g_bar_pooled)^2))

cat(sprintf("Pooled Var_s(mean_g): %.6f  (SD: %.6f)\n",
            var_s_pooled$var_s, sqrt(var_s_pooled$var_s)))

# Commodity contributions to shock variance (which commodities drive g-variation?)
comm_contrib <- pooled_shocks |>
  group_by(pink_var, grp) |>
  summarise(s_j_norm = first(s_j_norm),
            mean_g   = mean(g, na.rm = TRUE),
            sd_g     = sd(g, na.rm = TRUE),
            var_g    = var(g, na.rm = TRUE), .groups = "drop") |>
  mutate(
    g_bar_p  = sum(s_j_norm * mean_g),
    contrib  = s_j_norm * (mean_g - g_bar_p)^2,   # contribution to pooled Var_s
    share_var = s_j_norm * var_g                   # within-commodity variance (temporal)
  ) |>
  arrange(desc(abs(contrib)))

cat("\nCommodity contributions to pooled Var_s(mean_g):\n")
cat(sprintf("%-20s  %-12s  %8s  %8s  %8s  %10s\n",
            "Commodity", "Group", "s_j", "mean_g", "sd_g(t)", "contrib"))
cat(strrep("-", 70), "\n")
for (i in seq_len(nrow(comm_contrib))) {
  r <- comm_contrib[i, ]
  cat(sprintf("%-20s  %-12s  %8.4f  %8.4f  %8.4f  %10.6f\n",
              r$pink_var, r$grp, r$s_j_norm, r$mean_g, r$sd_g, r$contrib))
}

# =============================================================================
# Section 4: LaTeX table — annual Var_s(g) + Var(B) + HHI summary
# =============================================================================
message("\n=== Section 4: LaTeX table ===")

write_tex <- function(lines, filename) {
  path <- file.path(OUT_TABS, filename)
  writeLines(lines, path)
  message("  Saved: ", filename)
}

# Summary statistics block
tab <- c(
  "% BHJ (2022) shock-level variance decomposition",
  "",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{BHJ (2022) shock-level variance decomposition of the consumption Bartik shock}",
  "\\label{tab:bhj-variance}",
  "\\begin{adjustbox}{max width=\\textwidth}",
  "\\begin{threeparttable}\\small",
  "\\begin{tabular}{l rrrrr}",
  "\\toprule",
  "Year & $\\text{Var}_s(g_t)$ & $\\text{SD}_s(g_t)$ & $\\text{Var}(B_{it})$ & $\\text{SD}(B_{it})$ & $\\text{SD}_s(g)/\\text{SD}(B)$ \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(bhj_annual))) {
  r <- bhj_annual[i, ]
  tab <- c(tab, sprintf(
    "%d & %.5f & %.5f & %.5f & %.5f & %.1f \\\\",
    r$year, r$var_s_g, r$sd_s_g, r$var_B, r$sd_B, r$shock_to_bartik
  ))
}

tab <- c(tab,
  "\\midrule",
  sprintf("Mean & %.5f & %.5f & %.5f & %.5f & %.1f \\\\",
    mean(bhj_annual$var_s_g, na.rm=TRUE),
    mean(bhj_annual$sd_s_g,  na.rm=TRUE),
    mean(bhj_annual$var_B,   na.rm=TRUE),
    mean(bhj_annual$sd_B,    na.rm=TRUE),
    mean(bhj_annual$shock_to_bartik, na.rm=TRUE)),
  "\\midrule",
  sprintf("\\multicolumn{6}{l}{\\textit{Shock concentration: HHI} $= %.3f$, effective $N = %.1f$} \\\\",
    hhi, eff_n),
  "\\bottomrule",
  "\\end{tabular}",
  paste0(
    "\\begin{tablenotes}\\small",
    "\n\\item \\textit{Notes:} Following \\citet{borusyak2022quasi}. ",
    "Aggregate weight $s_j = N^{-1}\\sum_i z_{ij}$ where $z_{ij}$ is country $i$'s ",
    "time-invariant caloric consumption share for commodity $j$. ",
    "$\\text{Var}_s(g_t) = \\sum_j s_j(g_{jt} - \\bar{g}_t)^2$ is the share-weighted variance ",
    "of the $J=", length(commodities), "$ commodity log-price changes in year $t$. ",
    "$\\text{Var}(B_{it})$ is the cross-country variance of the Bartik shock. ",
    "HHI $= \\sum_j s_j^2$ measures share concentration; effective $N = 1/\\text{HHI}$ is the ",
    "equivalent number of independent shocks. ",
    "Note: this HHI is based on average consumption shares $s_j$; the Rotemberg identification-",
    "weight HHI reported elsewhere is higher because it weights by contribution to the LP estimate, ",
    "not by caloric exposure. ",
    "The SD ratio $\\text{SD}_s(g)/\\text{SD}(B)$ shows how much the macro shocks vary relative to ",
    "the cross-country Bartik variation; values $>1$ indicate share heterogeneity compresses variation.",
    "\n\\end{tablenotes}"
  ),
  "\\end{threeparttable}",
  "\\end{adjustbox}",
  "\\end{table}"
)

write_tex(tab, "tab_bhj_variance.tex")
message("\nDone.")
