# =============================================================================
# 07_panel.R
# Merge all processed datasets into the analysis panel
#
# Inputs:
#   fsi.rds            — FSI outcome (2006–2023)
#   shocks.rds         — Bartik shocks + lags (SAMPLE_START–SAMPLE_END)
#   shocks_group.rds   — commodity-group shocks + lags
#   income_groups.rds  — World Bank income classification (static)
#   covariates.rds     — GDP pc, rural share, trade/GDP, gov revenue
#
# Output:
#   panel.rds — iso3c × year, SAMPLE_START–SAMPLE_END
#
# LP outcome: dfsi_h{h} = FSI_{t+h} - FSI_{t-1}  for h = 0, 1, 2, 3, 4
# =============================================================================

source(here::here("R", "00_setup.R"))

# ---- Load inputs ------------------------------------------------------------
fsi          <- readRDS(file.path(DATA_PRO, "fsi.rds"))
shocks       <- readRDS(file.path(DATA_PRO, "shocks.rds"))
shocks_group <- readRDS(file.path(DATA_PRO, "shocks_group.rds"))
income_grps  <- readRDS(file.path(DATA_PRO, "income_groups.rds"))
covariates   <- readRDS(file.path(DATA_PRO, "covariates.rds"))

message("FSI:          ", nrow(fsi), " obs | ",
        n_distinct(fsi$iso3c), " countries | ",
        min(fsi$year), "–", max(fsi$year))
message("Shocks:       ", nrow(shocks), " obs | ",
        n_distinct(shocks$iso3c), " countries | ",
        min(shocks$year), "–", max(shocks$year))

# ---- Build LP outcome variables ---------------------------------------------
# dfsi_h{h} = FSI_{t+h} - FSI_{t-1}
#
# Use complete() to fill year gaps within each country before applying
# lead/lag — ensures lead(x, h) always gives the value at year t+h, not
# just the h-th next non-missing row.

fsi_lp <- fsi |>
  select(iso3c, year, fsi_total) |>
  complete(iso3c, year = full_seq(year, 1)) |>
  arrange(iso3c, year) |>
  group_by(iso3c) |>
  mutate(
    fsi_lag1 = lag(fsi_total,  1),
    dfsi_h0  = fsi_total           - fsi_lag1,
    dfsi_h1  = lead(fsi_total, 1)  - fsi_lag1,
    dfsi_h2  = lead(fsi_total, 2)  - fsi_lag1,
    dfsi_h3  = lead(fsi_total, 3)  - fsi_lag1,
    dfsi_h4  = lead(fsi_total, 4)  - fsi_lag1
  ) |>
  ungroup() |>
  filter(year >= SAMPLE_START, year <= SAMPLE_END) |>
  select(iso3c, year, fsi_total, fsi_lag1, dfsi_h0, dfsi_h1, dfsi_h2, dfsi_h3, dfsi_h4)

message("FSI LP vars:  ", nrow(fsi_lp), " obs | ",
        sum(!is.na(fsi_lp$dfsi_h0)), " non-NA at h=0 | ",
        sum(!is.na(fsi_lp$dfsi_h4)), " non-NA at h=4")

# ---- Merge ------------------------------------------------------------------
# Spine: shocks (all countries with FAO FBS + Pink Sheet coverage)
# Everything else is left-joined — missing rows become NA and are dropped
# naturally in regressions.

panel <- shocks |>
  left_join(fsi_lp,      by = c("iso3c", "year")) |>
  left_join(shocks_group, by = c("iso3c", "year")) |>
  left_join(income_grps,  by = "iso3c") |>
  left_join(covariates,   by = c("iso3c", "year")) |>
  arrange(iso3c, year)

# ---- Diagnostics ------------------------------------------------------------
message("\n=== Panel diagnostics ===")
message("Rows:      ", nrow(panel))
message("Countries: ", n_distinct(panel$iso3c))
message("Years:     ", min(panel$year), "–", max(panel$year))

cat("\nOutcome coverage:\n")
panel |>
  summarise(across(c(dfsi_h0, dfsi_h1, dfsi_h2, dfsi_h3, dfsi_h4),
                   \(x) sum(!is.na(x)))) |>
  print()

cat("\nShock coverage:\n")
panel |>
  summarise(across(c(shock_cons, shock_prod, shock_imp, shock_exp, shock_net,
                     shock_val_prod, shock_val_cons),
                   \(x) sum(!is.na(x)))) |>
  print()

cat("\nIncome group distribution:\n")
print(table(panel$income_group, useNA = "ifany"))

cat("\nCovariate coverage (non-NA obs):\n")
panel |>
  summarise(across(c(log_gdp_pc, rural_pct, trade_gdp, gov_rev),
                   \(x) sum(!is.na(x)))) |>
  print()

# ---- Save -------------------------------------------------------------------
saveRDS(panel, file.path(DATA_PRO, "panel.rds"))
message("\nSaved: data/processed/panel.rds")
message("Columns: ", ncol(panel))
