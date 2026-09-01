# =============================================================================
# 04_shocks.R
# Construct weighted food price shock variables
#
#     Shock_{i,t} = Σ_c share_{i,c,t} × Δlog p_{c,t}
#
# Five exposure types: consumption, production, import, export, net-import.
# Robustness: value-weighted production shock
#   s_val_prod = (prod_qty_ma × price_ma) / Σ_c (prod_qty_ma × price_ma)
#   Both quantity and price use a 5-yr backward MA so the share is
#   predetermined relative to the current-year price shock.
# =============================================================================

source(here::here("R", "00_setup.R"))

# ---- Load inputs ------------------------------------------------------------
prices <- readRDS(file.path(DATA_PRO, "pinksheet_annual.rds"))
shares <- readRDS(file.path(DATA_PRO, "fao_shares.rds"))

# ---- Helper: construct weighted shock ---------------------------------------
# shares_df: iso3c, year, pink_var, <share_col>
# prices_df: year, pink_var, <value_col>
# Joined on (pink_var, year) — works for both fixed and time-varying shares.

make_shock <- function(shares_df, prices_df, share_col, value_col,
                       shock_name = "shock") {
  shares_df |>
    select(iso3c, year, pink_var, share = all_of(share_col)) |>
    mutate(share = replace_na(share, 0)) |>
    inner_join(
      prices_df |> select(year, pink_var, value = all_of(value_col)),
      by = c("pink_var", "year")
    ) |>
    filter(!is.na(value)) |>
    group_by(iso3c, year) |>
    summarise(
      !!shock_name  := sum(share * value, na.rm = TRUE),
      n_commodities  = sum(!is.na(value)),
      .groups = "drop"
    )
}

# ---- Main shocks (quantity / kcal shares × Δlog p) -------------------------
shock_cons <- make_shock(shares, prices, "s_cons",    "dlog_p", "shock_cons")
message("Consumption: ", nrow(shock_cons), " obs, ",
        n_distinct(shock_cons$iso3c), " countries")

shock_prod <- make_shock(shares, prices, "s_prod_ma", "dlog_p", "shock_prod")
message("Production:  ", nrow(shock_prod), " obs, ",
        n_distinct(shock_prod$iso3c), " countries")

shock_imp  <- make_shock(shares, prices, "s_imp_ma",  "dlog_p", "shock_imp")
message("Import:      ", nrow(shock_imp),  " obs, ",
        n_distinct(shock_imp$iso3c),  " countries")

shock_exp  <- make_shock(shares, prices, "s_exp_ma",  "dlog_p", "shock_exp")
message("Export:      ", nrow(shock_exp),  " obs, ",
        n_distinct(shock_exp$iso3c),  " countries")

shock_net  <- make_shock(shares, prices, "s_net_ma",  "dlog_p", "shock_net")
message("Net-import:  ", nrow(shock_net),  " obs, ",
        n_distinct(shock_net$iso3c),  " countries")

# ---- Value-weighted production shock (robustness) ---------------------------
price_ma <- prices |>
  arrange(pink_var, year) |>
  group_by(pink_var) |>
  mutate(
    p_l1 = lag(price, 1), p_l2 = lag(price, 2),
    p_l3 = lag(price, 3), p_l4 = lag(price, 4), p_l5 = lag(price, 5),
    price_ma = rowMeans(cbind(p_l1, p_l2, p_l3, p_l4, p_l5), na.rm = TRUE),
    price_ma = if_else(is.nan(price_ma), price, price_ma)
  ) |>
  ungroup() |>
  select(year, pink_var, price_ma)

shares_val <- shares |>
  left_join(price_ma, by = c("pink_var", "year")) |>
  group_by(iso3c, year) |>
  mutate(
    val_prod_raw  = s_prod_ma * price_ma,
    val_cons_raw  = s_cons    * price_ma,
    val_prod_tot  = sum(val_prod_raw, na.rm = TRUE),
    val_cons_tot  = sum(val_cons_raw, na.rm = TRUE),
    s_val_prod_ma = if_else(val_prod_tot > 0, val_prod_raw / val_prod_tot, NA_real_),
    s_val_cons_ma = if_else(val_cons_tot > 0, val_cons_raw / val_cons_tot, NA_real_)
  ) |>
  ungroup()

shock_val_prod <- make_shock(shares_val, prices, "s_val_prod_ma", "dlog_p", "shock_val_prod")
message("Value-prod:  ", nrow(shock_val_prod), " obs, ",
        n_distinct(shock_val_prod$iso3c), " countries")

shock_val_cons <- make_shock(shares_val, prices, "s_val_cons_ma", "dlog_p", "shock_val_cons")
message("Value-cons:  ", nrow(shock_val_cons), " obs, ",
        n_distinct(shock_val_cons$iso3c), " countries")

# ---- Merge all variants and add lags ----------------------------------------
shocks <- shock_cons |>
  select(iso3c, year, shock_cons, n_commodities) |>
  left_join(shock_prod     |> select(iso3c, year, shock_prod),     by = c("iso3c", "year")) |>
  left_join(shock_imp      |> select(iso3c, year, shock_imp),      by = c("iso3c", "year")) |>
  left_join(shock_exp      |> select(iso3c, year, shock_exp),      by = c("iso3c", "year")) |>
  left_join(shock_net      |> select(iso3c, year, shock_net),      by = c("iso3c", "year")) |>
  left_join(shock_val_prod |> select(iso3c, year, shock_val_prod), by = c("iso3c", "year")) |>
  left_join(shock_val_cons |> select(iso3c, year, shock_val_cons), by = c("iso3c", "year")) |>
  arrange(iso3c, year) |>
  group_by(iso3c) |>
  mutate(
    shock_cons_l1     = lag(shock_cons,     1), shock_cons_l2     = lag(shock_cons,     2),
    shock_cons_l3     = lag(shock_cons,     3), shock_cons_l4     = lag(shock_cons,     4),
    shock_prod_l1     = lag(shock_prod,     1), shock_prod_l2     = lag(shock_prod,     2),
    shock_prod_l3     = lag(shock_prod,     3), shock_prod_l4     = lag(shock_prod,     4),
    shock_imp_l1      = lag(shock_imp,      1), shock_imp_l2      = lag(shock_imp,      2),
    shock_imp_l3      = lag(shock_imp,      3), shock_imp_l4      = lag(shock_imp,      4),
    shock_exp_l1      = lag(shock_exp,      1), shock_exp_l2      = lag(shock_exp,      2),
    shock_exp_l3      = lag(shock_exp,      3), shock_exp_l4      = lag(shock_exp,      4),
    shock_net_l1      = lag(shock_net,      1), shock_net_l2      = lag(shock_net,      2),
    shock_net_l3      = lag(shock_net,      3), shock_net_l4      = lag(shock_net,      4),
    shock_val_prod_l1 = lag(shock_val_prod, 1), shock_val_prod_l2 = lag(shock_val_prod, 2),
    shock_val_prod_l3 = lag(shock_val_prod, 3), shock_val_prod_l4 = lag(shock_val_prod, 4),
    shock_val_cons_l1 = lag(shock_val_cons, 1), shock_val_cons_l2 = lag(shock_val_cons, 2),
    shock_val_cons_l3 = lag(shock_val_cons, 3), shock_val_cons_l4 = lag(shock_val_cons, 4)
  ) |>
  ungroup() |>
  filter(year >= SAMPLE_START, year <= SAMPLE_END)

# ---- Diagnostics ------------------------------------------------------------
message("\n=== Shock diagnostics ===")
message("Countries: ", n_distinct(shocks$iso3c))
message("Years: ", min(shocks$year, na.rm = TRUE), "–",
        max(shocks$year, na.rm = TRUE))

cat("\nShock summaries:\n")
shocks |>
  select(shock_cons, shock_prod, shock_imp, shock_exp, shock_net,
         shock_val_prod, shock_val_cons) |>
  summary() |>
  print()

cat("\nCorrelation between shock variants:\n")
shocks |>
  select(shock_cons, shock_prod, shock_imp, shock_exp, shock_net,
         shock_val_prod, shock_val_cons) |>
  cor(use = "pairwise.complete.obs") |>
  round(3) |>
  print()

# ---- Save -------------------------------------------------------------------
saveRDS(shocks, file.path(DATA_PRO, "shocks.rds"))
message("\nSaved: data/processed/shocks.rds")
