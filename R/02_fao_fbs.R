# =============================================================================
# 02_fao_fbs.R
# Build FAO Food Balance Sheet consumption and production/trade shares
#
# Input:
#   data/raw/fao/FoodBalanceSheets_E_All_Data_(Normalized).csv
#   (New FBS methodology, 2010–2023; column "Element Code" 664 = kcal)
#
# Outputs:
#   fao_shares.rds    
#
# NOTE: fao_shares_cons_9909.rds requires the old FBSH file
#   FoodBalanceSheetsHistoric_E_All_Data_(Normalized).csv
#   Download from: https://www.fao.org/faostat/en/#data/FBSH
#   Place it in data/raw/fao/ before running.
# =============================================================================

source(here::here("R", "00_setup.R"))

FBS_CSV  <- file.path(DATA_RAW, "fao",
                      "FoodBalanceSheets_E_All_Data_(Normalized).csv")
FBSH_CSV <- file.path(DATA_RAW, "fao/Historic",
                      "FoodBalanceSheetsHistoric_E_All_Data_(Normalized).csv")

fbs_codes <- COMMODITY_MAP$fbs_code
MA_WINDOW <- 5L

# ---- Read FBS (2010–2023) ----------------------------------------------------
message("Reading FBS CSV (2010-2023)...")

fbs <- read_csv(FBS_CSV, show_col_types = FALSE,
                locale = locale(encoding = "latin1")) |>
  janitor::clean_names() |>
  filter(element_code %in% c(664L, 5511L, 5611L, 5911L, 5142L),
         item_code %in% fbs_codes) |>
  mutate(
    iso3c = countrycode::countrycode(area_code, "fao", "iso3c", warn = FALSE),
    val   = suppressWarnings(as.numeric(value))
  ) |>
  filter(!is.na(iso3c), !is.na(val), val >= 0) |>
  select(iso3c, item_code, element_code, year, val) |>
  left_join(COMMODITY_MAP |> select(pink_var, fbs_code, grp),
            by = c("item_code" = "fbs_code"))

message("  Rows: ", nrow(fbs),
        " | Years: ", min(fbs$year), "–", max(fbs$year),
        " | Countries: ", n_distinct(fbs$iso3c))

# ---- Consumption shares: 2010–2023 fixed average ----------------------------
message("\n=== Consumption shares (2010-2023) ===")

shares_cons <- fbs |>
  filter(element_code == 664L, !is.na(pink_var)) |>
  group_by(iso3c, pink_var, grp) |>
  summarise(kcal = mean(val, na.rm = TRUE), .groups = "drop") |>
  group_by(iso3c) |>
  mutate(
    kcal_total = sum(kcal, na.rm = TRUE),
    s_cons     = if_else(kcal_total > 0, kcal / kcal_total, NA_real_)
  ) |>
  ungroup() |>
  filter(!is.na(s_cons)) |>
  select(iso3c, pink_var, grp, s_cons)

message("  Countries: ", n_distinct(shares_cons$iso3c),
        " | Commodities: ", n_distinct(shares_cons$pink_var))
#saveRDS(shares_cons, file.path(DATA_PRO, "fao_shares_cons.rds"))
#message("  Saved: fao_shares_cons.rds")

# ---- Time-varying production/trade shares: 5-yr backward MA -----------------
message("\n=== Time-varying production/trade shares (", MA_WINDOW, "-yr MA) ===")

trade_annual <- fbs |>
  filter(element_code %in% c(5511L, 5611L, 5911L, 5142L), !is.na(pink_var)) |>
  pivot_wider(names_from = element_code, values_from = val, names_prefix = "el_") |>
  rename(prod        = el_5511,
         imports     = el_5611,
         exports     = el_5911,
         food_supply = el_5142)

totals_annual <- trade_annual |>
  group_by(iso3c, year) |>
  summarise(prod_total = sum(prod,        na.rm = TRUE),
            fs_total   = sum(food_supply, na.rm = TRUE),
            .groups = "drop")

shares_annual <- trade_annual |>
  left_join(totals_annual, by = c("iso3c", "year")) |>
  mutate(
    s_prod = if_else(prod_total > 0, prod / prod_total,               NA_real_),
    s_imp  = if_else(fs_total   > 0, imports / fs_total,             NA_real_),
    s_exp  = if_else(fs_total   > 0, exports / fs_total,             NA_real_),
    s_net  = if_else(fs_total   > 0, (imports - exports) / fs_total, NA_real_)
  ) |>
  select(iso3c, year, pink_var, grp, s_prod, s_imp, s_exp, s_net)

# 5-yr backward MA using explicit lags (zoo not needed)
shares_ma <- shares_annual |>
  arrange(iso3c, pink_var, year) |>
  group_by(iso3c, pink_var) |>
  mutate(across(c(s_prod, s_imp, s_exp, s_net),
                list(l1 = \(x) lag(x, 1), l2 = \(x) lag(x, 2),
                     l3 = \(x) lag(x, 3), l4 = \(x) lag(x, 4),
                     l5 = \(x) lag(x, 5)),
                .names = "{.col}_{.fn}")) |>
  mutate(
    s_prod_ma = rowMeans(cbind(s_prod_l1,s_prod_l2,s_prod_l3,s_prod_l4,s_prod_l5), na.rm=TRUE),
    s_imp_ma  = rowMeans(cbind(s_imp_l1, s_imp_l2, s_imp_l3, s_imp_l4, s_imp_l5),  na.rm=TRUE),
    s_exp_ma  = rowMeans(cbind(s_exp_l1, s_exp_l2, s_exp_l3, s_exp_l4, s_exp_l5),  na.rm=TRUE),
    s_net_ma  = rowMeans(cbind(s_net_l1, s_net_l2, s_net_l3, s_net_l4, s_net_l5),  na.rm=TRUE),
    # fall back to contemporaneous share when no lags exist (first years)
    s_prod_ma = if_else(is.nan(s_prod_ma), s_prod, s_prod_ma),
    s_imp_ma  = if_else(is.nan(s_imp_ma),  s_imp,  s_imp_ma),
    s_exp_ma  = if_else(is.nan(s_exp_ma),  s_exp,  s_exp_ma),
    s_net_ma  = if_else(is.nan(s_net_ma),  s_net,  s_net_ma)
  ) |>
  ungroup() |>
  select(iso3c, year, pink_var, grp, s_prod_ma, s_imp_ma, s_exp_ma, s_net_ma)

# Back-fill pre-2010 estimation-sample years using 2010 shares
shares_2010 <- shares_ma |> filter(year == 2010) |> select(-year)
shares_tv <- bind_rows(
  shares_ma,
  shares_2010 |>
    crossing(year = seq(SAMPLE_START, 2009L)) |>
    select(iso3c, year, pink_var, grp, s_prod_ma, s_imp_ma, s_exp_ma, s_net_ma)
) |>
  arrange(iso3c, pink_var, year)

message("  Countries: ", n_distinct(shares_tv$iso3c),
        " | Year range: ", min(shares_tv$year), "–", max(shares_tv$year))

shares = merge(shares_tv, shares_cons, by = c("iso3c", "pink_var", "grp"), all.x = TRUE) 

saveRDS(shares, file.path(DATA_PRO, "fao_shares.rds"))
message("  Saved: fao_shares.rds")

# ---- Consumption shares: 1999–2009 (requires old FBSH) ----------------------
message("\n=== Consumption shares (1999-2009, FBSH) ===")

if (!file.exists(FBSH_CSV)) {
  message("  SKIPPED — old FBSH file not found:")
  message("  ", FBSH_CSV)
  message("  Download from https://www.fao.org/faostat/en/#data/FBSH")
  message("  (select 'Food Balances (old methodology, 1961-2013)')")
} else {
  fbsh <- read_csv(FBSH_CSV, show_col_types = FALSE,
                   locale = locale(encoding = "latin1")) |>
    janitor::clean_names() |>
    filter(element_code == 664L, item_code %in% fbs_codes) |>
    mutate(
      iso3c = countrycode::countrycode(area_code, "fao", "iso3c", warn = FALSE),
      val   = suppressWarnings(as.numeric(value))
    ) |>
    filter(!is.na(iso3c), !is.na(val), val >= 0) |>
    select(iso3c, item_code, year, val) |>
    left_join(COMMODITY_MAP |> select(pink_var, fbs_code, grp),
              by = c("item_code" = "fbs_code"))
  
  shares_cons_9909 <- fbsh |>
    filter(year %in% 1999:2009, !is.na(pink_var)) |>
    group_by(iso3c, item_code, pink_var, grp) |>
    summarise(kcal = mean(val, na.rm = TRUE), n_yrs = n(), .groups = "drop") |>
    filter(n_yrs >= 3) |>
    group_by(iso3c) |>
    mutate(
      kcal_total = sum(kcal, na.rm = TRUE),
      s_cons     = if_else(kcal_total > 0, kcal / kcal_total, NA_real_)
    ) |>
    ungroup() |>
    filter(!is.na(s_cons)) |>
    select(iso3c, pink_var, grp, s_cons)
  
  message("  Countries: ", n_distinct(shares_cons_9909$iso3c),
          " | Commodities: ", n_distinct(shares_cons_9909$pink_var))
  saveRDS(shares_cons_9909, file.path(DATA_PRO, "fao_shares_cons_9909.rds"))
  message("  Saved: fao_shares_cons_9909.rds")
}

message("\nDone.")
