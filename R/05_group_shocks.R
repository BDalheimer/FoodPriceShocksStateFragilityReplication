# =============================================================================
# 05_group_shocks.R
# Commodity-group decomposition of food price shocks
#
#     shock_{i,g,t} = Σ_{c ∈ g} share_{i,c,t} × Δlog p_{c,t}
#
# Shares are NOT renormalized within groups, so group shocks are additive:
#     Σ_g shock_{i,g,t} = shock_{i,t}  (total shock from 04_shocks.R)
#
# Exposure types: cons, prod, imp, exp, net, val_prod, val_cons
# Groups: cereals, oils, animal, sugar_fruit, cash_crops
#
# Output: data/processed/shocks_group.rds
#   iso3c × year wide table with shock_{exposure}_{group} and four lags each
# =============================================================================

source(here::here("R", "00_setup.R"))

# ---- Load inputs ------------------------------------------------------------
prices <- readRDS(file.path(DATA_PRO, "pinksheet_annual.rds"))
shares <- readRDS(file.path(DATA_PRO, "fao_shares.rds"))

# ---- Price MA (5-yr backward, same as 04_shocks.R) --------------------------
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

# ---- Value shares (same construction as 04_shocks.R) -----------------------
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

# ---- Group shocks -----------------------------------------------------------
# Replace NA shares with 0: no production/trade in a commodity = zero exposure,
# not missing. Consistent with 04_shocks.R.
group_shocks <- shares_val |>
  inner_join(prices |> select(year, pink_var, dlog_p), by = c("year", "pink_var")) |>
  filter(!is.na(dlog_p), !is.na(grp)) |>
  mutate(across(c(s_cons, s_prod_ma, s_imp_ma, s_exp_ma, s_net_ma,
                  s_val_prod_ma, s_val_cons_ma),
                \(x) replace_na(x, 0))) |>
  group_by(iso3c, year, grp) |>
  summarise(
    shock_cons     = sum(s_cons         * dlog_p),
    shock_prod     = sum(s_prod_ma      * dlog_p),
    shock_imp      = sum(s_imp_ma       * dlog_p),
    shock_exp      = sum(s_exp_ma       * dlog_p),
    shock_net      = sum(s_net_ma       * dlog_p),
    shock_val_prod = sum(s_val_prod_ma  * dlog_p),
    shock_val_cons = sum(s_val_cons_ma  * dlog_p),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from  = grp,
    values_from = c(shock_cons, shock_prod, shock_imp, shock_exp, shock_net,
                    shock_val_prod, shock_val_cons),
    names_glue  = "{.value}_{grp}"
  )

message("Group shock columns: ",
        paste(grep("^shock_", names(group_shocks), value = TRUE), collapse = ", "))

# ---- Additive check ---------------------------------------------------------
# Sum of group shocks per exposure should equal the aggregate shock from 04.
groups     <- unique(COMMODITY_MAP$grp)
exp_types  <- c("cons", "prod", "imp", "exp", "net", "val_prod", "val_cons")

for (exp in exp_types) {
  grp_cols <- paste0("shock_", exp, "_", groups)
  grp_cols <- grp_cols[grp_cols %in% names(group_shocks)]
  group_shocks[[paste0("shock_", exp, "_sum")]] <-
    rowSums(group_shocks[grp_cols], na.rm = TRUE)
}

shocks_agg <- readRDS(file.path(DATA_PRO, "shocks.rds")) |>
  select(iso3c, year, shock_cons, shock_prod, shock_imp, shock_exp, shock_net,
         shock_val_prod, shock_val_cons)

check <- group_shocks |>
  inner_join(shocks_agg, by = c("iso3c", "year"))

for (exp in exp_types) {
  r <- cor(check[[paste0("shock_", exp, "_sum")]],
           check[[paste0("shock_", exp)]],
           use = "complete.obs")
  message("Additive check (", exp, "): r = ", round(r, 4), " (should be ≈ 1)")
}

# Remove helper sum columns
group_shocks <- group_shocks |>
  select(-ends_with("_sum"))

# ---- Add lags ---------------------------------------------------------------
shock_cols <- grep("^shock_", names(group_shocks), value = TRUE)

group_shocks <- group_shocks |>
  filter(year >= SAMPLE_START, year <= SAMPLE_END) |>
  arrange(iso3c, year) |>
  group_by(iso3c) |>
  mutate(across(all_of(shock_cols),
                list(l1 = \(x) lag(x, 1), l2 = \(x) lag(x, 2),
                     l3 = \(x) lag(x, 3), l4 = \(x) lag(x, 4)),
                .names = "{.col}_{.fn}")) |>
  ungroup()

# ---- Save -------------------------------------------------------------------
saveRDS(group_shocks, file.path(DATA_PRO, "shocks_group.rds"))
message("\nSaved: data/processed/shocks_group.rds")
message("Rows: ", nrow(group_shocks), " | Countries: ", n_distinct(group_shocks$iso3c))
message("Shock columns: ", length(shock_cols),
        " base + ", length(shock_cols) * 4, " lags")
