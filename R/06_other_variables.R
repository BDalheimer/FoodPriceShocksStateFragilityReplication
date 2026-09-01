# =============================================================================
# 06_other_variables.R
# Pull covariates from World Bank WDI and save for panel construction.
#
# Variables:
#   gdp_pc       — GDP per capita, constant 2015 USD       (NY.GDP.PCAP.KD)
#   rural_pct    — Rural population share, %               (SP.RUR.TOTL.ZS)
#   trade_gdp    — Trade (imports + exports) % of GDP      (NE.TRD.GNFS.ZS)
#   gov_rev      — Gov. revenue excl. grants, % GDP        (GC.REV.XGRT.GD.ZS)
#   trade_tax    — Taxes on international trade, % revenue (GC.TAX.INTT.RV.ZS)
#   gs_tax       — Goods & services tax revenue, % revenue (GC.TAX.GSRV.RV.ZS)
#   gov_exp      — General govt expenditure, % GDP         (GC.XPN.TOTL.GD.ZS)
#   (GC.BAL.CASH.GD.ZS and DT.DOD.DECT.GD.ZS not available via WDI API)
#
# Income group (World Bank classification) is handled separately in
# get_income_groups.R → data/processed/income_groups.rds
#
# Output: data/processed/covariates.rds
#   iso3c × year with: log_gdp_pc, rural_pct, trade_gdp, gov_rev,
#                      trade_tax, gs_tax, gov_exp
# =============================================================================

source(here::here("R", "00_setup.R"))

RAW_CSV <- file.path(DATA_RAW, "wdi", "wdi_covariates.csv")

# ---- Download (run once, then comment out) ----------------------------------
# library(WDI)
# wdi_raw <- WDI::WDI(
#   indicator = c(
#     gdp_pc     = "NY.GDP.PCAP.KD",
#     rural_pct  = "SP.RUR.TOTL.ZS",
#     trade_gdp  = "NE.TRD.GNFS.ZS",
#     gov_rev   = "GC.REV.XGRT.GD.ZS",
#     trade_tax = "GC.TAX.INTT.RV.ZS",
#     gs_tax    = "GC.TAX.GSRV.RV.ZS",
#     gov_exp   = "GC.XPN.TOTL.GD.ZS"
#   ),
#   country = "all",
#   start   = SAMPLE_START - 5,   # extra years for potential lags
#   end     = SAMPLE_END,
#   extra   = FALSE
# ) |>
#   filter(nchar(iso3c) == 3)
# dir.create(file.path(DATA_RAW, "wdi"), showWarnings = FALSE)
# write_csv(wdi_raw, RAW_CSV)
# message("Saved raw WDI download: ", RAW_CSV)

# ---- Load saved data --------------------------------------------------------
message("Reading: ", RAW_CSV)
wdi <- read_csv(RAW_CSV, show_col_types = FALSE)

FISCAL_VARS <- c("gov_rev", "trade_tax", "gs_tax", "gov_exp", "fiscal_bal", "ext_debt")

# ---- Process ----------------------------------------------------------------
covariates <- wdi |>
  filter(nchar(iso3c) == 3) |>
  mutate(
    iso3c      = as.character(iso3c),
    year       = as.integer(year),
    log_gdp_pc = log(gdp_pc)
  ) |>
  select(iso3c, year, log_gdp_pc, rural_pct, trade_gdp,
         any_of(FISCAL_VARS)) |>
  arrange(iso3c, year)

message("Covariates: ", nrow(covariates), " obs | ",
        n_distinct(covariates$iso3c), " countries | ",
        min(covariates$year), "–", max(covariates$year))

cat("\nCoverage (non-NA obs per variable):\n")
covariates |>
  summarise(across(where(is.numeric), \(x) sum(!is.na(x)))) |>
  print()

# ---- Save -------------------------------------------------------------------
saveRDS(covariates, file.path(DATA_PRO, "covariates.rds"))
message("Saved: data/processed/covariates.rds")

