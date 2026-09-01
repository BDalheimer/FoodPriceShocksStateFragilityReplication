# =============================================================================
# 00_setup.R
# Package installation, global paths, commodity-to-FAO mapping
# =============================================================================

# ---- Packages ---------------------------------------------------------------
pkgs <- c(
  "tidyverse", "readxl", "writexl",
  "fixest",          # fast FE & IV regressions
  "modelsummary",    # publication tables
  "countrycode",     # country identifier crosswalk
  "here",            # project-relative paths
  "zoo",             # rolling means
  "mFilter",         # HP filter
  "httr",            # HTTP downloads
  "janitor",         # clean_names()
  "haven",           # read .dta / .xls stata-style
  "patchwork",       # combine ggplot panels
  "scales",          # axis formatting
  "ggrepel"          # non-overlapping labels
)

new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs) > 0) install.packages(new_pkgs, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Paths ------------------------------------------------------------------
ROOT     <- here::here()
DATA_RAW <- file.path(ROOT, "data", "raw")
DATA_PRO <- file.path(ROOT, "data", "processed")
OUT_TABS <- file.path(ROOT, "output", "tables")
OUT_FIGS <- file.path(ROOT, "output", "figures")

# ---- Global settings --------------------------------------------------------

LP_HORIZONS   <- 0:4      # local projection horizons (0 = contemporaneous)
SAMPLE_START  <- 2005
SAMPLE_END    <- 2024


# ---- Commodity map: Pink Sheet series → FAO FBS item codes -----------------
#
# fbs_code: FAO Food Balance Sheet item code (kcal/capita/day element = 664)
# The caloric share of each item determines how much weight that commodity's
# global price change receives in country i's shock.
#
# For commodities with multiple Pink Sheet series (e.g. wheat HRW + SRW,
# beef AUS + US), we average the Pink Sheet prices before computing log changes.

COMMODITY_MAP <- tribble(
  ~pink_var,        ~fbs_code,  ~grp,          ~description,
  # Grains — primary staples, highest caloric shares, core food security channel
  "wheat",           2511,      "cereals",     "Wheat and products",
  "rice",            2805,      "cereals",     "Rice (Milled Equivalent)",
  "rice",            2807,      "cereals",     "Rice and products", # for new methodology
  "maize",           2514,      "cereals",     "Maize and products",
  "barley",          2513,      "cereals",     "Barley and products",
  "sorghum",         2520,      "cereals",     "Sorghum and products",
  # Oilseeds & oils — important calories, globally traded cooking oils
  "soybeans",        2555,      "oils",        "Soya Beans",
  "soybean_oil",     2571,      "oils",        "Soybean Oil",
  "palm_oil",        2573,      "oils",        "Palm Oil",
  "groundnut_oil",   2572,      "oils",        "Groundnut Oil",
  "sunflower_oil",   2574,      "oils",        "Sunflowerseed Oil",
  "coconut_oil",     2578,      "oils",        "Coconut Oil",
  "rapeseed_oil",    2576,      "oils",        "Rape and Mustard Oil",
  # Sugar & fruit — moderate calories, tropical production
  "sugar",           2537,      "sugar_fruit", "Sugar Raw Equivalent",
  "bananas",         2617,      "sugar_fruit", "Bananas",
  "oranges",         2612,      "sugar_fruit", "Oranges, Mandarines",
  # Animal products — income-elastic protein; opportunity cost channel for pastoralists
  "beef",            2731,      "animal",      "Bovine Meat",
  "chicken",         2734,      "animal",      "Poultry Meat",
  "lamb",            2732,      "animal",      "Mutton & Goat Meat",
  "shrimp",          2765,      "animal",      "Crustaceans",
  # Beverage crops — near-zero caloric share; shocks affect producer income, not consumer budgets
  "coffee",          2630,      "cash_crops",  "Coffee and products",
  "tea",             2635,      "cash_crops",  "Tea (including Mate)",
  "cocoa",           2633,      "cash_crops",  "Cocoa Beans and products"
)


