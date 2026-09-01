# =============================================================================
# 01_pinksheet.R
# Download and process World Bank Pink Sheet commodity prices
#
# Output: data/processed/pinksheet_annual.rds
#   - year, pink_var, price_usd, dlog_p (annual log price change)
# =============================================================================

source(here::here("R", "00_setup.R"))

# ---- Download ---------------------------------------------------------------
# The Pink Sheet monthly Excel is available from the World Bank commodity
# markets page. The URL below is for the historical monthly data.
# If this URL breaks, go to:
#   https://www.worldbank.org/en/research/commodity-markets
# and download "CMO-Historical-Data-Monthly.xlsx"

ps_url  <- "https://thedocs.worldbank.org/en/doc/5d903e848db1d1b83e0ec8f744e55570-0350012021/related/CMO-Historical-Data-Monthly.xlsx"
ps_file <- file.path(DATA_RAW, "pinksheet", "CMO-Historical-Data-Monthly.xlsx")

if (!file.exists(ps_file)) {
  message("Downloading Pink Sheet...")
  resp <- httr::GET(ps_url, httr::write_disk(ps_file, overwrite = TRUE),
                    httr::timeout(120))
  if (httr::status_code(resp) != 200) {
    stop("Download failed (", httr::status_code(resp), "). ",
         "Please download CMO-Historical-Data-Monthly.xlsx manually from ",
         "https://www.worldbank.org/en/research/commodity-markets ",
         "and place it in data/raw/pinksheet/")
  }
  message("Downloaded: ", ps_file)
} else {
  message("Pink Sheet already downloaded.")
}

# ---- Read -------------------------------------------------------------------
# The Pink Sheet has multiple sheets. "Monthly Prices" contains the time series.
# Rows 1-4 are headers/units; data starts at row 5.

sheets <- readxl::excel_sheets(ps_file)
message("Sheets found: ", paste(sheets, collapse = ", "))

# Find the monthly prices sheet (name varies slightly across versions)
monthly_sheet <- sheets[str_detect(tolower(sheets), "monthly")]
if (length(monthly_sheet) == 0) monthly_sheet <- sheets[1]
message("Using sheet: ", monthly_sheet[1])

raw <- readxl::read_excel(ps_file, sheet = monthly_sheet[1],
                          skip = 4, col_names = TRUE) |>
  janitor::clean_names()

# First column is the period (e.g. "1960M01")
names(raw)[1] <- "period"

# ---- Parse dates ------------------------------------------------------------
ps_long <- raw |>
  filter(!is.na(period), str_detect(period, "^\\d{4}M")) |>
  mutate(
    year  = as.integer(str_sub(period, 1, 4)),
    month = as.integer(str_sub(period, 6, 7))
  ) |>
  select(-period) |>
  pivot_longer(-c(year, month), names_to = "raw_var", values_to = "price") |>
  mutate(price = suppressWarnings(as.numeric(price))) |>
  filter(!is.na(price), price > 0)

# ---- Map Pink Sheet variables to our commodity names -----------------------
# We need to match column names (cleaned by janitor) to commodity groups.
# The cleaned names look like: "wheat_us_hrw", "rice_thai_5", "palm_oil", etc.
# We build a crosswalk by pattern matching.

ps_crosswalk <- tribble(
  ~raw_pattern,             ~pink_var,
  "wheat.*hrw",             "wheat",
  "wheat.*srw",             "wheat",        # will be averaged with HRW
  "rice.*thai",             "rice",
  "maize",                  "maize",
  "barley",                 "barley",
  "soybean.*oil|soy.*oil",  "soybean_oil",
  "soybean[^_]|soybeans",   "soybeans",
  "palm.*oil",              "palm_oil",
  "groundnut.*oil|peanut.*oil", "groundnut_oil",
  "sunflower.*oil",         "sunflower_oil",
  "coconut.*oil",           "coconut_oil",
  "rapeseed.*oil|canola",   "rapeseed_oil",
  "sugar.*world",           "sugar",
  "beef|bovine",            "beef",
  "chicken|poultry",        "chicken",
  "lamb|sheep",             "lamb",
  "shrimp",                 "shrimp",
  "banana.*us|banana.*eu|bananas", "bananas",
  "orange",                 "oranges",
  "coffee.*arabica",        "coffee",
  "coffee.*robusta",        "coffee",       # averaged with Arabica
  "tea.*mombasa|tea",       "tea",
  "cocoa",                  "cocoa",
  "sorghum",                "sorghum"
)

# Match each raw variable to a pink_var
assign_pink_var <- function(var_name) {
  for (i in seq_len(nrow(ps_crosswalk))) {
    if (str_detect(var_name, ps_crosswalk$raw_pattern[i])) {
      return(ps_crosswalk$pink_var[i])
    }
  }
  return(NA_character_)
}

var_map <- ps_long |>
  distinct(raw_var) |>
  mutate(pink_var = map_chr(raw_var, assign_pink_var))

message("\nVariable mapping:")
var_map |> filter(!is.na(pink_var)) |> print(n = 40)
unmapped <- var_map |> filter(is.na(pink_var)) |> pull(raw_var)
message("Unmapped (excluded): ", paste(head(unmapped, 10), collapse = ", "))

# ---- Annual averages --------------------------------------------------------
ps_annual <- ps_long |>
  left_join(var_map, by = "raw_var") |>
  filter(!is.na(pink_var)) |>
  # keep only our mapped commodities
  filter(pink_var %in% COMMODITY_MAP$pink_var) |>
  # average across months and across varieties within the same pink_var
  group_by(year, pink_var) |>
  summarise(price = mean(price, na.rm = TRUE), .groups = "drop")

# Check coverage
ps_annual |>
  group_by(pink_var) |>
  summarise(n = n(), first = min(year), last = max(year)) |>
  arrange(pink_var) |>
  print(n = 30)

# ---- Log price changes ------------------------------------------------------
ps_annual <- ps_annual |>
  arrange(pink_var, year) |>
  group_by(pink_var) |>
  mutate(
    log_p  = log(price),
    dlog_p = log_p - lag(log_p)         # annual log price change
  ) |>
  ungroup()

ps_annual <- ps_annual |>
  filter(year >= SAMPLE_START, year <= SAMPLE_END)

# ---- Save -------------------------------------------------------------------
saveRDS(ps_annual, file.path(DATA_PRO, "pinksheet_annual.rds"))
message("Saved: data/processed/pinksheet_annual.rds")
message("Commodities: ", paste(sort(unique(ps_annual$pink_var)), collapse = ", "))
message("Years: ", min(ps_annual$year, na.rm = TRUE), "-",
        max(ps_annual$year, na.rm = TRUE))

