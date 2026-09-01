# =============================================================================
# 03_fsi.R
# Process Fund for Peace Fragile States Index from raw xlsx files.
#
# Input:  data/raw/outcomes/fsi*_raw.xlsx  (2006–2023)
# Output: data/processed/fsi.rds
#   iso3c, year, fsi_total, [sub-indicators]
# =============================================================================

source(here::here("R", "00_setup.R"))

out_dir <- file.path(DATA_RAW, "outcomes")
fsi_rds <- file.path(DATA_PRO, "fsi.rds")

parse_fsi_xlsx <- function(path, year) {
  raw <- tryCatch(
    readxl::read_excel(path) |> janitor::clean_names(),
    error = function(e) { message("  Could not read: ", path); NULL }
  )
  if (is.null(raw)) return(NULL)
  
  ctry_col  <- intersect(c("country", "country_name", "nation"), names(raw))[1]
  total_col <- intersect(c("total", "total_score", "fsi_total", "score"), names(raw))[1]
  if (is.na(ctry_col) || is.na(total_col)) return(NULL)
  
  sub_cols <- names(raw)[grepl("^[cepsx][1-3](_|$)", names(raw)) &
                           names(raw) != total_col]
  
  raw |>
    select(country = all_of(ctry_col), fsi_total = all_of(total_col),
           all_of(sub_cols)) |>
    filter(!is.na(country)) |>
    mutate(year      = as.integer(year),
           fsi_total = suppressWarnings(as.numeric(fsi_total)),
           iso3c     = countrycode::countrycode(country, "country.name", "iso3c",
                                                warn = FALSE)) |>
    filter(!is.na(iso3c), !is.na(fsi_total))
}

xlsx_files <- list.files(out_dir,
                         pattern = "^[Ff][Ss][Ii][-_]?\\d{4}.*\\.xlsx$",
                         full.names = TRUE)
message("Found ", length(xlsx_files), " FSI xlsx files")

year_dfs <- map(xlsx_files, function(f) {
  yr <- as.integer(stringr::str_extract(basename(f), "\\d{4}"))
  parse_fsi_xlsx(f, yr)
}) |> setNames(stringr::str_extract(basename(xlsx_files), "\\d{4}"))

fsi <- bind_rows(year_dfs) |>
  arrange(iso3c, year) |>
  distinct(iso3c, year, .keep_all = TRUE)

message("FSI: ", nrow(fsi), " obs, ", n_distinct(fsi$iso3c), " countries, ",
        min(fsi$year), "–", max(fsi$year))

saveRDS(fsi, fsi_rds)
message("Saved: ", fsi_rds)
