# Replication Package: Food Price Shocks and State Fragility

Dalheimer, Bae, Foster, Balagtas (2026)  
*Nature Food* — "The Geopolitics of Food" focus issue

---

## Quick start — reproduce all paper figures and tables

All figures and tables can be reproduced directly from the pre-built processed
datasets in `data/processed/`. No data pipeline re-run is required.

**Step 1** — unzip the FAO Food Balance Sheet data:

```bash
unzip data/raw/fao/FBSH_normalized.zip -d data/raw/fao/
```

**Step 2** — open the R project file in RStudio (or set the working
directory to this folder) and run the script in the R folder sequentially


Figures are written to `output/figures/` and tables to `output/tables/`.
Running time is approximately 5–10 minutes.



## Software requirements

R ≥ 4.2.0. Required packages are installed automatically by `R/00_setup.R`:

```
tidyverse, readxl, writexl, fixest, modelsummary, countrycode,
here, zoo, mFilter, httr, janitor, haven, patchwork, scales, ggrepel, WDI
```

---
