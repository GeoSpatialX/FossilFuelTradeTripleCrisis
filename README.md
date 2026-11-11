# FossilFuelTradeTripleCrisis

R code for *Dynamic Evolution of Global Fossil Fuel Trade Networks Under Triple Crises: A Complex Network Analysis*.

The study builds monthly, directed, weighted trade networks for coal (HS 2701), crude petroleum (HS 2709),
refined petroleum (HS 2710) and petroleum gases (HS 2711) from 2016 to 2024, and compares how they responded to
the US–China trade tensions (onset July 2018), the COVID-19 pandemic (March 2020) and the Russia–Ukraine conflict
(February 2022). Running the scripts regenerates the figure panels (Figs. 1–5, A1–A8), Table A1 and the
numbers reported in the paper.

## Data

The raw data are **not** included in this repository. Please download them from the original sources below and
place them in the folders shown.

### UN Comtrade (main data)

Monthly bilateral trade records from the [UN Comtrade Database](https://comtradeplus.un.org/):

| Setting | Value |
|---|---|
| Type / frequency / classification | Goods (C), monthly (M), HS |
| Commodities | HS 2701–2716 (4-digit codes under Chapter 27) |
| Period | January 2016 – December 2024 |
| Reporters / partners | All reporters × all partners |
| Flows | Imports (M) and exports (X) |
| Variables used | `refYear`, `refMonth`, `reporterCode`, `partnerCode`, `flowCode`, `cmdCode`, `netWgt`, `isNetWgtEstimated`, `primaryValue` |

Save the files as CSV in `Data/UNComtrade/raw/`. `Scripts/1_data_process.R` reads every CSV in that folder, so
the way the download is split across files does not matter.

`Data/UNComtrade/download_hs27.py` is the Python script we used for the download. It calls the UN Comtrade
API and saves one file per reporter and year (`hs27m_<reporter>_<year>.csv`):

1. Get a free API subscription key from the [UN Comtrade API portal](https://comtradedeveloper.un.org/).
2. Save the key as a single line in `Data/UNComtrade/.comtrade_key`. This file is excluded by `.gitignore`;
   never commit it.
3. Run `python Data/UNComtrade/download_hs27.py`. It needs `pandas` and `requests`. Re-running it skips files
   that are already downloaded.

Records are revised over time, so a new download may differ slightly from ours (downloaded in August 2025).

### Validation data (Appendix, Fig. A1 only)

The second part of `1_data_process.R` compares the processed Comtrade panel with two other databases. They are not
needed for the main analysis (scripts 2–5).

- **CEPII-BACI**: annual HS 6-digit bilateral trade, 2016–2023, from
  [CEPII](https://www.cepii.fr/DATA_DOWNLOAD/baci/doc/baci_webpage.html) (we used the yearly files distributed by
  the [OEC](https://oec.world/en/resources/bulk-download/international), `trade_i_baci_a_07_<year>.csv`).
  Place them in `Data/BACI/`. The script extracts the HS 27 rows to `Data/BACI/baci_e.csv`.
- **Observatory of Economic Complexity (OEC)**: monthly bilateral records for HS 2709, 2710 and 2711 from
  [oec.world](https://oec.world/en/), saved in `Data/OEC/` as `<hs4>_<export|import>_<start>_<end>.csv`
  (e.g. `2709_export_201601201812.csv`), with the columns `Origin Country ID`, `Country ID`, `HS4 ID`, `Time`,
  `Trade Flow`, `Unit` and `Quantity`.

### Country boundaries

The scripts also read a country boundary shapefile, `Data/worldmap/worldmap.shp`. Each feature is a country
with the fields `label` (name), `iso3`, `m49` (UN M49 code) and `region` (the world regions used in Figs. 4–5).

## Workflow

Set the working directory to the repository root (the folder that holds `Data/` and `Scripts/`) and run the
scripts in order:

| Script | What it does | Main outputs |
|---|---|---|
| `1_data_process.R` | Reconciles import and export records into one bilateral panel using a unit-value plausibility screen; validates it against BACI and OEC | `Data/UNComtrade/comtrade_mirror.csv`, Fig. A1 |
| `2_build_nw.R` | Builds the monthly trade networks for HS 2701, 2709, 2710 and 2711 | `Data/Network/`, Table A1, Fig. 1 |
| `3_network.R` | Network properties, year-on-year change, link similarity (Jaccard, Bray–Curtis) and Infomap communities | Figs. A2–A7 |
| `4_crisis_impact.R` | PPML estimates of departures of bilateral trade volumes from pre-crisis trends, with dynamic and placebo specifications | Fig. 2 |
| `5_nodes.R` | Country-level resilience (shock depth, cumulative loss, end state), correlations with pre-crisis node properties, partner-composition case studies | Figs. 3–5, A8 |

`Scripts/functions.R` holds shared settings and helpers (crisis dates, HS codes, palettes, year-on-year and
similarity functions); scripts 3–5 load it automatically.

```r
setwd("path/to/FossilFuelTradeTripleCrisis")
for (s in c("1_data_process", "2_build_nw", "3_network", "4_crisis_impact", "5_nodes"))
  source(file.path("Scripts", paste0(s, ".R")), encoding = "UTF-8")
```

The full run takes about 6 minutes on a 12-core laptop.

## Outputs

- **Figures:** `Figures/fig_<n>/`, one folder per manuscript figure. File names start with the figure and panel,
  e.g. `fig_2j_...png` is Fig. 2j. Each folder also holds that figure's legend.
- **Tables:** `Data/*.csv`, e.g. `crisis_impact_estimates.csv` (Fig. 2), `resilience_correlation.csv` (Fig. 3),
  `partner_composition.csv` (Figs. 4–5), `validation_stats.csv` (Fig. A1) and `table_a1.csv`.

## Requirements

R 4.6.1 with dplyr, tidyr, igraph, tidygraph, fixest, clue, terra, tidyterra, ggplot2, ggforce, ggh4x, ggridges,
migest, viridis and scales. The Arial font must be installed.

The results are deterministic: Infomap (script 3) is seeded, and correlation inputs (script 5) are rounded to 10
significant digits so that tied values rank the same on every platform.

## License

The code is released under the [MIT License](LICENSE). The data remain subject to the terms of their providers.
