# Step 1. Builds the mirrored UN Comtrade bilateral panel (Data/UNComtrade/comtrade_mirror.csv;
# Section 2.1) and validates it against BACI and OEC (Appendix, Fig. A1).
library(dplyr)
library(terra)
library(tidyterra)
library(ggplot2)

countries <- vect('Data/worldmap/worldmap.shp') %>%
  select(label, iso3, m49) %>%
  as.data.frame()

# Process the UN Comtrade data.-----
# Recode legacy m49 codes, then drop domestic flows and non-country aggregates
# (e.g. "World", region groupings) to prevent double counting.
recodeM49 <- function(x) sprintf('%03d', case_when(
  x == 699 ~ 356, x == 251 ~ 250, x == 579 ~ 578, x == 842 ~ 840, x == 757 ~ 756, TRUE ~ x
))

comtrade <- lapply(list.files('Data/UNComtrade/raw', full.names = TRUE), function(f) {
  df <- tryCatch(read.csv(f), error = function(e) NULL)
  if (is.null(df)) return(NULL)
  df %>%
    transmute(
      year = refYear, month = refMonth,
      from_m49 = ifelse(flowCode == 'M', partnerCode, reporterCode),
      to_m49 = ifelse(flowCode == 'M', reporterCode, partnerCode),
      flow = flowCode, hs = cmdCode, weight = netWgt, value = primaryValue,
      source = ifelse(isNetWgtEstimated, "estimated", "reported")
    )
}) %>%
  bind_rows() %>%
  mutate(from_m49 = recodeM49(from_m49), to_m49 = recodeM49(to_m49)) %>%
  filter(from_m49 != to_m49, from_m49 %in% countries$m49, to_m49 %in% countries$m49)

# Median unit value per HS code, year and flow direction (UNSD methodology).
hs_year_median_uv <- comtrade %>%
  filter(weight > 0, source == "reported") %>%
  mutate(unit_value = value / weight) %>%
  group_by(hs, year, flow) %>%
  summarize(median_uv = median(unit_value, na.rm = TRUE), .groups = "drop")

# A record is plausible when its unit value is within 1/100-100x the benchmark.
comtrade <- comtrade %>%
  left_join(hs_year_median_uv, by = c('hs', 'year', 'flow')) %>%
  mutate(
    unit_value = ifelse(weight > 0, value / weight, NA_real_),
    plausible = !is.na(unit_value) & !is.na(median_uv) & median_uv > 0 &
      unit_value >= median_uv / 100 & unit_value <= median_uv * 100
  )

flowSide <- function(f, suffix) comtrade %>%
  filter(flow == f) %>%
  select(year, month, from_m49, to_m49, hs, weight, value, plausible, source, median_uv) %>%
  rename_with(~ paste0(.x, suffix), c(weight, value, plausible, source, median_uv))

# Use the plausible record if there is one (import first), otherwise impute the weight
# from value / median unit value. `source` is 'reported'/'estimated' when an observed
# weight is used, or 'imputed' when it is derived from value.
comtrade_mirror <- merge(
  flowSide('M', '_m'), flowSide('X', '_x'),
  by = c('year', 'month', 'from_m49', 'to_m49', 'hs'), all = TRUE
) %>%
  mutate(
    use_m = plausible_m %in% TRUE,
    use_x = !use_m & plausible_x %in% TRUE,
    weight = case_when(use_m ~ weight_m, use_x ~ weight_x,
                       !is.na(value_m) ~ value_m / median_uv_m, TRUE ~ value_x / median_uv_x),
    value = case_when(use_m ~ value_m, use_x ~ value_x, !is.na(value_m) ~ value_m, TRUE ~ value_x),
    source = case_when(use_m ~ source_m, use_x ~ source_x, TRUE ~ 'imputed')
  ) %>%
  filter(!is.na(weight), !is.na(value), weight > 0, value > 0) %>%
  left_join(countries, by = c('from_m49' = 'm49')) %>%
  rename(from_iso3 = iso3, from = label) %>%
  left_join(countries, by = c('to_m49' = 'm49')) %>%
  rename(to_iso3 = iso3, to = label) %>%
  select(year, month, from, from_iso3, to, to_iso3, hs, weight, value, source)

write.csv(comtrade_mirror, 'Data/UNComtrade/comtrade_mirror.csv', row.names = FALSE)

# =============================================================================
# Validation against BACI and OEC (Fig. A1)
# =============================================================================
# Each source is taken as it comes (no mirroring, imputation or plausibility
# filter) at its native resolution: BACI is annual, so Comtrade's months are
# summed to the year; OEC is monthly. Only flows present on both sides are
# compared. BACI reconciles its two reporters into one value (one row of
# panels); OEC keeps importer- and exporter-reported records apart (two rows).
validation_hs <- c(2701, 2709, 2710, 2711)

# kg -> tonnes.
comtrade_validation <- comtrade_mirror %>%
  filter(hs %in% validation_hs) %>%
  transmute(exporter = from_iso3, importer = to_iso3, year, month, hs, weight_t = weight / 1e3)

# BACI: Data/BACI/baci_e.csv holds the HS27 rows of the annual BACI files
# (2016-2023, CEPII-BACI as distributed by OEC). It is shipped; the block below
# rebuilds it only if it is missing and the original files are present.
if (!file.exists("Data/BACI/baci_e.csv")) {
  baci <- bind_rows(lapply(2016:2023, function(y) {
    read.csv(paste0("Data/BACI/trade_i_baci_a_07_", y, ".csv")) %>%
      select(year, exporter_iso3 = exporter_id, exporter_name, importer_iso3 = importer_id,
             importer_name, hs = hs_code, commodity = product_name, value, quantity) %>%
      mutate(hs = sprintf("%06s", hs)) %>%
      filter(grepl("^27\\d{4}$", hs)) %>%
      group_by(year, exporter_iso3, exporter_name, importer_iso3, importer_name, hs, commodity) %>%
      summarize(value = sum(value), quantity = sum(quantity), .groups = "drop")
  }))
  write.csv(baci, "Data/BACI/baci_e.csv", row.names = FALSE)
}

# HS6 truncated to HS4; `quantity` is already in tonnes.
baci_validation <- read.csv("Data/BACI/baci_e.csv") %>%
  mutate(hs = as.integer(substr(as.character(hs), 1, 4))) %>%
  filter(hs %in% validation_hs) %>%
  transmute(exporter = exporter_iso3, importer = importer_iso3, year, hs, weight_t = quantity)

# OEC: the per-HS downloads in Data/OEC/ carry a `Unit` on every record. Only
# MASS units are kept and converted exactly to kg; volumes, energy, counts and
# "Not Specified" are dropped rather than converted by an assumed density.
# There is no 2701 file.
# `Origin Country ID` is the reporter with a continent prefix ("afzaf"),
# `Country ID` the partner; `HS4 ID` carries a leading section digit
# ("52710"). An "Imports" row is reported by the importer, "Exports" by the
# exporter.-----
oec_mass_unit_kg <- c("Kilogram" = 1, "Liquid Kilogram" = 1, "Gram" = 1e-3, "Metric Tonne Air Dry" = 1e3)

oec_raw <- bind_rows(lapply(
  list.files("Data/OEC", pattern = "^27(09|10|11)_(export|import)_.*[.]csv$", full.names = TRUE),
  function(f) read.csv(f, check.names = FALSE) %>%
    transmute(
      reporter = substr(`Origin Country ID`, 3, 5),
      partner = substr(`Country ID`, 3, 5),
      hs = as.integer(substr(as.character(`HS4 ID`), 2, 5)),
      time = Time, flow = `Trade Flow`, unit = Unit, quantity = Quantity
    )
))

# OEC records by declared unit, and whether each unit is kept.
write.csv(oec_raw %>% count(unit) %>% mutate(pct = round(100 * n / sum(n), 1), kept = unit %in% names(oec_mass_unit_kg)) %>%
            arrange(desc(n)), "Data/validation_oec_units.csv", row.names = FALSE)

oec_mass <- oec_raw %>%
  filter(hs %in% validation_hs, flow %in% c("Imports", "Exports"),
         unit %in% names(oec_mass_unit_kg), quantity > 0) %>%
  mutate(
    exporter = ifelse(flow == "Imports", partner, reporter),
    importer = ifelse(flow == "Imports", reporter, partner),
    year = as.integer(substr(time, 1, 4)),
    month = as.integer(substr(time, 5, 6))
  )

# Two corrections, inferred from Data/validation_oec_unit_check.csv below:
# China's "Gram" records sit at exactly 1000x Comtrade, i.e. they are already
# kilograms and are read as such; the USA's "Kilogram" records scatter over
# five orders of magnitude against Comtrade, so the USA is excluded as a
# reporter. Every other reporter is taken at its declared unit.-----
oec_unit_override <- c(chn = "Gram")
oec_excluded_reporter <- c("usa")

oec_validation <- oec_mass %>%
  filter(!reporter %in% oec_excluded_reporter) %>%
  transmute(
    source = ifelse(flow == "Imports", "OEC (importer-reported)", "OEC (exporter-reported)"),
    exporter, importer, year, month, hs,
    weight_t = quantity * ifelse(
      reporter %in% names(oec_unit_override) & unit == oec_unit_override[reporter],
      1, oec_mass_unit_kg[unit]
    ) / 1e3
  )

# Audit table behind the two corrections: for every (reporter, declared unit),
# the share of matched monthly flows in each power-of-ten ratio to Comtrade,
# on the uncorrected declared-unit conversion.-----
oec_unit_check <- oec_mass %>%
  transmute(reporter, unit, exporter, importer, year, month, hs,
            oec_t = quantity * oec_mass_unit_kg[unit] / 1e3) %>%
  inner_join(comtrade_validation %>% filter(weight_t > 0) %>% rename(comtrade_t = weight_t),
             by = c("exporter", "importer", "hs", "year", "month")) %>%
  mutate(ratio_bin = round(log10(comtrade_t / oec_t))) %>%
  count(reporter, unit, ratio_bin) %>%
  group_by(reporter, unit) %>%
  mutate(n_flows = sum(n), share = round(n / n_flows, 3)) %>%
  ungroup() %>%
  select(reporter, unit, n_flows, ratio_bin, share) %>%
  arrange(desc(n_flows), reporter, unit, ratio_bin)

write.csv(oec_unit_check, "Data/validation_oec_unit_check.csv", row.names = FALSE)

# Sum to the keys given (HS6 -> HS4 for BACI; months -> year where annual).
aggregateValidation <- function(df, keys) df %>%
  group_by(across(any_of(c("source", keys)))) %>%
  summarise(weight_t = sum(weight_t, na.rm = TRUE), .groups = "drop") %>%
  filter(weight_t > 0)

# inner_join keeps only flows both sources record.
matchSource <- function(other, keys) inner_join(
  aggregateValidation(comtrade_validation, keys) %>% rename(weight_comtrade = weight_t),
  aggregateValidation(other, keys) %>% rename(weight_other = weight_t),
  by = keys
)

annual_keys <- c("exporter", "importer", "hs", "year")

validation_data <- bind_rows(
  matchSource(baci_validation, annual_keys) %>% mutate(source = "BACI"),
  matchSource(oec_validation, c(annual_keys, "month"))
) %>%
  mutate(
    hs = factor(hs, levels = validation_hs),
    source = factor(source, levels = c("BACI", "OEC (importer-reported)", "OEC (exporter-reported)"))
  )

write.csv(validation_data, "Data/validation_comtrade_vs_baci_oec.csv", row.names = FALSE)

# Agreement per panel: slope of a through-the-origin fit of Comtrade on the
# reference source (1 = agree on average), R^2, and the number of points
# (pair-years for BACI, pair-months for OEC).-----
validation_stats <- validation_data %>%
  group_by(source, hs) %>%
  summarise(
    n = n(),
    slope = sum(weight_other * weight_comtrade) / sum(weight_other^2),
    r2 = cor(weight_comtrade, weight_other)^2,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("slope = %.2f\nR² = %.2f\nn = %d", slope, r2, n),
    x_pos = max(validation_data$weight_other),
    y_pos = min(validation_data$weight_comtrade)
  )

write.csv(validation_stats %>% select(source, hs, n, slope, r2), "Data/validation_stats.csv", row.names = FALSE)

# Log-log axes; the dashed line is y = x (perfect agreement), not a fit.-----
p_validation <- ggplot(validation_data, aes(x = weight_other, y = weight_comtrade)) +
  geom_point(color = "#3D5A80", size = 0.25, alpha = 0.08, stroke = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#B23A48", linewidth = 0.4) +
  geom_text(
    data = validation_stats, aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE, hjust = 1, vjust = 0,
    family = "Arial", size = 2.3, color = "grey25", lineheight = 1.15
  ) +
  facet_grid(source ~ hs) +
  scale_x_log10(labels = scales::label_log()) +
  scale_y_log10(labels = scales::label_log()) +
  labs(x = "Reference source, unprocessed (tonnes)", y = "UN Comtrade, processed (tonnes)") +
  theme_minimal(base_size = 9, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    plot.title = element_blank(),
    strip.text = element_text(size = 8),
    strip.text.y = element_text(size = 8, angle = 0),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 11.25),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = unit(2, "pt")
  )

dir.create("Figures/fig_a1", showWarnings = FALSE, recursive = TRUE)
ggsave("Figures/fig_a1/fig_a1.png", p_validation, width = 190, height = 145, unit = "mm", dpi = 300)
