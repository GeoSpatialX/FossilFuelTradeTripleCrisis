# Shared settings and helpers, sourced by scripts 3-5. Run every script with the
# working directory set to the folder that holds Data/ and Scripts/.
library(dplyr)
library(tidyr)
library(tidygraph)
library(ggplot2)

fontstyle <- "Arial"

# Crisis onset months. Baseline = the 12 months before onset; observation =
# the 12 months starting at onset.
crises <- list(
  us_china_trade_tensions = c(year = 2018, month = 7),
  covid_19_pandemic = c(year = 2020, month = 3),
  russia_ukraine_conflict = c(year = 2022, month = 2)
)
# 2701 coal, 2709 crude, 2710 refined petroleum, 2711 petroleum gases.
hs_codes <- c(2701, 2709, 2710, 2711)

hs_color_palette <- c("2701" = "#B23A48", "2709" = "#3D5A80", "2710" = "#5F8B4C", "2711" = "#C97F2A")

region_palette <- c(
  "East Asia" = "#A45063", "Southeast Asia" = "#BB7D5C", "South Asia" = "#826800",
  "Central Asia" = "#819337", "West Asia & North Africa" = "#237A32", "Europe" = "#00A080",
  "North America" = "#007D81", "Latin America & Caribbean" = "#3196BC",
  "Sub-Saharan Africa" = "#5C66A6", "Oceania" = "#AC78BF", "Rest of the World" = "#A14C87",
  "Unknown" = "grey60"
)

region_lookup <- as.data.frame(terra::vect("Data/worldmap/worldmap.shp")) %>%
  transmute(iso3 = toupper(iso3), region)

# Output path for a file belonging to manuscript figure `fig` (e.g. "fig_2").
figPath <- function(fig, file)
{
  dir.create(file.path("Figures", fig), showWarnings = FALSE, recursive = TRUE)
  file.path("Figures", fig, file)
}

# n consecutive (year, month) pairs starting at (start_year, start_month).
getYM <- function(start_year, start_month, n = 12)
{
  idx <- (start_year * 12 + start_month - 1) + 0:(n - 1)
  data.frame(year = idx %/% 12, month = idx %% 12 + 1)
}

baselineYM <- function(ref) getYM(ref[["year"]] - 1, ref[["month"]])
observationYM <- function(ref) getYM(ref[["year"]], ref[["month"]])

# Row i pairs the same calendar month one year apart; `period_index` counts
# 1-12 from the onset month.
periodMap <- function(ref)
{
  b <- baselineYM(ref)
  data.frame(
    month = b$month, baseline_year = b$year,
    observation_year = observationYM(ref)$year, period_index = seq_len(nrow(b))
  )
}

ymLabel <- function(year, month) format(as.Date(sprintf("%d-%02d-01", year, month)), "%b %Y")

# Year-on-year comparison for every crisis x HS code. `data` holds
# (hs, year, month, <keys>, value); `gridFun(d, period_map)` returns the
# (month x keys) grid to fill, so a missing month stays as a row of NAs.
yoyByCrisis <- function(data, keys, gridFun)
{
  bind_rows(lapply(names(crises), function(crisis_name) {
    period_map <- periodMap(crises[[crisis_name]])
    bind_rows(lapply(hs_codes, function(hs_code) {
      d <- data %>% filter(hs == hs_code)
      v <- d %>% select(all_of(c(keys, "year", "month")), value)
      gridFun(d, period_map) %>%
        left_join(rename(v, baseline_year = year, baseline_value = value),
                  by = c("baseline_year", "month", keys)) %>%
        left_join(rename(v, observation_year = year, observation_value = value),
                  by = c("observation_year", "month", keys)) %>%
        mutate(
          crisis = crisis_name,
          hs = hs_code,
          abs_change = observation_value - baseline_value,
          pct_change = (observation_value - baseline_value) / baseline_value * 100,
          abs_change = ifelse(is.finite(abs_change), abs_change, NA_real_),
          pct_change = ifelse(is.finite(pct_change), pct_change, NA_real_)
        ) %>%
        select(crisis, hs, baseline_year, month, period_index, observation_year, all_of(keys),
               baseline_value, observation_value, abs_change, pct_change)
    }))
  }))
}

# Legend grob of a plot ("guide-box" or "guide-box-<position>", depending on
# the ggplot2 version), written to its own PDF. cairo_pdf() resolves system
# fonts such as `fontstyle`; base pdf() does not.
getLegend <- function(gg)
{
  g <- ggplotGrob(gg)
  candidates <- g$grobs[grepl("^guide-box", g$layout$name)]
  candidates[[which(!sapply(candidates, inherits, "zeroGrob"))[1]]]
}

saveLegend <- function(gg, filename, width = 2, height = 3)
{
  grDevices::cairo_pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(getLegend(gg))
  dev.off()
}

# One month's directed edges (exporter iso3 -> importer iso3, weight) from
# the networks built in "2_build_nw.R".
getEdges <- function(hs_code, year, month)
{
  nw <- readRDS(sprintf("Data/Network/hs%d/weight_y%d_m%d.rds", as.integer(hs_code), year, month))$nw
  iso3 <- as.data.frame(activate(nw, nodes))$iso3
  as.data.frame(activate(nw, edges)) %>%
    transmute(from_iso3 = iso3[from], to_iso3 = iso3[to], weight = weight)
}

safeDiv <- function(num, denom) if (denom == 0) NA_real_ else num / denom

# Similarity of two edge tables (`edges0` earlier, `edgest` later): Jaccard
# and its retention/formation/dissolution rates (which links exist),
# Bray-Curtis on weights, and Bray-Curtis on each edge's share of the
# period total (reallocation only).
nwSimilarity <- function(edges0, edgest)
{
  key <- function(d) paste(d$from_iso3, d$to_iso3, sep = "->")
  e0 <- unique(key(edges0))
  et <- unique(key(edgest))
  total <- function(w) { s <- sum(w, na.rm = TRUE); if (s == 0) NA_real_ else s }

  all_edges <- full_join(
    edges0 %>% transmute(from_iso3, to_iso3, w0 = weight),
    edgest %>% transmute(from_iso3, to_iso3, wt = weight),
    by = c("from_iso3", "to_iso3")
  ) %>%
    mutate(
      w0 = ifelse(is.na(w0), 0, w0), wt = ifelse(is.na(wt), 0, wt),
      p0 = w0 / total(edges0$weight), pt = wt / total(edgest$weight)
    )

  data.frame(
    n_edges_baseline = length(e0),
    n_edges_observation = length(et),
    jaccard = safeDiv(length(intersect(e0, et)), length(union(e0, et))),
    retention_rate = safeDiv(length(intersect(e0, et)), length(e0)),
    formation_rate = safeDiv(length(setdiff(et, e0)), length(e0)),
    dissolution_rate = safeDiv(length(setdiff(e0, et)), length(e0)),
    bray_curtis = safeDiv(
      sum(all_edges$wt + all_edges$w0) - sum(abs(all_edges$wt - all_edges$w0)),
      sum(all_edges$wt + all_edges$w0)
    ),
    bray_curtis_share = 1 - 0.5 * sum(abs(all_edges$pt - all_edges$p0))
  )
}
