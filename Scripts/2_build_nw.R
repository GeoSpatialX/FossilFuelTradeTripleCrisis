# Step 2. Builds the monthly directed, weighted trade networks for coal (2701),
# crude (2709), refined petroleum (2710) and petroleum gases (2711)
# (Data/Network/hs<code>/; Section 2.2), Table A1, and Fig. 1.
library(dplyr)
library(tidyr)
library(igraph)
library(tidygraph)
library(ggplot2)
library(terra)
library(tidyterra)
library(migest)
library(parallel)

hs_codes <- c(2701, 2709, 2710, 2711)

# Countries as nodes, located at their centroids.-----
worldmap <- vect('Data/worldmap/worldmap.shp') %>%
  select(label, iso3, region)

nodes <- worldmap %>%
  as.data.frame() %>%
  select(label, iso3) %>%
  mutate(iso3 = toupper(iso3)) %>%
  bind_cols(centroids(worldmap) %>% crds() %>% as.data.frame() %>% setNames(c('lon', 'lat')))

# Table A1: ISO3 codes, country names and regions.
write.csv(worldmap %>% as.data.frame() %>% rename(name = label) %>% mutate(iso3 = toupper(iso3)),
          "Data/table_a1.csv", row.names = FALSE)

# Trade data compiled in "1_data_process.R".-----
comtrade <- read.csv('Data/UNComtrade/comtrade_mirror.csv') %>%
  mutate(from_iso3 = toupper(from_iso3), to_iso3 = toupper(to_iso3))

# One month's network for one HS code, with the node metrics of Section 2.2
# (density, strength and ENP are computed from the edges in "3_network.R").-----
buildNw <- function(y, m, hs_code)
{
  edges <- comtrade %>%
    filter(year == y, month == m, hs == hs_code, from_iso3 %in% nodes$iso3, to_iso3 %in% nodes$iso3) %>%
    group_by(from = from_iso3, to = to_iso3) %>%
    summarise(weight = sum(weight, na.rm = TRUE), .groups = "drop") %>%
    na.omit()
  nw <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE, node_key = "iso3")
  w <- as.data.frame(activate(nw, edges))$weight
  for (d in c("in", "out")) {
    vertex_attr(nw, paste0("degree_", d)) <- igraph::degree(nw, mode = d)
    vertex_attr(nw, paste0("strength_", d)) <- igraph::strength(nw, mode = d, weights = w)
  }
  V(nw)$betweenness <- igraph::betweenness(nw, directed = TRUE, normalized = TRUE, weights = NA)
  V(nw)$transitivity <- igraph::transitivity(nw, type = "local")
  list(nw = nw, nw_metrics = c(global_efficiency = igraph::global_efficiency(nw, directed = TRUE, weights = NA)))
}

build_grid <- expand.grid(hs = hs_codes, y = 2016:2024, m = 1:12)
for (hs_code in hs_codes) dir.create(paste0("Data/Network/hs", hs_code), showWarnings = FALSE, recursive = TRUE)
cl <- makeCluster(max(1, min(10, detectCores() - 1)))
clusterExport(cl, c("build_grid", "nodes", "comtrade", "buildNw"))
invisible(clusterEvalQ(cl, { library(dplyr); library(igraph); library(tidygraph) }))
invisible(parLapply(cl, seq_len(nrow(build_grid)), function(i) {
  g <- build_grid[i, ]
  saveRDS(buildNw(g$y, g$m, g$hs), sprintf("Data/Network/hs%d/weight_y%d_m%d.rds", g$hs, g$y, g$m))
}))
stopCluster(cl)

# Fig. 1a: monthly trade weight by HS code (all of HS 2701-2716).-----
trade_weight_monthly <- comtrade %>%
  group_by(year, month, hs) %>%
  summarise(weight = sum(weight, na.rm = TRUE), .groups = "drop") %>%
  complete(year = 2016:2024, month = 1:12, hs, fill = list(weight = 0)) %>%
  mutate(period = paste0(year, "-", sprintf('%02d', month)), hs = as.factor(hs), weight = weight / 1e9)
write.csv(trade_weight_monthly, "Data/trade_weight_monthly.csv", row.names = FALSE)  # million tonnes

fig_1a <- ggplot(trade_weight_monthly) +
  geom_bar(aes(x = period, y = weight, fill = hs), stat = "identity", position = 'stack') +
  scale_fill_viridis_d(option = 'magma') +
  scale_x_discrete(breaks = unique(trade_weight_monthly$period)[seq(1, length(unique(trade_weight_monthly$period)), 3)]) +
  labs(x = "Period", y = "Trade volume (million tonnes)", fill = 'Category') +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.1),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 6, hjust = 1),
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
    axis.title = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.1),
    legend.title = element_text(size = 7),
    legend.key.size = unit(2, "mm"),
    legend.text = element_text(size = 7),
    plot.title = element_blank()
  )
dir.create("Figures/fig_1", showWarnings = FALSE, recursive = TRUE)
ggsave("Figures/fig_1/fig_1a_trade_weight.png", fig_1a, unit = 'mm', width = 150, height = 100, device = 'png', dpi = 300)

# Fig. 1b-e: chord diagrams of trade over 2016-2024. Colors and sector order
# are fixed across diagrams (by total HS 2701-2716 weight); only countries
# above `label_threshold` of the category's own import + export weight are
# labelled.-----
country_pal <- data.frame(iso3 = nodes$iso3, color = viridis::viridis(length(nodes$iso3), option = "D", direction = -1))

traders_seq <- bind_rows(comtrade %>% rename(iso3 = from_iso3), comtrade %>% rename(iso3 = to_iso3)) %>%
  group_by(iso3) %>%
  summarise(weight = sum(weight, na.rm = TRUE)) %>%
  arrange(desc(weight))

label_threshold <- 0.01
trader_shares <- list()

for (k in seq_along(hs_codes)) {
  hs_code <- hs_codes[k]
  edges <- comtrade %>%
    filter(hs == hs_code) %>%
    group_by(from_iso3, to_iso3) %>%
    summarise(weight = sum(weight, na.rm = TRUE), .groups = "drop") %>%
    filter(from_iso3 %in% traders_seq$iso3, to_iso3 %in% traders_seq$iso3) %>%
    na.omit()

  # Each country's share of the category's exports, imports, and import + export weight.
  trader_shares[[k]] <- full_join(
    edges %>% group_by(iso3 = from_iso3) %>% summarise(export = sum(weight)),
    edges %>% group_by(iso3 = to_iso3) %>% summarise(import = sum(weight)),
    by = "iso3"
  ) %>%
    mutate(across(c(export, import), ~ replace_na(.x, 0)), hs = hs_code,
           export_share = export / sum(export), import_share = import / sum(import),
           trade_share = (export + import) / sum(export + import))

  category_weight <- bind_rows(
    edges %>% select(iso3 = from_iso3, weight),
    edges %>% select(iso3 = to_iso3, weight)
  ) %>%
    group_by(iso3) %>%
    summarise(weight = sum(weight, na.rm = TRUE))
  share <- setNames(category_weight$weight / sum(category_weight$weight), category_weight$iso3)[traders_seq$iso3]
  share[is.na(share)] <- 0

  png(sprintf("Figures/fig_1/fig_1%s_chord_hs%d.png", letters[k + 1], hs_code), bg = "transparent",
      width = 150, height = 150, units = "mm", res = 300)
  par(mar = c(0, 0, 0, 0))
  mig_chord(
    x = edges,
    order = traders_seq$iso3,
    lab = setNames(ifelse(share > label_threshold, traders_seq$iso3, ""), traders_seq$iso3),
    gap.degree = 0.8,
    label_size = 0.6,
    no_axis = TRUE,
    zero_margin = TRUE,
    grid.col = country_pal$color[match(traders_seq$iso3, country_pal$iso3)]
  )
  title(main = paste0("HS", hs_code), cex.main = 1.5, line = -3)
  dev.off()
}

write.csv(bind_rows(trader_shares) %>% select(hs, iso3, export, import, export_share, import_share, trade_share) %>%
            arrange(hs, desc(trade_share)), "Data/trader_shares.csv", row.names = FALSE)
