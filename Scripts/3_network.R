# Step 3. Network-level response to the three crises (Sections 2.2 and 3.2):
# node/network property tables, year-on-year change (Fig. A2), link
# similarity (Fig. A3) and Infomap communities (Figs. A4-A7).
source("Scripts/functions.R")
library(terra)
library(tidyterra)
library(ggforce)

# Node- and network-level properties of every monthly network built in
# "2_build_nw.R", in chronological order (list.files() sorts m1, m10, m11, m2).-----
nw_files <- list.files("Data/Network", pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
file_ym <- do.call(rbind, regmatches(basename(nw_files), regexec("_y(\\d+)_m(\\d+)\\.rds$", basename(nw_files))))
file_hs <- sub("^hs", "", basename(dirname(nw_files)))
file_order <- order(file_hs, as.integer(file_ym[, 2]), as.integer(file_ym[, 3]))

# Effective number of partners: inverse Herfindahl index of a node's partner
# shares over `group_col` ("from" = ENP_out, "to" = ENP_in); NA without edges.-----
calcEnp <- function(edges_df, group_col, n_nodes)
{
  shares <- edges_df %>%
    group_by(.data[[group_col]]) %>%
    mutate(total = sum(weight, na.rm = TRUE)) %>%
    summarise(hhi = sum((weight / total)^2, na.rm = TRUE), total = first(total), .groups = "drop")
  out <- rep(NA_real_, n_nodes)
  out[shares[[group_col]]] <- ifelse(shares$total == 0, NA_real_, 1 / shares$hhi)
  out
}

properties <- lapply(file_order, function(i)
{
  hs <- file_hs[i]
  y <- as.integer(file_ym[i, 2])
  m <- as.integer(file_ym[i, 3])
  rds <- readRDS(nw_files[i])
  nodes_df <- as.data.frame(activate(rds$nw, nodes))
  edges_df <- as.data.frame(activate(rds$nw, edges)) %>% select(from, to, weight)
  n_nodes <- nrow(nodes_df)

  list(
    nw = data.frame(
      hs = hs, year = y, month = m,
      metric = c(names(rds$nw_metrics), "global_strength", "global_density"),
      value = c(as.numeric(rds$nw_metrics), sum(edges_df$weight, na.rm = TRUE),
                igraph::edge_density(rds$nw, loops = FALSE))
    ),
    node = nodes_df %>%
      select(-lon, -lat) %>%
      pivot_longer(-c(label, iso3), names_to = "metric", values_to = "value") %>%
      mutate(hs = hs, year = y, month = m) %>%
      select(hs, year, month, iso3, label, metric, value) %>%
      bind_rows(data.frame(
        hs = hs, year = y, month = m, iso3 = nodes_df$iso3, label = nodes_df$label,
        metric = rep(c("enp_out", "enp_in"), each = n_nodes),
        value = c(calcEnp(edges_df, "from", n_nodes), calcEnp(edges_df, "to", n_nodes))
      ))
  )
})

nw_properties <- bind_rows(lapply(properties, `[[`, "nw"))
write.csv(nw_properties, "Data/nw_properties.csv", row.names = FALSE)
write.csv(bind_rows(lapply(properties, `[[`, "node")), "Data/node_properties.csv", row.names = FALSE)
rm(properties)

# Year-on-year change of network-level metrics: each observation month vs the
# same calendar month a year earlier.-----
nw_metric <- c("global_strength", "global_density", "global_efficiency")

crisis_response_nw <- yoyByCrisis(
  nw_properties %>% filter(hs %in% hs_codes, metric %in% nw_metric), "metric",
  function(d, period_map) tidyr::crossing(period_map, metric = nw_metric)
)
write.csv(crisis_response_nw, "Data/crisis_response_nw.csv", row.names = FALSE)

# Year-on-year link similarity, month by month (see nwSimilarity()).-----
crisis_similarity <- bind_rows(lapply(names(crises), function(crisis_name) {
  b <- baselineYM(crises[[crisis_name]])
  o <- observationYM(crises[[crisis_name]])
  bind_rows(lapply(hs_codes, function(hs_code) bind_rows(lapply(seq_len(nrow(b)), function(i) data.frame(
    crisis = crisis_name, hs = hs_code, period_index = i, month = b$month[i],
    baseline_year = b$year[i], observation_year = o$year[i],
    nwSimilarity(getEdges(hs_code, b$year[i], b$month[i]), getEdges(hs_code, o$year[i], o$month[i]))
  )))))
})) %>%
  select(crisis, hs, period_index, month, baseline_year, observation_year, n_edges_baseline,
         n_edges_observation, jaccard, retention_rate, formation_rate, dissolution_rate,
         bray_curtis, bray_curtis_share)
write.csv(crisis_similarity, "Data/crisis_similarity.csv", row.names = FALSE)

# Infomap communities on each period's 12 months of edges summed into one
# network, for the baseline and the observation period.-----
aggEdges <- function(hs_code, ym) bind_rows(lapply(seq_len(nrow(ym)), function(i) getEdges(hs_code, ym$year[i], ym$month[i]))) %>%
  group_by(from_iso3, to_iso3) %>%
  summarise(weight = sum(weight, na.rm = TRUE), .groups = "drop")

# `nodes_ref` fixes the 239-country vertex set, so isolated countries still
# get a community.
detectCommunity <- function(edges, nodes_ref)
{
  g <- tbl_graph(nodes = nodes_ref, edges = rename(edges, from = from_iso3, to = to_iso3), directed = TRUE, node_key = "iso3")
  data.frame(iso3 = nodes_ref$iso3, community = igraph::cluster_infomap(g, e.weights = as.data.frame(activate(g, edges))$weight)$membership)
}

# Relabels `target_comm`'s communities to match `original_comm`'s by the
# one-to-one mapping that maximises shared members (Hungarian algorithm);
# unmatched target communities get fresh labels. No country changes community.-----
alignCommunity <- function(original_comm, target_comm)
{
  original_labels <- sort(unique(original_comm$community))
  target_labels <- sort(unique(target_comm$community))
  joined <- inner_join(rename(original_comm, original = community), rename(target_comm, target = community), by = "iso3")
  overlap <- matrix(
    table(factor(joined$original, levels = original_labels), factor(joined$target, levels = target_labels)),
    nrow = length(original_labels), ncol = length(target_labels),
    dimnames = list(as.character(original_labels), as.character(target_labels))
  )
  # solve_LSAP needs ncol >= nrow; zero-overlap dummies leave a baseline-only community unmatched.
  pad <- nrow(overlap) - ncol(overlap)
  cost <- if (pad <= 0) overlap else
    cbind(overlap, matrix(0, nrow(overlap), pad, dimnames = list(NULL, paste0("dummy", seq_len(pad)))))

  pairs <- data.frame(original = rownames(cost), target = colnames(cost)[clue::solve_LSAP(cost, maximum = TRUE)]) %>%
    filter(target %in% colnames(overlap))
  pairs <- pairs[overlap[cbind(pairs$original, pairs$target)] > 0, ]

  label_map <- setNames(as.numeric(pairs$original), pairs$target)
  unmatched <- setdiff(as.character(target_labels), names(label_map))
  label_map[unmatched] <- max(as.numeric(original_labels)) + seq_along(unmatched)
  target_comm %>% mutate(community = as.numeric(label_map[as.character(community)]))
}

community_nodes_ref <- readRDS("Data/Network/hs2701/weight_y2016_m1.rds")$nw %>%
  activate(nodes) %>%
  as.data.frame() %>%
  select(label, iso3)

# Infomap is stochastic; the seed makes the partitions reproducible. The
# community changes reported in Section 3.2 held for every one of 31 seeds tested.
set.seed(1)
community_detection <- bind_rows(lapply(names(crises), function(crisis_name) {
  ref <- crises[[crisis_name]]
  bind_rows(lapply(hs_codes, function(hs_code) {
    original_comm <- detectCommunity(aggEdges(hs_code, baselineYM(ref)), community_nodes_ref)
    target_comm <- alignCommunity(original_comm, detectCommunity(aggEdges(hs_code, observationYM(ref)), community_nodes_ref))
    bind_rows(
      original_comm %>% mutate(crisis = crisis_name, hs = hs_code, period = "baseline"),
      target_comm %>% mutate(crisis = crisis_name, hs = hs_code, period = "observation")
    )
  }))
})) %>%
  select(crisis, hs, period, iso3, community)
write.csv(community_detection, "Data/community_detection.csv", row.names = FALSE)

# Figs. A4-A7: community maps.-----
worldmap_community <- vect("Data/worldmap/worldmap.shp") %>%
  select(label, iso3, region) %>%
  mutate(iso3 = toupper(iso3))

centroids_df <- terra::centroids(worldmap_community) %>%
  terra::crds() %>%
  as.data.frame() %>%
  setNames(c("lon", "lat")) %>%
  mutate(iso3 = worldmap_community$iso3)

# Community colors for one period. Five opaque anchors: China's community red,
# the USA's blue, Russia's gold (unless already taken), and Europe's two
# largest communities by European-country count (excluding the first three).
# Other communities with > 5 members take `extra_palette` by size at alpha
# 0.4; the rest stay grey. Content, not label identity, decides the color.-----
communityColors <- function(period_comm)
{
  community_of <- function(iso) { v <- period_comm$community[period_comm$iso3 == iso]; if (length(v) == 0) NA else v[1] }
  china_comm <- community_of("CHN")
  usa_comm <- community_of("USA")
  rus_comm <- community_of("RUS")

  europe_counts <- period_comm %>%
    filter(!(community %in% c(china_comm, usa_comm, rus_comm))) %>%
    left_join(region_lookup, by = "iso3") %>%
    filter(region == "Europe") %>%
    count(community, sort = TRUE)
  europe_comm <- if (nrow(europe_counts) >= 1) europe_counts$community[1] else NA
  europe2_comm <- if (nrow(europe_counts) >= 2) europe_counts$community[2] else NA

  levels_ <- sort(unique(period_comm$community))
  colors <- setNames(rep("grey75", length(levels_)), levels_)
  colors[as.character(europe_comm)] <- "#4C9F70"
  colors[as.character(europe2_comm)] <- "#3C8F99"
  colors[as.character(usa_comm)] <- "#3B6FA0"
  colors[as.character(china_comm)] <- "#C1443D"
  if (!is.na(rus_comm) && !(rus_comm %in% c(china_comm, usa_comm, europe_comm, europe2_comm))) {
    colors[as.character(rus_comm)] <- "#E0A526"
  }

  extra_palette <- c("#756300", "#83932C", "#337100", "#009BB6", "#585DA3",
                     "#9B7DC8", "#884A9A", "#B773BD", "#9A4288", "#C86EA0", "#9F4466")
  extra_ids <- period_comm %>%
    filter(!(community %in% c(china_comm, usa_comm, europe_comm, europe2_comm, rus_comm))) %>%
    count(community, sort = TRUE) %>%
    filter(n > 5) %>%
    pull(community) %>%
    head(length(extra_palette)) %>%
    as.character()
  dot_alpha <- setNames(rep(1, length(levels_)), levels_)
  colors[extra_ids] <- extra_palette[seq_along(extra_ids)]
  dot_alpha[extra_ids] <- 0.4
  list(color = colors, alpha = dot_alpha)
}

# All 24 maps' data first, so the dot-size scale is shared across them.
# Panels: rows = crises, columns = baseline (left) / observation (right).-----
community_fig <- c("2701" = "fig_a4", "2709" = "fig_a5", "2710" = "fig_a6", "2711" = "fig_a7")
community_map_data <- list()
for (k in seq_along(crises)) for (hs_code in hs_codes) for (p in 1:2) {
  crisis_name <- names(crises)[k]
  period_name <- c("baseline", "observation")[p]
  period_comm <- community_detection %>% filter(crisis == crisis_name, hs == hs_code, period == period_name)
  period_colors <- communityColors(period_comm)
  period_edges <- aggEdges(hs_code, list(baselineYM, observationYM)[[p]](crises[[crisis_name]]))

  # Total import + export weight per country sizes the dots.
  weight_lookup <- bind_rows(
    period_edges %>% select(iso3 = from_iso3, weight),
    period_edges %>% select(iso3 = to_iso3, weight)
  ) %>%
    group_by(iso3) %>%
    summarise(total_weight = sum(weight, na.rm = TRUE), .groups = "drop")

  comm_plot <- period_comm %>%
    left_join(centroids_df, by = "iso3") %>%
    left_join(weight_lookup, by = "iso3") %>%
    mutate(
      color = period_colors$color[as.character(community)],
      dot_alpha = period_colors$alpha[as.character(community)],
      total_weight = replace_na(total_weight, 0)
    )

  # Edges take their community's color only when both ends share it.
  comm_lookup <- setNames(comm_plot$community, comm_plot$iso3)
  edges_plot <- period_edges %>%
    left_join(rename(centroids_df, lon_from = lon, lat_from = lat), by = c("from_iso3" = "iso3")) %>%
    left_join(rename(centroids_df, lon_to = lon, lat_to = lat), by = c("to_iso3" = "iso3")) %>%
    mutate(
      comm_from = comm_lookup[from_iso3],
      comm_to = comm_lookup[to_iso3],
      within_community = !is.na(comm_from) & comm_from == comm_to,
      edge_color = ifelse(within_community, period_colors$color[as.character(comm_from)], "grey55"),
      edge_alpha = ifelse(within_community, 0.15, 0.03)
    )

  fig <- community_fig[[as.character(hs_code)]]
  community_map_data[[length(community_map_data) + 1]] <- list(
    file = figPath(fig, sprintf("%s%s_%s_%s.png", fig, letters[(k - 1) * 2 + p], crisis_name, period_name)),
    comm_plot = comm_plot, edges_plot = edges_plot
  )
}

community_weight_max <- max(sapply(community_map_data, function(d) max(d$comm_plot$total_weight)))

# Trade weight -> dot radius (lon/lat degrees), sqrt-scaled with a floor.
weightToRadius <- function(w, r_range = c(0.15, 6))
{
  r_range[1] + (r_range[2] - r_range[1]) * sqrt(pmin(w, community_weight_max) / community_weight_max)
}

# geom_circle(), not geom_point(): point glyphs leave a ring at alpha < 1.
# ylim crops the empty polar margins (all centroids lie within -54.5 to 78.9).
for (d in community_map_data) {
  p <- ggplot() +
    geom_spatvector(data = worldmap_community, fill = "grey88", color = "white", linewidth = 0.05) +
    geom_segment(
      data = d$edges_plot,
      aes(x = lon_from, y = lat_from, xend = lon_to, yend = lat_to, color = edge_color, alpha = edge_alpha),
      linewidth = 0.15
    ) +
    geom_circle(data = d$comm_plot %>% mutate(r = weightToRadius(total_weight)),
                aes(x0 = lon, y0 = lat, r = r, fill = color, alpha = dot_alpha), color = NA) +
    scale_color_identity() +
    scale_fill_identity() +
    scale_alpha_identity() +
    coord_sf(ylim = c(-60, 80), expand = FALSE) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  ggsave(d$file, p, width = 200, height = 78, unit = "mm", dpi = 300)
}

# Node-size legend shared by Figs. A4-A7: stacked reference circles (a
# geom_circle() radius has no automatic legend), breaks within the data range.-----
community_legend_breaks <- pretty(c(0, community_weight_max), n = 4)
community_legend_data <- data.frame(value = community_legend_breaks[community_legend_breaks <= community_weight_max]) %>%
  mutate(r = weightToRadius(value)) %>%
  arrange(desc(r)) %>%
  mutate(x = 0)
community_legend_data$y <- Reduce(
  function(y_prev, i) y_prev - community_legend_data$r[i - 1] - 0.3 - community_legend_data$r[i],
  seq(2, nrow(community_legend_data)), init = -community_legend_data$r[1], accumulate = TRUE
)
label_x <- max(community_legend_data$r) + 0.6

community_legend <- ggplot(community_legend_data) +
  geom_circle(aes(x0 = x, y0 = y, r = r), fill = "grey45", color = NA) +
  geom_text(aes(x = label_x, y = y, label = format(value, big.mark = ",", scientific = FALSE)),
            hjust = 0, size = 3, family = fontstyle) +
  labs(title = "Trade weight") +
  coord_fixed(clip = "off") +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(family = fontstyle, size = 10, hjust = 0), plot.margin = margin(5, 30, 5, 5))

for (fig in community_fig) {
  grDevices::cairo_pdf(figPath(fig, paste0(fig, "_legend.pdf")), width = 2.2, height = 3)
  print(community_legend)
  dev.off()
}

# Figs. A2 and A3: crises down the rows, metrics across the columns, each row
# on its own calendar dates (facet_grid2 `independent = "x"`; plain
# facet_grid would share x per column). Panels are tagged a, b, c...-----
crisis_label <- c(
  us_china_trade_tensions = "US-China trade tensions",
  covid_19_pandemic = "COVID-19 pandemic",
  russia_ukraine_conflict = "Russia-Ukraine conflict"
)
nw_metric_label <- c(global_strength = "Trade volume", global_density = "Density", global_efficiency = "Efficiency")
nw_similarity_label <- c(jaccard = "Jaccard", retention_rate = "Retention rate",
                         formation_rate = "Formation rate", bray_curtis_share = "Bray-Curtis share")
yoy_label <- "Year-on-year change (%)"

nwGrid <- function(df, metric_levels, y_label, y_scale)
{
  df <- df %>%
    mutate(
      date = as.Date(sprintf("%d-%02d-01", observation_year, month)),
      crisis = factor(crisis, levels = names(crises), labels = crisis_label),
      hs = factor(hs)
    )
  tags <- tidyr::crossing(crisis = factor(crisis_label, levels = crisis_label),
                          metric = factor(metric_levels, levels = metric_levels)) %>%
    arrange(crisis, metric) %>%
    mutate(label = letters[row_number()])

  ggplot(df, aes(x = date, y = value, color = hs)) +
    (if (y_label == yoy_label) geom_hline(yintercept = 0, linetype = "dashed", color = "grey40")) +
    geom_line(linewidth = 0.4) +
    geom_point(size = 0.7) +
    geom_text(data = tags, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
              hjust = -0.6, vjust = 1.4, family = fontstyle, size = 3.2, fontface = "bold") +
    ggh4x::facet_grid2(rows = vars(crisis), cols = vars(.data[["metric"]]),
                       scales = "free_x", independent = "x", switch = "y") +
    scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
    y_scale +
    scale_color_manual(values = hs_color_palette) +
    labs(x = "Month", y = y_label, color = "HS code") +
    theme_minimal(base_size = 9, base_family = fontstyle) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_blank(),
      legend.position = "bottom",
      strip.placement = "outside",
      strip.text.x = element_text(size = 9),
      strip.text.y.left = element_text(size = 8, angle = 90),
      axis.text = element_text(size = 7),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(size = 11.25),
      axis.line = element_line(colour = "grey20", linewidth = 0.3),
      axis.ticks = element_line(colour = "grey20", linewidth = 0.3),
      axis.ticks.length = grid::unit(2, "pt")
    )
}

# One symmetric y-range for Fig. A2, sized to the largest |change|.
nw_y_max <- max(abs(crisis_response_nw$pct_change), na.rm = TRUE)
nw_y_breaks <- c(-30, -20, -10, 0, 10, 20, 30)
if (nw_y_max >= 40) nw_y_breaks <- sort(c(nw_y_breaks, -40, 40))

fig_a2 <- crisis_response_nw %>%
  mutate(value = pct_change, metric = factor(metric, levels = nw_metric, labels = nw_metric_label[nw_metric])) %>%
  nwGrid(unname(nw_metric_label[nw_metric]), yoy_label,
         scale_y_continuous(limits = c(-nw_y_max, nw_y_max), breaks = nw_y_breaks))
ggsave(figPath("fig_a2", "fig_a2.png"), fig_a2, width = 190, height = 185, unit = "mm", dpi = 300)

fig_a3 <- crisis_similarity %>%
  pivot_longer(all_of(names(nw_similarity_label)), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = names(nw_similarity_label), labels = nw_similarity_label)) %>%
  nwGrid(unname(nw_similarity_label), "Index (0-1)", scale_y_continuous(limits = c(0, 1)))
ggsave(figPath("fig_a3", "fig_a3.png"), fig_a3, width = 190, height = 185, unit = "mm", dpi = 300)
