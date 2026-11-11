# Step 5. Node-level response to the three crises: the three resilience
# measures (Section 2.4), their correlation with pre-crisis node properties
# (Fig. 3, Fig. A8), and the partner-composition case studies (Figs. 4-5).
# Input: Data/node_properties.csv from "3_network.R".
source("Scripts/functions.R")

node_metric <- c("degree_in", "degree_out", "strength_in", "strength_out", "enp_in", "enp_out",
                 "betweenness", "transitivity")
strength_metric <- c("strength_in", "strength_out")

# =============================================================================
# 1. Year-on-year node metrics over every country, and the three resilience
# measures (strength) plus pre-crisis baselines (all node metrics).
# =============================================================================
crisis_response_node <- yoyByCrisis(
  read.csv("Data/node_properties.csv") %>% filter(hs %in% hs_codes, metric %in% node_metric),
  c("iso3", "metric"),
  function(d, period_map) tidyr::crossing(period_map, iso3 = unique(d$iso3[d$metric %in% strength_metric]), metric = node_metric)
)

# One node-direction's 12-month series of proportional YoY change `y`,
# clipped to [-1, cap] (downside only), summarised so that higher = more
# resilient: shock depth D = min(0, y); cumulative loss C = 1 + sum(min(y, 0)) / n;
# end state E = mean of the last three months.-----
nodeResilience <- function(y, cap = 0)
{
  if (all(is.na(y))) return(data.frame(shock_depth = NA_real_, cumulative_loss = NA_real_, end_state = NA_real_))
  y <- pmin(pmax(y, -1), cap)
  end_state <- mean(y[max(1, length(y) - 2):length(y)], na.rm = TRUE)
  data.frame(
    shock_depth = min(c(0, y), na.rm = TRUE),
    cumulative_loss = 1 + sum(pmin(y, 0), na.rm = TRUE) / length(y),
    end_state = if (is.nan(end_state)) NA_real_ else end_state
  )
}
resilience_measure_label <- c(shock_depth = "Shock depth", cumulative_loss = "Cumulative loss", end_state = "End state")

node_resilience <- crisis_response_node %>%
  filter(metric %in% strength_metric) %>%
  mutate(direction = ifelse(metric == "strength_in", "import", "export"), y = pct_change / 100) %>%
  group_by(crisis, hs, iso3, direction) %>%
  arrange(period_index, .by_group = TRUE) %>%
  summarise(nodeResilience(y), .groups = "drop")

# Pre-crisis position = mean of the 12 baseline months per metric; degree,
# strength and ENP are direction-matched, betweenness and transitivity not.-----
node_baseline_avg <- crisis_response_node %>%
  group_by(crisis, hs, iso3, metric) %>%
  summarise(baseline_mean = mean(baseline_value, na.rm = TRUE), .groups = "drop")

baseline_metric <- list(
  import = c(degree = "degree_in", strength = "strength_in", enp = "enp_in", betweenness = "betweenness", transitivity = "transitivity"),
  export = c(degree = "degree_out", strength = "strength_out", enp = "enp_out", betweenness = "betweenness", transitivity = "transitivity")
)

resilience_data <- bind_rows(lapply(names(baseline_metric), function(dir_name) {
  m <- baseline_metric[[dir_name]]
  Reduce(function(out, as_name) left_join(
    out,
    node_baseline_avg %>% filter(metric == m[[as_name]]) %>% select(crisis, hs, iso3, !!as_name := baseline_mean),
    by = c("crisis", "hs", "iso3")
  ), names(m), node_resilience %>% filter(direction == dir_name))
})) %>%
  left_join(region_lookup, by = "iso3") %>%
  mutate(region = ifelse(is.na(region), "Unknown", region))
write.csv(resilience_data, "Data/resilience_data.csv", row.names = FALSE)

# =============================================================================
# 2. Correlation of each resilience measure with each baseline node metric,
# per crisis x HS code x direction (Fig. 3, Fig. A8).
# =============================================================================
correlation_metric_label <- c(degree = "Degree", strength = "Strength", enp = "Effective partners",
                              betweenness = "Betweenness", transitivity = "Transitivity")

# Pearson if both vectors pass a Shapiro-Wilk test, Spearman otherwise (all
# 360 tests turn out Spearman); NA for a group too small or degenerate to
# test. Inputs are rounded to 10 significant digits so that values tied in
# exact arithmetic stay tied in ranks, whatever floating-point noise the
# platform or a CSV round trip adds.-----
corTest <- function(x, y, alpha = 0.05)
{
  ok <- is.finite(x) & is.finite(y)
  x <- signif(x[ok], 10)
  y <- signif(y[ok], 10)
  n <- length(x)
  if (n < 4 || sd(x) == 0 || sd(y) == 0) {
    return(data.frame(n = n, method = NA_character_, shapiro_p_x = NA_real_,
                      shapiro_p_y = NA_real_, estimate = NA_real_, p_value = NA_real_))
  }
  shapiroP <- function(v) stats::shapiro.test(if (n > 5000) sample(v, 5000) else v)$p.value
  p_x <- shapiroP(x)
  p_y <- shapiroP(y)
  method <- if (p_x < alpha || p_y < alpha) "spearman" else "pearson"
  ct <- suppressWarnings(stats::cor.test(x, y, method = method, exact = FALSE))
  data.frame(n = n, method = method, shapiro_p_x = p_x, shapiro_p_y = p_y,
             estimate = unname(ct$estimate), p_value = ct$p.value)
}

# Benjamini-Hochberg adjustment over all 360 coefficients.-----
set.seed(1)
correlation_grid <- expand_grid(
  crisis = names(crises), hs = hs_codes, direction = c("import", "export"),
  resilience_measure = names(resilience_measure_label), node_metric = names(correlation_metric_label)
)
resilience_correlation <- bind_rows(lapply(seq_len(nrow(correlation_grid)), function(i) {
  g <- correlation_grid[i, ]
  d <- resilience_data %>% filter(crisis == g$crisis, hs == g$hs, direction == g$direction)
  data.frame(g, corTest(d[[g$resilience_measure]], d[[g$node_metric]]))
})) %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significance = cut(p_adjusted, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", "ns"))
  ) %>%
  select(crisis, hs, direction, resilience_measure, node_metric, n, method,
         shapiro_p_x, shapiro_p_y, estimate, p_value, p_adjusted, significance)
write.csv(resilience_correlation, "Data/resilience_correlation.csv", row.names = FALSE)

# Node metrics ordered by median coefficient.-----
crisis_short_label <- c(us_china_trade_tensions = "U.S.–China", covid_19_pandemic = "COVID-19",
                        russia_ukraine_conflict = "Russia–Ukraine")
metric_order <- resilience_correlation %>%
  group_by(node_metric) %>%
  summarise(med = median(estimate, na.rm = TRUE), .groups = "drop") %>%
  arrange(med) %>%
  pull(node_metric)

correlation_plot_data <- resilience_correlation %>%
  filter(!is.na(estimate)) %>%
  mutate(
    hs = factor(hs, levels = hs_codes),
    crisis = factor(crisis, levels = names(crises), labels = crisis_short_label),
    direction = factor(direction, levels = c("import", "export"), labels = c("Import", "Export")),
    resilience_measure = factor(resilience_measure, levels = names(resilience_measure_label), labels = resilience_measure_label),
    node_metric = factor(node_metric, levels = metric_order, labels = correlation_metric_label[metric_order])
  )

# Fig. 3: every coefficient, diverging fill centred on 0, labelled with its
# value and significance stars; import blocks above export blocks.-----
correlation_limit <- max(abs(correlation_plot_data$estimate))
fig_3 <- ggplot(correlation_plot_data, aes(x = hs, y = node_metric, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = paste0(sprintf("%.2f", estimate), ifelse(significance == "ns", "", as.character(significance))),
        color = abs(estimate) > 0.45),
    size = 2.6, family = fontstyle, show.legend = FALSE
  ) +
  facet_grid(direction + resilience_measure ~ crisis) +
  scale_fill_gradient2(low = "#B23A48", mid = "#F0F0F0", high = "#3D5A80", midpoint = 0,
                       limits = c(-correlation_limit, correlation_limit), name = "Spearman correlation") +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  labs(x = "HS code", y = NULL) +
  theme_minimal(base_size = 9, base_family = fontstyle) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 8),
    strip.text.y = element_text(size = 8, angle = 0),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11.25)
  )
ggsave(figPath("fig_3", "fig_3_heatmap.png"), fig_3, width = 230, height = 175, unit = "mm", dpi = 300)
saveLegend(
  fig_3 +
    guides(fill = guide_colorbar(direction = "horizontal", barwidth = unit(6, "cm"), barheight = unit(0.4, "cm"))) +
    theme(legend.position = "bottom"),
  figPath("fig_3", "fig_3_legend.pdf"), width = 7, height = 1.2
)

# Fig. A8: the coefficients collapsed by node metric, crisis and HS code, as
# ridgelines per resilience measure (import blue, export red, overlaid; the
# rule in each ridge is its mean).-----
correlation_box_data <- bind_rows(lapply(c(`Node metric` = "node_metric", Crisis = "crisis", `HS code` = "hs"), function(col)
  correlation_plot_data %>% transmute(level = as.character(.data[[col]]), direction, resilience_measure, estimate)
), .id = "grouping") %>%
  mutate(
    grouping = factor(grouping, levels = c("Node metric", "Crisis", "HS code")),
    level = factor(level, levels = c(rev(as.character(correlation_metric_label[metric_order])),
                                     as.character(crisis_short_label), as.character(hs_codes)))
  )

fig_a8 <- ggplot(
  correlation_box_data %>% mutate(level = factor(level, levels = rev(levels(level)))),
  aes(x = estimate, y = level, fill = direction, vline_color = direction)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  ggridges::stat_density_ridges(
    quantile_lines = TRUE, quantile_fun = function(x, ...) mean(x, na.rm = TRUE), vline_linewidth = 0.45,
    scale = 0.95, rel_min_height = 0.01, colour = "white", linewidth = 0.25, alpha = 0.5
  ) +
  scale_discrete_manual(aesthetics = "vline_color", values = c(Import = "#1E3050", Export = "#6E1F2B"), guide = "none") +
  # facet_grid2 frees y per panel (not per row), dropping each panel's unused levels.
  ggh4x::facet_grid2(resilience_measure ~ grouping, scales = "free_y", independent = "y", switch = NULL) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(mult = c(0.02, 0.22))) +
  scale_fill_manual(values = c("#3D5A80", "#B23A48")) +
  labs(x = "Correlation with resilience measure", y = NULL, fill = NULL) +
  theme_minimal(base_size = 9, base_family = fontstyle) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 9),
    strip.text.y = element_text(size = 9, angle = -90),
    panel.spacing.x = grid::unit(6, "pt"),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11.25),
    axis.line = element_line(colour = "grey20", linewidth = 0.3),
    axis.ticks = element_line(colour = "grey20", linewidth = 0.3),
    axis.ticks.length = grid::unit(2, "pt")
  )
ggsave(figPath("fig_a8", "fig_a8_ridge.png"), fig_a8, width = 190, height = 165, unit = "mm", dpi = 300)
saveLegend(fig_a8 + theme(legend.position = "right"), figPath("fig_a8", "fig_a8_legend.pdf"), width = 2, height = 1.2)

# The numbers behind each ridge; `mean` is the rule drawn.-----
correlation_summary <- correlation_box_data %>%
  group_by(grouping, level, direction, resilience_measure) %>%
  summarise(n = n(), median = median(estimate), q1 = quantile(estimate, 0.25), q3 = quantile(estimate, 0.75),
            min = min(estimate), max = max(estimate), mean = mean(estimate), .groups = "drop") %>%
  mutate(across(c(median, q1, q3, min, max, mean), ~ round(.x, 3)))
write.csv(correlation_summary, "Data/resilience_correlation_summary.csv", row.names = FALSE)

# =============================================================================
# 3. Figs. 4-5: partner-region composition of selected countries' monthly
# import/export weight (all four HS codes) over the 24-month window.
# =============================================================================
# Panels a, b, ... = (country, direction) in the order listed.
case_studies <- list(
  fig_4 = list(crisis = "us_china_trade_tensions", iso3 = c("CHN", "USA")),
  fig_5 = list(crisis = "russia_ukraine_conflict", iso3 = c("CHN", "USA", "RUS"))
)
case_tick_index <- seq(1, 24, by = 3)

# One country's monthly weight by partner region.-----
partnerComposition <- function(month_edges, iso3, direction)
{
  own_edges <- if (direction == "import") {
    month_edges %>% filter(to_iso3 == iso3) %>% rename(partner = from_iso3)
  } else {
    month_edges %>% filter(from_iso3 == iso3) %>% rename(partner = to_iso3)
  }
  own_edges %>%
    left_join(region_lookup, by = c("partner" = "iso3")) %>%
    mutate(region = ifelse(is.na(region), "Unknown", region)) %>%
    group_by(period_index, region) %>%
    summarise(weight = sum(weight, na.rm = TRUE), .groups = "drop") %>%
    complete(period_index = unique(month_edges$period_index), region = names(region_palette), fill = list(weight = 0))
}

# 100%-stacked areas with a dashed line at onset. Tick labels are drawn as
# text layers because the ggplot2 4.0 axis guide ignores axis.text's margin
# and justification.-----
y_tick_labels <- data.frame(y = c(0, 0.25, 0.5, 0.75, 1)) %>% mutate(label = scales::percent(y))

plotPartnerComposition <- function(comp, tick_labels, direction)
{
  ggplot(comp %>% mutate(region = factor(region, levels = names(region_palette))),
         aes(x = period_index, y = weight, fill = region)) +
    geom_area(position = "fill", color = "white", linewidth = 0.15) +
    geom_vline(xintercept = 12.5, linetype = "dashed", color = "white", linewidth = 0.4) +
    geom_text(data = data.frame(x = case_tick_index, label = tick_labels[case_tick_index]),
              aes(x = x, y = -0.03, label = label),
              inherit.aes = FALSE, family = fontstyle, size = 9 / .pt, vjust = 1) +
    geom_text(data = y_tick_labels, aes(x = 0.8, y = y, label = label),
              inherit.aes = FALSE, family = fontstyle, size = 9 / .pt, hjust = 1) +
    scale_fill_manual(values = region_palette, drop = FALSE) +
    labs(x = "Month", y = sprintf("Share of %s weight", direction), fill = "Region") +
    theme_minimal(base_size = 9, base_family = fontstyle) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_text(size = 11.25),
      axis.ticks.length = unit(0, "pt"),
      legend.position = "none"
    )
}

region_legend <- ggplot(
  data.frame(x = 1, y = 1, region = factor(names(region_palette), levels = names(region_palette))),
  aes(x = x, y = y, fill = region)
) +
  geom_tile() +
  scale_fill_manual(values = region_palette, drop = FALSE) +
  labs(fill = "Region") +
  theme_minimal(base_size = 9, base_family = fontstyle) +
  theme(legend.position = "right")

partner_composition <- list()
for (fig in names(case_studies)) {
  cs <- case_studies[[fig]]
  ref <- crises[[cs$crisis]]
  ym <- bind_rows(baselineYM(ref), observationYM(ref)) %>% mutate(period_index = row_number())
  # Monthly edges summed over the four HS codes.
  month_edges <- bind_rows(lapply(seq_len(nrow(ym)), function(i) {
    bind_rows(lapply(hs_codes, function(hs_code) getEdges(hs_code, ym$year[i], ym$month[i]))) %>%
      group_by(from_iso3, to_iso3) %>%
      summarise(weight = sum(weight, na.rm = TRUE), .groups = "drop") %>%
      mutate(period_index = ym$period_index[i])
  }))
  panels <- expand_grid(iso3 = cs$iso3, direction = c("import", "export"))
  for (i in seq_len(nrow(panels))) {
    comp <- partnerComposition(month_edges, panels$iso3[i], panels$direction[i])
    partner_composition[[length(partner_composition) + 1]] <- comp %>%
      mutate(crisis = cs$crisis, iso3 = panels$iso3[i], direction = panels$direction[i])
    ggsave(figPath(fig, sprintf("%s%s_%s_%s.png", fig, letters[i], panels$iso3[i], panels$direction[i])),
           plotPartnerComposition(comp, ymLabel(ym$year, ym$month), panels$direction[i]),
           width = 180, height = 100, unit = "mm", dpi = 300)
  }
  saveLegend(region_legend, figPath(fig, paste0(fig, "_legend.pdf")))
}

# The shares plotted in Figs. 4-5 (period_index 13 = onset month).
write.csv(
  bind_rows(partner_composition) %>%
    group_by(crisis, iso3, direction, period_index) %>%
    mutate(share = weight / sum(weight)) %>%
    ungroup() %>%
    select(crisis, iso3, direction, period_index, region, weight, share),
  "Data/partner_composition.csv", row.names = FALSE
)
