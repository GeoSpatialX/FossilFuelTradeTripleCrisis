# Step 4. Crisis-related departures of bilateral energy trade from pre-crisis trends
# (Section 2.3, Fig. 2). PPML on the exporter -> importer pair x month panel
# from "1_data_process.R", one per crisis and HS code:
#
#   E[weight_ijt] = exp(gamma * After + phi * After * MonthsSince
#                       + delta * Time + MonthFE + Pair)
#
# MonthFE = calendar-month fixed effects, Pair = exporter-importer fixed
# effects; standard errors clustered by pair (fixest::fepois).-----
source("Scripts/functions.R")

comtrade <- read.csv("Data/UNComtrade/comtrade_mirror.csv") %>%
  filter(hs %in% hs_codes, weight > 0) %>%
  transmute(
    ym = year * 12 + month - 1,  # month index, continuous across years
    year, month,
    pair = paste(from_iso3, to_iso3, sep = "->"),
    hs, weight
  )

# Panel for one crisis x HS code: every month of the window for each pair
# active in >= `min_active_months` baseline months, missing months as 0.
# `onset_ym` is the (possibly placebo) onset, `end_ym` the last month kept.-----
buildPanel <- function(df, onset_ym, obs_len = 12, baseline_len = 24, min_active_months = 6,
                       end_ym = onset_ym + obs_len - 1)
{
  start_ym <- onset_ym - baseline_len
  active_pairs <- df %>%
    filter(ym >= start_ym, ym < onset_ym) %>%
    count(pair) %>%
    filter(n >= min_active_months) %>%
    pull(pair)

  tidyr::crossing(pair = active_pairs, ym = start_ym:end_ym) %>%
    left_join(df %>% filter(ym >= start_ym, ym <= end_ym) %>% select(pair, ym, weight), by = c("pair", "ym")) %>%
    mutate(
      weight = ifelse(is.na(weight), 0, weight),
      month = ym %% 12 + 1,
      time = ym - start_ym + 1,
      after = as.integer(ym >= onset_ym),
      months_since = pmax(ym - onset_ym, 0) * after,
      post_month = factor(ifelse(after == 1, as.character(months_since), "baseline"),
                          levels = c("baseline", as.character(0:(obs_len - 1)))),
      month_fe = factor(month, levels = 1:12)
    )
}

# PPML with absorbed pair and calendar-month fixed effects, pair-clustered SEs
# with a G/(G-1) correction only. All-zero pairs and collinear terms are
# dropped by fepois.-----
fitPairPPML <- function(panel, rhs)
{
  fit <- fixest::fepois(
    as.formula(paste("weight ~", rhs, "| pair + month_fe")),
    data = panel, cluster = ~pair,
    ssc = fixest::ssc(adj = FALSE, cluster.adj = TRUE), notes = FALSE
  )
  data.frame(
    term = sub("^post_month", "post_month_", names(coef(fit))),
    estimate = unname(coef(fit)),
    std_error = unname(fixest::se(fit)),
    n_obs = fit$nobs,
    n_pairs = unname(fit$fixef_sizes["pair"])
  ) %>%
    mutate(t_value = estimate / std_error, p_value = 2 * pnorm(-abs(t_value)),
           pct_change = (exp(estimate) - 1) * 100)
}

impact <- bind_rows(lapply(names(crises), function(crisis_name) {
  onset_ym <- crises[[crisis_name]][["year"]] * 12 + crises[[crisis_name]][["month"]] - 1
  bind_rows(lapply(hs_codes, function(hs_code) {
    df <- comtrade %>% filter(hs == hs_code)
    panel <- buildPanel(df, onset_ym)
    # Placebo onset 8 months earlier, inside the real baseline (16 + 8 = 24 months).
    panel_placebo <- buildPanel(df, onset_ym - 8, end_ym = onset_ym - 1, baseline_len = 16)
    bind_rows(
      fitPairPPML(panel, "after + months_since + time") %>% mutate(model = "main"),
      fitPairPPML(panel, "post_month + time") %>% mutate(model = "dynamic"),
      fitPairPPML(panel_placebo, "after + months_since + time") %>% mutate(model = "placebo")
    ) %>%
      mutate(crisis = crisis_name, hs = hs_code)
  }))
})) %>%
  mutate(term = case_when(
    term == "after" ~ "gamma (level shift at onset)",
    term == "months_since" ~ "phi (monthly slope change after onset)",
    term == "time" ~ "delta (pre-crisis monthly trend)",
    grepl("^post_month_", term) ~ paste("month", sub("^post_month_", "", term), "after onset"),
    TRUE ~ term
  )) %>%
  select(crisis, hs, model, term, estimate, std_error, t_value, p_value, pct_change, n_obs, n_pairs)

write.csv(impact, "Data/crisis_impact_estimates.csv", row.names = FALSE)

# Fitted main-model path, gamma + phi * months_since, as % of the counterfactual.-----
impact_path <- impact %>%
  filter(model == "main") %>%
  select(crisis, hs, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(gamma = `gamma (level shift at onset)`, phi = `phi (monthly slope change after onset)`) %>%
  tidyr::crossing(months_since = 0:11) %>%
  mutate(log_effect = gamma + phi * months_since, pct_effect = (exp(log_effect) - 1) * 100)

write.csv(impact_path, "Data/crisis_impact_fitted_path.csv", row.names = FALSE)

# Fig. 2 panels, all one size and without titles, axis titles or legends
# (added in the figure layout; legends saved once as PDFs).-----
panel_theme <- theme_minimal(base_size = 9, base_family = fontstyle) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 9),
    axis.title = element_blank(),
    axis.line = element_line(colour = "grey20", linewidth = 0.3),
    axis.ticks = element_line(colour = "grey20", linewidth = 0.3),
    axis.ticks.length = grid::unit(2, "pt")
  )
savePanel <- function(p, file) ggsave(figPath("fig_2", file), p, width = 105, height = 80, unit = "mm", dpi = 300)
ci <- function(df) mutate(df, lower = (exp(estimate - 1.96 * std_error) - 1) * 100,
                          upper = (exp(estimate + 1.96 * std_error) - 1) * 100,
                          hs = factor(hs, levels = hs_codes))

# Fig. 2d-o: monthly estimates with 95% bands (solid) over the fitted main
# path (dashed). x = months since onset (0-11), y = effect on bilateral trade
# weight (%). Rows = HS code (d-f coal ... m-o petroleum gases), columns = crisis.-----
effectPanel <- function(crisis_name, hs_code)
{
  dyn <- impact %>%
    filter(crisis == crisis_name, model == "dynamic", hs == hs_code, grepl("^month ", term)) %>%
    mutate(months_since = as.integer(sub("month (\\d+) after onset", "\\1", term))) %>%
    ci()
  path <- impact_path %>%
    filter(crisis == crisis_name, hs == hs_code) %>%
    mutate(hs = factor(hs, levels = hs_codes))

  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_ribbon(data = dyn, aes(x = months_since, ymin = lower, ymax = upper, fill = hs), alpha = 0.15) +
    geom_line(data = dyn, aes(x = months_since, y = pct_change, color = hs), linewidth = 0.4) +
    geom_point(data = dyn, aes(x = months_since, y = pct_change, color = hs), size = 0.8) +
    geom_line(data = path, aes(x = months_since, y = pct_effect, color = hs), linetype = "dashed", linewidth = 0.4) +
    scale_x_continuous(breaks = path$months_since) +
    scale_color_manual(values = hs_color_palette, drop = FALSE) +
    scale_fill_manual(values = hs_color_palette, drop = FALSE) +
    labs(color = "HS code") +
    panel_theme
}

# Fig. 2a-c: level shift at onset (gamma), real vs placebo onset. x = HS
# code, y = level shift at onset (%).-----
gamma_plot <- impact %>%
  filter(model %in% c("main", "placebo"), grepl("^gamma", term)) %>%
  ci() %>%
  mutate(model = factor(model, levels = c("main", "placebo"),
                        labels = c("Real onset", "Placebo onset (8 months earlier)")))

gammaPanel <- function(crisis_name)
{
  ggplot(gamma_plot %>% filter(crisis == crisis_name), aes(x = hs, y = pct_change, color = hs, shape = model)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, position = position_dodge(width = 0.6)) +
    geom_point(position = position_dodge(width = 0.6), size = 2) +
    scale_color_manual(values = hs_color_palette, guide = "none", drop = FALSE) +
    scale_shape_manual(values = c(16, 1), drop = FALSE) +
    labs(shape = NULL) +
    panel_theme
}

for (k in seq_along(crises)) {
  crisis_name <- names(crises)[k]
  savePanel(gammaPanel(crisis_name), sprintf("fig_2%s_%s_level_shift.png", letters[k], crisis_name))
  for (h in seq_along(hs_codes)) {
    savePanel(effectPanel(crisis_name, hs_codes[h]), sprintf("fig_2%s_%s_hs%d.png", letters[3 * h + k], crisis_name, hs_codes[h]))
  }
}

# Shared legends: HS-code colour key (line + point) and real/placebo shape key.-----
saveLegend(
  ggplot(data.frame(x = c(1, 2), y = 1, hs = factor(rep(hs_codes, each = 2), levels = hs_codes)),
         aes(x = x, y = y, color = hs)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = hs_color_palette) +
    labs(color = "HS code") +
    theme_minimal(base_size = 9, base_family = fontstyle),
  figPath("fig_2", "fig_2_legend_hs.pdf")
)
saveLegend(
  gammaPanel(names(crises)[1]) +
    theme(legend.position = "right", legend.text = element_text(family = fontstyle, size = 9)),
  figPath("fig_2", "fig_2_legend_onset.pdf"), width = 3, height = 1
)

