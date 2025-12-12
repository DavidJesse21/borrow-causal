options(box.path = "R")

box::use(
  fs,
  data.table[...],
  ggplot2[...],
  scales[label_percent],
  patchwork[wrap_plots]
)

box::use(
  simfuns/metrics[...],
  simfuns/algos[dt_algos],
  simfuns/gen_data[make_scenarios]
)

pal_oi = function(which = NULL, unname = TRUE) {
  li_pal_oi = c(
    orange = "#E69F00",
    light_blue = "#56B4E9",
    green = "#009E73",
    yellow = "#F0E442",
    blue = "#0072B2",
    red = "#D55E00",
    purple = "#CC79A7",
    grey = "#999999",
    black = "#000000",
    sky_blue = "#56B4E9",
    dark_yellow = "#F5C710"
  )
  
  out = if (is.null(which)) li_pal_oi else li_pal_oi[which]
  out = if (unname) unname(out) else out
  
  return(out)
}

theme_set(theme_bw())



# Data ----

dt_algos[, algo.name := factor(
  c("Current data", "Pooled data", "PS-MEM", "PS-Commensurate", "BART"),
  levels = c("Current data", "Pooled data", "PS-MEM", "PS-Commensurate", "BART")
)]

dt_scenarios = make_scenarios()

dt_res = readRDS(fs$path("simulation", "results", "simres_2025-08-08", ext = "rds"))
dt_res = merge(dt_res, dt_scenarios[, .(scenario.id, b_trt)], by = "scenario.id") |>
  _[, .(bias = calc_bias(est, b_trt),
        emp_se = calc_emp_se(est),
        mse = calc_mse(est, b_trt),
        ci_coverage = calc_ci_coverage(ci_lower, ci_upper, b_trt),
        ci_mean_width = calc_ci_mean_width(ci_lower, ci_upper),
        rejection = calc_rejection_rate(pval, level = 0.025),
        num_na = calc_num_na(est)),
    by = .(scenario.id, algo.id)]



# Plots ----

## Supporting objects ----

# Reference lines for type 1 error and power
ref_lines = data.table(
  b_trt = c(0, 0, 4.5),
  yintercept = c(0.015, 0.035, 0.6)
)

# Space between plot and axis titles
axis_title_margins = theme(
  axis.title.x = element_text(margin = margin(t = 8, unit = "pt")),
  axis.title.y = element_text(margin = margin(r = 8, unit = "pt"))
)

# Opaque color for the legends
legend_color = guides(color = guide_legend(override.aes = list(alpha = 1)))


## Power / Type 1 error big picture ----

dt_scenarios[n_hist %in% c(60, 120) & b_trt %in% c(0, 4.5)]

# "Global" plot of type 1 error / power ~ shift of u
p_global_u = dt_res |>
  # Prepare data set for plotting
  _[algo.id != 2] |>
  merge(
    dt_scenarios[n_hist %in% c(60, 120) & b_trt %in% c(0, 4.5)],
    by = "scenario.id"
  ) |>
  merge(dt_algos, by = "algo.id") |>
  # Start plot
  ggplot(aes(x = factor(shift_u), y = rejection, color = algo.name)) +
  geom_point(position = position_jitterdodge(), alpha = 0.5, size = 2) +
  facet_grid(
    rows = vars(b_trt),
    cols = vars(n_hist),
    scales = "free_y",
    labeller = labeller(
      b_trt = \(x) paste0("Delta == ", x),
      n_hist = \(x) paste0("N[hist] == ", x),
      .default = label_parsed
    )
  ) +
  # Reference lines for nominal level and power for two-sample t-test
  geom_hline(
    data = ref_lines,
    aes(yintercept = yintercept),
    linetype = "dashed",
    alpha = 0.5
  ) +
  # Axes
  scale_y_continuous(
    name = "Type 1 error / Power (in %)",
    labels = scales::label_percent(suffix = NULL)
  ) +
  scale_x_discrete(
    name = expression("Mean shift of " * u * " (unobserved heterogeneity)")
  ) +
  axis_title_margins +
  # Legend
  theme(legend.position = "top") +
  scale_color_manual(
    name = "Method",
    values = pal_oi(c("blue", "orange", "red", "purple"))
  ) +
  legend_color

print(p_global_u)

ggsave(
  fs$path("simulation", "analysis", "plot_reject_u_slides", ext = "png"),
  plot = p_global_u,
  dpi = 600, width = 12, height = 7
)


## Type 1 error: Impact comparison ----

# In the following plot, we compare the impact of "confounding" that is either caused by
# prognostic baseline covariates or by latent/unobserved factors.
# To isolate this from other factors, we set the respective other parameter of the 
# data-generating model to 0.

p_t1e_u = dt_res |>
  # Prepare data set for plotting
  _[algo.id != 2] |>
  merge(dt_scenarios[b_trt == 0 & shift_x2 == 0], by = "scenario.id") |>
  merge(dt_algos, by = "algo.id") |>
  # Start plot
  ggplot(aes(x = factor(shift_u), y = rejection, color = algo.name)) +
  geom_point(position = position_jitterdodge(), alpha = 0.8, size = 2.5) +
  # Reference lines for nominal level and power for two-sample t-test
  geom_hline(
    data = ref_lines[1:2],
    aes(yintercept = yintercept),
    linetype = "dashed",
    alpha = 0.5
  ) +
  # Axes
  scale_y_continuous(
    name = "Type 1 error (in %)",
    labels = scales::label_percent(suffix = NULL),
    limits = c(0, 0.08)
  ) +
  scale_x_discrete(
    name = expression("Mean shift of " * u * " (unobserved heterogeneity)")
  ) +
  axis_title_margins +
  # Legend
  theme(legend.position = "top") +
  scale_color_manual(
    name = "Method",
    values = pal_oi(c("blue", "orange", "red", "purple"))
  ) +
  legend_color

p_t1e_x2 = dt_res |>
  # Prepare data set for plotting
  _[algo.id != 2] |>
  merge(dt_scenarios[b_trt == 0 & shift_u == 0], by = "scenario.id") |>
  merge(dt_algos, by = "algo.id") |>
  # Start plot
  ggplot(aes(x = factor(shift_x2), y = rejection, color = algo.name)) +
  geom_point(position = position_jitterdodge(), alpha = 0.8, size = 2.5) +
  # Reference lines for nominal level and power for two-sample t-test
  geom_hline(
    data = ref_lines[1:2],
    aes(yintercept = yintercept),
    linetype = "dashed",
    alpha = 0.5
  ) +
  # Axes
  scale_y_continuous(
    name = "Type 1 error (in %)",
    labels = scales::label_percent(suffix = NULL),
    limits = c(0, 0.08)
  ) +
  scale_x_discrete(
    name = expression("Mean shift of " * x[2] * " (observed prognostic covariate)")
  ) +
  axis_title_margins +
  # Legend
  theme(legend.position = "top") +
  scale_color_manual(
    name = "Method",
    values = pal_oi(c("blue", "orange", "red", "purple"))
  ) +
  legend_color

p_t1e_compare_u_x2 = wrap_plots(
  p_t1e_u, p_t1e_x2,
  nrow = 1, axes = "collect_y", guides = "collect"
) &
  theme(legend.position = "top")

print(p_t1e_compare_u_x2)

ggsave(
  fs$path("simulation", "analysis", "plot_t1e_u_x2_slides", ext = "png"),
  plot = p_t1e_compare_u_x2,
  dpi = 600, width = 10, height = 5
)



## Impact of historical data sample size (MSE) ----

p_mse_nhist = dt_res |>
  # Prepare data set for plotting
  _[algo.id %in% c(1, 3:5)] |>
  merge(dt_scenarios[shift_u == 0], by = "scenario.id") |>
  merge(dt_algos, by = "algo.id") |>
  # Start plot
  ggplot(aes(x = factor(n_hist), y = mse, color = algo.name)) +
  geom_point(position = position_jitterdodge(), alpha = 0.8, size = 2.5) +
  # Axes
  scale_y_continuous(
    name = "Mean squared error"
  ) +
  scale_x_discrete(
    name = expression(N[hist])
  ) +
  axis_title_margins +
  # Legend
  theme(legend.position = "top") +
  scale_color_manual(
    name = "Method",
    values = pal_oi(c("blue", "orange", "red", "purple"))
  ) +
  legend_color

print(p_mse_nhist)

ggsave(
  fs$path("simulation", "analysis", "plot_mse_nhist_slides", ext = "png"),
  plot = p_mse_nhist,
  dpi = 600, width = 10, height = 6
)
