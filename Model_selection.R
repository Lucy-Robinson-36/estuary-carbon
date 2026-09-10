####----------------------------CARBON-MAP 4--------------------------------
### Comparison of performance metrics of initial models

###--------------------------SET UP WORKSPACE-------------------------------
## Clear the environment
rm(list=ls()) 

## Set working directory
setwd("path/Carbon_Map")

## Load packages
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(ggokabeito)
library(ggtext)

## Set ggplot defaults
theme_set(
  theme_bw(base_size = 12, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "right"
    )
)

## Define file path
path <- "outputs/initial_models/"

## Read in files
files_list <- list.files(path, pattern = "_results\\.csv$", full.names = TRUE)
combined_metrics <- map_dfr(files_list, read.csv)


###-------------------------MODEL ASSESSMENT------------------------------
# using holdout error metrics as a primary measure of performance

## Summarise error metrics for each model
holdout_summary <- combined_metrics %>%
  filter(source == "holdout", metric %in% c("RMSE", "MAE", "R2")) %>%
  group_by(model, metric) %>%
  summarise(
    mean_value = mean(value),
    sd_value   = sd(value),
    min_value  = min(value),
    max_value  = max(value),
    n_runs     = n(),
    .groups = "drop"
  ) %>%
  arrange(metric, mean_value)

holdout_summary <- holdout_summary %>%
  group_by(metric) %>%
  mutate(is_best = case_when(
    metric == "R2" & mean_value == max(mean_value) ~ TRUE,
    metric %in% c("RMSE", "MAE") & mean_value == min(mean_value) ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  ungroup()

holdout_summary <- holdout_summary %>%
  mutate(model = recode(model, "ensemble_weighted" = "ensemble"))

## Error metrics for each model run per scenario
holdout_runs <- combined_metrics %>%
  filter(source == "holdout", metric %in% c("RMSE", "MAE", "R2")) %>%
  mutate(model = recode(model, "ensemble_weighted" = "ensemble"))

holdout_runs$model <- factor(holdout_runs$model, levels = rev(c("rf", "cubist", "gbm", "glmnet", "ensemble")))

## Labels for scenarios
labels = c("fullsurfaceOCl" = "Full Surface %OC<sub>L</sub>",
           "marshsurfaceOCl" = "Marsh Surface %OC<sub>L</sub>",
           "marshhumifiedOCl" = "Marsh Humified %OC<sub>L</sub>",
           "fullsurfaceTOC" = "Full Surface TOC",
           "marshsurfaceTOC" = "Marsh Surface TOC",
           "marshhumifiedTOC" = "Marsh Humified TOC"
)

## Create and save scatter plot of holdout metrics with mean and sd
holdoutmetrics_plot <- ggplot(holdout_summary, aes(x = model, y = mean_value)) +
  geom_point(data = holdout_runs, aes(x = model, y = value, fill = run_id),
             shape = 21, alpha = 0.7, size = 2, colour = "gray40", stroke = 0.7,
             position = position_jitter(width = 0.1, height = 0)) +
  scale_fill_okabe_ito(labels = labels) +
  geom_pointrange(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
                  size = 0.8, fatten = 2) +
  facet_wrap(~ factor(metric, c("R2", "RMSE", "MAE")), scales = "free_x") +
  coord_flip() +
  labs(x = NULL, y = "Value (mean ± SD)", fill = "Scenario") +
  theme(legend.position = "bottom", legend.text = element_markdown(), panel.spacing.x = unit(1, "cm"))

holdoutmetrics_plot

ggsave(
  filename = file.path(paste0(path, "holdoutmetrics_plot.png")),
  plot = holdoutmetrics_plot,
  device = ragg::agg_png,
  width = 8, height = 6, units = "in",
  dpi = 300
)

## Compare area of applicability (as % of AoI) for each model type
combined_metrics %>%
  filter(source == "aoa", metric == "pct_of_AoI") %>%
  group_by(model) %>%
  summarise(mean_pct_AoI = mean(value), .groups = "drop")


###-----------------------SELECTED PREDICTORS-------------------------------
## Read in files
preds_list <- list.files(path, pattern = "sel_preds\\.csv$", full.names = TRUE)

selpreds <- preds_list %>%
  set_names(basename(preds_list)) %>%
  map_dfr(read.csv, .id = "file_name")

selpreds <- selpreds %>%
  mutate(run_id = sub("sel_preds\\.csv$", "", file_name)) %>%
  extract(run_id, into = c("Area", "Depth", "OC"),
          regex = "^(full|marsh)(surface|humified)(TOC|OCl)$") %>%
  select(!(file_name))

## Count predictor use
selpreds_count <- selpreds %>%
  count(Area, Depth, OC, variable)

selpreds_count <- selpreds_count %>%
  mutate(
    AreaDepth = factor(
      paste(Area, Depth, sep = "_"),
      levels = c(
        "full_surface",
        "marsh_surface",
        "marsh_humified"
      )
    )
  )

## Plot of predictor use frequency
ggplot(selpreds_count, aes(x = variable, y = n, fill = variable)) +
  geom_col() +
  facet_grid(rows = vars(OC), cols = vars(AreaDepth)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Count")

