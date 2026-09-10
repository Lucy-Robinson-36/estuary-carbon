####----------------------------CARBON-MAP 3--------------------------------
### Generating initial models for surface %OCl across the full AoI

###--------------------------SET UP WORKSPACE-------------------------------
## Clear the environment
rm(list=ls()) 

## Set working directory
setwd("path/Carbon_Map")

## Load packages
library(terra)
library(dplyr)
library(caret)
library(ggplot2)
library(sf)
library(CAST)
library(doParallel)
library(gbm)
library(glmnet)
library(tidyterra)
library(viridis)
library(tidyr)
library(systemfonts)
library(ragg)
library(tibble)

## Set ggplot defaults for uniform figures
theme_set(
  theme_bw(base_size = 12, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "right"
    )
)

## Define output directory
path <- "outputs/initial_models/"

## Read in files
surface_OCl <- read_sf("outputs/model_inputs/surface_OCl.shp")
predictors <- rast("outputs/model_inputs/predictors.tif")
AoI <- read_sf("outputs/model_inputs/full_AoI.shp")

## Define response type and run ID
run_id <- "fullsurfaceOCl"
resp_type <- "OCl"
resp_unit <- "% of TOC"


###-------------------------DATA PREPARATION--------------------------------
# Advise an initial separate investigation of response and predictor variables to ensure
# the model methodology used here is appropriate
## Log transformation of response
surface_OCl <- surface_OCl %>%
  mutate(log_OCl = log(mean_perc_OCl))

## View histograms
hist(surface_OCl$mean_perc_OCl, breaks = 40, main = "", xlab = paste0(resp_type, " (", resp_unit, ")"))
hist(surface_OCl$log_OCl, breaks = 40, main = "", xlab = paste0("log ",resp_type))

## Create test-train split
surface_OCl_train <- surface_OCl %>%
  filter(!Site %in% c("Site1", "Site2"))

surface_OCl_test <- surface_OCl %>%
  filter(Site %in% c("Site1", "Site2"))

## Extract predictor values
surface_OCl_train_sf <- extract(predictors, surface_OCl_train, bind = TRUE)
surface_OCl_train_sf <- st_as_sf(surface_OCl_train_sf)
surface_OCl_train <- st_drop_geometry(surface_OCl_train_sf)

surface_OCl_test_sf <- extract(predictors, surface_OCl_test, bind = TRUE)
surface_OCl_test_sf <- st_as_sf(surface_OCl_test_sf)
surface_OCl_test <- st_drop_geometry(surface_OCl_test_sf)

## Predictors
predictor_names <- names(predictors)


###--------------------------CREATE SITE FOLDS-------------------------------
## Define site folds (one fold per site) and set up cross-validation
site_folds <- CreateSpacetimeFolds(surface_OCl_train_sf,
                                   spacevar = "Site",
                                   k = length(unique(surface_OCl_train_sf$Site)))

sapply(site_folds$indexOut, function(idx) unique(surface_OCl_train_sf$Site[idx]))

ctrl <- trainControl(method = "cv",
                     index = site_folds$index,
                     indexOut = site_folds$indexOut,
                     savePredictions = "final")


###--------------------------FFS MODEL TRAINING------------------------------
## Run forward feature selection for each of the model methods (ran parallel for speed)
models_to_test <- c("rf", "cubist", "gbm", "glmnet")
model_list <- list()

nCores <- detectCores()
cl <- makePSOCKcluster(nCores - 1)
registerDoParallel(cl)

for (m in models_to_test) {
  set.seed(123)
  t0 <- Sys.time()
  model_list[[m]] <- ffs(surface_OCl_train[, predictor_names],
                         surface_OCl_train$log_OCl,
                         method = m,
                         trControl = ctrl,
                         tuneLength = 4)
  cat("Finished ffs:", m, "-- took",
      round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
  cat("  selected vars:", paste(model_list[[m]]$selectedvars, collapse = ", "), "\n")
}

stopCluster(cl)
registerDoSEQ()

## Extract selected variables
selected_vars_df <- data.frame(
  model = rep(models_to_test, sapply(models_to_test, function(m) length(model_list[[m]]$selectedvars))),
  variable = unlist(lapply(models_to_test, function(m) model_list[[m]]$selectedvars))
)

print(selected_vars_df)

## Get best tune per model
get_best_pred <- function(model) {
  merge(model$pred, model$bestTune)
}

best_preds <- lapply(model_list, get_best_pred)

## Duan's smear estimate for back-transform
get_smear <- function(pred_df) {
  resid_log <- pred_df$obs - pred_df$pred
  mean(exp(resid_log))
}

smear_factors <- sapply(best_preds, get_smear)
print(smear_factors)

## Performance metrics on original scale
get_orig_scale_metrics <- function(pred_df, smear) {
  obs_orig  <- exp(pred_df$obs)
  pred_orig <- exp(pred_df$pred) * smear
  c(RMSE = sqrt(mean((obs_orig - pred_orig)^2)),
    MAE  = mean(abs(obs_orig - pred_orig)),
    R2   = cor(obs_orig, pred_orig)^2)
}

comparison <- t(mapply(get_orig_scale_metrics, best_preds, smear_factors))
print(comparison)


###-------------------------ENSEMBLE MODEL--------------------------------
## Weight ensemble by inverse of RMSE (higher RMSE = lower weighting)
ensemble_weights <- (1 / comparison[, "RMSE"]) / sum(1 / comparison[, "RMSE"])
cat("\nEnsemble weights (inverse-RMSE):\n")
print(round(ensemble_weights, 3))

## Wide table of predictions from each model
obs_ref <- best_preds[[1]][, c("rowIndex", "obs")]
pred_wide <- obs_ref

for (m in models_to_test) {
  df <- best_preds[[m]][, c("rowIndex", "pred")]
  names(df)[2] <- paste0("pred_", m)
  pred_wide <- merge(pred_wide, df, by = "rowIndex")
}

## Back transform all predictions using smear factors
for (m in models_to_test) {
  pred_wide[[paste0("orig_", m)]] <- exp(pred_wide[[paste0("pred_", m)]]) * smear_factors[[m]]
}

## Build ensemble predictions as simple and weighted averages
pred_wide$ensemble_simple_avg <- rowMeans(pred_wide[, paste0("orig_", models_to_test)])
pred_wide$ensemble_weighted   <- as.matrix(pred_wide[, paste0("orig_", models_to_test)]) %*%
  ensemble_weights[models_to_test]

obs_orig <- exp(pred_wide$obs)

## Compute prediction metrics on original scale
get_ensemble_metrics <- function(pred_orig) {
  c(RMSE = sqrt(mean((obs_orig - pred_orig)^2)),
    MAE  = mean(abs(obs_orig - pred_orig)),
    R2   = cor(obs_orig, pred_orig)^2)
}

## Combine with comparison from individual models
ensemble_cv_comparison <- rbind(
  comparison,
  simple_avg = get_ensemble_metrics(pred_wide$ensemble_simple_avg),
  weighted_avg = get_ensemble_metrics(pred_wide$ensemble_weighted)
)

cat("\nCV comparison including ensemble (original scale):\n")
print(ensemble_cv_comparison)


###------------------------MODEL PREDICTIONS------------------------------
## Define function to generate model predictions
predict_positive_raster <- function(model, raster_stack, smear) {
  vars_needed <- model$selectedvars
  raster_subset <- raster_stack[[vars_needed]]
  pred_log <- predict(raster_subset, model, na.rm = TRUE)
  exp(pred_log) * smear
}

## Apply function to the predictor set for each of the models
pred_rasters <- list()
for (m in models_to_test) {
  pred_rasters[[m]] <- predict_positive_raster(model_list[[m]],
                                               predictors,
                                               smear_factors[[m]])
  cat("Finished predicting:", m, "\n")
}

pred_stack <- rast(pred_rasters)
names(pred_stack) <- models_to_test

## Generate prediction maps
pred_maps <- ggplot() +
  geom_spatraster(data = pred_stack) +
  facet_wrap(~ lyr) +
  scale_fill_viridis_c(name = paste0(resp_type, "\n(", resp_unit, ")"), na.value = "transparent") +
  labs(x = NULL, y = NULL)

## Weighted ensemble predictions
ensemble_raster <- sum(pred_stack * ensemble_weights[models_to_test])

names(ensemble_raster) <- "ensemble_weighted"

ensemble_map <- ggplot() +
  geom_spatraster(data = ensemble_raster) +
  scale_fill_viridis_c(name = paste0(resp_type, "\n(", resp_unit, ")"), na.value = "transparent") +
  labs(title = "Ensemble (weighted average) prediction", x = NULL, y = NULL)


###-----------------------AREA OF APPLICABILITY--------------------------
## Use aoa function to determine area of applicability for each model
aoa_list <- list()
for (m in models_to_test) {
  imp_check <- tryCatch(varImp(model_list[[m]])$importance,
                        error = function(e) NA)
  use_weight <- !anyNA(imp_check) && !any(is.nan(unlist(imp_check)))
  
  aoa_list[[m]] <- aoa(newdata = predictors, model = model_list[[m]],
                       useWeight = use_weight)
  cat("Finished AOA:", m, "-- weighted:", use_weight, "\n")
}

## Generate AOA maps for each of the models
AOA_stack <- rast(lapply(aoa_list, function(x) x$AOA))
names(AOA_stack) <- models_to_test

AOA_stack_f <- AOA_stack %>%
  mutate(across(everything(), ~ factor(.x, levels = c(0, 1),
                                       labels = c("Outside AOA", "Inside AOA"))))

aoa_maps <- ggplot() +
  geom_spatraster(data = AOA_stack_f) +
  facet_wrap(~ lyr) +
  scale_fill_manual(values = c("Outside AOA" = "grey85", "Inside AOA" = "steelblue"),
                    name = NULL, na.value = "transparent") +
  labs(title = "Area of Applicability", x = NULL, y = NULL) +
  theme_bw(base_size = 12)

## AoA summary for each model (area in km2 and percent coverage of AoI)
aoa_area_summary <- do.call(rbind, lapply(models_to_test, function(m) {
  area_df <- expanse(AOA_stack[[m]], byValue = TRUE, unit = "km")
  inside  <- area_df$area[area_df$value == 1]
  total   <- sum(area_df$area)
  data.frame(model = m,
             area_inside_AOA_km2 = ifelse(length(inside) == 0, 0, inside),
             pct_of_AoI = ifelse(length(inside) == 0, 0, 100 * inside / total))
}))

cat("\nArea of Applicability summary:\n")
print(aoa_area_summary) 

## Compute maps for dissimilarity index
DI_stack <- rast(lapply(aoa_list, function(x) x$DI))
names(DI_stack) <- models_to_test

di_maps <- ggplot() +
  geom_spatraster(data = DI_stack) +
  facet_wrap(~ lyr) +
  scale_fill_viridis_c(name = "DI", option = "viridis", na.value = "transparent") +
  labs(title = "Dissimilarity Index", x = NULL, y = NULL)

## Mask model predictions to the area of applicability
pred_masked_list <- lapply(models_to_test, function(m) {
  mask(pred_rasters[[m]], aoa_list[[m]]$AOA, maskvalues = 0)  # AOA == 0 = outside AOA
})
names(pred_masked_list) <- models_to_test

pred_masked_stack <- rast(pred_masked_list)
names(pred_masked_stack) <- models_to_test

masked_maps <- ggplot() +
  geom_spatraster(data = pred_masked_stack) +
  facet_wrap(~ lyr) +
  scale_fill_viridis_c(name = paste0(resp_type, "\n(", resp_unit, ")"), na.value = "transparent") +
  labs(title = "Prediction (masked to AOA)", x = NULL, y = NULL)


###-----------------------EVALUATE TEST SITES----------------------------
## Generate model predictions for the test dataset and get metrics on original scale
holdout_predictors <- surface_OCl_test[, predictor_names]

get_holdout_metrics <- function(model, smear) {
  vars_needed <- model$selectedvars
  pred_log  <- predict(model, newdata = holdout_predictors[, vars_needed])
  pred_orig <- exp(pred_log) * smear
  obs_orig  <- surface_OCl_test[["mean_perc_OCl"]]
  c(RMSE = sqrt(mean((obs_orig - pred_orig)^2)),
    MAE  = mean(abs(obs_orig - pred_orig)),
    R2   = cor(obs_orig, pred_orig)^2)
}

holdout_comparison <- t(mapply(get_holdout_metrics, model_list, smear_factors))

holdout_preds_by_model <- sapply(models_to_test, function(m) {
  vars_needed <- model_list[[m]]$selectedvars
  pred_log <- predict(model_list[[m]], newdata = holdout_predictors[, vars_needed])
  as.numeric(exp(pred_log) * smear_factors[[m]])
})

## Predictions and metrics for weighted ensemble model
holdout_ensemble_pred <- as.numeric(holdout_preds_by_model %*% ensemble_weights[models_to_test])
holdout_obs <- as.numeric(surface_OCl_test[["mean_perc_OCl"]])

holdout_ensemble_metrics <- c(
  RMSE = sqrt(mean((holdout_obs - holdout_ensemble_pred)^2)),
  MAE  = mean(abs(holdout_obs - holdout_ensemble_pred)),
  R2   = cor(holdout_obs, holdout_ensemble_pred)^2
)

## Add ensemble metrics to comparison
holdout_comparison <- rbind(holdout_comparison, ensemble_weighted = holdout_ensemble_metrics)

cat("\nIndependent holdout performance (", nrow(surface_OCl_test),
    "points from AG and CS", "):\n")
print(holdout_comparison)

## Check whether the holdout sites fall inside each model's AOA -
# if they fall outside, weaker holdout performance is expected/explained
# rather than a sign the model is simply bad.
holdout_aoa_check <- function(aoa_result, holdout_sf_pts) {
  terra::extract(aoa_result$AOA, vect(holdout_sf_pts))
}

for (m in models_to_test) {
  aoa_vals <- holdout_aoa_check(aoa_list[[m]], surfaceOC_test_sf)
  cat(m, "-- holdout points inside AOA:",
      sum(aoa_vals[, 2] == 1, na.rm = TRUE), "of", nrow(surface_OCl_test_sf), "\n")
}

## Plot observed vs. predicted for the holdout set, per model
holdout_plot_df <- do.call(rbind, lapply(models_to_test, function(m) {
  data.frame(
    model    = m,
    Observed = surface_OCl_test[["mean_perc_OCl"]],
    Predicted = holdout_preds_by_model[, m]
  )
}))

holdout_plots <- ggplot(holdout_plot_df, aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  facet_wrap(~ model) +
  coord_equal() +
  labs(x = paste0("Observed ", resp_type, " (", resp_unit, ")"),
       y = paste0("Predicted ", resp_type, " (", resp_unit, ")"))

# 
###----------------------------WRITE OUTPUTS-------------------------------
## save plots
plots_to_save <- list(
  pred_maps = pred_maps,
  ensemble_map = ensemble_map,
  aoa_maps = aoa_maps,
  di_maps = di_maps,
  masked_maps = masked_maps,
  holdout_plots = holdout_plots
)

for (p in names(plots_to_save)) {
  ggsave(
    filename = file.path(paste0(path, run_id, "_", p, ".png")),
    plot = plots_to_save[[p]],
    device = ragg::agg_png,
    width = 8, height = 6, units = "in",
    dpi = 300
  )
}

## save metrics
# function to pivot into long format
tidy_matrix_metrics <- function(mat, source_label) {
  as.data.frame(mat) %>%
    rownames_to_column("model") %>%
    pivot_longer(cols = -model, names_to = "metric", values_to = "value") %>%
    mutate(source = source_label, run_id = run_id, .before = 1)
}

# prepare results dfs
cv_results <- tidy_matrix_metrics(comparison, "cv")
ensemble_results <- tidy_matrix_metrics(ensemble_cv_comparison, "cv_ensemble")
holdout_results <- tidy_matrix_metrics(holdout_comparison, "holdout")

aoa_results <- aoa_area_summary %>%
  pivot_longer(cols = c(area_inside_AOA_km2, pct_of_AoI),
               names_to = "metric", values_to = "value") %>%
  rename(model = model) %>%
  mutate(source = "aoa", run_id = run_id, .before = 1)

smear_results <- data.frame(
  run_id = run_id, source = "smear", model = names(smear_factors),
  metric = "smear_factor", value = as.numeric(smear_factors)
)

ensemble_weight_results <- data.frame(
  run_id = run_id, source = "ensemble_weight",
  model = names(ensemble_weights), metric = "weight",
  value = as.numeric(ensemble_weights)
)

# combine all results and write csv
all_results <- bind_rows(cv_results, ensemble_results, holdout_results,
                         aoa_results, smear_results, ensemble_weight_results) %>%
  mutate(predictor_set = paste(predictor_names, collapse = ", "),
         timestamp = Sys.time())

write.csv(all_results, file.path(paste0(path, run_id, "_results.csv")), row.names = FALSE)
write.csv(selected_vars_df, file.path(paste0(path, run_id, "sel_preds.csv")), row.names = FALSE)



