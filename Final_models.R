####----------------------------CARBON-MAP 5--------------------------------
### Generating final models using identified method

###--------------------------SET UP WORKSPACE-------------------------------
## Clear the environment
rm(list=ls()) 

## Set working directory
setwd("path/Carbon_Map")

## Read in packages
library(dplyr)
library(sf)
library(terra)
library(caret)
library(CAST)
library(doParallel)
library(gbm)
library(tibble)

## Set ggplot defaults
theme_set(
  theme_bw(base_size = 12, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "right"
    )
)

## Set output data path
path <- "outputs/final_models/"

## Read in data
# Response variables
full_surface_TOC <- read_sf("outputs/model_inputs/surface_TOC.shp")
full_surface_OCl <- read_sf("outputs/model_inputs/surface_OCl.shp")
marsh_depth_TOC <- read_sf("outputs/model_inputs/humified_TOC.shp")
marsh_depth_OCl <- read_sf("outputs/model_inputs/depth_OCl.shp")

# Predictor variables
full_predictors <- rast("outputs/model_inputs/predictors.tif")
marsh_predictors <- rast("outputs/model_inputs/marsh_predictors.tif")

# Area of interest
full_AoI <- read_sf("outputs/model_inputs/full_AoI.shp")
marsh_AoI <- read_sf("outputs/model_inputs/marsh_AoI.shp")


###-------------------------PREPARE DATA----------------------------------
## Prepare marsh surface data
marsh_surface_TOC <- full_surface_TOC %>%
  semi_join(st_drop_geometry(marsh_depth_TOC), by = c("Site", "Transect", "Core"))

marsh_surface_OCl <- full_surface_OCl %>%
  semi_join(st_drop_geometry(marsh_depth_OCl), by = c("Site", "Transect", "Core"))

## Rename response columns to "resp"
full_surface_TOC <- full_surface_TOC %>%
  dplyr::rename(resp = mean_MgC.ha)

marsh_surface_TOC <- marsh_surface_TOC %>%
  dplyr::rename(resp = mean_MgC.ha)

marsh_depth_TOC <- marsh_depth_TOC %>%
  dplyr::rename(resp = mean_MgC.ha)

full_surface_OCl <- full_surface_OCl %>%
  dplyr::rename(resp = mean_perc_OCl)

marsh_surface_OCl <- marsh_surface_OCl %>%
  dplyr::rename(resp = mean_perc_OCl)

marsh_depth_OCl <- marsh_depth_OCl %>%
  dplyr::rename(resp = mean_perc_OCl)

## Configure modelling scenarios
scenario_config <- list(
  fsTOC = list(response = full_surface_TOC, predictors = full_predictors, aoi = full_AoI),
  fsOCl = list(response = full_surface_OCl, predictors = full_predictors, aoi = full_AoI),
  msTOC = list(response = marsh_surface_TOC, predictors = marsh_predictors, aoi = marsh_AoI),
  msOCl = list(response = marsh_surface_OCl, predictors = marsh_predictors, aoi = marsh_AoI),
  mdTOC = list(response = marsh_depth_TOC, predictors = marsh_predictors, aoi = marsh_AoI),
  mdOCl = list(response = marsh_depth_OCl, predictors = marsh_predictors, aoi = marsh_AoI)
)


###---------------------------MODEL FUNCTION--------------------------------
run_gbm_scenario <- function(scenario_name, response_df, predictor_stack, aoi) {
  
  cat("\n=================================\n")
  cat("Running scenario:", scenario_name, "\n")
  cat("=================================\n")
  
  predictor_names <- names(predictor_stack)
  
  ## Prepare full dataset (train + holdout combined) 
  data_sf <- response_df %>%
    mutate(log_resp = log(resp))
  
  data_sf <- terra::extract(predictor_stack, data_sf, bind = TRUE)
  data_sf <- sf::st_as_sf(data_sf)
  data    <- st_drop_geometry(data_sf)
  
  ## Site-based CV folds (all sites now, since holdout is included)
  site_folds <- CreateSpacetimeFolds(data_sf,
                                     spacevar = "Site",
                                     k = length(unique(data_sf$Site)))
  
  ctrl <- trainControl(method = "cv",
                       index = site_folds$index,
                       indexOut = site_folds$indexOut,
                       savePredictions = "final")
  
  ## Fit gbm via ffs 
  nCores <- detectCores()
  cl <- makePSOCKcluster(nCores - 1)
  registerDoParallel(cl)
  
  set.seed(123)
  t0 <- Sys.time()
  gbm_model <- ffs(data[, predictor_names],
                   data$log_resp,
                   method = "gbm",
                   trControl = ctrl,
                   tuneLength = 4,
                   verbose = FALSE)
  
  stopCluster(cl)
  registerDoSEQ()
  
  cat("  Finished ffs -- took",
      round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
  cat("  Selected vars:", paste(gbm_model$selectedvars, collapse = ", "), "\n")
  
  ## Duan's smear estimate 
  best_pred <- merge(gbm_model$pred, gbm_model$bestTune)
  resid_log <- best_pred$obs - best_pred$pred
  smear <- mean(exp(resid_log))
  
  ## CV performance metrics (original scale)
  obs_orig  <- exp(best_pred$obs)
  pred_orig <- exp(best_pred$pred) * smear
  
  cv_metrics <- data.frame(
    scenario = scenario_name,
    metric   = c("RMSE", "MAE", "R2"),
    value    = c(sqrt(mean((obs_orig - pred_orig)^2)),
                 mean(abs(obs_orig - pred_orig)),
                 cor(obs_orig, pred_orig)^2)
  )
  
  ## Predict to raster, masked to AoI 
  vars_needed   <- gbm_model$selectedvars
  raster_subset <- predictor_stack[[vars_needed]]
  
  pred_log    <- terra::predict(raster_subset, gbm_model, na.rm = TRUE)
  pred_raster <- exp(pred_log) * smear
  pred_raster <- terra::mask(pred_raster, aoi)
  pred_raster <- terra::crop(pred_raster, aoi)
  names(pred_raster) <- paste0(scenario_name, "_pred")
  
  ## AOA / DI 
  imp_check  <- tryCatch(varImp(gbm_model)$importance, error = function(e) NA)
  use_weight <- !anyNA(imp_check) && !any(is.nan(unlist(imp_check)))
  
  aoa_result <- aoa(newdata = predictor_stack, model = gbm_model, useWeight = use_weight)
  
  aoa_masked <- terra::mask(terra::crop(aoa_result$AOA, aoi), aoi)
  di_masked  <- terra::mask(terra::crop(aoa_result$DI, aoi), aoi)
  
  area_df <- expanse(aoa_masked, byValue = TRUE, unit = "km")
  inside  <- area_df$area[area_df$value == 1]
  total   <- sum(area_df$area)
  
  aoa_metrics <- data.frame(
    scenario = scenario_name,
    metric   = c("area_inside_AOA_km2", "pct_of_AoI"),
    value    = c(ifelse(length(inside) == 0, 0, inside),
                 ifelse(length(inside) == 0, 0, 100 * inside / total))
  )
  
  ## Prediction masked to AOA
  pred_masked <- terra::mask(pred_raster, aoa_masked, maskvalues = 0)
  names(pred_masked) <- paste0(scenario_name, "_pred_masked")
  
  ##Variable importance 
  varimp_df <- tryCatch({
    vi <- varImp(gbm_model, scale = TRUE)$importance
    vi_df <- vi %>%
      rownames_to_column("variable") %>%
      rename(importance = Overall)
    
    # reindex against this scenario's full predictor set so unselected vars = 0
    data.frame(variable = predictor_names) %>%
      left_join(vi_df, by = "variable") %>%
      mutate(importance = ifelse(is.na(importance), 0, importance))
  }, error = function(e) {
    cat("  No varImp available for:", scenario_name, "\n")
    data.frame(variable = predictor_names, importance = NA)
  })
  
  varimp_df <- varimp_df %>%
    mutate(scenario = scenario_name, model = "gbm")
  
  write.csv(varimp_df,
            file.path(path, paste0(scenario_name, "_varimp.csv")),
            row.names = FALSE)
  
  ## Save rasters 
  terra::writeRaster(pred_raster, file.path(path, paste0(scenario_name, "_pred.tif")),
                     overwrite = TRUE)
  terra::writeRaster(pred_masked, file.path(path, paste0(scenario_name, "_pred_masked.tif")),
                     overwrite = TRUE)
  terra::writeRaster(aoa_masked, file.path(path, paste0(scenario_name, "_AOA.tif")),
                     overwrite = TRUE)
  terra::writeRaster(di_masked, file.path(path, paste0(scenario_name, "_DI.tif")),
                     overwrite = TRUE)
  
  ## Combine numeric results
  scenario_results <- bind_rows(
    cv_metrics %>% mutate(source = "cv"),
    aoa_metrics %>% mutate(source = "aoa"),
    data.frame(scenario = scenario_name, metric = "smear_factor",
               value = smear, source = "smear")
  ) %>%
    mutate(model = "gbm",
           selected_vars = paste(vars_needed, collapse = ", "),
           timestamp = Sys.time())
  
  write.csv(scenario_results,
            file.path(path, paste0(scenario_name, "_gbm_results.csv")),
            row.names = FALSE)

  cat("  Done:", scenario_name, "\n")
  
  return(list(model = gbm_model, results = scenario_results,
              pred_raster = pred_raster, pred_masked = pred_masked,
              aoa = aoa_masked, di = di_masked, varimp = varimp_df))
}


###-----------------------------RUN MODELS----------------------------------
all_scenario_outputs <- list()

for (scen in names(scenario_config)) {
  cfg <- scenario_config[[scen]]
  all_scenario_outputs[[scen]] <- run_gbm_scenario(
    scenario_name   = scen,
    response_df     = cfg$response,
    predictor_stack = cfg$predictors,
    aoi             = cfg$aoi
  )
}


###-----------------------VARIABLE IMPORTANCE SUMMARY------------------------
all_varimp <- bind_rows(lapply(all_scenario_outputs, function(x) x$varimp))

write.csv(all_varimp, file.path(path, "all_scenarios_varimp.csv"), row.names = FALSE)
