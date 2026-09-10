####----------------------------CARBON-MAP 2--------------------------------
### Preparing data for models

###--------------------------SET UP WORKSPACE-------------------------------
## Clear the environment
rm(list=ls()) 

## Set working directory
setwd("path/Carbon_Map")

## Load packages
library(sf)
library(dplyr)
library(terra)
library(stars)
library(stringr)


###----------------------------READ IN DATA---------------------------------
## Geochemical data
surface_TOC_data <- read.csv("outputs/core_results/Surface_TOC.csv")
depth_TOC_data <- read.csv("outputs/core_results/Depth_TOC.csv")
surface_OCl_data <- read.csv("outputs/core_results/Surface_OCl.csv")
depth_OCl_data <- read.csv("outputs/core_results/Depth_OCl.csv")

## Predictor rasters
NDVI <- rast("predictors/NDVI.tif")
elevation <- rast("predictors/Elevation.tif")
riverinf <- rast("predictors/Riverinf.tif")
slope <- rast("predictors/Slope.tif")
fetch <- rast("predictors/Fetch.tif")

## Coordinates
coords <- read.csv("samples/Coords.csv")

## Area of interest polygons
full_AoI <- vect("predictors/Full_AoI.shp")
marsh_AoI <- vect("predictors/Marsh_AoI.shp")


###-------------------------PREPARE RASTERS-------------------------------
## Check extent, projection and resolution of rasters
crs(NDVI); crs(elevation); crs(riverinf); crs(slope); crs(fetch)
res(NDVI); res(elevation); res(riverinf); res(slope); res(fetch)
ext(NDVI); ext(elevation); ext(riverinf); ext(slope); ext(fetch)

## Align rasters
# create a template cropped to AoI
template <- rast(elevation)
template <- crop(template, AoI)
# resample all layers using template
NDVI      <- resample(NDVI, template)
riverinf  <- resample(riverinf, template)
slope     <- resample(slope, template)
elevation <- resample(elevation, template)
fetch     <- resample(fetch, template)

## Check extent, projection and resolution of rasters
crs(NDVI); crs(elevation); crs(riverinf); crs(slope); crs(fetch)
res(NDVI); res(elevation); res(riverinf); res(slope); res(fetch)
ext(NDVI); ext(elevation); ext(riverinf); ext(slope); ext(fetch)

## Create raster stack
predictors <- c(NDVI, elevation, riverinf, slope, fetch)
names(predictors) <- c("NDVI", "elevation", "riverinf", "slope", "fetch")


###----------------------CROP RASTERS TO MARSH----------------------------
## Crop each raster to marsh AoI
NDVI      <- mask(NDVI, marsh_AoI)
riverinf  <- mask(riverinf, marsh_AoI)
slope     <- mask(slope, marsh_AoI)
elevation <- mask(elevation, marsh_AoI)
fetch     <- mask(fetch, marsh_AoI)

## Check extent, projection and resolution of rasters
crs(NDVI); crs(elevation); crs(riverinf); crs(slope); crs(fetch)
res(NDVI); res(elevation); res(riverinf); res(slope); res(fetch)
ext(NDVI); ext(elevation); ext(riverinf); ext(slope); ext(fetch)

## Create raster stack
marsh_predictors <- c(NDVI, elevation, riverinf, slope, fetch)
names(marsh_predictors) <- c("NDVI", "elevation", "riverinf", "slope", "fetch")


###-----------------------------FUNCTIONS----------------------------------
## Convert units to MgC/ha
convert_units <- function(data_df) {
  data_df$mean_MgC.ha <- data_df$mean_gC.cm2*100
  data_df
}

## Combine response and coordinates to make simple feature layer
make_sf_points <- function(stock_df, coords, from_crs = 4326, to_crs = 27700) {
  stock_df <- left_join(stock_df, coords, by = "Core_ID", keep = FALSE)
  pts <- st_as_sf(stock_df, coords = c("Longitude", "Latitude"))
  pts <- st_set_crs(pts, from_crs)
  st_transform(pts, to_crs)
}

## Find grid cell centroid means
grid_centroid_means <- function(fishnet, sf_points, value_col) {
  joined <- st_join(fishnet, sf_points)
  joined <- na.omit(joined)
  summarized <- joined %>%
    group_by(Core_ID) %>%
    summarize(!!value_col := mean(.data[[value_col]], na.rm = TRUE))
  st_point_on_surface(summarized)
}

## Fix site codes
fix_site_codes <- function(df, threshold = 35) {
  df <- separate(df, Core_ID, into = c("Site", "Transect", "Core"), sep = "_")
}


###-----------------------PREPARE RESPONSE VARS-----------------------------
## Convert units of TOC dataframe
surface_TOC_data <- convert_units(surface_TOC_data)
depth_TOC_data <- convert_units(depth_TOC_data)

## Create sf dataframes for response variables
surface_TOC <- make_sf_points(surface_TOC_data)
depth_TOC <- make_sf_points(depth_TOC_data)
surface_OCl <- make_sf_points(surface_OCl_data)
depth_OCl <- make_sf_points(depth_OCl_data)

## Check projections
crs(surface_TOC); crs(depth_TOC); crs(surface_OCl); crs(depth_OCl)


###----------------------GRID CELL CENTROID MEANS----------------------------
## Create fishnet
min_ext <- sum(predictors)
min_ext[min_ext == 0] <- NA
fishnet <- st_as_sf(stars::st_as_stars(min_ext), as_points = FALSE, merge = FALSE)
fishnet$ID <- 1:nrow(fishnet)


## Find means
surface_TOC_grid   <- grid_centroid_means(fishnet, surface_TOC, "mean_MgC.ha")
depth_TOC_grid  <- grid_centroid_means(fishnet, depth_TOC, "mean_MgC.ha")
surface_OCl_grid <- grid_centroid_means(fishnet, surface_OCl, "mean_perc_OCl")
depth_OCl_grid <- grid_centroid_means(fishnet, depth_OCl, "mean_perc_OCl")


###-----------------------------SAVE FILES----------------------------------
## Fix site names
surface_TOC_grid  <- fix_site_codes(surface_TOC_grid)
depth_TOC_grid <- fix_site_codes(depth_TOC_grid)
surface_OCl_grid <- fix_site_codes(surface_OCl_grid)
depth_OCl_grid <- fix_site_codes(depth_OCl_grid)

## Save files
writeRaster(predictors, "outputs/model_inputs/predictors.tif", overwrite = TRUE)
writeRaster(marsh_predictors, "outputs/model_inputs/marsh_predictors.tif", overwrite = TRUE)
writeVector(full_AoI, "outputs/model_inputs/full_AoI.shp", overwrite = TRUE)
writeVector(marsh_AoI, "outputs/model_inputs/marsh_AoI.shp", overwrite = TRUE)

write_sf(surface_TOC_grid, dsn = "outputs/model_inputs", layer = "surface_TOC", driver = "ESRI Shapefile")
write_sf(depth_TOC_grid, dsn = "outputs/model_inputs", layer = "depth_TOC", driver = "ESRI Shapefile")
write_sf(surface_OCl_grid, dsn = "outputs/model_inputs", layer = "surface_OCl", driver = "ESRI Shapefile")
write_sf(depth_OCl_grid, dsn = "outputs/model_inputs", layer = "depth_OCl", driver = "ESRI Shapefile")

