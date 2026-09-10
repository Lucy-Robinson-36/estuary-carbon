####----------------------------CARBON-MAP 6--------------------------------
### Extracting predictions from the final models

###--------------------------SET UP WORKSPACE-------------------------------
## Clear the environment
rm(list=ls()) 

## Set working directory
setwd("path/Carbon_Map")

## Read in packages
library(dplyr)
library(sf)
library(terra)
library(stars)
library(ggplot2)
library(ggtext)
library(biscale)
library(cowplot)
library(purrr)
library(tidyterra)

## Read in files
# Model predictions
full_surface_OCl <- rast("outputs/final_models/fsOCl_pred.tif")
full_surface_TOC <- rast("outputs/final_models/fsTOC_pred.tif")

marsh_surface_OCl <- rast("outputs/final_models/msOCl_pred.tif")
marsh_surface_TOC <- rast("outputs/final_models/msTOC_pred.tif")

marsh_depth_OCl <- rast("outputs/final_models/mdOCl_pred.tif")
marsh_depth_TOC <- rast("outputs/final_models/mdTOC_pred.tif")

# Polygons
saltmarsh <- vect("predictors/saltmarsh.shp")
tidal_flat <- vect("predictors/tidal_flat.shp")
ind_marshes <- vect("predictors/individual_marshes.shp")


###---------------------------SINGLE-POLYGON----------------------------------
## Create lists
rasters <- list(full_surface_OCl, full_surface_TOC, marsh_surface_OCl, 
                marsh_surface_TOC, marsh_depth_OCl, marsh_depth_TOC)
raster_names <- c("fs_OCl", "fs_TOC", "ms_OCl", "ms_TOC", "md_OCl", "md_TOC")

vectors <- list(saltmarsh, tidal_flat)  
vector_names <- c("saltmarsh", "tidal flat")

TOC_rasters <- list(full_surfaceTOC, marsh_surfaceTOC, marsh_depthTOC)
TOC_names <- c("fs_TOC", "ms_TOC", "md_TOC")

## Build results tables
# Averages
results <- data.frame()

for (i in seq_along(rasters)) {
  for (j in seq_along(vectors)) {
    vals <- terra::extract(rasters[[i]], vectors[[j]], ID = FALSE)[[1]]
    results <- rbind(results, data.frame(
      raster = raster_names[i],
      vector = vector_names[j],
      mean = mean(vals, na.rm = TRUE),
      sd = sd(vals, na.rm = TRUE)
    ))
  }
}

results

fullarea_averages <- results

# Total TOC
results <- data.frame()

cell_area_ha <- 0.01

for (i in seq_along(TOC_rasters)) {
  for (j in seq_along(vectors)) {
    vals <- terra::extract(TOC_rasters[[i]], vectors[[j]], ID = FALSE)[[1]]
    total_stock <- sum(vals * cell_area_ha, na.rm = TRUE)
    results <- rbind(results, data.frame(
      raster = TOC_names[i],
      vector = vector_names[j],
      total_Mg = total_stock
    ))
  }
}

results

fullarea_totals <- results



###-------------------------------MULTI-FEATURE-------------------------------
## Averages
marsh_averages <- data.frame()

for (i in seq_along(rasters)) {
  means <- terra::extract(rasters[[i]], ind_marshes, fun = mean, na.rm = TRUE, ID = FALSE)[[1]]
  sds   <- terra::extract(rasters[[i]], ind_marshes, fun = sd,   na.rm = TRUE, ID = FALSE)[[1]]
  
  marsh_averages <- rbind(marsh_averages, data.frame(
    raster  = raster_names[i],
    vector = "ind_marshes",
    feature = if (!is.null(ind_marshes$Name)) ind_marshes$Name else seq_along(means),
    mean    = means,
    sd      = sds
  ))
}


## Totals
marsh_totals <- data.frame()

for (i in seq_along(TOC_rasters)) {
  TOC_Mg <- TOC_rasters[[i]] * 0.01
  sums <- terra::extract(TOC_Mg, ind_marshes, fun = sum, na.rm = TRUE, ID = FALSE)[[1]]
  
  marsh_totals <- rbind(marsh_totals, data.frame(
    raster  = TOC_names[i],
    vector  = "ind_marshes",
    feature = if (!is.null(ind_marshes$Name)) ind_marshes$Name else seq_along(sums),
    total_Mg = sums
  ))
}


###-------------------------------FORMATTING---------------------------------
## Marsh totals
marsh_total_TOC <- marsh_totals %>%
  separate(raster, c("Model_Code", "TOC")) %>%
  rename(Marsh = feature,
         `TOC_stock  [Mg]` = total_Mg) %>%
  separate(Model_Code, c("Area", "Sed_layer"), sep=1) %>%
  select(Area, Sed_layer, Marsh, `TOC_stock  [Mg]`)


marsh_totalOC$Area <- recode_values(marsh_total_TOC$Area, 
                                    "f" ~ "Full", "m" ~ "Saltmarsh")
marsh_totalOC$Sed_layer <- recode_values(marsh_total_TOC$Sed_layer, 
                                     "s" ~ "Surface", "d" ~ "Depth")

## Marsh Averages
marsh_averages <- marsh_averages %>%
  separate(raster, c("Model_Code", "Variable")) %>%
  rename(Marsh = feature) %>%
  separate(Model_Code, c("Area", "Depth"), sep=1) %>%
  select(Area, Depth, Marsh, Variable, mean, sd)

marsh_averages$Area <- recode_values(marsh_averages$Area, 
                                     "f" ~ "Full", "m" ~ "Saltmarsh")
marsh_averages$Sed_layer <- recode_values(marsh_averages$Sed_layer, 
                                      "s" ~ "Surface", "d" ~ "Depth")  

marsh_averages <- pivot_wider(marsh_averages,
                              names_from = Variable,
                              values_from = c(mean, sd))

marsh_averages <- marsh_averages %>%
  rename(`TOC_avg  [Mg.ha]` = mean_TOC,
         `TOC_sd  [Mg.ha]` = sd_TOC,
         `OCl_avg  [%]` = mean_OCl,
         `OCl_sd  [%]` = sd_OCl)

## Individual marshes combined df
marsh_values <- merge(marsh_averages, marsh_total_TOC, by = c("Area", "Depth", "Marsh"))


## Combined area averages
region_averages <- fullarea_averages %>%
  separate(raster, c("Model_Code", "Variable")) %>%
  rename(Region = vector) %>%
  separate(Model_Code, c("Area", "Sed_layer"), sep=1) %>%
  select(Area, Sed_layer, Region, Variable, mean, sd)

region_averages$Area <- recode_values(region_averages$Area, 
                                      "f" ~ "Full", "m" ~ "Saltmarsh")
region_averages$Sed_layer <- recode_values(region_averages$Sed_layer, 
                                       "s" ~ "Surface", "d" ~ "Depth")

region_averages <- pivot_wider(region_averages,
                               names_from = Variable,
                               values_from = c(mean, sd))

region_averages <- region_averages %>%
  rename(`TOC_avg  [Mg.ha]` = mean_TOC,
         `TOC_sd  [Mg.ha]` = sd_TOC,
         `OCl_avg  [%]` = mean_OCl,
         `OCl_sd  [%]` = sd_OCl)


## Combined area totals
region_totals <- fullarea_totals %>%
  separate(raster, c("Model_Code", "Variable")) %>%
  rename(Region = vector,
         `TOC_stock  [Mg]` = total_Mg) %>%
  separate(Model_Code, c("Area", "Sed_layer"), sep=1) %>%
  select(Area, Sed_layer, Region, `TOC_stock  [Mg]`)

region_totals$Area <- recode_values(region_totals$Area, 
                                    "f" ~ "Full", "m" ~ "Saltmarsh")
region_totals$Sed_layer <- recode_values(region_totals$Sed_layer, 
                                     "s" ~ "Surface", "d" ~ "Depth")

## Region combined df
region_values <- merge(region_averages, region_totals, by = c("Area", "Sed_layer", "Region"))


###---------------------------MODEL COMPARISON-------------------------------
# Comparison of surface models for full AoI and marsh AoI
## Combined dataframe for averages
marsh_averages <- marsh_averages %>%
  rename(Region = Marsh)

combined_averages <- rbind(region_averages, marsh_averages)


## Surface marsh only
surface_averages <- combined_averages %>%
  filter(Sed_layer == "Surface",
         Region != ("tidal flat"))

surface_averages$Region <- factor(surface_averages$Region, 
                                  levels = c("Marsh1", "Marsh2", "Marsh3", "etc"))

## Plot
TOC <- ggplot(surface_averages, aes(x = Region, y = `TOC_avg  [Mg.ha]`, fill = Area)) +
  geom_col(position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#CC79A7", "#00D199")) +
  geom_errorbar(aes(ymin = `TOC_avg  [Mg.ha]` - `TOC_sd  [Mg.ha]`,
                    ymax = `TOC_avg  [Mg.ha]` + `TOC_sd  [Mg.ha]`),
                width = 0.2, position = position_dodge(width = 0.9)) +
  labs(y = "mean TOC Mg ha<sup>-1</sup>") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_markdown())

OCl <- ggplot(surface_averages, aes(x = Region, y = `OCl_avg  [%]`, fill = Area)) +
  geom_col(position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#CC79A7", "#00D199")) +
  geom_errorbar(aes(ymin = `OCl_avg  [%]` - `OCl_sd  [%]`,
                    ymax = `OCl_avg  [%]` + `OCl_sd  [%]`),
                width = 0.2, position = position_dodge(width = 0.9)) +
  labs(y = "mean %OC <sub>L</sub>") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_markdown())

combined <- plot_grid(TOC, OCl, ncol = 1)

ggsave(
  filename = "outputs/final_models/area_comparison.png",
  plot = combined,
  device = ragg::agg_png,
  width = 8, height = 8, units = "in",
  dpi = 300,
  bg = "white"
)

###----------------------------EXPORT TABLES--------------------------------
write.csv(marsh_values, "outputs/final_models/marsh_values.csv")
write.csv(region_values, "outputs/final_models/region_values.csv")


###-------------------------SELECTED VARIABLES-------------------------------
## Read in files
# list output metric files
path <- "outputs/final_models/"
results_list <- list.files(path, pattern = "_results\\.csv$", full.names = TRUE)

# read in and combine to one long df
combined_results <- map_dfr(results_list, read.csv)

# filter df for just one row of each scenario
selpreds <- combined_results %>%
  filter(metric == "RMSE") %>%
  select(scenario, selected_vars)


###-----------------------------BIVARIATE PLOT-------------------------------
## Define bilinear palette
bilinear_mix <- function(c00, c10, c01, c11, n = 3) {
  # c00 = low-low, c10 = high x/low y, c01 = low x/high y, c11 = high-high
  rgb00 <- grDevices::col2rgb(c00) / 255
  rgb10 <- grDevices::col2rgb(c10) / 255
  rgb01 <- grDevices::col2rgb(c01) / 255
  rgb11 <- grDevices::col2rgb(c11) / 255
  
  steps <- seq(0, 1, length.out = n)
  pal <- character(0)
  labels <- character(0)
  
  for (yi in seq_along(steps)) {
    for (xi in seq_along(steps)) {
      x <- steps[xi]; y <- steps[yi]
      col <- (1 - x) * (1 - y) * rgb00 +
        x       * (1 - y) * rgb10 +
        (1 - x) * y       * rgb01 +
        x       * y       * rgb11
      pal <- c(pal, grDevices::rgb(col[1], col[2], col[3]))
      labels <- c(labels, paste(xi, yi, sep = "-"))
    }
  }
  setNames(pal, labels)
}

## Define colours for palette
custom_pal <- bilinear_mix(
  c00 = "#E8E8E8",     # low TOC, low OCl
  c10 = "#83C1E5",   # high TOC, low OCl
  c01 = "#E7C474",   # low TOC, high OCl
  c11 = "#829D71",   # high TOC, high OCl
  n = 3
)

## Combine rasters for surface TOC and OCl and check distribution
combined_OC <- as.data.frame(c(full_surface_TOC, full_surface_OCl), xy = TRUE)
names(combined_OC)[3:4] <- c("TOC", "OCl")

hist(combined_OC$OCl)
hist(combined_OC$TOC)

## Create bivariate dataframe
# ammend style to reflect data distribution
bivar_df <- bi_class(combined_OC, x = TOC, y = OCl, style = "fisher", dim = 3)

bivar_df <- bivar_df %>%
  select(x, y, bi_class)

bivar_df$bi_class <- as.factor(bivar_df$bi_class)
bivar_df$class_id <- as.integer(bivar_df$bi_class)

bivar_OC <- rast(bivar_df[, c("x", "y", "class_id")], type = "xyz")

lvls <- data.frame(
  id = seq_along(levels(bivar_df$bi_class)),
  category = levels(bivar_df$bi_class)
)
levels(bivar_OC) <- lvls

crs(bivar_OC) <- crs(full_surfaceTOC)


###---------------------------COMBINED PLOT-----------------------------------
## Create maps for TOC, OCl, and bivariate vilnerability map
fs_OCl_map <- ggplot() +
  geom_spatraster(data = full_surfaceOCl) +
  scale_fill_stepsn(
    colours = c("#f7f0be", "#cf8f09"), breaks = seq(50, 100, by = 10),
    limits = c(50, 100), name = "OC<sub>L</sub> %", na.value = "transparent") +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = "none")

fs_TOC_map <- ggplot() +
  geom_spatraster(data = full_surfaceTOC) +
  scale_fill_stepsn(
    colours = c("#c0e3f7", "#0072b2"), breaks = seq(0, 300, by = 60),
    limits = c(0, 300), name = "TOC Mg ha<sup>-1</sup>", na.value = "transparent") +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = "none")

bivar_map <- ggplot() +
  geom_spatraster(data = bivar_OC, aes(fill = category), show.legend = FALSE) +
  bi_scale_fill(pal = custom_pal, dim = 3, na.value = "transparent") +
  coord_sf() +
  theme_bw()

## Extract/create legends
fs_OCl_legend <- get_legend(fs_OCl_map + theme(legend.position = "right", legend.title = element_markdown()))

fs_TOC_legend <- get_legend(fs_TOC_map + theme(legend.position = "right", legend.title = element_markdown()))

bivar_legend <- bi_legend(pal = custom_pal, dim = 3, 
                          xlab = "Higher TOC", ylab = "Higher OCl", size = 8) +
  labs(y = expression(paste("Higher OC "[L], "%" %->% "")))

## Combine plots and legends using plot_grid
row_TOC   <- plot_grid(fs_TOC_map,   fs_TOC_legend,   ncol = 2, rel_widths = c(1, 0.25),
                       align = "h", axis = "t")
row_OCl   <- plot_grid(fs_OCl_map,   fs_OCl_legend,   ncol = 2, rel_widths = c(1, 0.25),
                       align = "h", axis = "t")
row_bivar <- plot_grid(bivar_map, bivar_legend, ncol = 2, rel_widths = c(1, 0.25),
                       align = "h", axis = "t")

final_plot <- plot_grid(
  row_TOC, row_OCl, row_bivar, ncol = 1, labels = c("A", "B", "C"), label_size = 12,
  align = "v", axis = "lr")

## View and save final figure
final_plot

ggsave(
  filename = "outputs/final_models/bivariate_figure.png",
  plot = final_plot,
  device = ragg::agg_png,
  width = 8, height = 8, units = "in",
  dpi = 300,
  bg = "white"
)
