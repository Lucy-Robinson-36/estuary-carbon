####----------------------------CARBON-MAP 1--------------------------------
### Calculating TOC stock and average %OCl per core

###--------------------------SET UP WORKSPACE-------------------------------
## Clear the environment
rm(list=ls()) 

## Set working directory
setwd("path/Carbon_Map")

## Load in packages
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(rjags)

## Read in files
soli_data <- read.csv("samples/geochemical_results.csv")
layers <- read.csv("samples/marsh_horizons.csv")
DBD <- read.csv("samples/bulk_density.csv")

###--------------------------PREPARE DATA-----------------------------------
## Filter soli data
soli_data_clean <- soli_data %>%
  # remove blanks and standards
  filter(!(Name %in% c("blnk", "silty soil (toc/roc/tic)"))) %>%
  # calculate percentage of labile OC relative to total OC
  mutate(perc_OCl = 100 * `TOC400  [%]`/(`TOC  [%]`)) %>%
  # expand Name column
  separate(Name, c("Site", "Transect", "Core", "Depth"), sep = "_", remove = TRUE) %>%
  # expand Depth column
  separate(Depth, c("Depth_1", "Depth_2"), sep = "-", remove = TRUE) %>%
  # select columns
  select("Site", "Transect", "Core", "Depth_1", "Depth_2", "TOC  [%]", "perc_OCl")

## Add bulk density data
soli_DBD <- merge(soli_data_clean, DBD, by = c("Site", "Transect", "Core",  "Depth_1", "Depth_2"), all.x = FALSE)

## Depth of saltmarsh soils
marsh_depth <- layers %>%
  mutate(Marsh_depth_cm = Fibrous_peat_thickness_cm + Humified_peat_thickness_cm  + Transitional_thickness_cm) %>%
  select("Site", "Transect", "Core", "Marsh_depth_cm")

## Select surface samples
surface_data <- soli_DBD %>%
  filter("Depth_2" <= 10) %>%
  # add weighting
  mutate(Weighting = case_when(Depth_1 == 0  ~ 4,
                               Depth_1 == 4  ~ 6))


## Select saltmarsh soil samples
# First merge soli data and marsh depth dataframes
marsh_depth_data <- merge(soli_DBD, marsh_depth, by = c("Site", "Transect", "Core"), all.x = FALSE)

saltmarsh_data <- marsh_depth_data %>%
  # add and "end_depth" column - this is the lower depth that the sample represents
  # e.g. the 0-1cm sample represents the sediments from 0 to 4cm, where the next sample was taken
  mutate(end_depth = case_when(
    Depth_1 == 0 ~ 4,
    Depth_1 == 4 ~ 10,
    Depth_1 == 10 ~ 20,
    Depth_1 == 20 ~ 30,
    Depth_1 == 30 ~ 40,
    Depth_1 == 40 ~ 50,
    Depth_1 == 50 ~ 60,
    Depth_1 == 60 ~ 70,
    Depth_1 == 70 ~ 80,
    Depth_1 == 80 ~ 90,
    Depth_1 == 90 ~ 100)) %>%
  # filter out samples that are only basal sediments
  filter(Depth_1 <= Marsh_depth_cm) %>%
  # add weighting to reflect depth interval
  mutate(Weighting = case_when(
    # samples containing basal sediments
      end_depth > Marsh_depth_cm ~ Marsh_depth_cm - Depth_1,
    # samples representing just saltmarsh sediments
      Depth_1 == 0  ~ 4,
      Depth_1 == 4  ~ 6,
      Depth_1 >= 10 ~ 10)) %>%
  # drop redundant columns
  select(-Marsh_depth_cm, -end_depth)

###---------------------------SURFACE TOC----------------------------------
## Create ID and index columns
surface_data$Core_ID <- paste(surface_data$Site, surface_data$Transect, surface_data$Core, sep = "_")
surface_data$core_index <- as.numeric(as.factor(surface_data$Core_ID))

n_cores <- length(unique(surface_data$core_index))

## Baseline stock calculations
test_surfaceTOC <- surface_data %>%
  group_by(Core_ID) %>%
  summarise(Stock_g.cm2 = sum((`TOC  [%]`/100) * DBD_g.cm3 * Weighting))

## Create model string
# Adjust priors to reflect data
model_string <- "
model {

  # ----- Interval level -----
  for (i in 1:N) {

    # Observations
    OC_obs[i] ~ dnorm(OC_true[i], tau_OC[i])
    BD_obs[i] ~ dnorm(BD_true[i], tau_BD[i])

    # True latent values
    OC_true[i] ~ dnorm(OC_mean, tau_OC_mean)
    BD_true[i] ~ dnorm(BD_mean, tau_BD_mean)

    # Stock per interval
    stock[i] <- OC_true[i] * BD_true[i] * weight[i]
  }

  # ----- Sum per core -----
  for (j in 1:n_cores) {
    surface_TOC[j] <- sum(stock[] * equals(core_index[], j))
  }

  # ----- Priors -----
  OC_mean ~ dnorm(0.05, 100)
  tau_OC_mean ~ dgamma(0.001, 0.001)

  BD_mean ~ dnorm(0.3, 10)
  tau_BD_mean ~ dgamma(0.001, 0.001)
}
"

## Prepare data for JAGS
data_jags <- list(
  N = nrow(surface_data),
  n_cores = n_cores,
  core_index = surface_data$core_index,
  OC_obs = surface_data$`TOC  [%]` / 100,
  BD_obs = surface_data$DBD_g.cm3,
  weight = surface_data$Weighting,
  tau_OC = rep(100, nrow(surface_data)),
  tau_BD = rep(100, nrow(surface_data))
)

## Initialise jags model
tau_OC <- 1 / (0.04^2)
tau_BD <- 1 / (0.05^2)

## Run JAGS model
model <- jags.model(
  textConnection(model_string),
  data = data_jags,
  n.chains = 3
)

update(model, 5000)

samples <- coda.samples(
  model,
  variable.names = c("surface_TOC"),
  n.iter = 50000
)

## Extract results
surfacesamples_mat <- as.matrix(samples)
stock_cols <- grep("^surface_TOC", colnames(surfacesamples_mat))

surfaceTOC_results <- data.frame(
  Core_ID = levels(as.factor(surface_data$Core_ID)),
  `mean_gC.cm2` = colMeans(surfacesamples_mat[, stock_cols]),
  `sd_gC.cm2` = apply(surfacesamples_mat[, stock_cols], 2, sd)
)

## Compare deterministic and jags results
compare_vals <- test_surfaceTOC %>%  left_join(surfaceTOC_results, by = "Core_ID")

compare_vals <- compare_vals %>%
  separate(Core_ID, c('Site', 'Transect', "Core"), sep = "_", remove = FALSE)

ggplot(compare_vals, aes(x = Stock_g.cm2, y = mean_gC.cm2, colour = Site)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = (mean_gC.cm2 - sd_gC.cm2), ymax = (mean_gC.cm2 + sd_gC.cm2)), width = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Deterministic",
    y = "JAGS (mean ± sd)",
    title = "Surface TOC: deterministic vs JAGS"
  ) +
  theme_minimal()

## Outlier analysis
compare_vals$residual <- compare_vals$mean_gC.cm2 - compare_vals$Stock_g.cm2

Q1 <- quantile(compare_vals$residual, 0.25)
Q3 <- quantile(compare_vals$residual, 0.75)
IQR_val <- Q3 - Q1

compare_vals$outlier <- compare_vals$residual < (Q1 - 1.5*IQR_val) |
  compare_vals$residual > (Q3 + 1.5*IQR_val)

subset(compare_vals, outlier == TRUE)

ggplot(compare_vals, aes(x = Stock_g.cm2, y = mean_gC.cm2)) +
  geom_point(aes(color = outlier), size = 3) +
  geom_abline(linetype = "dashed") +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal()

# With all cores
mean_all <- mean(compare_vals$mean_gC.cm2)
print(mean_all)

# Without outliers
mean_no_outliers <- mean(compare_vals$mean_gC.cm2[!compare_vals$outlier])
print(mean_no_outliers) 

## Export to csv
write.csv(surfaceTOC_results, "outputs/core_results/Surface_TOC.csv")


###-----------------------------SURFACE OCl--------------------------------
## Baseline stock calculations
test_surfaceOCl <- surface_data %>%
  group_by(Core_ID) %>%
  summarise(det_perc_OCl = sum(perc_OCl*Weighting)/sum(Weighting))

## Create model string
# Adjust priors to reflect data
## Create model string
model_string <- "
model {

  # ----- Interval level -----
  for (i in 1:N) {

    # Observations
    OCl_obs[i] ~ dnorm(OCl_true[i], tau_OCl[i])

    # True latent values
    OCl_true[i] ~ dnorm(OCl_mean, tau_OCl_mean)

    # Weight per interval
    weighted_OCl[i] <- OCl_true[i] * weight[i]
  }

  # ----- Sum per core -----
  
  for (j in 1:n_cores) {
    surface_OCl[j] <-
      sum(weighted_OCl[] * equals(core_index[], j)) /
      sum(weight[] * equals(core_index[], j))
  }


  # ----- Priors -----
  OCl_mean ~ dnorm(0.05, 100)
  tau_OCl_mean ~ dgamma(0.001, 0.001)
}
"

## Prepare data for JAGS
data_jags <- list(
  N = nrow(surface_data),
  n_cores = n_cores,
  core_index = surface_data$core_index,
  OCl_obs = surface_data$perc_OCl,
  weight = surface_data$Weighting,
  tau_OCl = rep(100, nrow(surface_data))
)


## Run JAGS model
model <- jags.model(
  textConnection(model_string),
  data = data_jags,
  n.chains = 3
)

update(model, 5000)

samples <- coda.samples(
  model,
  variable.names = c("surface_OCl"),
  n.iter = 50000
)

## Extract results
surfacesamples_mat <- as.matrix(samples)
OCl_cols <- grep("^surface_OCl", colnames(surfacesamples_mat))

surfaceOCl_results <- data.frame(
  Core_ID = levels(as.factor(surface_data$Core_ID)),
  mean_perc_OCl = colMeans(surfacesamples_mat[, OCl_cols]),
  sd_perc_OCl = apply(surfacesamples_mat[, OCl_cols], 2, sd)
)

## Compare deterministic and jags results
compare_vals <- test_surfaceOCl %>%  left_join(surfaceOCl_results, by = "Core_ID")

compare_vals <- compare_vals %>%
  separate(Core_ID, c('Site', 'Transect', "Core"), sep = "_", remove = FALSE)

ggplot(compare_vals, aes(x = det_perc_OCl, y = mean_perc_OCl, colour = Site)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = (mean_perc_OCl - sd_perc_OCl), ymax = (mean_perc_OCl + sd_perc_OCl)), width = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Deterministic",
    y = "JAGS (mean ± sd)",
    title = "Surface OCl: deterministic vs JAGS"
  ) +
  theme_minimal()

## Outlier analysis
compare_vals$residual <- compare_vals$mean_perc_OCl - compare_vals$det_perc_OCl

Q1 <- quantile(compare_vals$residual, 0.25)
Q3 <- quantile(compare_vals$residual, 0.75)
IQR_val <- Q3 - Q1

compare_vals$outlier <- compare_vals$residual < (Q1 - 1.5*IQR_val) |
  compare_vals$residual > (Q3 + 1.5*IQR_val)

subset(compare_vals, outlier == TRUE)

ggplot(compare_vals, aes(x = det_perc_OCl, y = mean_perc_OCl)) +
  geom_point(aes(color = outlier), size = 3) +
  geom_abline(linetype = "dashed") +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal()

# With all cores
mean_all <- mean(compare_vals$mean_perc_OCl)
print(mean_all)

# Without outliers
mean_no_outliers <- mean(compare_vals$mean_perc_OCl[!compare_vals$outlier])
print(mean_no_outliers) 

## Export to csv
write.csv(surfaceOCl_results, "outputs/core_results/Surface_OCl.csv")
  

###---------------------------DEPTH-INTEGRATED TOC----------------------------------
## Create ID and index columns
saltmarsh_data$Core_ID <- paste(saltmarsh_data$Site, saltmarsh_data$Transect, saltmarsh_data$Core, sep = "_")
saltmarsh_data$core_index <- as.numeric(as.factor(saltmarsh_data$Core_ID))

n_cores <- length(unique(saltmarsh_data$core_index))

## Baseline stock calculations
test_depthTOC <- saltmarsh_data %>%
  group_by(Core_ID) %>%
  summarise(Stock_g.cm2 = sum((`TOC  [%]`/100) * DBD_g.cm3 * Weighting))

## Create model string
# Adjust priors to reflect data
model_string <- "
model {

  # ----- Interval level -----
  for (i in 1:N) {

    # Observations
    OC_obs[i] ~ dnorm(OC_true[i], tau_OC[i])
    BD_obs[i] ~ dnorm(BD_true[i], tau_BD[i])

    # True latent values
    OC_true[i] ~ dnorm(OC_mean, tau_OC_mean)
    BD_true[i] ~ dnorm(BD_mean, tau_BD_mean)

    # Stock per interval
    stock[i] <- OC_true[i] * BD_true[i] * weight[i]
  }

  # ----- Sum per core -----
  for (j in 1:n_cores) {
    depth_TOC[j] <- sum(stock[] * equals(core_index[], j))
  }

  # ----- Priors -----
  OC_mean ~ dnorm(0.05, 100)
  tau_OC_mean ~ dgamma(0.001, 0.001)

  BD_mean ~ dnorm(0.3, 10)
  tau_BD_mean ~ dgamma(0.001, 0.001)
}
"

## Prepare data for JAGS
data_jags <- list(
  N = nrow(saltmarsh_data),
  n_cores = n_cores,
  core_index = saltmarsh_data$core_index,
  OC_obs = saltmarsh_data$`TOC  [%]` / 100,
  BD_obs = saltmarsh_data$DBD_g.cm3,
  weight = saltmarsh_data$Weighting,
  tau_OC = rep(100, nrow(saltmarsh_data)),
  tau_BD = rep(100, nrow(saltmarsh_data))
)

## Initialise jags model
tau_OC <- 1 / (0.04^2)
tau_BD <- 1 / (0.05^2)

## Run JAGS model
model <- jags.model(
  textConnection(model_string),
  data = data_jags,
  n.chains = 3
)

update(model, 5000)

samples <- coda.samples(
  model,
  variable.names = c("depth_TOC"),
  n.iter = 50000
)

## Extract results
depthsamples_mat <- as.matrix(samples)
stock_cols <- grep("^depth_TOC", colnames(depthsamples_mat))

depthTOC_results <- data.frame(
  Core_ID = levels(as.factor(depth_data$Core_ID)),
  `mean_gC.cm2` = colMeans(depthsamples_mat[, stock_cols]),
  `sd_gC.cm2` = apply(depthsamples_mat[, stock_cols], 2, sd)
)

## Compare deterministic and jags results
compare_vals <- test_depthTOC %>%  left_join(depthTOC_results, by = "Core_ID")

compare_vals <- compare_vals %>%
  separate(Core_ID, c('Site', 'Transect', "Core"), sep = "_", remove = FALSE)

ggplot(compare_vals, aes(x = Stock_g.cm2, y = mean_gC.cm2, colour = Site)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = (mean_gC.cm2 - sd_gC.cm2), ymax = (mean_gC.cm2 + sd_gC.cm2)), width = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Deterministic",
    y = "JAGS (mean ± sd)",
    title = "Surface TOC: deterministic vs JAGS"
  ) +
  theme_minimal()

## Outlier analysis
compare_vals$residual <- compare_vals$mean_gC.cm2 - compare_vals$Stock_g.cm2

Q1 <- quantile(compare_vals$residual, 0.25)
Q3 <- quantile(compare_vals$residual, 0.75)
IQR_val <- Q3 - Q1

compare_vals$outlier <- compare_vals$residual < (Q1 - 1.5*IQR_val) |
  compare_vals$residual > (Q3 + 1.5*IQR_val)

subset(compare_vals, outlier == TRUE)

ggplot(compare_vals, aes(x = Stock_g.cm2, y = mean_gC.cm2)) +
  geom_point(aes(color = outlier), size = 3) +
  geom_abline(linetype = "dashed") +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal()

# With all cores
mean_all <- mean(compare_vals$mean_gC.cm2)
print(mean_all)

# Without outliers
mean_no_outliers <- mean(compare_vals$mean_gC.cm2[!compare_vals$outlier])
print(mean_no_outliers) 

## Export to csv
write.csv(depthTOC_results, "outputs/core_results/Depth_TOC.csv")


###-----------------------------DEPTH-INTEGRATED OCl--------------------------------
## Baseline stock calculations
test_depthOCl <- saltmarsh_data %>%
  group_by(Core_ID) %>%
  summarise(det_perc_OCl = sum(perc_OCl*Weighting)/sum(Weighting))

## Create model string
# Adjust priors to reflect data
## Create model string
model_string <- "
model {

  # ----- Interval level -----
  for (i in 1:N) {

    # Observations
    OCl_obs[i] ~ dnorm(OCl_true[i], tau_OCl[i])

    # True latent values
    OCl_true[i] ~ dnorm(OCl_mean, tau_OCl_mean)

    # Weight per interval
    weighted_OCl[i] <- OCl_true[i] * weight[i]
  }

  # ----- Sum per core -----
  
  for (j in 1:n_cores) {
    depth_OCl[j] <-
      sum(weighted_OCl[] * equals(core_index[], j)) /
      sum(weight[] * equals(core_index[], j))
  }


  # ----- Priors -----
  OCl_mean ~ dnorm(0.05, 100)
  tau_OCl_mean ~ dgamma(0.001, 0.001)
}
"

## Prepare data for JAGS
data_jags <- list(
  N = nrow(depth_data),
  n_cores = n_cores,
  core_index = depth_data$core_index,
  OCl_obs = depth_data$perc_OCl,
  weight = depth_data$Weighting,
  tau_OCl = rep(100, nrow(depth_data))
)

## Run JAGS model
model <- jags.model(
  textConnection(model_string),
  data = data_jags,
  n.chains = 3
)

update(model, 5000)

samples <- coda.samples(
  model,
  variable.names = c("depth_OCl"),
  n.iter = 50000
)

## Extract results
depthsamples_mat <- as.matrix(samples)
OCl_cols <- grep("^depth_OCl", colnames(depthsamples_mat))

depthOCl_results <- data.frame(
  Core_ID = levels(as.factor(saltmarsh_data$Core_ID)),
  mean_perc_OCl = colMeans(depthsamples_mat[, OCl_cols]),
  sd_perc_OCl = apply(depthsamples_mat[, OCl_cols], 2, sd)
)

## Compare deterministic and jags results
compare_vals <- test_depthOCl %>%  left_join(depthOCl_results, by = "Core_ID")

compare_vals <- compare_vals %>%
  separate(Core_ID, c('Site', 'Transect', "Core"), sep = "_", remove = FALSE)

ggplot(compare_vals, aes(x = det_perc_OCl, y = mean_perc_OCl, colour = Site)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = (mean_perc_OCl - sd_perc_OCl), ymax = (mean_perc_OCl + sd_perc_OCl)), width = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Deterministic",
    y = "JAGS (mean ± sd)",
    title = "Depth OCl: deterministic vs JAGS"
  ) +
  theme_minimal()

## Outlier analysis
compare_vals$residual <- compare_vals$mean_perc_OCl - compare_vals$det_perc_OCl

Q1 <- quantile(compare_vals$residual, 0.25)
Q3 <- quantile(compare_vals$residual, 0.75)
IQR_val <- Q3 - Q1

compare_vals$outlier <- compare_vals$residual < (Q1 - 1.5*IQR_val) |
  compare_vals$residual > (Q3 + 1.5*IQR_val)

subset(compare_vals, outlier == TRUE)

ggplot(compare_vals, aes(x = det_perc_OCl, y = mean_perc_OCl)) +
  geom_point(aes(color = outlier), size = 3) +
  geom_abline(linetype = "dashed") +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal()

# With all cores
mean_all <- mean(compare_vals$mean_perc_OCl)
print(mean_all)

# Without outliers
mean_no_outliers <- mean(compare_vals$mean_perc_OCl[!compare_vals$outlier])
print(mean_no_outliers) 

## Export to csv
write.csv(depthOCl_results, "outputs/core_results/Depth_OCl.csv")

