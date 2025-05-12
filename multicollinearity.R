library(terra)
library(usdm)
library(jsonlite)
source("functions.R")  # For the function "data_rename()"


## Load data

selected_bioclim_vars <- c(1, 4, 7, 10, 11, 12)
selected_landuse_vars <- c(2:3, 5:19)

original_data <- data_rename(c(
    rast("data/climate/baseline.tif"),
    rast("data/landuse/baseline.tif"),
    rast("data/population/baseline.tif"),
    rast("data/topography/slope.tif"),
    rast("data/topography/roughness.tif")
))

selected_data <- data_rename(c(
    rast("data/climate/baseline.tif")[[selected_bioclim_vars]],
    rast("data/landuse/baseline.tif")[[selected_landuse_vars]],
    rast("data/population/baseline.tif"),
    rast("data/topography/slope.tif")
))


## Run VIF analysis

original_vif <- vif(original_data)
selected_vif <- vif(selected_data)


# Save results

save(original_vif, selected_vif, file="output/multicollinearity/multicollinearity.Rdata")
write_json(selected_vif$Variables, path="output/multicollinearity/selected_vars.json", pretty=TRUE)
