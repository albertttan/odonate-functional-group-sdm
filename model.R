library(terra)
library(dismo)
library(predicts)
library(jsonlite)
library(geosphere)
source("functions.R")


## Set parameters

initial_run = FALSE
selected_bioclim_vars <- c(1, 4, 7, 10, 11, 12)
selected_landuse_vars <- c(2:3, 5:19)
selected_species_index <- c(1)
selected_years <- c("2050")
selected_ssps <- c("ssp1")


## Load data from files

odonata_data <- read.csv("output/cleaning/presence.csv")
selected_species <- fromJSON("output/cleaning/selected_species.json")

baseline_data <- data_rename(c(
    rast("data/climate/baseline.tif")[[selected_bioclim_vars]],
    rast("data/landuse/baseline.tif")[[selected_landuse_vars]],
    rast("data/population/baseline.tif"),
    rast("data/topography/slope.tif")
))

future_data <- list()
for (year in selected_years) {
    for (ssp in selected_ssps) {
        future_data[[year]][[ssp]] <- data_rename(c(
            rast(paste0("data/climate/", year, "/", ssp, ".tif"))[[selected_bioclim_vars]],
            rast(paste0("data/landuse/", year, "/", ssp, ".tif"))[[selected_landuse_vars]],
            rast(paste0("data/population/", year, "/", ssp, ".tif")),
            rast("data/topography/slope.tif")
        ))
    }
}


## Load background points

if (file.exists("output/model/background.RData")) {
    load("output/model/background.RData")
} else {
    set.seed(0)
    background_data <- as.data.frame(spatSample(x=baseline_data, size=8192, values=FALSE, na.rm=TRUE, xy=TRUE))
    colnames(background_data) <- c("decimalLongitude", "decimalLatitude")
    save(background_data, file="output/model/background.RData")
}


## Create output files if necessary

if (!file.exists("output/model/baseline.csv")) {
    write(paste0("species,suitability,range,centroid_x,centroid_y,runtime,threshold,training_auc,testing_auc,", paste(names(baseline_data), collapse=",")), file="output/model/baseline.csv")
    write("species,year,ssp,suitability,range,centroid_x,centroid_y,runtime,suitability_change,range_change,distance_change", file="output/model/future.csv")
}


## Main loop

for (species in selected_species[selected_species_index]) {
    model_baseline(species, odonata_data, background_data, baseline_data, initial_run, FALSE, FALSE)
    if (!initial_run) {
        for (year in selected_years) {
            for (ssp in selected_ssps) {
                model_future(species, year, ssp, future_data)
            }
        }
    }
}
