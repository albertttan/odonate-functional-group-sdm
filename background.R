library(terra)
library(dismo)
library(predicts)
library(jsonlite)
library(geosphere)
source("functions.R")


## Load data from files

odonata_data <- read.csv("output/cleaning/presence.csv")
selected_species <- "Argia sedula"
# selected_species <- fromJSON("output/cleaning/selected_species.json")

baseline_data <- data_rename(c(
    rast("data/climate/baseline.tif"),
    rast("data/landuse/baseline.tif"),
    rast("data/population/baseline.tif"),
    rast("data/topography/slope.tif"),
    rast("data/topography/roughness.tif")
))


## Create output file if necessary

if (!file.exists("output/background/background.csv")) {
    write("points,runtime,training_auc", file="output/background/background.csv")
}


## Load background points and run models

for (power in 14:16) {
    set.seed(0)
    size = 2 ^ power
    timestamp <- Sys.time()
    background_data <- as.data.frame(spatSample(x=baseline_data, size=size, values=FALSE, na.rm=TRUE, xy=TRUE))
    colnames(background_data) <- c("decimalLongitude", "decimalLatitude")
    model_baseline(selected_species, odonata_data, background_data, baseline_data, TRUE, size, timestamp)
}
