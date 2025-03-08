library(terra)
library(usdm)
library(jsonlite)

bioclim_data <- rast("data/climate/2020.tif") [[c(1, 4, 7, 10, 11, 12)]]
topography_data_1 <- rast("data/topography/slope.tif")
topography_data_2 <- rast("data/topography/roughness.tif")
population_data <- rast("data/population/2020.tif")
landuse_data <- rast("data/landuse/2020.tif")[[c(2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)]]
environment_data <- c(bioclim_data, topography_data_1, population_data, landuse_data)
vif_result <- vif(environment_data)
print(vif_result)
save(vif_result, file="output/multicollinearity/multicollinearity_selected.Rdata")
write_json(vif_result$Variables, path="output/multicollinearity/multicollinearity.json", pretty=TRUE)
