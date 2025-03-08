library(terra)
library(geodata)
library(dismo)
library(predicts)
library(ncdf4)
library(jsonlite)


# Loading Data

selected_bioclim_vars <- c(1, 4, 7, 10, 11, 12)
selected_landuse_vars <- c(1:3, 7:19)
odonata_data <<- read.delim("data/odonata/odonata.csv")
bioclim_data <- rast("data/climate/2020.tif")[[selected_bioclim_vars]]
topography_data <- rast("data/topography/slope.tif")
population_data <- rast("data/population/2020.tif")
landuse_data <- rast("data/landuse/2020.tif")[[selected_landuse_vars]]
environment_data <<- c(bioclim_data, topography_data, population_data, landuse_data)
world_map <- world(resolution=3, path="data/")
cropped_map <<- crop(x=world_map, y=ext(x = c(-180, -50, 5, 85)))

selected_species <- fromJSON("data/odonata/selected_species.json")
background_points <- c(16000, 32000, 64000)

process_species <- function(species, size) {
    start_time <- Sys.time()
    print(paste("Processing", species, "with", size, "background points"))
    
    # Presence-Absence Data Points
    set.seed(1)
    obs_data <- odonata_data[odonata_data$"species" == species, c("decimalLatitude", "decimalLongitude")]
    background <- spatSample(x=environment_data, size=size, values=FALSE, na.rm=TRUE, xy=TRUE)
    
    plot(cropped_map, axes=TRUE, col="grey95")
    points(background, col="grey30", pch=1, cex=0.5)
    points(x=obs_data$decimalLongitude, y=obs_data$decimalLatitude, col="salmon", pch=20, cex=0.5)
    
    presence <- obs_data
    presence$pa <- 1
    absence <- as.data.frame(background)
    colnames(absence) <- c("decimalLongitude", "decimalLatitude")
    absence$pa <- 0
    
    points_extract <- rbind(presence, absence)
    print("Presence and absence points extracted...")
    print(Sys.time() - start_time)
    
    # Adding Environment Data
    environment_extract <- extract(x=environment_data, y=points_extract[, c("decimalLongitude", "decimalLatitude")], ID=FALSE)
    points_environment_extract <- cbind(points_extract, environment_extract)
    drop_cols <- which(colnames(points_environment_extract) %in% c("decimalLongitude", "decimalLatitude"))
    points_environment_extract <- points_environment_extract[, -drop_cols]
    
    fold <- folds(x=points_environment_extract, k=5, by=points_environment_extract$pa)
    testing <- points_environment_extract[fold == 1, ]
    training <- points_environment_extract[fold != 1, ]
    print("Environment data included...")
    print(Sys.time() - start_time)
    
    # MaxEnt Model Building
    maxent_model <- maxent(x=training[, colnames(training)[colnames(training) != 'pa']], p=training$pa, removeDuplicates=FALSE, path="output/maxent_files")
    maxent_predict <- predict(environment_data, maxent_model, na.rm=TRUE, type="response")
    print("MaxEnt model constructed...")
    print(Sys.time() - start_time)
    
    # MaxEnt Model Evaluation
    maxent_eval <- pa_evaluate(p=testing[testing$pa == 1, ], a=testing[testing$pa == 0, ], model=maxent_model, type="response")
    maxent_threshold <- maxent_eval@thresholds$max_spec_sens
    binary_prediction <- maxent_predict > maxent_threshold
    raw_results <<- c(model = maxent_model, predict = maxent_predict, eval = maxent_eval)
    
    # Variables and Evaluations
    end_time <- Sys.time()
    eval_training_auc <- maxent_model@results[5]
    eval_testing_auc <- maxent_eval@stats$auc
    eval_runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
    print(paste("Finished processing", species))
   
    return(list(
        raw_results = c(model=maxent_model, predict=maxent_predict, eval=maxent_eval),
        eval_training_auc = eval_training_auc,
        eval_testing_auc = eval_testing_auc,
        eval_runtime = eval_runtime
    ))
}

training_auc <- c()
testing_auc <- c()
runtime <- c()
for (points in background_points) {
    result <- process_species("Aeshna interrupta", points)
    training_auc <- c(training_auc, result$eval_training_auc)
    testing_auc <- c(testing_auc, result$eval_testing_auc)
    runtime <- c(runtime, result$eval_runtime)
}

output <- data.frame(Size = background_points, Training_AUC = training_auc, Testing_AUC = testing_auc, Time = runtime)
write_json(output, path="output/background_points.json", pretty=TRUE)
