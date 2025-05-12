data_rename <- function(data) {

    ## Rename bioclim variables to "BIO" + numbers
    bio_cols <- grep("wc2.1_2.5m_bioc", names(data), value = TRUE)
    name_mapping <- setNames(paste0("BIO", sub(".*_(\\d+)$", "\\1", bio_cols)), bio_cols)
    names(data)[match(names(name_mapping), names(data))] <- unname(name_mapping)

    ## Rename population density variable to "Population"
    pop_col <- grep("^(?i)ssp\\d+(?:_\\d+)?$", names(data), value = TRUE)
    names(data)[names(data) == pop_col] <- "Human Population Density"

    ## Capitalize "Slope" and "Roughness"
    slope_col <- grep("slope", names(data), value = TRUE)
    names(data)[names(data) == slope_col] <- "Slope"
    roughness_col <- grep("roughness", names(data), value = TRUE)
    names(data)[names(data) == roughness_col] <- "Roughness"
    
    return(data)
}


model_baseline <- function(species, odonata_data, background_data, baseline_data, initial_run, background_run, background_start) {

    print(paste(species, "|", ifelse(initial_run, "initial", "baseline")))
    timestamp <- Sys.time()


    ## Load presence and environmental data; extract background points

    data <- baseline_data
    if (!initial_run) {
        selected_vars <- fromJSON(paste0("output/model/", species, "/selected_vars.json"))
        data <- data[[selected_vars]]
    }
    
    presence <- odonata_data[odonata_data$"species" == species, c("decimalLatitude", "decimalLongitude")]
    presence$pa <- 1
    
    absence <- background_data
    absence$pa <- 0

    points_extract <- rbind(presence, absence)


    ## Append environmental variables to training & testing sets

    environment_extract <- extract(x=data, y=points_extract[, c("decimalLongitude", "decimalLatitude")], ID=FALSE)
    points_environment_extract <- cbind(points_extract, environment_extract)

    drop_cols <- which(colnames(points_environment_extract) %in% c("decimalLongitude", "decimalLatitude"))
    points_environment_extract <- points_environment_extract[, -drop_cols]

    set.seed(0)
    fold <- folds(x=points_environment_extract, k=5, by=points_environment_extract$pa)
    testing <- points_environment_extract[fold == 1, ]
    training <- points_environment_extract[fold != 1, ]


    ## MaxEnt model fitting

    set.seed(0)
    maxent_model <- maxent(
        x=training[, colnames(training)[colnames(training) != 'pa']], 
        p=training$pa, removeDuplicates=FALSE
    )


    ## Output selected variables and end function if background run or initial run

    if (background_run) {
        # Return runtime and training AUC
        runtime <- as.numeric(difftime(Sys.time(), background_start, units = "secs"))
        training_auc <- maxent_model@results[5]
        print("Completed!")
        write(paste(background_run, runtime, training_auc, sep=","), file="output/background/background.csv", append=TRUE)
        return()
    }

    if (initial_run) {
        # Write selected variables to JSON
        dir.create(paste0("output/model/", species))
        write_json(
            gsub("\\.", " ", gsub("\\.permutation\\.importance$", "",
                rownames(maxent_model@results[grepl(
                    "permutation.importance", rownames(maxent_model@results)
                ) & maxent_model@results >= 2.5, , drop = FALSE
            ]))),
            path=paste0("output/model/", species, "/selected_vars.json")
        )
        print("Completed!")
        return()
    }
    
    print("Checkpoint: MaxEnt model fitted...")


    ## Continue to prediction and evaluation if final run

    save(maxent_model, file=paste0("output/model/", species, "/model.RData"))  # Save model for future predictions
    maxent_predict <- predict(data, maxent_model, na.rm=TRUE, type="response")
    maxent_eval <- pa_evaluate(
        p=testing[testing$pa == 1, ], a=testing[testing$pa == 0, ], 
        model=maxent_model, type="response"
    )
    maxent_threshold <- maxent_eval@thresholds$max_spec_sens  # Save threshold for future predictions

    print("Completed!")
    runtime <- as.numeric(difftime(Sys.time(), timestamp, units = "secs"))


    ## Save results as RData, raster, and CSV

    writeRaster(maxent_predict, paste0("output/model/", species, "/baseline.tif"), overwrite=TRUE)
    binary_prediction <- maxent_predict > maxent_threshold
    
    baseline_values <<- list()
    baseline_values[["suitability"]] <<- mean(maxent_predict[], na.rm=TRUE)  # Habitat suitability
    baseline_values[["range"]] <<- sum(binary_prediction[], na.rm=TRUE)  # Range size
    baseline_values[["centroid"]] <<- colMeans(xyFromCell(binary_prediction, which(binary_prediction[], arr.ind=FALSE)))  # Range centroid
    baseline_values[["threshold"]] <<- maxent_threshold

    var_importance <- maxent_model@results[grepl("permutation.importance", rownames(maxent_model@results)), ]
    var_importance_output <- ""
    for (var in names(baseline_data)) {
        if (var %in% selected_vars) {
            var_importance_output <- paste0(var_importance_output, ",", var_importance[paste0(var, ".permutation.importance")])
        } else {
            var_importance_output <- paste0(var_importance_output, ",")
        }
    }

    write(paste(species,
        baseline_values[["suitability"]],
        baseline_values[["range"]],
        baseline_values[["centroid"]]["x"],
        baseline_values[["centroid"]]["y"],
        runtime, baseline_values[["threshold"]],
        maxent_model@results[5], maxent_eval@stats$auc,  # Training & testing AUC
        substring(var_importance_output, 2),  # Permutation importance of each variable
    sep=","), file="output/model/baseline.csv", append=TRUE)
}


model_future <- function(species, year, ssp, future_data) {

    print(paste(species, "|", year, "|", ssp))
    timestamp <- Sys.time()


    ## Load environmental data; MaxEnt prediction and evaluation

    selected_vars <- fromJSON(paste0("output/model/", species, "/selected_vars.json"))
    data <- future_data[[year]][[ssp]][[selected_vars]]

    load(paste0("output/model/", species, "/model.RData"))
    maxent_predict <- predict(data, maxent_model, na.rm=TRUE, type="response")  # Need current model from previous baseline runs

    print("Completed!")
    runtime <- as.numeric(difftime(Sys.time(), timestamp, units = "secs"))


    ## Save results as raster and CSV

    writeRaster(maxent_predict, paste0("output/model/", species, "/", year, "_", ssp, ".tif"), overwrite=TRUE)  # Prediction

    binary_prediction <- maxent_predict > baseline_values[["threshold"]]
    centroid <- colMeans(xyFromCell(binary_prediction, which(binary_prediction[], arr.ind=FALSE)))

    write(paste(species, year, ssp,
        mean(maxent_predict[], na.rm=TRUE),  # Habitat suitability
        sum(binary_prediction[], na.rm=TRUE),  # Range size
        centroid["x"], centroid["y"], runtime,
        mean(maxent_predict[], na.rm=TRUE) / baseline_values[["suitability"]],
        sum(binary_prediction[], na.rm=TRUE)  / baseline_values[["range"]],
        distm(centroid, baseline_values[["centroid"]], fun = distGeo) / 1000,  # Centroid shift distance in km
    sep=","), file="output/model/future.csv", append=TRUE)
}
