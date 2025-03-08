options(java.parameters = "-Xmx6g")
library(terra)
library(geodata)
library(dismo)
library(predicts)
library(ncdf4)
library(jsonlite)


# Settings

initial_run = FALSE
selected_bioclim_vars <- c(1, 4, 7, 10, 11, 12)
selected_landuse_vars <- c(2:3, 5:19)
selected_species_index <- c(29, 30, 26, 14) # Need finals for 14+26+29+30; 14 will be quick
selected_years <- c("2020", "2050", "2075", "2100")
selected_ssps <- c("ssp1", "ssp2", "ssp5")


# Loading Data

odonata_data <- read.delim("data/odonata/odonata.csv")
odonata_data <- odonata_data[odonata_data$year >= 2020, ]
world_map <- world(resolution=3, path="data/")
cropped_map <- crop(x=world_map, y=ext(x = c(-180, -50, 5, 85)))
selected_species <- fromJSON("data/odonata/selected_species.json")

environment_data <- c()
for (year in selected_years) {
    if (year == "2020") {
        bioclim_data <- rast("data/climate/2020.tif")[[selected_bioclim_vars]]
        population_data <- rast("data/population/2020.tif")
        topography_data <- rast("data/topography/slope.tif")
        landuse_data <- rast("data/landuse/2020.tif")[[selected_landuse_vars]]
        environment_data[["2020"]] <- c(bioclim_data, population_data, topography_data, landuse_data)
    } else {
        temp_data <- list()
        for (ssp in selected_ssps) {
            bioclim_data <- rast(paste0("data/climate/", year, "/", ssp, ".tif"))[[selected_bioclim_vars]]
            population_data <- rast(paste0("data/population/", year, "/", ssp, ".tif"))
            topography_data <- rast("data/topography/slope.tif")
            landuse_data <- rast(paste0("data/landuse/", year, "/", ssp, ".tif"))[[selected_landuse_vars]]
            temp_data[[ssp]] <- c(bioclim_data, population_data, topography_data, landuse_data)
        }
        environment_data[[year]] <- temp_data
    }
}

eval_plots <- function(tr_stats, species) {
    pdf(paste0("output/maxent_files/", species, ifelse(initial_run, "/initial_output/plots.pdf", "/final_output/plots.pdf")))
    plot_data <- data.frame(
        threshold = tr_stats$treshold,
        TPR = as.numeric(tr_stats$TPR),    # Recall
        FPR = as.numeric(tr_stats$FPR),
        PPP = as.numeric(tr_stats$PPP)     # Precision
    )
    plot_data <- plot_data[is.finite(plot_data$TPR) & 
                               is.finite(plot_data$FPR) & 
                               is.finite(plot_data$PPP), ]
    plot_data <- plot_data[order(plot_data$threshold), ]
    par(mfrow = c(1, 2), pty = "s")
    plot(plot_data$FPR, plot_data$TPR, 
         type = "l", 
         col = "blue",
         lwd = 2,
         xlab = "1 – Specificity",
         ylab = "Sensitivity",
         main = "ROC Curve",
         xlim = c(0, 1),
         ylim = c(0, 1),
         panel.first = grid())
    abline(0, 1, lty = 2, col = "gray")
    box()
    plot(plot_data$TPR, plot_data$PPP,
         type = "l",
         col = "blue",
         lwd = 2,
         xlab = "Recall",
         ylab = "Precision",
         main = "Precision-Recall",
         xlim = c(0, 1),
         ylim = c(0, 1),
         panel.first = grid())
    box()
    par(mfrow = c(1, 1), pty = "m")
    dev.off()
}

env_data_rename <- function(env_data) {
    bio_cols <- grep("wc2.1_2.5m_bioc", names(env_data), value = TRUE)
    bio_numbers <- sub(".*_(\\d+)$", "\\1", bio_cols)
    new_names <- paste0("Bioclim", bio_numbers)
    name_mapping <- setNames(new_names, bio_cols)
    names(env_data)[match(names(name_mapping), names(env_data))] <- unname(name_mapping)
    pop_col <- grep("^(?i)ssp\\d+(?:_\\d+)?$", names(env_data), value = TRUE)
    names(env_data)[names(env_data) == pop_col] <- "Population"
    return(env_data)
}

clean_output_json <- function(output) {
    if (output$year == "2020") {
        selected_outputs <- c("species", "year", "ssp", "range_mean_suitability", "range_size", "range_centroid_x", "range_centroid_y", "eval_runtime", "eval_training_auc", "eval_testing_auc", "output_threshold")
        cleaned_output <- output[selected_outputs]
        for (var_name in names(env_data_rename(environment_data$"2020")[[current_selected_vars]])) {
            cleaned_output[var_name] <- output$eval_var_importance[paste0(var_name, ".permutation.importance")]
        }
    } else {
        selected_outputs <- c("species", "year", "ssp", "range_mean_suitability", "range_size", "range_centroid_x", "range_centroid_y", "eval_runtime")
        cleaned_output <- output[selected_outputs]
    }
    return(cleaned_output)
}

process_species <- function(species, year, ssp = NULL) {
    start_time <- Sys.time()
    set.seed(0)
    if (year == "2020") {
        print(paste("2020 |", species))
        env_data <- env_data_rename(environment_data$"2020")[[current_selected_vars]]
        # Proper model
        obs_data <- odonata_data[odonata_data$"species" == species, c("decimalLatitude", "decimalLongitude")]
        background <- spatSample(x=env_data, size=8000, values=FALSE, na.rm=TRUE, xy=TRUE)
        # Quick tests with smaller sample size
        # obs_data <- head(odonata_data[odonata_data$"species" == species, c("decimalLatitude", "decimalLongitude")], 10)
        # background <- spatSample(x=env_data, size=30, values=FALSE, na.rm=TRUE, xy=TRUE)
        
        presence <- obs_data
        presence$pa <- 1
        absence <- as.data.frame(background)
        colnames(absence) <- c("decimalLongitude", "decimalLatitude")
        absence$pa <- 0
        points_extract <- rbind(presence, absence)
        print("Presence and absence points extracted...")
        print(Sys.time() - start_time)
        
        environment_extract <- extract(x=env_data, y=points_extract[, c("decimalLongitude", "decimalLatitude")], ID=FALSE)
        points_environment_extract <- cbind(points_extract, environment_extract)
        drop_cols <- which(colnames(points_environment_extract) %in% c("decimalLongitude", "decimalLatitude"))
        points_environment_extract <- points_environment_extract[, -drop_cols]
        fold <- folds(x=points_environment_extract, k=5, by=points_environment_extract$pa)
        testing <- points_environment_extract[fold == 1, ]
        training <- points_environment_extract[fold != 1, ]
        print("Environment data added...")
        print(Sys.time() - start_time)
        
        maxent_model <- maxent(x=training[, colnames(training)[colnames(training) != 'pa']], p=training$pa, removeDuplicates=FALSE, path=paste0("output/maxent_files/", species, ifelse(initial_run, "/initial_output", "/final_output")))
        current_model <<- maxent_model
        print("MaxEnt model constructed...")
        print(Sys.time() - start_time)
        
        if (initial_run) {
            output_model <- maxent_model
            output_predict <- NULL
            output_eval <- NULL
            output_threshold <- NULL
            eval_var_importance <- maxent_model@results[grepl("permutation.importance", rownames(maxent_model@results)), ]
            eval_var_selected <- gsub("\\.permutation\\.importance$", "", rownames(maxent_model@results[grepl("permutation.importance", rownames(maxent_model@results)) & maxent_model@results >= 2.5, , drop = FALSE]))
            eval_training_auc <- NULL
            eval_testing_auc <- NULL
            binary_prediction <- NULL
            range_mean_suitability <- NULL
            range_size <- NULL
            range_centroid_x <- NULL
            range_centroid_y <- NULL
        } else {
            maxent_predict <- predict(env_data, maxent_model, na.rm=TRUE, type="response")
            maxent_eval <- pa_evaluate(p=testing[testing$pa == 1, ], a=testing[testing$pa == 0, ], model=maxent_model, type="response")
            current_threshold <<- maxent_eval@thresholds$max_spec_sens
            print("MaxEnt model fitted...")
            print(Sys.time() - start_time)
            
            output_model <- maxent_model
            output_predict <- maxent_predict
            output_eval <- maxent_eval
            output_threshold <- current_threshold
            eval_var_importance <- maxent_model@results[grepl("permutation.importance", rownames(maxent_model@results)), ]
            eval_var_selected <- NULL
            eval_training_auc <- maxent_model@results[5]
            eval_testing_auc <- maxent_eval@stats$auc
            eval_plots(maxent_eval@tr_stats, species)
            
            binary_prediction <- maxent_predict > current_threshold
            range_mean_suitability <- mean(maxent_predict[], na.rm=TRUE)
            range_size <- sum(binary_prediction[], na.rm=TRUE)
            range_centroid <- colMeans(xyFromCell(binary_prediction, which(binary_prediction[], arr.ind=FALSE)))
            range_centroid_x <- range_centroid["x"]
            range_centroid_y <- range_centroid["y"]
        }
    } else {
        print(paste(year, "|", ssp, "|", species))
        env_data <- env_data_rename(environment_data[[year]][[ssp]])[[current_selected_vars]]
        maxent_predict <- predict(env_data, current_model, na.rm=TRUE, type="response")
        print("MaxEnt model fitted...")
        print(Sys.time() - start_time)
        
        output_model <- current_model
        output_predict <- maxent_predict
        output_eval <- NULL
        output_threshold <- current_threshold
        eval_var_importance <- NULL
        eval_var_selected <- NULL
        eval_training_auc <- NULL
        eval_testing_auc <- NULL
        
        binary_prediction <- maxent_predict > current_threshold
        range_mean_suitability <- mean(maxent_predict[], na.rm=TRUE)
        range_size <- sum(binary_prediction[], na.rm=TRUE)
        range_centroid <- colMeans(xyFromCell(binary_prediction, which(binary_prediction[], arr.ind=FALSE)))
        range_centroid_x <- range_centroid["x"]
        range_centroid_y <- range_centroid["y"]
    }

    end_time <- Sys.time()
    eval_runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
    print(paste("Completed!"))
   
    return(list(
        species = species,
        year = year,
        ssp = ssp,
        output_model = output_model,
        output_predict = output_predict,
        output_eval = output_eval,
        output_threshold = output_threshold,
        range_mean_suitability = range_mean_suitability,
        range_size = range_size,
        range_centroid_x = range_centroid_x,
        range_centroid_y = range_centroid_y,
        eval_var_importance = eval_var_importance,
        eval_var_selected = eval_var_selected,
        eval_training_auc = eval_training_auc,
        eval_testing_auc = eval_testing_auc,
        eval_runtime = eval_runtime
    ))
}

for (species in selected_species[selected_species_index]) {
    if (initial_run) {
        current_selected_vars <- names(env_data_rename(environment_data$"2020"))
        selected_years <- c("2020")
    } else {
        current_selected_vars <- fromJSON(paste0("output/selected_vars/", species, ".json"))
    }
    for (year in selected_years) {
        if (year == "2020") {
            result <- process_species(species, "2020")
            if (initial_run) {
                # dir.create(path=paste0("output/selected_vars"), recursive=TRUE)
                dir.create(paste0("output/maxent_files/", species, "/results"), recursive=TRUE)
                save(current_model, file=paste0("output/maxent_files/", species, "/initial_output/model.RData"))
                write_json(result$eval_var_selected, path=paste0("output/selected_vars/", species, ".json"))
            } else {
                save(current_model, file=paste0("output/maxent_files/", species, "/final_output/model.RData"))
                write_json(clean_output_json(result), path=paste0("output/maxent_files/", species, "/results/2020_final.json"), null="null", auto_unbox=TRUE, pretty=TRUE)
                writeRaster(result$output_predict, paste0("output/maxent_files/", species, "/results/2020_final.tif"), overwrite=TRUE)
            }
        } else {
            for (ssp in selected_ssps) {
                result <- process_species(species, year, ssp)
                write_json(clean_output_json(result), path=paste0("output/maxent_files/", species, "/results/", year, "_", ssp, ".json"), null="null", auto_unbox=TRUE, pretty=TRUE)
                writeRaster(result$output_predict, paste0("output/maxent_files/", species, "/results/", year, "_", ssp, ".tif"), overwrite=TRUE)
            }
        }
    }
}
