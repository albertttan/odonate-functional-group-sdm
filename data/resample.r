# Crop & resample all raster files to North America at 1km resolution
# All input rasters should be placed in the folder `input`

library(terra)
library(fs)
library(tools)

process_file <- function(file_path, output_base_path) {
  cat("Processing file:", file_path, "\n")
  data <- rast(file_path)
  background <- rast("blank.tif")

  cat("Cropping data...\n")
  data <- crop(x = data, y = ext(x = c(-180, -50, 5, 85)))
  cat("Resampling data...\n")
  data <- resample(data, background)

  relative_path <- path_rel(file_path, start = "input")
  output_path <- file.path(output_base_path, path_ext_set(relative_path, "tif"))
  dir_create(path_dir(output_path))

  cat("Saving file to:", output_path, "\n")
  writeRaster(data, output_path)
  cat("File processed successfully!\n")
}

process_files_recursively <- function(input_dir, output_dir) {
  cat("Scanning for .tif and .nc files in:", input_dir, "\n")
  files <- dir_ls(input_dir, recurse = TRUE, regexp = "\\.(tif|nc)$")
  for (file in files) {
    process_file(file, output_dir)
  }
}

input_dir <- "input"
output_dir <- "output"

cat("Starting processing...\n")
process_files_recursively(input_dir, output_dir)
cat("Processing completed!\n")
