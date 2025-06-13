# Forecasting divergence: Climate-driven habitat shifts in North American Odonates depend on functional groups

Author: Yunchao Tan

Code by Yunchao Tan

This repository contains all data, scripts, and documentation necessary to reproduce the analyses and figures presented in our manuscript. The manuscript has been submitted to PeerJ. Upon acceptance, this README will be updated with the manuscript DOI.


## Table of Contents

- [Project Overview](#project-overview)
- [Data](#data)
- [Scripts](#scripts)
- [Dependencies](#dependencies)
- [Reproducibility](#reproducibility)


## Project Overview

My study predicts the impacts of different climate change scenarios on the future habitat and distribution of Odonates. I retrieved data from Odonata Central, WorldClim, EarthEnv, and other published data. I used MaxEnt to construct species distribution models (SDMs) for 30 North American Odonate species across seven functional groups. Each model was applied to three future years and three different Shared Socio-economic Pathways (SSPs).


## Data

All data is contained within the `Data` folder. After the original data is retrieved from its respective sources, the raster files are cropped to the spatial extent of North America and resampled into $1 \mathrm{km^2}$ grids.

### Odonates Presence Data

- `presence.csv`: Presence-only entries of North American Odonates recorded between 2015–2025 (Source: [GBIF](https://doi.org/10.15468/dl.rj2qh2))

### Environmental Variables

- `climate/`: Future predictions of 19 bioclimatic variables under the HadGEM3-GC31-LL model ([WorldClim](https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html))
- `landuse/`: Future predictions of land use (plant functional types) under the HadGEM2-ES model ([Chen et al.](https://doi.org/10.1038/s41597-020-00669-x))
- `population/`: Future predictions of human population density ([Wang et al.](https://doi.org/10.1038/s41597-022-01675-x))
- `topography/`: Terrain slope and roughness raster data ([EarthEnv](https://www.earthenv.org/topography))


## Scripts

All analysis scripts are contained within the main directory. Please execute shell script `workflow` to run the scripts in order.

- `cleaning.py`: Clean the presence-only entries data to select the relevant rows and columns
- `background.R`: Test different background points numbers to balance performance and prediction accuracy
- `multicollinearity.R`: Select variables to ensure low multicollinearity
- `model_initial.R`: Execute initial variable screening models
- `model_final.R`: Execute final species distribution models
- `analysis.R`: Data analysis on outputs from species distribution models
- `functions.R`: Underlying functions referenced by `background.R`, `multicollinearity.R`, `model_initial.R`, and `model_final.R`


## Dependencies

The analyses were conducted in Python (3.13.3) and R (version 4.4.3).

### Data Manipulation and Processing
- `numpy` (Python, version 2.2.6): Numerical computing
- `pandas` (Python, version 2.2.3): Data manipulation
- `stringr` (R, version 1.5.1): String manipulation
- `reshape2` (R, version 1.4.4): Data reshaping

### Data, Import, Export, and Visualization
- `json` (Python, version 2.0.9): JSON parsing and generation
- `jsonlite` (R, version 1.8.9): JSON parsing and generation
- `ggplot2` (R, version 3.5.2): Data visualization
- `matplotlib` (Python, version 3.10.0): Data visualization

### Spatial Analysis and Modeling
- `usdm` (R, version 2.1.7): Uncertainty analysis for species distribution models
- `dismo` (R, version 1.3.16): Species distribution modeling
- `terra` (R, version 1.8.5): Spatial data analysis
- `predicts` (R, version 0.1.17): Spatial distribution prediction
- `geosphere` (R, version 1.5.20): Spherical trignometry

### Statistical Analysis
- `car` (R, version 3.1.3): Companion to Applied Regression
- `boot` (R, version 1.3.31): Bootstrap methods
- `broom` (R, version 1.0.7): Converting statistical objects into tidy dataframes
- `rlang` (R, version 1.1.4): Functions for Base Types and Core R and 'Tidyverse' Features
- `DHARMa` (R, version 0.4.7): Residual diagnostics for hierarchical models
- `emmeans` (R, version 1.10.7): Estimated marginal means


## Reproducibility

To reproduce the analysis and figures:

1. Clone the repository using `git clone [repository URL]`;
2. Set your working directory to the cloned repository;
3. Execute the shell script `workflow`.

Results will be saved in the folder `output`.


## Contact

For questions about the code or data, please [contact Yunchao Tan](mailto:albert.yunchao.tan@gmail.com).
