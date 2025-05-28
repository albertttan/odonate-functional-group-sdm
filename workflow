#!/bin/sh

set -e

pip3 install numpy pandas matplotlib
Rscript -e "install.packages(c('stringr', 'reshape2', 'jsonlite', 'ggplot2', 'usdm', 'terra', 'dismo', 'predicts', 'geosphere', 'car', 'boot', 'broom', 'rlang', 'DHARMa', 'emmeans'), repos='https://cran.rstudio.com/')"

mkdir output

mkdir output/cleaning
python3 cleaning.py

mkdir output/background
Rscript background.R

mkdir output/multicollinearity
Rscript multicollinearity.R

mkdir output/model
Rscript model_initial.R
Rscript model_final.R

mkdir output/analysis
mkdir output/analysis/images
Rscript analysis.R
