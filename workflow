#!/bin/sh

mkdir output

mkdir output/cleaning
python3 cleaning.py

mkdir output/background
Rscript background.R

mkdir output/multicollinearity
Rscript multicollinearity.R

mkdir output/model
Rscript model.R

mkdir output/analysis
mkdir output/analysis/images
Rscript analysis.R
