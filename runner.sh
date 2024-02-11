#!/bin/bash

# Define the input parameters
ITER_SIZE=100
X="..."
BETA="..."
SIGMA_U="..."
GROUP_SIZES="..."
FITFAM="..."
FORMULA="..."
FAMILY="..."
# Add other parameters as needed

# Run the R script
Rscript simmulate.R "$ITER_SIZE" "$X" "$BETA" "$SIGMA_U" "$GROUP_SIZES" "$FITFAM" "$FORMULA" "$FAMILY"
