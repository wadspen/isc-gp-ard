#!/bin/bash

#SBATCH --partition=general
#SBATCH --time=11:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=120
#SBATCH --mem=500G

module purge
module load gsl/2.8 udunits/2.2.28-gcc14.2 cuda/11.6 freetype/2.12.1 gdal/3.9.2 r/4.4.2 proj geos cmake
#
#Rscript ./cyc_shuff_test.R
Rscript ./four_scramble_test.R
