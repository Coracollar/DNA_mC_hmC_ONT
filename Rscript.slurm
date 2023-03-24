#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --time=0-40:00:00
#SBATCH --mem=500G

source ~/anaconda3/etc/profile.d/conda.sh

conda activate r_env

Rscript DSS.R
