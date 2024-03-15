#!/bin/bash

#SBATCH -J gobai-create
#SBATCH -A hindcasts
#SBATCH -q batch
#SBATCH -p hercules
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 40
#SBATCH --exclusive
#SBATCH -t 8:00:00
#SBATCH -o output/%x-%j.out
#SBATCH -e errors/%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jonathan.sharp@noaa.gov

source /home/sharp/shell_scripts/dir.txt
cd $GOBAI_DIR
module load matlab
matlab -nodisplay -r "gobai_o2_initiate; load_standard_config_files; load('Config/base_config_RFROM.mat'); load('Config/predict_years_config_19.mat'); gobai_o2_predict_ffnn;"
