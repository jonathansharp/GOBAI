#!/bin/bash

#SBATCH -J gobai-kfold
#SBATCH -A hindcasts
#SBATCH -q batch
#SBATCH -p hercules
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 80
#SBATCH --exclusive
#SBATCH -t 8:00:00
#SBATCH -o output/%x-%j.out
#SBATCH -e errors/%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jonathan.sharp@noaa.gov

source /home/sharp/shell_scripts/dir.txt
cd $GOBAI_DIR
module load matlab
matlab -nodisplay -r "gobai_o2_initiate; load('Config/base_config_RFROM.mat'); load('Config/load_data_config_D.mat'); load('Config/cluster_config_15.mat'); load('Config/kfold_config_5.mat'); load('Config/ffnn_config_70_15_15.mat'); gobai_o2_kfold_train_ffnn;"
