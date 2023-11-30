#!/bin/bash

#SBATCH -J gobai-kfold
#SBATCH -A hindcasts
#SBATCH -q batch
#SBATCH -p hercules
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 40
#SBATCH -t 8:00:00
#SBATCH -o output/gobai_o2_kfold_train_ffnn.out
#SBATCH -e errors/gobai_o2_kfold_train_ffnn.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jonathan.sharp@noaa.gov

module load matlab
cd /work2/noaa/hindcasts/GOBAI-O2-New_JS
srun matlab -nodisplay -r "create_config_files; load('Config/load_data_config_A_D.mat'); load('Config/cluster_config_10.mat'); load('Config/kFold_config_10_5.mat'); load('Config/ffnn_config_70_15_15.mat'); gobai_o2_kfold_train_ffnn;"
