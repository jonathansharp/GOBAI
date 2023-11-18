#!/bin/bash

#SBATCH -J gobai
#SBATCH -A hindcasts
#SBATCH -q batch
#SBATCH -p hercules
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 40
#SBATCH -t 8:00:00
#SBATCH -o gobai_o2_kfold_train_gbm.out
#SBATCH -e gobai_o2_kfold_train_gbm.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jonathan.sharp@noaa.gov

module load matlab
cd /work2/noaa/hindcasts/GOBAI-O2-New_JS
srun matlab -nodisplay -r "create_config_files; load('Config/load_data_config_D.mat'); load('Config/cluster_config_20.mat'); load('Config/kFold_config_20_5.mat'); load('Config/gbm_config_400.mat'); gobai_o2_kfold_train_gbm;"
