#!/bin/bash

#SBATCH -J gobai-load
#SBATCH -A hindcasts
#SBATCH -q batch
#SBATCH -p hercules
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 40
#SBATCH -t 8:00:00
#SBATCH -o output/%x-%j.out
#SBATCH -e errors/%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jonathan.sharp@noaa.gov

module load matlab
cd /work2/noaa/hindcasts/GOBAI-O2-New_JS
srun matlab -nodisplay -r "gobai_o2_initiate; load('Config/load_data_config_A_D.mat'); gobai_o2_load;"

