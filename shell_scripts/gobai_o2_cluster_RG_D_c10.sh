#!/bin/bash

#SBATCH -J gobai
#SBATCH -A hindcasts
#SBATCH -q batch
#SBATCH -p hercules
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 40
#SBATCH -t 8:00:00
#SBATCH -o gobai_o2_cluster.out
#SBATCH -e gobai_o2_cluster.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jonathan.sharp@noaa.gov

module load matlab
cd /work2/noaa/hindcasts/GOBAI-O2-New_JS
srun matlab -nodisplay -r "create_config_files; load('Config/load_data_config_D.mat'); load('Config/cluster_config_10.mat'); gobai_o2_cluster_rg;"