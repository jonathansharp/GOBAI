#!/bin/bash

#SBATCH -J gobai-kfold
#SBATCH -A hindcasts
#SBATCH -q batch
#SBATCH -p hercules
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 40
#SBATCH -t 8:00:00
#SBATCH -o output/$(basename "$0" .sh).out
#SBATCH -e errors/$(basename "$0" .sh).err