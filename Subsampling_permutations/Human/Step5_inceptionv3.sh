#!/bin/bash
# run_extract_H.slurm
# SLURM script to extract InceptionV3 features for Human permutations (GPU jobs)
#SBATCH --job-name=extractH
#SBATCH --array=1-10
#SBATCH --gres=gpu:2            # request one GPU per task
#SBATCH --mem=32G               # memory per job
#SBATCH --partition=shared-gpu
#SBATCH --time=04:00:00
#SBATCH --output=logs/extractH_%a.out
#SBATCH --error=logs/extractH_%a.err

#micromamba activate inception
cd /home/users/j/javed/humous_reviews/
python Step5_inceptionv3.py

