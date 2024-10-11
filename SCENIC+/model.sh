#!/bin/bash

#SBATCH --job-name=pycisTopic_job     # Job name
#SBATCH --output=pycisTopic_%j.log      # Output log file
#SBATCH --error=pycisTopic_%j.err       # Error log file
#SBATCH --ntasks=1                       # Number of tasks (default is 1)
#SBATCH --cpus-per-task=10               # Number of CPU cores per task
#SBATCH --mem=250G                       # Total memory limit
#SBATCH --time=3-00:00:00
#SBATCH --partition=medium

# Load necessary modules (if applicable)
module load Python
cd /home/twyang2/
# Run your Python script
python model.py