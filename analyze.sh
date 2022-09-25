#!/bin/bash

#SBATCH -J analyze_sse
#SBATCH --mail-user=rczhang@u.northwestern.edu
#SBATCH --error=analyzeerr.out
#SBATCH --output=analyzeoutput.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem=50G
#SBATCH --time=10:00:00
#SBATCH --account=b1095
#SBATCH --partition=grail-std

source activate cosmic
python analyze_sse.py

