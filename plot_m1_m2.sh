#!/bin/bash

#SBATCH -J plotm1m2
#SBATCH --mail-user=rczhang@u.northwestern.edu
#SBATCH --error=ploterr.out
#SBATCH --output=plotoutput.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem=50G
#SBATCH --time=10:00:00
#SBATCH --account=b1094
#SBATCH --partition=ciera-std

source activate cosmic
python plot_m1_m2.py
