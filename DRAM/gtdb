#!/bin/bash
#SBATCH --nodes=1                  # Number of requested nodes
#SBATCH --ntasks=10                 # Cores per node requested
#SBATCH --time=7-50:00:00        	# Timelimit
#SBATCH --mem=100gb                # Job memory (1024gb=1tb)
#SBATCH --mail-type=BEGIN,END,FAIL         
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

cd /home/projects/Sporulation_AMG/hmp/VirHostMatcher/

gtdbtk classify_wf --extension fna --genome_dir /home/projects/Sporulation_AMG/hmp/VirHostMatcher/host --out_dir /home/projects/Sporulation_AMG/hmp/VirHostMatcher/GTDB_1.5.0_jrr
