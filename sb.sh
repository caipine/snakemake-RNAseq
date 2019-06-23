#!/bin/bash
#SBATCH -J R61            # job name
#SBATCH -o R61.o%j        # output and error file name (%j expands to jobID)
#SBATCH -N 1                # number of nodes requested
#SBATCH -n 48               # total number of mpi tasks requested
#SBATCH -p normal      # queue (partition) -- normal, development, etc.
#SBATCH -t 48:00:00         # run time (hh:mm:ss) - 1.5 hours

# Slurm email notifications are now working on Lonestar 5 
#SBATCH --mail-user=qcai1@mdanderson.org
#SBATCH --mail-type=begin   # email me when the job starts
#SBATCH --mail-type=end     # email me when the job finishes
#SBATCH -A PGK-inhibitors

cd /scratch/03988/qscai/RNAseq/RNAseq147/rna-seq-star-deseq2-master
bash
snakemake  --unlock
snakemake --cores 48 --rerun-incomplete          
