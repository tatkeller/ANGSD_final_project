#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=multiqc
#SBATCH --mem=10G

# spack load -r py-multiqc

spack load singularity@2.6.0

singularity exec /athena/angsd/scratch/simg/multiqc-1.8.simg multiqc .

