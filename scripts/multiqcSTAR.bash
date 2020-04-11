#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=multiqc
#SBATCH --mem=10G

spack load -r py-multiqc

multiqc -o ../multiqc/starResults/  ../alignments/

