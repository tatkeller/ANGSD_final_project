#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=multiqc
#SBATCH --mem=10G

# Place this script in the fastq folder to generate a multiqc report

spack load -r py-multiqc

multiqc . 



