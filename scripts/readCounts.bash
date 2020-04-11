#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=features
#SBATCH --mem=20G

spack load subread

featureCounts -a /home/tak76/tak76/finalProject/extra/hg38/hg38.refGene.gtf \
              -o /home/tak76/tak76/finalProject/featureCounts/featureCounts.txt\
              /home/tak76/tak76/finalProject/alignments/*.bam       

