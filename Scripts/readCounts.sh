#! /bin/bash

spack load subread

featureCounts -a /home/tak76/tak76/finalProject/hg38/hg38.refGene.gtf \
              -o featureCounts.txt\
              /home/tak76/tak76/finalProject/alignments/*.bam       

