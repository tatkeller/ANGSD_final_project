#! /bin/bash

spack load star@2.7.0e

mkdir -p /scratchLocal/tak76

STAR --runMode genomeGenerate \
     --runThreadN 1 \
     --genomeDir  /home/tak76/tak76/finalProject/hg38index \
     --genomeFastaFiles /home/tak76/tak76/finalProject/hg38/hg38.fa \
     --sjdbGTFfile  /home/tak76/tak76/finalProject/hg38/hg38.refGene.gtf \
     --sjdbOverhang 124 \
     --outTmpDir /scratchLocal/tak76/hg38_genome

rm -rf /scratchLocal/tak76
