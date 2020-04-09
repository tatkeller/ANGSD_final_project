#!/bin/bash




spack load star@2.7.0e

FILES="SRR4888615
SRR4888616
SRR4888618"

for file in $FILES; do
 STAR --runMode alignReads \
      --runThreadN 4 \
      --genomeDir /home/tak76/tak76/finalProject/hg38index \
      --readFilesIn /home/tak76/tak76/finalProject/fastq/${file}/* \
      --readFilesCommand zcat \
      --outFileNamePrefix /home/tak76/tak76/finalProject/alignments/${file}. \
      --outFilterMultimapNmax 1 \
      --outSAMtype BAM SortedByCoordinate \
      --twopassMode Basic \
      --alignIntronMin 1 \
      --alignIntronMax 3000
done
