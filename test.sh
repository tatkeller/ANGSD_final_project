#! /bin/bash -l

spack load star@2.7.0e

FILENAME=`ls rawFastq`
echo $FILENAME

for file in $FILENAME; do
 STAR --runMode alignReads \
      --runThreadN 4 \
      --genomeDir /home/tak76/tak76/finalProject/hg38index \
      --readFilesIn /home/tak76/tak76/finalProject/rawFastq/${file}/*.fastq.gz \
      --readFilesCommand zcat \
      --outFileNamePrefix /home/tak76/tak76/finalProject/alignments/${file}. \
      --outFilterMultimapNmax 1 \
      --outSAMtype BAM SortedByCoordinate \
      --twopassMode Basic \
      --alignIntronMin 1 \
      --alignIntronMax 100000
done
