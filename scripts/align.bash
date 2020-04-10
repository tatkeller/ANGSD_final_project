#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --job-name=align
#SBATCH --mem=40G

spack load star@2.7.0e

FOLDERNAME=`ls -d ../fastq/*`
FILESDONE=`ls ../alignments/*`

for file in $FOLDERNAME; do
 if [[ ! -d ${file} ]]
 then
  continue
 fi
 name=`grep -oP '(?<=\/)[Acute_Survior|Acute_Fatal]*[0-9]*' <<< ${file}` 
 if [[ -f ../alignments/${name}.Aligned.sortedByCoord.out.bam ]]
 then
  continue
 fi
 if [[ ${file} == "../fastq/multiqc_data"  ]]
 then
  continue
 fi
 echo "${name}"
 STAR --runMode alignReads \
      --runThreadN 8 \
      --genomeDir /home/tak76/tak76/finalProject/extra/hg38index \
      --readFilesIn /home/tak76/tak76/finalProject/fastq/${name}/*.fastq.gz \
      --readFilesCommand zcat \
      --outFileNamePrefix /home/tak76/tak76/finalProject/alignments/${name}. \
      --outFilterMultimapNmax 1 \
      --outSAMtype BAM SortedByCoordinate \
      --twopassMode Basic \
      --alignIntronMin 1 \
      --alignIntronMax 100000
done
