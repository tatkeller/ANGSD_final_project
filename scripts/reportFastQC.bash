#!/bin/bash

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=fastq
#SBATCH --mem=20G

spack load fastqc

for file in ../fastq/*; do
	#if it is a folder
	if [[ -d $file ]]
	then 
		#get the files inside the folder
		insidefiles=`ls -d $PWD/${file}/*`
		count=`ls -1 ${file}/*.zip 2>/dev/null | wc -l`
		#ignore the multiqc
		if [[ $file == "../fastq/multiqc_data" ]]
		then
			continue
		fi
		#ignore ones that have been done already
		if [[ $count != 0 ]]
		then
			continue
		fi
		fastqc ${insidefiles}
	fi
done

