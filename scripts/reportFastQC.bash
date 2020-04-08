#!/bin/bash

# SBATCH --partition=angsd_class
# SBATCH --nodes=1
# SBATCH --ntasks=1
# SBATCH --job-name=fastqRun
# SBATCH --mem=40G

spack load fastqc

for folders in ../rawFastq/*/; do
	files=`ls -d $PWD/${folders}/*`
	#filenames=`ls ${folders}`
	#echo $filesnames
	#echo $files	
	fastqc ${files}
done

