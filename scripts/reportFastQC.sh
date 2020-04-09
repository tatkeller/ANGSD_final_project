#!/bin/bash

spack load fastqc

for folders in ../rawFastq/*/; do
	files=`ls -d $PWD/${folders}/*`
	fastqc ${files}
done

