#!/bin/bash

reportFile=report.json

batches=$(grep -P '{' ${reportFile})

IFS=$'\n'

surviveCount=0
fatalCount=0

for batch in ${batches}; do
	title=$(grep -oP '((?<=: )|(?<=:))[Acute] ?.*(?=",)' <<< "$batch")
	class=$(grep -oP '(.*(?=_[1-9]+))' <<< "$title")
        file1=$(grep -oP '(?<=")ftp.*(?=;)' <<< "$batch")
        file2=$(grep -oP '(?<=;)ftp.*(?=")' <<< "$batch")
	if [[ -d ../fastq/$title ]]
	then
		continue
	fi
	if [[ $class == "Acute_Survivor" ]]
	then
		if [[ $surviveCount -lt 20  ]]
		then
			mkdir ../fastq/${title}
			wget -O ../fastq/${title}/${title}_1.fastq.gz ftp://${file1}
                        wget -O ../fastq/${title}/${title}_2.fastq.gz ftp://${file2}
			surviveCount=$((surviveCount+1))
		fi
	fi
	if [[ $class == "Acute_Fatal" ]]
	then
                if [[ $fatalCount -lt 20  ]]
                then
			mkdir ../fastq/${title}
                        wget -O ../fastq/${title}/${title}_1.fastq.gz ftp://${file1}
			wget -O ../fastq/${title}/${title}_2.fastq.gz ftp://${file2}
                        fatalCount=$((fatalCount+1))
                fi
	fi
done
