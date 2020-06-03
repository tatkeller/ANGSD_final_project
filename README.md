# Determining a genetic difference between survivors and victims of Ebola virus disease from RNA-seq data

## Overview

This is a complete genomics project encompassing genomic alignment, feature counting, QC reporting, and differential expression analysis. The goal was to determine what genetic differences survivors and victims of the Ebola virus disease may have. The dataset used contained the RNA-seq data from survivors and victims of the 2014 West African Ebola outbreak. Fastq files were aligned to hg38. Differential gene analysis was conducted to uncover the genes that had the most significant differences in expression. This analysis revealed victims to have had higher expression levels of genes related to immune response, and ultimately shows how the dataset is limited in both size and temporal detail.    

### Tools used:

  #### Bash:
  - STAR
  - subread
  - FastQC
  - MultiQC
  #### R:
  - DESeq2

## Comments

This project was completed as part of the graduate course "Analysis of Next Generation Sequencing Data" at Weill Cornell Medical College. 

Course Website:
https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/schedule_2020/

Instructions:
https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/project_description_and_ideas.html

Dataset: 
https://www.ebi.ac.uk/ena/browser/view/PRJNA352396
