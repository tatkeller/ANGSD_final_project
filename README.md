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

## Report

# Introduction

Ebola Virus Disease (EBOV) is a disease caused by the ebola virus, causing fatalities in an average of 60% of patients (Kadanali, 2015). Recent outbreaks in West Africa and Congo have resulted in great stress on the health care systems in these regions, and led to international concern over the virus's potential to spread globally. Previous EBOV studies have revealed a correlation between host immune response and survival, and genomic analyses of RNA expression in EBOV patients have tracked viral evolution and revealed a correlation between co-infections and mortality in EBOV patients (Reynard et al., 2019; Holmes et al., 2016; Carroll et al, 2017). Further analyses of RNA expression in EBOV patients could reveal genomic advantages in survivors of the disease, and can be used to develop treatments to increase the chance of survival amongst EBOV patients. Here, RNA sequence data from victims and survivors of EBOV was examined to determine if there exists significant gene expression differences between the two groups. 

# Results

In 2017, a dataset containing RNA-seq data from survivors and victims of the 2014 West African Ebola outbreak was published (ENA: PRJNA352396). Here, 110 random samples from the data set, 55 from survivors and 55 from victims, were used to analyze the difference in transcriptome between survivors and victims. Results show 8 genes that have significant differences in expression levels between survivors and victims. 

Overall, survivors have higher expressions of genes that are not obvious for their role in treatment against ebola, and these genes (PRSS8, CKB, UCA1) should be studied more in the context of ebola infection to conclude their role. Victims have higher expressions of genes involved in immune response, which is very obvious for the role in treatment against ebola. The outstanding question is whether or not this increase in expression is: 1) a sign of a poor immune response, and patients who express these genes less are more likely survive, 2) a sign of a proportional response to a particularly critical infection of ebola, or 3) a result of the ebola virus interferring with the immune system function. Prior research has suggested that a measured immune response to ebola results in a higher chance of survival, so it is possible that the genes identified here (ORM1, SAA1, CCN1, C7, IL6) are good targets for possible immune supression treatment (Reynard et al., 2019). 

Ultimately this dataset is limited mostly in its size. In order to find a strong correlation between survival and gene expression, millions of samples would need to be collected - but hopefully there will never be a scenario in which that many patients exist. A better way to track this sort of data would be to include other biometrics, co-morbidities, symptoms, days since infection, and treatments that the patient recieved while infected with ebola. Future studies should focus on the genes identified here in the context of other external data.

# Methods and Figures

A complete repository of workflow and code can be found at: https://github.com/tatkeller/ANGSD_final_project.

The dataset was downloaded from the European Nucleotide Archive, under project PRJNA352396 (https://www.ebi.ac.uk/ena/browser/view/PRJNA352396). Fastq files were downloaded to the Weill Cornell server. The files were downloaded using a shell script (https://github.com/tatkeller/ANGSD_final_project/blob/master/scripts/downloadFastq.sh). Next, fastqc files were created for each sample using FastQC version 6.3.0 (Andrews, 2010). Next, each sample was aligned to the Human Genome 38 (http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/). STAR was first used to index the Human Genome 38 data, and then used to align each sample (STAR version 2.7.0e, Dobin et al, 2013). Finally, the transcription level of each gene was counted using the featureCounts function from the subread package (version 6.3.0, Liao et al, 2014). Finally, to analyze the quality of the reads, a MultiQC report was created for all samples' FastQC, Alignments, and Feature counts (version 1.8, Philip et al, 2016. Each of these functions were called using various shell scripts which can be found at: https://github.com/tatkeller/ANGSD_final_project/tree/master/scripts. The multiqc report can be found on the github repository. Some samples have bacterial content that result in low quality scores in the multiqc report, but all samples were left in since it is not certain which samples have the bacterial content.

The featureCounts were loaded into R-markdown.

```{r message=FALSE}
library(magrittr)
library(data.table)
genesFull <- data.table::fread(file = '../featureCounts/featureCounts.txt', header = TRUE)%>% as.data.frame
originalNames <- names(genesFull)
names(genesFull) <- gsub(".*(Acute_Fatal|Acute_Survivor)(_[0-9]+).*", "\\1\\2", originalNames)
newNames = names(genesFull)
```

DESeq2 library was used to analyze gene expression for each sample.

```{r message=FALSE}
library(DESeq2)
row.names(genesFull) <- make.names(genesFull$Geneid)
genesFull <- genesFull[ , -c(1:6)]
sample_info <- DataFrame(condition = gsub("_[0-9]+", "", names(genesFull)),row.names = names(genesFull) )
DESeq.ds <- DESeqDataSetFromMatrix(countData = genesFull,colData = sample_info,design = ~ condition)
```

The DESeq object had genes with zero counts removed.

```{r}
keep_genes <- rowSums(counts(DESeq.ds)) > 0
DESeq.ds <- DESeq.ds[ keep_genes, ]
```

Then the object was normalized, and a heatmap of correlations between the classes was created:

```{r fig.width=6, fig.height=6}
DESeq.vst <- vst(DESeq.ds, blind = TRUE)
vst.norm.counts <- assay(DESeq.vst)
vst.norm.counts = vst.norm.counts[ , order(colnames(vst.norm.counts))]
corr_coeff <- cor(vst.norm.counts, method = "pearson")
as.dist(1-corr_coeff, upper = TRUE) %>% as.matrix %>% pheatmap::pheatmap(., main = "Pearson correlation")
```

There was no correlation found between the two groups for the whole genome. This is mostly to be expected, since they are random patients, and aside from genes that have some role in outcome of survival, no difference in genome should be found.

The magnitude and significance of differential gene expression between survivors and victims was calculated for each gene. DESeq calculates p-values of analyses by a Wald statistic, and then uses a correction tool, the Benjamini-Hochberg formula, to adjust the pvalues. 

```{r message=FALSE}
#Perform DGE tests
DESeq.ds$condition <- relevel(DESeq.ds$condition, ref="Acute_Survivor")
DESeq.ds <- DESeq(DESeq.ds)
DESeq.ds <- estimateSizeFactors(DESeq.ds)
DESeq.ds <- estimateDispersions(DESeq.ds)
DESeq.ds <- nbinomWaldTest(DESeq.ds)
DGE.results <- results(DESeq.ds, independentFiltering = TRUE, alpha = 0.05)

#Plot histogram of Adjusted P values
hist(DGE.results$padj, col="grey", border="white", xlab="", ylab="", main="frequencies of padj\n(all genes)", cex = 0.4)
DGE.results.sorted <- DGE.results[order(DGE.results$padj),]
```

Few genes score with a low adjusted p value, meaning that few genes are significantly differentially spread (beyond chance). 

An enhanced volcano plot with the threshold for adjusted p-value equal to 0.05 and threshold for log fold change equal to 2.5, reveals 8 genes of interest.

```{r}
library(EnhancedVolcano)
vp1 <- EnhancedVolcano(DGE.results, lab = rownames(DGE.results), x = 'log2FoldChange', y = 'padj', pCutoff = 0.05, FCcutoff = 2.5, title = "Fatal / Survivor", ylim = c(0,7), xlim = c(-5,5))
print(vp1)
```

A heatmap of z-score (standard score) for each gene per sample shows a less random grouping of patients. The bottom left of the heatmap shows a small group of survivors who have extremely high expression of UCA1, PRSS8, IL6, and CKB genes. The right side of the heatmap shows a group of majority victims who have extremely high expression of ORM1, SAA1, CCN1, and C7 genes. This does not show the whole story for both classes though so each gene was further analyzed.

```{r fig.width=7, fig.height=3}
library(pheatmap)
DGEgenes <- rownames(subset(subset(DGE.results.sorted, padj < 0.05), abs(log2FoldChange) > 2.5 ))
vst.dge <- DESeq.vst[DGEgenes,] %>% assay
pheatmap(vst.dge, scale="row", show_rownames = TRUE, main = "DGE (row-based z-score)")
```

A closer look at the distribution of counts for each gene makes it easier to discover possible explanations.

```{r}
# More Survival
plotCounts(DESeq.ds, gene="PRSS8", normalized = TRUE) 

```

For PRSS8, only 3 of the 55 survivors showed a relatively high expression level. PRSS8 plays a role in sodium channel regulation, and it is possible that very high expression of this gene can increase chances of survival, but more examples would need to support this before finding a possible explanation (Koda, 2009).


```{r}
plotCounts(DESeq.ds, gene="CKB", normalized = TRUE) 
```

Similarly, only 3 survivors of 55 had a relatively high expression level. CKB is a creatine kinase that is responsible for energy catalysis in muscles, the brain, and the heart (Mariman, 1989). It is possible that the increased ability to replenish ATP in important organs and muscles could increase chances of survival, but more testing would need to be done to prove this.

```{r}
plotCounts(DESeq.ds, gene="UCA1", normalized = TRUE)
```

UCA1 is a gene that is involved in cell proliferation, and is often upregulated in bladder cancer (Xiao-Song, 2006). Interestingly, survivors had higher expression of this gene than victims did.  

```{r}
plotCounts(DESeq.ds, gene="ORM1", normalized = TRUE) 
```

ORM1 functions in modulating the immune system response to acute diseases (NCBI, 2020). Victims of EBOV show a much higher increase in expression of this gene as compared to survivors. However, this does not necesarily suggest that lower expression of the gene is correlated with survival, since the immune system may be responding dynamically to an infection, and patients who are struggling to combat the infection may require a higher immune response. That being said, EBOV is known to attack the immune system in humans, and it is possible that this is one gene that increases in expression directly as a result of ebola infection rather than as a response to ebola infection. 

```{r}
plotCounts(DESeq.ds, gene="IL6", normalized = TRUE) 
```

IL6 is another gene involved in the immune system, and codes for cytokines (Tanaka, 2014). Cytokine storms are often one of the causes of death amongst ebola victims, and depending on when these samples were taken, cytokine levels may be significantly higher in patients in critical condition (Reynard, 2019). As with ORM1, it is hard to conclude whether or not the increase in IL6 amongst fatal patients is a poor immune system strategy in reaction to an ebola infection, or if it is a natural result of a particularly critical case of an ebola infection. 

```{r}
plotCounts(DESeq.ds, gene="C7", normalized = TRUE) 
```

C7 is responsible for cell lysis attacks from the host immune system, typically against bacterial infections (Manuel, 1991). Previous studies have linked undiagnosed bacterial infections as a co-morbidity for EBOV, so it is possible that patients with an increased level of C7 are also fighting a bacterial infection (Reynard, 2019). Like IL6 and ORM1, it is hard to know whether this increase in expression is a poor reaction to ebola that should be supressed, or if it is a result of a critical case of the disease.

```{r}
plotCounts(DESeq.ds, gene="SAA1", normalized = TRUE)
```

SAA1 is another protein involved in immune response, and for the same reasoning as the above proteins, more data and studies must be conducted (Husebekk, 1985).


```{r}
plotCounts(DESeq.ds, gene="CCN1", normalized = TRUE) 
```

CCN1 has a wide range of functions in the human body including proliferation and apoptosis, and could be increased in expression as a result of the body reacting to rapid cell death (Lau, 2011).

# Discussion of Limitations

Results show that on average, survivors and victims do contain slight differences in gene expressions, especially amongst genes responsible for immune response. However, much more data is needed to understand the nature of these findings. Additionally, the dataset is limited in that it only reveals the patient outcome, and does not reveal other potentially important outcomes for survival such as age, gender, treatments, etc. 

Some limitations of this study include resources for alignment of human genome samples. Fastq files for one sample were as large as 10GB in some cases, and for greater than 100 samples, that is not practical. Samples were aligned in batches to avoid the space issue. Given more time, the entire dataset of 174 samples would have been analyzed. 

Further downstream analyses could be done to further the support of the significance of the gene expressions identified here. A GO analysis was attempted, but HG38 is not supported by R packages such as goseq. 

# Key Data Sets

The key data sets generated are:

Alignments: all *.bam files saved to the alignments folder on the Weill Cornell server. 

Feature Counts: The feature counts table detailing the counts of frequency of counts of each gene. This is featureCounts.txt. "genesFull" is and R data table of this file.

DESeq Object: DESeq.ds is the DESeq object, and DESeq.vst is the normalized counterpart.

Genes of interest: DGEgenes is the object of gene expression analysis with an adjusted pvalue of less than 0.05 and absolute value log2 Fold Change of greater than 2.5.

# References

Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc

Carroll, Miles W et al. (2017) “Deep Sequencing of RNA from Blood and Oral Swab Samples Reveals the Presence of Nucleic Acid from a Number of Pathogens in Patients with Acute Ebola Virus Disease and Is Consistent with Bacterial Translocation across the Gut.” mSphere vol. 2,4 e00325-17. doi:10.1128/mSphereDirect.00325-17

Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., & Gingeras, T. R. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics (Oxford, England), 29(1), 15–21. https://doi.org/10.1093/bioinformatics/bts635

Holmes, Edward C et al. (2016) “The evolution of Ebola virus: Insights from the 2013-2016 epidemic.” Nature vol. 538,7624: 193-200. doi:10.1038/nature19790

Husebekk, A., SKOGEN, B., HUSBY, G. and MARHAUG, G. (1985), Transformation of Amyloid Precursor SAA to Protein AA and Incorporation in Amyloid Fibrils in Vivo. Scandinavian Journal of Immunology, 21: 283-287. doi:10.1111/j.1365-3083.1985.tb01431.x

Kadanali, A., & Karagoz, G. (2015). An overview of Ebola virus disease. Northern clinics of Istanbul, 2(1), 81–86. https://doi.org/10.14744/nci.2015.97269

Koda, A., Wakida, N., Toriyama, K. et al. Urinary prostasin in humans: relationships among prostasin, aldosterone and epithelial sodium channel activity. Hypertens Res 32, 276–281 (2009). https://doi.org/10.1038/hr.2009.6

Lau L. F. (2011). CCN1/CYR61: the very model of a modern matricellular protein. Cellular and molecular life sciences : CMLS, 68(19), 3149–3163. https://doi.org/10.1007/s00018-011-0778-3

Liao Y, Smyth GK and Shi W (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-30. http://www.ncbi.nlm.nih.gov/pubmed/24227677

Manuel C. Peitsch, Jürg Tschopp, (1991) Assembly of macromolecular pores by immune defense systems, Current Opinion in Cell Biology, 3(4) 710-716, https://doi.org/10.1016/0955-0674(91)90045-Z.

Mariman, E. C., Schepens, J. T., & Wieringa, B. (1989). Complete nucleotide sequence of the human creatine kinase B gene. Nucleic acids research, 17(15), 6385. https://doi.org/10.1093/nar/17.15.6385

NCBI (2020). ORM1. URL: https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=ShowDetailView&TermToSearch=5004

Philip Ewels, Måns Magnusson, Sverker Lundin and Max Käller (2016) MultiQC: Summarize analysis results for multiple tools and samples in a single report, Bioinformatics doi: 10.1093/bioinformatics/btw354 PMID: 27312411

Reynard, S., Journeaux, A., Gloaguen, E., Schaeffer, J., Varet, H., Pietrosemoli, N., Mateo, M., Baillet, N., Laouenan, C., Raoul, H., Mullaert, J., & Baize, S. (2019). Immune parameters and outcomes during Ebola virus disease. JCI insight, 4(1), e125106. Advance online publication. https://doi.org/10.1172/jci.insight.125106

Tanaka, Toshio et al. (2014) “IL-6 in inflammation, immunity, and disease.” Cold Spring Harbor perspectives in biology vol. 6,10 a016295. doi:10.1101/cshperspect.a016295

Xiao-Song Wang, Zheng Zhang, Hong-Cheng Wang, Jian-Liang Cai, Qing-Wen Xu, Meng-Qiang Li, Yi-Cheng Chen, Xiao-Ping Qian, Tian-Jing Lu, Li-Zhang Yu, Yu Zhang, Dian-Qi Xin, Yan-Qun Na and Wei-Feng Chen (2006) Rapid Identification of UCA1 as a Very Sensitive and Specific Unique Marker for Human Bladder Carcinoma.Clin Cancer Res (12) (16) 4851-4858; DOI: 10.1158/1078-0432.CCR-06-0134
