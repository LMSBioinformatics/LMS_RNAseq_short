


Analysis of RNAseq data
========================================================
author:MRC LMS Bioinformatics Core
date: 27 September 2018
width: 1440
height: 1100
autosize: true
font-import: <link href='http://fonts.googleapis.com/css?family=Slabo+27px' rel='stylesheet' type='text/css'>
font-family: 'Slabo 27px', serif;
css:style.css


RNA-seq Introduction 
========================================================

RNA Sequencing (**RNA-seq**) is a well established technique which allows for studying expression profiles of thousands of genes/transcripts simultaneously. 

<div align="center">
<img src="figures/nrg2484-f1.jpg" alt="gene" height="700" width="600">
</div>

* Nature Reviews Genetics 10, 57-63 (January 2009)


RNA-seq analysis pipeline
====================================

<div align="center">
<img src="figures/RNAseqpipeline.png" alt="pipe" height="800" width="1027">
</div>


Contents
====================================
* [Quality Assessment](#/quality).
* [Read alignment](#/alignment).
* [Read feature counting](#/counting).
* [Differential gene expression analysis](#/de).
* [Exploring and saving results](#/explore).
* [Functional analysis - Gene Ontology](#/go).



Quality Assessment - FastQC (1/2)
====================================
id: quality

<div align="center">
<img src="figures/fastqc_1.png" alt="pipe" height="700" width="1100">
</div>


Quality Assessment - FastQC (2/2)
====================================

<div align="center">
<img src="figures/fastqc_2.png" alt="pipe" height="700" width="1100">
</div>


Read alignment - gapped aligner (1/2)
====================================
id: alignment

First step in RNA-Seq data analysis is to align the raw reads to genome/transcriptome. 

For RNASeq data alignment, aligner need to be able to align across the exon-exon junction.

<div align="center">
<img src="figures/aligner.png" alt="pipe" height="450" width="1100">
</div>

* Garber et al., 2011. Nat Methods. 8(6):469-77. Table 1


Read alignment - gapped aligner (2/2)
====================================

<div align="center">
<img src="figures/aligner2.png" alt="pipe" height="600" width="1100">
</div>

* Garber et al., 2011. Nat Methods. 8(6):469-77. Figure 1


Check before you carry on ...
========================================================

* What genome is involved?
   - Is it a well annotated organism or a poorly annotated one?
   - A different analysis strategy depending on the genome?

* Single end or Paired end sequencing experiment?

* Strand specific or non stand specific
  - Useful to distinguish sense and antisense transcripts  


unstranded RNAseq
========================================================

<div align="center">
<img src="figures/Stranded_non.png" height="700" width="950">
</div>

Masp2 (sense) and Tardbp (antisense)

fr-firststrand RNAseq
========================================================

<div align="center">
<img src="figures/Stranded_one.png" height="700" width="950">
</div>

Masp2 (sense) and Tardbp (antisense)

========================================================

<div align="center">
<img src="figures/Stranded_two.png" height="700" width="950">
</div>

Masp2 (sense) and Tardbp (antisense)



Read Counting (1/2)
====================================
id: counting

After alignment, next step is to assign reads to features such as genes, transcripts
or exons.

* htseq-count

	Anders et al., 2014. Bioinformatics. 31(2):166-9
	
* featurecount function in Rsubread R package

	Liao et al., 2014. Bioinformatics. 30(7):923-30



Read Counting (2/2)
====================================

<div align="center">
<img src="figures/readcount.png" alt="pipe" height="700" width="1100">
</div>
	
http://htseq.readthedocs.io/en/master/count.html



Check before you carry on again ...
========================================================

* How many factors involved in experiment (genotype, treatment, time course etc.)

   - Single factor vs multi factor analysis? 
   
   - Interactions?
   
* In this session, we will only discuss the single factor analysis
   
   

RNA-seq analysis pipeline
====================================

<div align="center">
<img src="figures/RNAseqpipeline.png" alt="pipe" height="800" width="1027">
</div>



Set working directory
========================================================

Before running any of the code, we need to set the working directory to the folder we unarchived. 
You may navigate to the unarchived LMS_RNAseq_short/course folder in the Rstudio menu

**Session -> Set Working Directory -> Choose Directory**

or in the console, 

Use the getwd() to see where your current directory is


```r
getwd()
```

Use setwd() to set up your directory in the console


```r
setwd("/PathToMyDownload/LMS_RNAseq_short/course")
```


Material
====================================

Dataset
* [GSE68360] (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68360)

We have generated the following files for this session.

- **target** file that contain the sample information.

targets.txt

- **count** matrix  containing count table

AllCounts.csv






Read sample information
========================================================
id: de

Read in **targets.txt** file by using read.table function


```r
targets <- read.table("targets.txt",sep="\t",header=TRUE)  

targets
```

```
  Sample Group Batch        InputFile       InputFile2 OutputFile
1   Viv1   Viv     a SRR2001243.fastq SRR2001244.fastq   Viv1.bam
2   Viv2   Viv     b SRR2001245.fastq SRR2001246.fastq   Viv2.bam
3   Viv3   Viv     c SRR2001247.fastq SRR2001248.fastq   Viv3.bam
4   Hfd1   Hfd     a SRR2001249.fastq SRR2001250.fastq   Hfd1.bam
5   Hfd2   Hfd     b SRR2001251.fastq SRR2001252.fastq   Hfd2.bam
6   Hfd3   Hfd     c SRR2001253.fastq SRR2001254.fastq   Hfd3.bam
```

Read count data
========================================================

Load count data by using read.csv function


```r
AllCounts<-read.csv(file="AllCounts.csv",row.names = 1)
```

 Show first few rows of "AllCounts" object and the class of AllCounts object


```r
head(AllCounts)
```

```
                   Viv1 Viv2 Viv3 Hfd1 Hfd2 Hfd3
ENSMUSG00000090025    0    0    0    0    0    0
ENSMUSG00000064842    0    0    0    0    0    0
ENSMUSG00000051951    0    1    1    3    0    0
ENSMUSG00000089699    0    0    0    0    0    0
ENSMUSG00000088390    0    0    0    0    0    0
ENSMUSG00000089420    0    0    0    0    0    0
```

```r
class(AllCounts)
```

```
[1] "data.frame"
```


Prepare deseqdataset object 
========================================================

 Collect sample information


```r
cData<-data.frame(name=targets$Sample,                                                                          Group=targets$Group,Batch=targets$Batch)

rownames(cData)<-cData[,1]
```



Prepare deseqdataset object (Continued)
========================================================

 Construct deseqdataset object

The class used by the DESeq2 package to store the read counts is **DESeqDataSet**.



```r
library(DESeq2)

dds<-DESeqDataSetFromMatrix(
    countData= AllCounts,colData=cData,
    design=~Group)
```

Differential expression analysis
========================================================

The standard differential expression analysis steps are wrapped into a single function, **DESeq**,
which performs normalization, fitting to the model and statistical testing.



```r
dds<-DESeq(dds)

?DESeq
```

The function DESeq runs the following functions in order,

**estimateSizeFactors()**

**estimateDispersions()**

**nbinomWaldTest()**


DESeq function - estimateSizeFactors()
========================================================


1 - **Estimation of size factors**

The sizeFactors vector assigns to each column of the count matrix a value, the size factor, such that  count  values  in  the  columns  can  be  brought  to  a  common  scale  by  dividing  by  the  corresponding size factor.


```r
#  estimateSizeFactors()

sizeFactors(dds)
```

```
     Viv1      Viv2      Viv3      Hfd1      Hfd2      Hfd3 
1.2430187 0.7755226 1.0501449 0.9457439 1.0124687 1.0515602 
```

Calculating  sizeFactors
========================================================

<div align="center">
<img src="sizefactor.png" alt="gene" height="768" width="924">
</div>

DESeq2 calculation of the sizeFactors 
========================================================

Please refer to Equation (5):

**Anders S and Huber W. (2010). Differential expression analysis for sequence count data. *Genome Biol.*:11(10):R106.**

$$
\hat{s_{j}}=\begin{matrix}
  median \\
  i
 \end{matrix}\frac{k_{ij}}{(\prod_{v=1}^{m}k_{iv})^{1/m}}
 
 \\
 \hat{s_{j}}:\text{estimated size factor for the jth sample}
 \\
 
$$
The denominator of this expression can be interpreted as a pseudo-reference sample obtained by taking the geometric mean across samples. Thus, each estimated size factor is computed as the median of the ratios of the j-th sample's counts to those of the pseudo-reference.


DESeq function - estimateDispersions()
========================================================

2 - **Estimation of dispersion**

This function obtains gene-wide dispersion estimates. Then, a curve is fit to the estimates to capture the overall trend of dispersion-mean dependence. 


```r
#  estimateDispersions()

head(dispersions(dds))
```

```
[1]       NA       NA 2.495857       NA       NA       NA
```

```r
plotDispEsts(dds)
```

![plot of chunk unnamed-chunk-11](Practical_ShortRNAseq-figure/unnamed-chunk-11-1.png)


DESeq function - nbinomWaldTest()
========================================================
3- **Hypothesis test for differential expression**


For significance testing, DESeq2 by default uses a Wald test, where
the function tests whether each model coefficient differs significantly from zero, using previously calculated sizeFactors and dispersion estimates.

The Wald test P values  are adjusted for multiple testing using the procedure of Benjamini and Hochberg.



```r
nbinomWaldTest()
```



Getting results
========================================================

Results tables are generated using the function **results()**, which extracts a results table with log2 fold changes, p  values and adjusted p values.

 Use the function "results()"


```r
res<-results(dds) 
```

 Order results by adjusted p value


```r
resOrdered<-res[order(res$padj),]

?results      
```


Getting results
========================================================


```r
head(resOrdered)
```

```
log2 fold change (MLE): Group Viv vs Hfd 
Wald test p-value: Group Viv vs Hfd 
DataFrame with 6 rows and 6 columns
                           baseMean    log2FoldChange             lfcSE
                          <numeric>         <numeric>         <numeric>
ENSMUSG00000032080  14291.817820761 -5.04940305923945 0.251952097635319
ENSMUSG00000024526  465.87904057095 -6.89875038071526  0.34659786372613
ENSMUSG00000069170 343.928323594871 -4.05808344768568 0.253507932847091
ENSMUSG00000042041 651.002326876615 -3.18879830769903  0.20366517802941
ENSMUSG00000032231 883.324337619713 -2.58604293458338 0.176729903620997
ENSMUSG00000026043 833.867936897615 -2.46575350121176 0.171552106859386
                                stat               pvalue
                           <numeric>            <numeric>
ENSMUSG00000032080 -20.0411233191956 2.41261013993144e-89
ENSMUSG00000024526 -19.9041918682062 3.74278923935366e-88
ENSMUSG00000069170 -16.0077177945095 1.12874721852708e-57
ENSMUSG00000042041 -15.6570619413327 2.97318548586217e-55
ENSMUSG00000032231 -14.6327411581078 1.73647573696414e-48
ENSMUSG00000026043 -14.3732044237313 7.62213368291837e-47
                                   padj
                              <numeric>
ENSMUSG00000032080 4.02471623543363e-85
ENSMUSG00000024526 3.12186050454489e-84
ENSMUSG00000069170 6.27658703315625e-54
ENSMUSG00000042041 1.23996700687882e-51
ENSMUSG00000032231 5.79357764880717e-45
ENSMUSG00000026043 2.11920723497407e-43
```

Add Gene symbol
========================================================


```r
# Add Gene Symbols

library(biomaRt)

mart=useMart('ENSEMBL_MART_ENSEMBL',
dataset='mmusculus_gene_ensembl',
host="may2012.archive.ensembl.org")

bm<-getBM(attributes=c('ensembl_gene_id','mgi_symbol'),
filters ='ensembl_gene_id',
values=rownames(resOrdered), mart=mart)

# see the first few rows of "bm" object
head(bm)      
```

```
     ensembl_gene_id mgi_symbol
1 ENSMUSG00000000078       Klf6
2 ENSMUSG00000000325      Arvcf
3 ENSMUSG00000000340        Dbt
4 ENSMUSG00000000915      Hip1r
5 ENSMUSG00000001095    Slc13a2
6 ENSMUSG00000001119     Col6a1
```

Add Gene symbol (Continued)
========================================================


```r
# merge the Gene_symbol to our DE dataset
resAnnotated <- merge(as.data.frame(resOrdered),bm,by.x=0,by.y=1)
head(resAnnotated)
```

```
           Row.names     baseMean log2FoldChange     lfcSE        stat
1 ENSMUSG00000000001 2.438865e+03    -0.01014313 0.1128504 -0.08988119
2 ENSMUSG00000000003 0.000000e+00             NA        NA          NA
3 ENSMUSG00000000028 3.803191e+01    -0.27265098 0.3757113 -0.72569288
4 ENSMUSG00000000031 1.419411e+01    -2.87882064 0.6654750 -4.32596354
5 ENSMUSG00000000037 4.993369e-01    -2.44601888 2.6091665 -0.93747136
6 ENSMUSG00000000049 3.907640e+04     0.11127605 0.1400379  0.79461373
        pvalue         padj mgi_symbol
1 9.283816e-01 0.9763669677      Gnai3
2           NA           NA       Pbsn
3 4.680271e-01 0.7711998763      Cdc45
4 1.518665e-05 0.0003671649        H19
5 3.485161e-01           NA      Scml2
6 4.268382e-01 0.7439677545       Apoh
```

```r
# change the column name
colnames(resAnnotated)[1]<-"ensembl_gene_id"
```

Add Gene symbol (Continued)
========================================================


```r
# Order results by adjusted p value
resAnnotated<-resAnnotated[order(resAnnotated$pvalue,decreasing=F),]

# show the result with gene symbol annotation
head(resAnnotated)
```

```
         ensembl_gene_id   baseMean log2FoldChange     lfcSE      stat
8820  ENSMUSG00000032080 14291.8178      -5.049403 0.2519521 -20.04112
4582  ENSMUSG00000024526   465.8790      -6.898750 0.3465979 -19.90419
20354 ENSMUSG00000069170   343.9283      -4.058083 0.2535079 -16.00772
12526 ENSMUSG00000042041   651.0023      -3.188798 0.2036652 -15.65706
8894  ENSMUSG00000032231   883.3243      -2.586043 0.1767299 -14.63274
5372  ENSMUSG00000026043   833.8679      -2.465754 0.1715521 -14.37320
            pvalue         padj    mgi_symbol
8820  2.412610e-89 4.024716e-85         Apoa4
4582  3.742789e-88 3.121861e-84         Cidea
20354 1.128747e-57 6.276587e-54         Gpr98
12526 2.973185e-55 1.239967e-51 2010003K11Rik
8894  1.736476e-48 5.793578e-45         Anxa2
5372  7.622134e-47 2.119207e-43        Col3a1
```

Saving DEseq2 results
========================================================

Save the result in ".txt" or ".csv" format


```r
write.table(resAnnotated,file="DESeq_result.txt",sep="\t")
write.csv(resAnnotated,file="DESeq_result.csv",  row.names=F)
```

<div align="center">
<img src="figures/DESeq2_output.png" alt="gene" height="500" width="827">
</div>

Exploring results
========================================================

We can summarize some basic tallies using the summary function.


```r
summary(res)
```

```

out of 22605 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 992, 4.4%
LFC < 0 (down)     : 1459, 6.5%
outliers [1]       : 7, 0.031%
low counts [2]     : 5916, 26%
(mean count < 3)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

```r
# How many adjusted p-values were less than 0.05?
sum(res$padj < 0.05, na.rm=TRUE)
```

```
[1] 1970
```

How DESeq2 defines outliers - cooksCutoff
========================================================

theshold on Cook's distance, such that if one or more samples for a row have a distance higher, the p-value for the row is set to NA.


```r
summary(res)
```

```

out of 22605 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 992, 4.4%
LFC < 0 (down)     : 1459, 6.5%
outliers [1]       : 7, 0.031%
low counts [2]     : 5916, 26%
(mean count < 3)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

How DESeq2 filters low counts - independentFiltering
========================================================

A threshold on the filter statistic is found which optimizes the number of adjusted p values lower than a significance level alpha (we use the standard variable name for significance level, though it is unrelated to the dispersion parameter α). The theory behind independent filtering is discussed in greater detail below. The adjusted p values for the genes which do not pass the filter threshold are set to NA.


```r
summary(res)
```

```

out of 22605 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 992, 4.4%
LFC < 0 (down)     : 1459, 6.5%
outliers [1]       : 7, 0.031%
low counts [2]     : 5916, 26%
(mean count < 3)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

MA plot
========================================================

The  function **plotMA()** shows  the  log2  fold  changes  attributable  to  a  given  variable  over  the  mean of normalized counts.  Points will be colored red if the adjusted p value is less than 0.1.  Points which fall out of the window are plotted as open triangles pointing either up or down.


```r
plotMA(res, main="DESeq2", ylim=c(-4,4))
```

![plot of chunk unnamed-chunk-23](Practical_ShortRNAseq-figure/unnamed-chunk-23-1.png)

After calling plotMA, we can use the function identify to interactively detect the row number of individual
genes by clicking on the plot. 
Recover the gene identifiers by saving the resulting indices:


```r
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]
```


Plot counts
========================================================
 **Plot of normalized counts for a single gene on log scale**


```r
plotCounts(dds,gene=which.min(res$padj),intgroup="Group")
```

![plot of chunk unnamed-chunk-25](Practical_ShortRNAseq-figure/unnamed-chunk-25-1.png)

Coffee break and Exercise
=========================================================


* [RNAseq DE Exercises](exercise/Exercise_BasicDifferentialAnalysis.html)

Solutions
=========================================================

* [RNAseq DE Solutions](exercise/Answers_BasicDifferentialAnalysis.html)



Transformation of count data
========================================================
id: explore

In order to test for differential expression, we operate on raw counts. However for other downstream analyses ( e.g. visualization or clustering) it is useful to work with transformed versions of the count data.

Aim of these transformations,is  to  remove  the  dependence  of  the  variance on  the  mean, particularly  the  high  variance  of  the  logarithm  of  count  data  when  the  mean  is  low. 

    y = log2(n + n0)
            
    n represents the count values and n0 is a positive constant.
    
 There are two alternative approaches of choosing the parameter equivalent to n0 above.

**1. The regularized logarithm or rlog**

**2. variance stabilizing transformation or vst**


Both transforms the original count data to the log2 scale normalized to library size.



```r
rld<-rlog(dds)  

vsd <- varianceStabilizingTransformation(dds)
```

Transformation of count data - raw count
========================================================

<div align="left">
<img src="figures/gene.png" alt="gene" height="400" width="600">
</div>


Transformation of count data - rlog
========================================================

<div align="topright">
<img src="figures/rlog.png" alt="gene" height="400" width="600">
</div>

Transformation of count data - vst
========================================================

<div align="center">
<img src="figures/vst.png" alt="gene" height="400" width="600">
</div>


Data quality assessment 
========================================================

Data quality assessment and quality control  are essential steps
of any data analysis. These steps should typically be performed very early in the analysis of a new data set,
preceding or in parallel to the differential expression testing.

We will use following visualization tools to assess the data quality.

* Heatmap of count matrix

* Heatmap of sample to sample distances

* Principal component analysis plot


Heatmap of the count matrix
========================================================

To explore a counts matrix, it is often useful to look it as heatmap.




```r
library("pheatmap")
select <-order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]


pheatmap(assay(rld)[select,])
```

![plot of chunk unnamed-chunk-27](Practical_ShortRNAseq-figure/unnamed-chunk-27-1.png)



Heatmap of the count matrix
========================================================


```r
pheatmap(assay(vsd)[select,])
```

![plot of chunk unnamed-chunk-28](Practical_ShortRNAseq-figure/unnamed-chunk-28-1.png)



Heatmap of sample to sample distances
========================================================

Another use of the transformed data is sample clustering. Here, we apply the dist function to the transpose
of the transformed count matrix to get sample-to-sample distances. 

The assay function is used to extract the matrix of normalized values.


```r
rlogcount <- assay(rld)    
sampleDists <- as.matrix(dist(t(rlogcount)))
```

========================================================

 Save the plot as png file


```r
library(RColorBrewer)
showcols <- brewer.pal(8, "Set1")[1:length(unique(colData(dds)$Group))]

library(gplots)
png(file="sample_dis_map.png")
  heatmap.2(as.matrix(sampleDists), key=F, trace="none",
    col=colorpanel(100, "black", "white"),
    ColSideColors=showcols[colData(dds)$Group], 
    RowSideColors=showcols[colData(dds)$Group],
    margin=c(10, 10), main="Sample Distance Matrix")
dev.off()
```

```
quartz_off_screen 
                2 
```


=======================================================

<div align="center">
<img src="sample_dis_map.png" alt="gene" height="900" width="900">
</div>


Principal component plot of the samples
========================================================

**PCA plot** is useful to spot individual sample outliers. 


```r
plotPCA(rld, intgroup="Group")
```

![plot of chunk unnamed-chunk-31](Practical_ShortRNAseq-figure/unnamed-chunk-31-1.png)

```r
# save the plot

library(ggplot2)
ggsave(file="PCA_plot_version1.png")
```

Variations to standard workflow
========================================================


**Factor levels**

By default, R will choose a reference level for factors based on alphabetical order. 
The comparisons will be based on the alphabetical order of the levels. We can either explicitly tell results which comparison to make using the contrast argument (this will be shown later), or
we can explicitly set the factors levels.

**Setting the factor levels**

cData$Group <- factor(cData$Group, levels=c("Viv","Hfd"))

or using relevel, just specifying the reference level:


```r
cData$Group<-relevel(cData$Group,ref="Viv")
```


=========================================================

**Using Contrasts**

Contrasts enable the user to generate results for all  possible
comparisons: 

Consider an experimental design containing a factor with three levels, say A, B and C.
We can use contrasts to compare B vs A, of C vs A, and C vs B.


```r
    res_contrast<-results(dds,contrast=c("Group","Hfd","Viv")) 

    
    summary(res_contrast)
```

```

out of 22605 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 1459, 6.5%
LFC < 0 (down)     : 992, 4.4%
outliers [1]       : 7, 0.031%
low counts [2]     : 5916, 26%
(mean count < 3)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```



Exercises
=========================================================

* [RNAseq Visualization Exercises](exercise/Exercises_Visualization.html)

Solutions
=========================================================

* [RNAseq Visualization Solutions](exercise/Answers_Visualization.html)




Gene Ontology and Pathway Enrichment Analysis
========================================================
id: go

We will perform GO analysis using goseq package.

 In order to perform a GO analysis of your RNA-seq data,goseq only requires a simple named vector, which contains two pieces of information.

1 **Measured genes**
   
 all genes for which RNA-seq data was gathered for your experiment.  Each element of your vector should be named by a unique gene identifier.

2 **Differentially expressed genes**
   
 each element of your vector should be either a 1 or 0, where 1 indicates that the gene is differentially expressed and 0 that it is not.


========================================================


```r
library(KEGG.db)
library(goseq)

# remove the NAs

resdat<- res[complete.cases(res$padj),]
 
degenes<-as.integer(resdat$padj<0.05)
names(degenes)<-rownames(resdat)
    
# remove duplicate gene names
degenes<-degenes[match(unique(names(degenes)),                                                                                              names(degenes))]

table(degenes)
```

```
degenes
    0     1 
14712  1970 
```


========================================================

**Fitting the probability weighting function (PWF)**

We first need to obtain a weighting for each gene, 
depending on its length, given by the PWF



```r
pwf=nullp(degenes,genome="mm9",'ensGene', plot.fit=FALSE)
  
  head(pwf)
```

```
                   DEgenes bias.data       pwf
ENSMUSG00000025902       0    3190.5 0.1316076
ENSMUSG00000033845       1     830.0 0.1052353
ENSMUSG00000025903       0     938.0 0.1103089
ENSMUSG00000033813       0    2565.0 0.1316076
ENSMUSG00000062588       0     604.0 0.0932679
ENSMUSG00000033793       0    1907.0 0.1312570
```

========================================================


```r
   plotPWF(pwf)
```

![plot of chunk unnamed-chunk-36](Practical_ShortRNAseq-figure/unnamed-chunk-36-1.png)

========================================================

Change the Keggpath id to name in the goseq output


```r
xx <- as.list(KEGGPATHID2NAME)
temp <- cbind(names(xx),unlist(xx))
    
addKeggTogoseq <- function(JX,temp){
  for(l in 1:nrow(JX)){
      if(JX[l,1] %in% temp[,1]){
          JX[l,"term"] <- temp[temp[,1] %in% JX[l,1],2]
          JX[l,"ontology"] <- "KEGG"
      }
  }
  return(JX)
}
```

========================================================

 Calculate  the  over  and  under  expressed  GO
categories among DE genes


```r
go<-goseq(pwf,genome="mm9",'ensGene', test.cats=c("GO:BP","GO:MF","KEGG"))
```

========================================================

```r
head(go)
```

```
        category over_represented_pvalue under_represented_pvalue
10138 GO:0044281            8.068650e-26                        1
2663  GO:0006082            8.278639e-26                        1
9891  GO:0043436            4.079637e-25                        1
6089  GO:0019752            6.418361e-24                        1
7717  GO:0032787            2.213142e-22                        1
3021  GO:0006629            1.574901e-20                        1
      numDEInCat numInCat                                  term ontology
10138        327     1509      small molecule metabolic process       BP
2663         214      832        organic acid metabolic process       BP
9891         209      814             oxoacid metabolic process       BP
6089         200      781     carboxylic acid metabolic process       BP
7717         140      481 monocarboxylic acid metabolic process       BP
3021         229     1001               lipid metabolic process       BP
```

========================================================

```r
restemp<-addKeggTogoseq(go,temp)   

head(restemp)
```

```
        category over_represented_pvalue under_represented_pvalue
10138 GO:0044281            8.068650e-26                        1
2663  GO:0006082            8.278639e-26                        1
9891  GO:0043436            4.079637e-25                        1
6089  GO:0019752            6.418361e-24                        1
7717  GO:0032787            2.213142e-22                        1
3021  GO:0006629            1.574901e-20                        1
      numDEInCat numInCat                                  term ontology
10138        327     1509      small molecule metabolic process       BP
2663         214      832        organic acid metabolic process       BP
9891         209      814             oxoacid metabolic process       BP
6089         200      781     carboxylic acid metabolic process       BP
7717         140      481 monocarboxylic acid metabolic process       BP
3021         229     1001               lipid metabolic process       BP
```


========================================================

```r
write.table(restemp,file="GO_Kegg_Wallenius.txt", row.names=F,sep="\t")

write.csv(restemp,file="GO_Kegg_Wallenius.csv", row.names=F)
```


Exercises
=========================================================

* [RNAseq GO Exercises](exercise/Exercise_FunctionalAnalysis.html)


Solutions
=========================================================

* [RNAseq GO Solutions](exercise/Answers_FunctionalAnalysis.html)

Session Information
========================================================

```r
    sessionInfo()
```

```
R version 3.5.1 (2018-07-02)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS High Sierra 10.13.6

Matrix products: default
BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] gplots_3.0.1                biomaRt_2.36.1             
 [3] org.Mm.eg.db_3.6.0          KEGG.db_3.2.3              
 [5] AnnotationDbi_1.42.1        pheatmap_1.0.10            
 [7] ggplot2_3.0.0               RColorBrewer_1.1-2         
 [9] goseq_1.32.0                geneLenDataBase_1.16.0     
[11] BiasedUrn_1.07              DESeq2_1.20.0              
[13] SummarizedExperiment_1.10.1 DelayedArray_0.6.5         
[15] BiocParallel_1.14.2         matrixStats_0.54.0         
[17] Biobase_2.40.0              GenomicRanges_1.32.6       
[19] GenomeInfoDb_1.16.0         IRanges_2.14.11            
[21] S4Vectors_0.18.3            BiocGenerics_0.26.0        
[23] edgeR_3.22.3                limma_3.36.3               
[25] knitr_1.20                 

loaded via a namespace (and not attached):
 [1] nlme_3.1-137             bitops_1.0-6            
 [3] bit64_0.9-7              progress_1.2.0          
 [5] httr_1.3.1               tools_3.5.1             
 [7] backports_1.1.2          R6_2.2.2                
 [9] KernSmooth_2.23-15       rpart_4.1-13            
[11] Hmisc_4.1-1              DBI_1.0.0               
[13] lazyeval_0.2.1           mgcv_1.8-24             
[15] colorspace_1.3-2         nnet_7.3-12             
[17] withr_2.1.2              tidyselect_0.2.4        
[19] gridExtra_2.3            prettyunits_1.0.2       
[21] bit_1.1-14               compiler_3.5.1          
[23] htmlTable_1.12           rtracklayer_1.40.6      
[25] caTools_1.17.1.1         scales_1.0.0            
[27] checkmate_1.8.5          genefilter_1.62.0       
[29] stringr_1.3.1            digest_0.6.16           
[31] Rsamtools_1.32.3         foreign_0.8-71          
[33] XVector_0.20.0           base64enc_0.1-3         
[35] pkgconfig_2.0.2          htmltools_0.3.6         
[37] highr_0.7                htmlwidgets_1.2         
[39] rlang_0.2.2              rstudioapi_0.7          
[41] RSQLite_2.1.1            bindr_0.1.1             
[43] gtools_3.8.1             acepack_1.4.1           
[45] dplyr_0.7.6              RCurl_1.95-4.11         
[47] magrittr_1.5             GO.db_3.6.0             
[49] GenomeInfoDbData_1.1.0   Formula_1.2-3           
[51] Matrix_1.2-14            Rcpp_0.12.18            
[53] munsell_0.5.0            stringi_1.2.4           
[55] zlibbioc_1.26.0          plyr_1.8.4              
[57] grid_3.5.1               blob_1.1.1              
[59] gdata_2.18.0             crayon_1.3.4            
[61] lattice_0.20-35          Biostrings_2.48.0       
[63] splines_3.5.1            GenomicFeatures_1.32.2  
[65] annotate_1.58.0          hms_0.4.2               
[67] locfit_1.5-9.1           pillar_1.3.0            
[69] codetools_0.2-15         geneplotter_1.58.0      
[71] XML_3.98-1.16            glue_1.3.0              
[73] evaluate_0.11            latticeExtra_0.6-28     
[75] data.table_1.11.4        gtable_0.2.0            
[77] purrr_0.2.5              assertthat_0.2.0        
[79] xtable_1.8-3             survival_2.42-6         
[81] tibble_1.4.2             GenomicAlignments_1.16.0
[83] memoise_1.1.0            bindrcpp_0.2.2          
[85] cluster_2.0.7-1         
```
