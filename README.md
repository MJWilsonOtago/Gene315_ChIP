# Gene315_ChIP


Week 1 Using GenomicRanges to find overlapping and subsets of peak regions

## Using GenomicRanges to find common and unique peak regions

If you are working on your own computer and do not have the required packages installed, you can download them from BiocManager using the code below. **Student Desktop already has these packaged installed - move onto the next step.

```{r message=FALSE, results='hide'}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicRanges")
BiocManager::install("rtracklayer")
```
Note: if you are using an older version of R (earlier than 3.6), install the packages using Bioconductor
```{r message=FALSE, results='hide'}
source("http://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
biocLite("rtracklayer")
```

1. Load the required packages into R
```{r message =FALSE, results='hide'}
library("rtracklayer")
library("GenomicRanges")
```

2. Read the peak region bed files into R for each data set. Remember- you need to change the .bed file name to one specific for your research question (eg the folders you downloaded earlier (task2) will contain the appropriate bed files eg limb H3K4me1 bed file is LImbH3K4me1E11.bed)

```{r}
pk1=read.table("H3K4me1_all_peaks.bed")
pk2=read.table("H3K27Ac_all_peaks.bed")
```

3. Convert these into an appropriate data frame

```{r}
pk1.gr=makeGRangesFromDataFrame(pk1, seqnames.field=c("V1"),start.field=c("V2"), end.field=c("V3"))
pk2.gr=makeGRangesFromDataFrame(pk2, seqnames.field=c("V1"),start.field=c("V2"), end.field=c("V3"))
```
4. Generate a subset of peak regions: peak regions that a common to both data sets. In this example, this function will give me a list of the genome regions enriched for both the H3Kme1 and K27me3 modifications in ES cells.   If there is an error with exporting subset.bed -  check if rtracklayer package is installed correctly (ie if you type library("rtracklayer") from above it shouldn’t have displayed an error).  

```{r}
subset <- subsetByOverlaps(pk1.gr,pk2.gr)
#How many peak regions were common to both datasets?
length(subset)
#save the results to a new bed file
export.bed(subset, "subset.bed")
```

5. How many peak regions are unique to the first data set (in my case, H3K4me1, loaded into R as pk1 (above))

```{r}
setdiff1 <- setdiff(pk1.gr, pk2.gr)
length(setdiff1)
export.bed(setdiff1, "setdiff1.bed")
```

6. What about for the second dataset (in this case pk2=H3K27ac)?

```{r}
setdiff2 <- setdiff(pk2.gr, pk1.gr)
length(setdiff2)
export.bed(setdiff2, "setdiff2.bed")
```

Now you will have 3 new bed files in your working directory. One for peak regions common to both (ie. both histone marks are present at these genomic co-ordinates) and two files with the unique peak regions (ie. one of the histone marks is present).

Next we will annotate where these regions are in relation to gene features (eg do they overlap with gene promoter regions?)

7. If you are working on your own computer and do not have the required packages installed, you can download them from BiocManager using the code below. **StudentDesktop already has these packaged installed - move onto the next step.
```{r}
BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("clusterProfiler")
```
8. Load this packages into R

```{r}
library(clusterProfiler)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
```
9. Assign the annotation database for the mouse genome


```{r}
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
```
10. Make a group of your files and give them a name appropriate for what is in each bed file (for example, the subset file contains regions that overlap, for my analysis the setdiff1 file had the regions unique to the H3K27ac dataset).

```{r}
files <- list(overlap = "subset.bed",H3K27ac = "setdiff1.bed", H3K4me3 = "setdiff2.bed")
print(files)
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-2000, 2000), verbose=FALSE)
```
11. Have a look at what is stored in peakAnno
```{r}
peakAnnoList
```
12. And to visulise this, we will make a bar graph to more easily compare the datasets.

```{r}
plotAnnoBar(peakAnnoList)
```
Save this figure for your poster (Export...)

13. We can also looks the peak locations over the whole genome using covplot function - this calculates the coverage of peaks over each chromosome. Save this figure for your poster (Export...)
```{r}
covplot(subset, title = "ChIP Peaks over Chromosomes",xlab = "Chromosome Size (bp)")
```
14. Create a profile of ChIP peaks around the transcriptional start site regions. Save this figure for your poster (Export...)
```{r}
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(subset, windows=promoter)

plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Peak Count Frequency")
```

Now go back to the Week1 webpage and continue to the next task. 


