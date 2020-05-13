# Gene315_ChIP


Week 1 Using GenomicRanges to find overlapping and subsets of peak regions

## Using GenomicRanges to find common and unique peak regions

First you need to make sure you have the required packages installed. If not, you can download them from bioconductor using the code below:

```{r message=FALSE, results='hide'}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicRanges")
BiocManager::install("RIPSeeker")
```
Note: if you are using an older version of R (earlier than 3.6), install the packages using Bioconductor
```{r message=FALSE, results='hide'}
source("http://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
biocLite("RIPSeeker")
```

1. Load the required packages into R
```{r message =FALSE, results='hide'}
library("RIPSeeker")
library(GenomicRanges)
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
4. Generate a subset of peak regions: peak regions that a common to both data sets. In this example, this function will give me a list of the genome regions enriched for both the H3Kme1 and K27me3 modifications in ES cells.   If there is an error with exporting subset.bed -  check if RIPSeeker package is installed correctly (ie if you type library("RIPSeeker") from above it shouldn’t have displayed an error).  

```{r}
subset <- subsetByOverlaps(pk1.gr,pk2.gr)
#How many peak regions were common to both datasets?
length(subset)
#save the results to a new bed file
export.bed(subset, "subset.bed")
```

5. How many peak regions are unique to the first data set (in my case, H3K4me1, loaded into R as pk1 (above))

```{r}
setdiff <- setdiff(pk1.gr, pk2.gr)
length(setdiff)
export.bed(setdiff, "setdiff1.bed")
```

6. What about for the second dataset (in this case pk2=H3K27ac)?

```{r}
setdiff <- setdiff(pk2.gr, pk1.gr)
length(setdiff)
export.bed(setdiff, "setdiff2.bed")
```

Now you will have 3 new bed files in your working directory. One for peak regions common to both (ie. both histone marks are present at these genomic co-ordinates) and two files with the unique peak regions (ie. one of the histone marks is present).

Now go back to the Week1 webpage and continue to the next task. 

### Week 2 Motif Analysis: find transcription factor sites overrepresented in your peak regions of interest

## How common are the top DNA binding motifs?

Using RSAT we identified candidate transcription factor (TF) DNA binding motifs, enriched within the subset of ChIP peak regions. 

1. Download the required packages from bioconductor (if required)
```{r message =FALSE, results='hide'}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
BiocManager::install("JASPAR2018")
BiocManager::install("TFBSTools")
BiocManager::install("ggplot2")
```
Note: if you are using an older version of R (earlier than 3.6), install the packages using Bioconductor
```{r message =FALSE, results='hide'}
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Mmusculus.UCSC.mm10")
biocLite("JASPAR2018")
biocLite("GO.db")
biocLite("TFBSTools")
biocLite("ggplot2")
```


2. Load the required packages into Rstudio

```{r message =FALSE, results='hide'}
require(rtracklayer)
#load the JASPAR2018 package
require(JASPAR2018)
#load the TFBSTool package
require(TFBSTools)
#load the rtracklayer package
require(rtracklayer)
#load the mouse genome assembly mm10
require(BSgenome.Mmusculus.UCSC.mm10)
```

3. Import into Rstudo your peak bed file (shared peak regions)

```{r}
peaks <- import.bed("subset.bed")
```

4. Extract the DNA sequences for each peak region from the mouse genome

```{r}
genome <- BSgenome.Mmusculus.UCSC.mm10
peaks.seq <- getSeq(genome, peaks)
```

5. Now you will use the motif IDs (from RSAT) to find the motifs in each sequence. In my example, my two top TF binding sites were for Errb and Nkx2-8.
+ ID number:MA0141 	
+ ID for my other factor of interest was MA0673

*remember the ID number you use will be different to the example here* Only use the base ID (eg here I use MA0141 not MA0141.2")

6. First I will load the DNA motif details for Errb. 

```{r} 
pfm <- getMatrixByID(JASPAR2018, ID="MA0141")
pfm
```

7. Now find all the sequences with this motif present

```{r}
hitsErrb <- lapply(peaks.seq, function(s) matchPWM(as.matrix(pfm), s, min.score="75%"))
# how many peak sequences contained this motif?
countErrb <- sapply(hitsErrb, length)
sum(countErrb >= 1)
```

8. Repeat for the second motif. 

```{r} 
pfm <- getMatrixByID(JASPAR2018, ID="MA0673")
pfm
hitsNkx2 <- lapply(peaks.seq, function(s) matchPWM(as.matrix(pfm), s, min.score="75%"))
# how many peak sequences contained this second motif?
countNkx2 <- sapply(hitsNkx2, length)
sum(countNkx2 >= 1)
```

9. Now we are ready to plot the frequency of each motif in our peak sequences
```{r}
#create a data.frame with the number of motif hits and the corresponding peak source group as columns
plotDF <- data.frame(motif_hits = c(countErrb, countNkx2), group = rep(c("Errb peaks", "Nkx2 peaks"), c(length(countErrb), length(countNkx2))))
require(ggplot2)
qplot(x=motif_hits, fill=group, data=plotDF, geom="bar", facets=group~.)

#You can change the limits of the Y and X axes using the code below to improve the graph.
qplot(x=motif_hits, fill=group, data=plotDF, geom="bar", facets=group~.) +xlim(0,5)+ylim(0,10000)
```
Around 20,000 peak regions (areas of the genome with both the H3K4me1 and H3K27ac modifications) do not have a binding site for Erbb. Over 10,000 peak regions do have a single Errb site, and almost 5,000 have two Errb binding sites.Potentially, these may represent active enhancer elements bound by Errb in ES cells (the source of the chromatin used in the original ChIP experiment), as both of these histone modifications are associated with active enhancers in the genome. 

However, I need to consider that some/many of these sites are false positives. How could I validate these analyses?

