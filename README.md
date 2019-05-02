# Gene315_ChIP


Week 1 Using GenomicRanges to find overlapping and subsets of peak regions

https://rawgit.com/MJWilsonOtago/Gene315_ChIP/master/GenomicRanges.html


Week 2 Motif Analysis: find transcription factor sites overrepresented in your peak regions of interest

## How common are the top DNA binding motifs?

Using RSAT we identified candidate transcription factor (TF) DNA binding motifs, enriched within the subset of ChIP peak regions. 

1. Download the required packages from bioconductor (if required)

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

5. Now you need to obtain the TF motif details in order to use TBStools to search for this motif in your peak sequences. In my example, my two top TF binding sites were for Errb and Nkx2-8.

To get the Transcription factor (TF) motif information from the JASPAR website. Go to http://jaspar.genereg.net/ 
Select JASPAR CORE vertebrata
+ I searched for Errb (name of the transcription factor)
+ Went to the page with the  information regarding the DNA binding motif eg how it was identified and has a sequence logo you can download for your poster. 
+ Noted down the matrix ID number:MA0141.2 	
+ The matrix ID for my other factor of interest was MA0673.1

*remember the ID number you use will be different to the example here*

6. First I will load the DNA motif details for Errb

```{r} 
pfm <- getMatrixByID(JASPAR2018, ID="MA0141.2")
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
pfm <- getMatrixByID(JASPAR2018, ID="MA0673.1")
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
```


Around 20,000 peak regions (areas of the genome with both the H3K4me1 and H3K27ac modifications) do not have a binding site for Erbb. Over 10,000 peak regions do have a single Errb site, and almost 5,000 have two Errb binding sites.Potentially, these may represent active enhancer elements bound by Errb in ES cells (the source of the chromatin used in the original ChIP experiment), as both of these histone modifications are associated with active enhancers in the genome. 

However, I need to consider that some/many of these sites are false positives. How could I validate these analyses?

https://rawgit.com/MJWilsonOtago/Gene315_ChIP/master/Motifs.html
