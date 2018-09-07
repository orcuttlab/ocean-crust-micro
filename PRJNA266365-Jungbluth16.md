``` r
library(ShortRead)
```

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind,
    ##     colMeans, colnames, colSums, dirname, do.call, duplicated,
    ##     eval, evalq, Filter, Find, get, grep, grepl, intersect,
    ##     is.unsorted, lapply, lengths, Map, mapply, match, mget, order,
    ##     paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind,
    ##     Reduce, rowMeans, rownames, rowSums, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which, which.max,
    ##     which.min

    ## Loading required package: BiocParallel

    ## Loading required package: Biostrings

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: Rsamtools

    ## Loading required package: GenomeInfoDb

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomicAlignments

    ## Loading required package: SummarizedExperiment

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     type

    ## The following objects are masked from 'package:base':
    ## 
    ##     aperm, apply

``` r
library(dada2)
```

    ## Loading required package: Rcpp

``` r
library(phyloseq)
```

    ## 
    ## Attaching package: 'phyloseq'

    ## The following object is masked from 'package:SummarizedExperiment':
    ## 
    ##     distance

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     sampleNames

    ## The following object is masked from 'package:GenomicRanges':
    ## 
    ##     distance

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     distance

Set to directory with reads

``` r
setwd("~/Desktop/Ocean-Crust-Synthesis/PRJNA266365//")
```

``` r
path <- "~/Desktop/Ocean-Crust-Synthesis/PRJNA266365//"
list.files(path)
```

    ##  [1] "74701.cfe1.ER"            "74701.cfe1.OU"           
    ##  [3] "filtered"                 "jf16_taxa"               
    ##  [5] "jf16_trim_seqtab"         "jungbluth16-forwards.Rmd"
    ##  [7] "pbs_sra1.txt"             "ps_j16"                  
    ##  [9] "SRR1646879_.fastq.gz"     "SRR1646880_.fastq.gz"    
    ## [11] "SRR1646881_.fastq.gz"     "SRR1646882_.fastq.gz"    
    ## [13] "SRR1646883_.fastq.gz"     "SRR1646884_.fastq.gz"    
    ## [15] "SRR1646885_.fastq.gz"     "SRR1646886_.fastq.gz"    
    ## [17] "SRR1646887_.fastq.gz"     "SRR1646888_.fastq.gz"    
    ## [19] "SRR1646889_.fastq.gz"     "SRR1646890_.fastq.gz"    
    ## [21] "SRR1646891_.fastq.gz"     "SRR1646892_.fastq.gz"    
    ## [23] "SRR1646893_.fastq.gz"     "SRR1646894_.fastq.gz"    
    ## [25] "SRR1646895_.fastq.gz"     "SRR1646896_.fastq.gz"    
    ## [27] "SRR1646897_.fastq.gz"     "SRR1646898_.fastq.gz"    
    ## [29] "SRR1646901_.fastq.gz"

``` r
fnFs <- sort(list.files(path, pattern=".fastq.gz"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
```

``` r
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_filt.fastq.gz"))
```

``` r
plotQualityProfile(fnFs[1:2])
```

![](jungbluth16-forwards_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
out <- filterAndTrim(fnFs, filtFs, trimLeft = 20, maxN = 0, maxEE = 2, truncQ = 2)
head(out)
```

    ##                      reads.in reads.out
    ## SRR1646879_.fastq.gz    19813     19319
    ## SRR1646880_.fastq.gz    23917     23518
    ## SRR1646881_.fastq.gz    25639     25129
    ## SRR1646882_.fastq.gz    39527     38820
    ## SRR1646883_.fastq.gz     6123      5994
    ## SRR1646884_.fastq.gz    38153     37484

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## 61932987 total bases in 472771 reads from 21 samples will be used for learning the error rates.
    ## Initializing error rates to maximum possible estimate.
    ## selfConsist step 1 .....................
    ##    selfConsist step 2
    ##    selfConsist step 3
    ##    selfConsist step 4
    ##    selfConsist step 5
    ##    selfConsist step 6
    ## Convergence after  6  rounds.

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](jungbluth16-forwards_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
derepFs <- derepFastq(filtFs, verbose=TRUE)
```

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646879_filt.fastq.gz

    ## Encountered 4145 unique sequences from 19319 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646880_filt.fastq.gz

    ## Encountered 4961 unique sequences from 23518 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646881_filt.fastq.gz

    ## Encountered 6540 unique sequences from 25129 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646882_filt.fastq.gz

    ## Encountered 9754 unique sequences from 38820 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646883_filt.fastq.gz

    ## Encountered 1490 unique sequences from 5994 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646884_filt.fastq.gz

    ## Encountered 9889 unique sequences from 37484 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646885_filt.fastq.gz

    ## Encountered 7346 unique sequences from 30578 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646886_filt.fastq.gz

    ## Encountered 9946 unique sequences from 38626 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646887_filt.fastq.gz

    ## Encountered 3854 unique sequences from 16911 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646888_filt.fastq.gz

    ## Encountered 1909 unique sequences from 8202 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646889_filt.fastq.gz

    ## Encountered 2795 unique sequences from 13712 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646890_filt.fastq.gz

    ## Encountered 4832 unique sequences from 20439 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646891_filt.fastq.gz

    ## Encountered 3868 unique sequences from 20813 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646892_filt.fastq.gz

    ## Encountered 3740 unique sequences from 26616 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646893_filt.fastq.gz

    ## Encountered 2381 unique sequences from 23235 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646894_filt.fastq.gz

    ## Encountered 3781 unique sequences from 22817 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646895_filt.fastq.gz

    ## Encountered 5098 unique sequences from 28782 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646896_filt.fastq.gz

    ## Encountered 4157 unique sequences from 17562 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646897_filt.fastq.gz

    ## Encountered 3937 unique sequences from 19027 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646898_filt.fastq.gz

    ## Encountered 2896 unique sequences from 11002 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA266365///filtered/SRR1646901_filt.fastq.gz

    ## Encountered 6223 unique sequences from 24185 total sequences read.

``` r
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
```

``` r
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 19319 reads in 4145 unique sequences.
    ## Sample 2 - 23518 reads in 4961 unique sequences.
    ## Sample 3 - 25129 reads in 6540 unique sequences.
    ## Sample 4 - 38820 reads in 9754 unique sequences.
    ## Sample 5 - 5994 reads in 1490 unique sequences.
    ## Sample 6 - 37484 reads in 9889 unique sequences.
    ## Sample 7 - 30578 reads in 7346 unique sequences.
    ## Sample 8 - 38626 reads in 9946 unique sequences.
    ## Sample 9 - 16911 reads in 3854 unique sequences.
    ## Sample 10 - 8202 reads in 1909 unique sequences.
    ## Sample 11 - 13712 reads in 2795 unique sequences.
    ## Sample 12 - 20439 reads in 4832 unique sequences.
    ## Sample 13 - 20813 reads in 3868 unique sequences.
    ## Sample 14 - 26616 reads in 3740 unique sequences.
    ## Sample 15 - 23235 reads in 2381 unique sequences.
    ## Sample 16 - 22817 reads in 3781 unique sequences.
    ## Sample 17 - 28782 reads in 5098 unique sequences.
    ## Sample 18 - 17562 reads in 4157 unique sequences.
    ## Sample 19 - 19027 reads in 3937 unique sequences.
    ## Sample 20 - 11002 reads in 2896 unique sequences.
    ## Sample 21 - 24185 reads in 6223 unique sequences.

``` r
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
```

    ## [1]   21 2522

``` r
table(nchar(getSequences(seqtab)))
```

    ## 
    ##  131 
    ## 2522

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 82 bimeras out of 2522 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]   21 2440

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.991409

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
```

``` r
saveRDS(taxa, file = "jf16_taxa")
```

``` r
saveRDS(seqtab.nochim, file = "jf16_trim_seqtab")
```

``` r
ps_j16 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))
ps_j16
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 2440 taxa and 21 samples ]
    ## tax_table()   Taxonomy Table:    [ 2440 taxa by 6 taxonomic ranks ]

``` r
colnames(tax_table(ps_j16)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
```

``` r
saveRDS(ps_j16, file = "~/Desktop/Ocean-Crust-Synthesis/PRJNA266365/ps_j16")
```
