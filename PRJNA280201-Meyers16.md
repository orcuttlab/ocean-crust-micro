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

Directory with reads

``` r
setwd("~/Desktop/Ocean-Crust-Synthesis/PRJNA280201/")
```

``` r
path <- "~/Desktop/Ocean-Crust-Synthesis/PRJNA280201/"
list.files(path)
```

    ##  [1] "filtered"              "meyer-dada2.Rmd"      
    ##  [3] "meyers_readtracks.csv" "meyers_taxa"          
    ##  [5] "meyers_trim_seqtab"    "meyers1_taxa"         
    ##  [7] "meyers1_trim_seqtab"   "pbs_sra.txt"          
    ##  [9] "PRJNA280201_files"     "PRJNA280201.Rmd"      
    ## [11] "ps_m12"                "SRR1973000_1.fastq"   
    ## [13] "SRR1973000_2.fastq"    "SRR1973001_1.fastq"   
    ## [15] "SRR1973001_2.fastq"    "SRR1973003_1.fastq"   
    ## [17] "SRR1973003_2.fastq"    "SRR1973004_1.fastq"   
    ## [19] "SRR1973004_2.fastq"    "SRR1973005_1.fastq"   
    ## [21] "SRR1973005_2.fastq"    "SRR1973006_1.fastq"   
    ## [23] "SRR1973006_2.fastq"    "SRR1973010_1.fastq"   
    ## [25] "SRR1973010_2.fastq"    "SRR1973012_1.fastq"   
    ## [27] "SRR1973012_2.fastq"

``` r
# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern="_1.fastq"))
fnRs <- sort(list.files(path, pattern="_2.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
```

``` r
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
```

trim off primed region, trimleft=19.

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(75,75), trimLeft = 19,
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
```

    ##                    reads.in reads.out
    ## SRR1973000_1.fastq   241340    227997
    ## SRR1973001_1.fastq    84082     78110
    ## SRR1973003_1.fastq    48306     44929
    ## SRR1973004_1.fastq    96866     90907
    ## SRR1973005_1.fastq    95053     89018
    ## SRR1973006_1.fastq   312491    294664

``` r
errF1 <- learnErrors(filtFs, multithread=TRUE)
```

    ## 60953200 total bases in 1088450 reads from 8 samples will be used for learning the error rates.
    ## Initializing error rates to maximum possible estimate.
    ## selfConsist step 1 ........
    ##    selfConsist step 2
    ##    selfConsist step 3
    ##    selfConsist step 4
    ##    selfConsist step 5
    ##    selfConsist step 6
    ##    selfConsist step 7
    ##    selfConsist step 8
    ##    selfConsist step 9
    ## Convergence after  9  rounds.

``` r
errR1 <- learnErrors(filtRs, multithread=TRUE)
```

    ## 60953200 total bases in 1088450 reads from 8 samples will be used for learning the error rates.
    ## Initializing error rates to maximum possible estimate.
    ## selfConsist step 1 ........
    ##    selfConsist step 2
    ##    selfConsist step 3
    ##    selfConsist step 4
    ##    selfConsist step 5
    ## Convergence after  5  rounds.

``` r
plotErrors(errF1, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](PRJNA280201_files/figure-markdown_github/unnamed-chunk-8-1.png)

``` r
plotErrors(errR1, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](PRJNA280201_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
derepFs <- derepFastq(filtFs, verbose=TRUE)
```

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA280201//filtered/SRR1973000_F_filt.fastq.gz

    ## Encountered 15218 unique sequences from 227997 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA280201//filtered/SRR1973001_F_filt.fastq.gz

    ## Encountered 6755 unique sequences from 78110 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA280201//filtered/SRR1973003_F_filt.fastq.gz

    ## Encountered 4199 unique sequences from 44929 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA280201//filtered/SRR1973004_F_filt.fastq.gz

    ## Encountered 8258 unique sequences from 90907 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA280201//filtered/SRR1973005_F_filt.fastq.gz

    ## Encountered 7705 unique sequences from 89018 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA280201//filtered/SRR1973006_F_filt.fastq.gz

    ## Encountered 19473 unique sequences from 294664 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA280201//filtered/SRR1973010_F_filt.fastq.gz

    ## Encountered 21820 unique sequences from 225070 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA280201//filtered/SRR1973012_F_filt.fastq.gz

    ## Encountered 6305 unique sequences from 37755 total sequences read.

``` r
derepRs <- derepFastq(filtRs, verbose=TRUE)
```

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA280201//filtered/SRR1973000_R_filt.fastq.gz

    ## Encountered 12469 unique sequences from 227997 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA280201//filtered/SRR1973001_R_filt.fastq.gz

    ## Encountered 8108 unique sequences from 78110 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA280201//filtered/SRR1973003_R_filt.fastq.gz

    ## Encountered 3512 unique sequences from 44929 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA280201//filtered/SRR1973004_R_filt.fastq.gz

    ## Encountered 6842 unique sequences from 90907 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA280201//filtered/SRR1973005_R_filt.fastq.gz

    ## Encountered 6268 unique sequences from 89018 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA280201//filtered/SRR1973006_R_filt.fastq.gz

    ## Encountered 15747 unique sequences from 294664 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA280201//filtered/SRR1973010_R_filt.fastq.gz

    ## Encountered 17953 unique sequences from 225070 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA280201//filtered/SRR1973012_R_filt.fastq.gz

    ## Encountered 5642 unique sequences from 37755 total sequences read.

``` r
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```

``` r
dadaFs <- dada(derepFs, err=errF1, multithread=TRUE)
```

    ## Sample 1 - 227997 reads in 15218 unique sequences.
    ## Sample 2 - 78110 reads in 6755 unique sequences.
    ## Sample 3 - 44929 reads in 4199 unique sequences.
    ## Sample 4 - 90907 reads in 8258 unique sequences.
    ## Sample 5 - 89018 reads in 7705 unique sequences.
    ## Sample 6 - 294664 reads in 19473 unique sequences.
    ## Sample 7 - 225070 reads in 21820 unique sequences.
    ## Sample 8 - 37755 reads in 6305 unique sequences.

``` r
dadaRs <- dada(derepRs, err=errR1, multithread=TRUE)
```

    ## Sample 1 - 227997 reads in 12469 unique sequences.
    ## Sample 2 - 78110 reads in 8108 unique sequences.
    ## Sample 3 - 44929 reads in 3512 unique sequences.
    ## Sample 4 - 90907 reads in 6842 unique sequences.
    ## Sample 5 - 89018 reads in 6268 unique sequences.
    ## Sample 6 - 294664 reads in 15747 unique sequences.
    ## Sample 7 - 225070 reads in 17953 unique sequences.
    ## Sample 8 - 37755 reads in 5642 unique sequences.

``` r
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
```

    ## 226309 paired-reads (in 474 unique pairings) successfully merged out of 226999 (in 594 pairings) input.

    ## 77136 paired-reads (in 323 unique pairings) successfully merged out of 77534 (in 381 pairings) input.

    ## 44404 paired-reads (in 159 unique pairings) successfully merged out of 44490 (in 175 pairings) input.

    ## 89854 paired-reads (in 296 unique pairings) successfully merged out of 90120 (in 344 pairings) input.

    ## 87788 paired-reads (in 278 unique pairings) successfully merged out of 88178 (in 319 pairings) input.

    ## 292606 paired-reads (in 582 unique pairings) successfully merged out of 293680 (in 723 pairings) input.

    ## 221981 paired-reads (in 724 unique pairings) successfully merged out of 223144 (in 963 pairings) input.

    ## 35889 paired-reads (in 416 unique pairings) successfully merged out of 37066 (in 518 pairings) input.

``` r
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```

    ##                                                                  sequence
    ## 1   AACCTTACCATCCCTTGACATCCAGAGAAGAGACTAGAGATAGACTTGTGCCTTCGGGAACTCTGTGAC
    ## 2 AACCTTACCTAGCCTTGACATTGATTGAACAGAGTAGAGATACACTGGTGCCCTTCGGGGAACATGAAAAC
    ## 3   AACCTTACCAGGCCTTGACATCCAATGAACTTTCTAGAGATAGATTGGTGCCTTCGGGAACATTGAGAC
    ## 4   AACCTTACCTGGGTTTGACATCCTTGGACCGCTACAGAGATGTAGTTTTCTCTTCGGAGACCAAGTGAC
    ## 5   AACCTTAGCATCCCTTGACATCCAGAGAAGAGACTAGAGATAGACTTGTGCCTTCGGGAACTCTGTGAC
    ## 6   AACCTCACCATCCCTTGACATCCAGAGAAGAGACTAGAGATAGACTTGTGCCTTCGGGAACTCTGTGAC
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1     50794       1       1     43         0      0      2   TRUE
    ## 2     15958       2       2     41         0      0      2   TRUE
    ## 3      9685       3       3     43         0      0      2   TRUE
    ## 4      9158       4       4     43         0      0      2   TRUE
    ## 5      5772       5       1     43         0      0      2   TRUE
    ## 6      4898       6       1     43         0      0      2   TRUE

``` r
seqtab <- makeSequenceTable(mergers)
```

    ## The sequences being tabled vary in length.

``` r
dim(seqtab)
```

    ## [1]    8 2497

``` r
table(nchar(getSequences(seqtab)))
```

    ## 
    ##  57  60  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76  77 
    ##   1   2  11   3  11   8 209 116 160 804 111 336  82 388  39  38   7  21 
    ##  78  80  81  82  86  88  91  93  94  95  96  98  99 
    ##  42  79   2   2   2   1   1   1   4   2   5   3   6

``` r
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(65,80)]
dim(seqtab2)
```

    ## [1]    8 2440

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 959 bimeras out of 2440 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]    8 1481

``` r
sum(seqtab.nochim)/sum(seqtab2)
```

    ## [1] 0.8390693

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
```

``` r
saveRDS(taxa, file = "meyers_taxa")
```

``` r
saveRDS(seqtab.nochim, file = "meyers_trim_seqtab")
```

``` r
tab_m12 <- readRDS(file="~/Desktop/Ocean-Crust-Synthesis/PRJNA280201/meyers_trim_seqtab")
taxam12 <-readRDS(file="~/Desktop/Ocean-Crust-Synthesis/PRJNA280201/meyers_taxa")
```

``` r
ps_m12 <- phyloseq(otu_table(tab_m12, taxa_are_rows=FALSE), 
               tax_table(taxam12))
ps_m12
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 1481 taxa and 8 samples ]
    ## tax_table()   Taxonomy Table:    [ 1481 taxa by 6 taxonomic ranks ]

``` r
colnames(tax_table(ps_m12)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
```

``` r
saveRDS(ps_m12, file = "~/Desktop/Ocean-Crust-Synthesis/PRJNA280201/ps_m12")
```
