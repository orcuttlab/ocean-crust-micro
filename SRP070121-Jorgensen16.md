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

``` r
setwd("~/Desktop/Ocean-Crust-Synthesis/SRP070121/")
```

``` r
path <- "~/Desktop/Ocean-Crust-Synthesis/SRP070121"
list.files(path)
```

    ##  [1] "74710.cfe1.ER"               "74710.cfe1.OU"              
    ##  [3] "filtered"                    "jorgenson_taxa"             
    ##  [5] "jorgenson_trim_seqtab"       "jorgeson_iontorrent.Rmd"    
    ##  [7] "pbs_sra.txt"                 "ps_j16"                     
    ##  [9] "SRP070121-Jorgensen16_files" "SRP070121-Jorgensen16.Rmd"  
    ## [11] "SRR3169157_.fastq.gz"        "SRR3169184_.fastq.gz"       
    ## [13] "SRR3169185_.fastq.gz"        "SRR3169186_.fastq.gz"       
    ## [15] "SRR3169187_.fastq.gz"        "SRR3169188_.fastq.gz"       
    ## [17] "SRR3169189_.fastq.gz"        "SRR3169191_.fastq.gz"       
    ## [19] "SRR3169192_.fastq.gz"        "SRR3169264_.fastq.gz"       
    ## [21] "SRR3169265_.fastq.gz"        "SRR3169266_.fastq.gz"       
    ## [23] "SRR3169275_.fastq.gz"        "SRR3169278_.fastq.gz"       
    ## [25] "SRR3169287_.fastq.gz"        "SRR3169288_.fastq.gz"       
    ## [27] "SRR3169295_.fastq.gz"        "SRR3169318_.fastq.gz"       
    ## [29] "SRR3169345_.fastq.gz"        "SRR3169347_.fastq.gz"       
    ## [31] "SRR3169348_.fastq.gz"        "SRR3169351_.fastq.gz"       
    ## [33] "SRR3169353_.fastq.gz"        "SRR3169354_.fastq.gz"       
    ## [35] "SRR3169355_.fastq.gz"        "SRR3169356_.fastq.gz"       
    ## [37] "SRR3169357_.fastq.gz"        "SRR3169415_.fastq.gz"       
    ## [39] "SRR3169416_.fastq.gz"        "SRR3169417_.fastq.gz"       
    ## [41] "SRR3169418_.fastq.gz"        "SRR3169419_.fastq.gz"       
    ## [43] "SRR3169420_.fastq.gz"        "Untitled.Rmd"

``` r
# Sort ensures forward/reverse reads are in same order
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

    ## Scale for 'y' is already present. Adding another scale for 'y', which
    ## will replace the existing scale.

![](SRP070121-Jorgensen16_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
out <- filterAndTrim(fnFs, filtFs, trimLeft = 20, maxN = 0, maxEE = 2, truncQ = 2, truncLen = 250, n = 1e+06)
head(out)
```

    ##                      reads.in reads.out
    ## SRR3169157_.fastq.gz    29247     22521
    ## SRR3169184_.fastq.gz    62259     47932
    ## SRR3169185_.fastq.gz    74866     55907
    ## SRR3169186_.fastq.gz    29247     22521
    ## SRR3169187_.fastq.gz    62259     47932
    ## SRR3169188_.fastq.gz    46865     35509

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 111816570 total bases in 486159 reads from 11 samples will be used for learning the error rates.
    ## Initializing error rates to maximum possible estimate.
    ## selfConsist step 1 ...........
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

![](SRP070121-Jorgensen16_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
derepFs <- derepFastq(filtFs, verbose=TRUE)
```

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169157_filt.fastq.gz

    ## Encountered 7640 unique sequences from 22521 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169184_filt.fastq.gz

    ## Encountered 12539 unique sequences from 47932 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169185_filt.fastq.gz

    ## Encountered 13952 unique sequences from 55907 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169186_filt.fastq.gz

    ## Encountered 7640 unique sequences from 22521 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169187_filt.fastq.gz

    ## Encountered 12539 unique sequences from 47932 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169188_filt.fastq.gz

    ## Encountered 8233 unique sequences from 35509 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169189_filt.fastq.gz

    ## Encountered 10866 unique sequences from 44180 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169191_filt.fastq.gz

    ## Encountered 13229 unique sequences from 55436 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169192_filt.fastq.gz

    ## Encountered 15626 unique sequences from 61640 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169264_filt.fastq.gz

    ## Encountered 9802 unique sequences from 38124 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169265_filt.fastq.gz

    ## Encountered 9866 unique sequences from 54457 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169266_filt.fastq.gz

    ## Encountered 7214 unique sequences from 54021 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169275_filt.fastq.gz

    ## Encountered 5622 unique sequences from 43173 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169278_filt.fastq.gz

    ## Encountered 12740 unique sequences from 87484 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169287_filt.fastq.gz

    ## Encountered 14611 unique sequences from 46407 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169288_filt.fastq.gz

    ## Encountered 16210 unique sequences from 58116 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169295_filt.fastq.gz

    ## Encountered 15125 unique sequences from 56011 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169318_filt.fastq.gz

    ## Encountered 13492 unique sequences from 46343 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169345_filt.fastq.gz

    ## Encountered 14412 unique sequences from 50347 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169347_filt.fastq.gz

    ## Encountered 14973 unique sequences from 53456 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169348_filt.fastq.gz

    ## Encountered 11747 unique sequences from 39773 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169351_filt.fastq.gz

    ## Encountered 10974 unique sequences from 46654 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169353_filt.fastq.gz

    ## Encountered 14085 unique sequences from 54527 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169354_filt.fastq.gz

    ## Encountered 14516 unique sequences from 55192 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169355_filt.fastq.gz

    ## Encountered 10382 unique sequences from 40103 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169356_filt.fastq.gz

    ## Encountered 14285 unique sequences from 60535 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169357_filt.fastq.gz

    ## Encountered 11378 unique sequences from 45057 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169415_filt.fastq.gz

    ## Encountered 24772 unique sequences from 68668 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169416_filt.fastq.gz

    ## Encountered 13090 unique sequences from 56653 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169417_filt.fastq.gz

    ## Encountered 12877 unique sequences from 48015 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169418_filt.fastq.gz

    ## Encountered 12708 unique sequences from 57234 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169419_filt.fastq.gz

    ## Encountered 14501 unique sequences from 73394 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/SRP070121/filtered/SRR3169420_filt.fastq.gz

    ## Encountered 1924 unique sequences from 5192 total sequences read.

``` r
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
```

``` r
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 22521 reads in 7640 unique sequences.
    ## Sample 2 - 47932 reads in 12539 unique sequences.
    ## Sample 3 - 55907 reads in 13952 unique sequences.
    ## Sample 4 - 22521 reads in 7640 unique sequences.
    ## Sample 5 - 47932 reads in 12539 unique sequences.
    ## Sample 6 - 35509 reads in 8233 unique sequences.
    ## Sample 7 - 44180 reads in 10866 unique sequences.
    ## Sample 8 - 55436 reads in 13229 unique sequences.
    ## Sample 9 - 61640 reads in 15626 unique sequences.
    ## Sample 10 - 38124 reads in 9802 unique sequences.
    ## Sample 11 - 54457 reads in 9866 unique sequences.
    ## Sample 12 - 54021 reads in 7214 unique sequences.
    ## Sample 13 - 43173 reads in 5622 unique sequences.
    ## Sample 14 - 87484 reads in 12740 unique sequences.
    ## Sample 15 - 46407 reads in 14611 unique sequences.
    ## Sample 16 - 58116 reads in 16210 unique sequences.
    ## Sample 17 - 56011 reads in 15125 unique sequences.
    ## Sample 18 - 46343 reads in 13492 unique sequences.
    ## Sample 19 - 50347 reads in 14412 unique sequences.
    ## Sample 20 - 53456 reads in 14973 unique sequences.
    ## Sample 21 - 39773 reads in 11747 unique sequences.
    ## Sample 22 - 46654 reads in 10974 unique sequences.
    ## Sample 23 - 54527 reads in 14085 unique sequences.
    ## Sample 24 - 55192 reads in 14516 unique sequences.
    ## Sample 25 - 40103 reads in 10382 unique sequences.
    ## Sample 26 - 60535 reads in 14285 unique sequences.
    ## Sample 27 - 45057 reads in 11378 unique sequences.
    ## Sample 28 - 68668 reads in 24772 unique sequences.
    ## Sample 29 - 56653 reads in 13090 unique sequences.
    ## Sample 30 - 48015 reads in 12877 unique sequences.
    ## Sample 31 - 57234 reads in 12708 unique sequences.
    ## Sample 32 - 73394 reads in 14501 unique sequences.
    ## Sample 33 - 5192 reads in 1924 unique sequences.

``` r
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
```

    ## [1]   33 5268

``` r
table(nchar(getSequences(seqtab)))
```

    ## 
    ##  230 
    ## 5268

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 154 bimeras out of 5268 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]   33 5114

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.9944626

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
```

``` r
saveRDS(taxa, file = "jorgenson_taxa")
```

``` r
saveRDS(seqtab.nochim, file = "jorgenson_trim_seqtab")
```

``` r
tab_j16 <- readRDS(file="~/Desktop/Ocean-Crust-Synthesis/SRP070121/jorgenson_trim_seqtab")
taxaj16 <-readRDS(file="~/Desktop/Ocean-Crust-Synthesis/SRP070121/jorgenson_taxa")
```

``` r
ps_j16 <- phyloseq(otu_table(tab_j16, taxa_are_rows=FALSE), 
               tax_table(taxaj16))
ps_j16
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 5114 taxa and 33 samples ]
    ## tax_table()   Taxonomy Table:    [ 5114 taxa by 6 taxonomic ranks ]

``` r
colnames(tax_table(ps_j16)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
```

``` r
saveRDS(ps_j16, file = "~/Desktop/Ocean-Crust-Synthesis/SRP070121/ps_j16")
```
