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
library(vegan)
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.5-2

``` r
library(ggplot2)
library(plyr)
```

    ## 
    ## Attaching package: 'plyr'

    ## The following object is masked from 'package:ShortRead':
    ## 
    ##     id

    ## The following object is masked from 'package:matrixStats':
    ## 
    ##     count

    ## The following object is masked from 'package:XVector':
    ## 
    ##     compact

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     desc

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     rename

``` r
setwd("~/Desktop/Ocean-Crust-Synthesis/PRJNA308331/")
```

``` r
path <- "~/Desktop/Ocean-Crust-Synthesis/PRJNA308331/"
list.files(path)
```

    ##  [1] "filtered"              "pbs_sra.txt"          
    ##  [3] "ps_r17"                "reese-454_taxa"       
    ##  [5] "reese-454_trim_seqtab" "reese-454.Rmd"        
    ##  [7] "SRR3110000_.fastq"     "SRR3110001_.fastq"    
    ##  [9] "SRR3110002_.fastq"     "SRR3110003_.fastq"    
    ## [11] "SRR3110004_.fastq"     "SRR3110005_.fastq"

``` r
# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern=".fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
```

``` r
print(sample.names)
```

    ## [1] "SRR3110000" "SRR3110001" "SRR3110002" "SRR3110003" "SRR3110004"
    ## [6] "SRR3110005"

``` r
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_filt.fastq"))
```

``` r
plotQualityProfile(fnFs[1:2])
```

    ## Scale for 'y' is already present. Adding another scale for 'y', which
    ## will replace the existing scale.

![](reese-454_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
out <- filterAndTrim(fnFs, filtFs, trimLeft = 20, maxN = 0, maxEE = 2, truncQ = 2, truncLen = 300)
head(out)
```

    ##                   reads.in reads.out
    ## SRR3110000_.fastq     9068      7594
    ## SRR3110001_.fastq     9342      8226
    ## SRR3110002_.fastq    11757     10510
    ## SRR3110003_.fastq    10090      8065
    ## SRR3110004_.fastq     8486      7167
    ## SRR3110005_.fastq     5235      3562

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 12634720 total bases in 45124 reads from 6 samples will be used for learning the error rates.
    ## Initializing error rates to maximum possible estimate.
    ## selfConsist step 1 ......
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

![](reese-454_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
derepFs <- derepFastq(filtFs, verbose=TRUE)
```

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA308331//filtered/SRR3110000_filt.fastq

    ## Encountered 2685 unique sequences from 7594 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA308331//filtered/SRR3110001_filt.fastq

    ## Encountered 2907 unique sequences from 8226 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA308331//filtered/SRR3110002_filt.fastq

    ## Encountered 3693 unique sequences from 10510 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA308331//filtered/SRR3110003_filt.fastq

    ## Encountered 2873 unique sequences from 8065 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA308331//filtered/SRR3110004_filt.fastq

    ## Encountered 2293 unique sequences from 7167 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA308331//filtered/SRR3110005_filt.fastq

    ## Encountered 1444 unique sequences from 3562 total sequences read.

``` r
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
```

``` r
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
```

    ## Sample 1 - 7594 reads in 2685 unique sequences.
    ## Sample 2 - 8226 reads in 2907 unique sequences.
    ## Sample 3 - 10510 reads in 3693 unique sequences.
    ## Sample 4 - 8065 reads in 2873 unique sequences.
    ## Sample 5 - 7167 reads in 2293 unique sequences.
    ## Sample 6 - 3562 reads in 1444 unique sequences.

``` r
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
```

    ## [1]   6 594

``` r
table(nchar(getSequences(seqtab)))
```

    ## 
    ## 280 
    ## 594

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 125 bimeras out of 594 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]   6 469

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.8276244

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
```

``` r
saveRDS(taxa, file = "reese-454_taxa")
```

``` r
saveRDS(seqtab.nochim, file = "reese-454_trim_seqtab")
```

``` r
tab_r17 <-  readRDS(file="~/Desktop/Ocean-Crust-Synthesis/PRJNA308331/reese-454_trim_seqtab")
taxar17 <-readRDS(file="~/Desktop/Ocean-Crust-Synthesis/PRJNA308331/reese-454_taxa")
```

``` r
ps_r17 <- phyloseq(otu_table(tab_r17, taxa_are_rows=FALSE), 
               tax_table(taxar17))
ps_r17
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 469 taxa and 6 samples ]
    ## tax_table()   Taxonomy Table:    [ 469 taxa by 6 taxonomic ranks ]

``` r
colnames(tax_table(ps_r17)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
```

``` r
saveRDS(ps_r17, file = "~/Desktop/Ocean-Crust-Synthesis/PRJNA308331/ps_r17")
```
