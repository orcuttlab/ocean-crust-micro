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
setwd("~/Desktop/Ocean-Crust-Synthesis/PRJNA433243/")
```

``` r
path <- "~/Desktop/Ocean-Crust-Synthesis/PRJNA433243/"
list.files(path)
```

    ##  [1] "filtered"                  "names.csv"                
    ##  [3] "pbs_sra.txt"               "PRJNA433243-Zinke18_files"
    ##  [5] "PRJNA433243-Zinke18.Rmd"   "SRR6715661_1.fastq"       
    ##  [7] "SRR6715662_1.fastq"        "SRR6715663_1.fastq"       
    ##  [9] "SRR6715664_1.fastq"        "SRR6715665_1.fastq"       
    ## [11] "SRR6715666_1.fastq"        "SRR6715667_1.fastq"       
    ## [13] "SRR6715668_1.fastq"        "SRR6715669_1.fastq"       
    ## [15] "SRR6715670_1.fastq"        "SRR6715671_1.fastq"       
    ## [17] "SRR6715672_1.fastq"        "SRR6715673_1.fastq"       
    ## [19] "SRR6715674_1.fastq"        "SRR6715675_1.fastq"       
    ## [21] "SRR6715676_1.fastq"        "SRR6715677_1.fastq"       
    ## [23] "SRR6715678_1.fastq"        "SRR6715679_1.fastq"       
    ## [25] "SRR6715679.fastq"          "zinke-final"              
    ## [27] "zinke-rmd.Rmd"

``` r
# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern="_1.fastq"))
fnRs <- sort(list.files(path, pattern="_2.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
sample.names
```

    ##  [1] "SRR6715661" "SRR6715662" "SRR6715663" "SRR6715664" "SRR6715665"
    ##  [6] "SRR6715666" "SRR6715667" "SRR6715668" "SRR6715669" "SRR6715670"
    ## [11] "SRR6715671" "SRR6715672" "SRR6715673" "SRR6715674" "SRR6715675"
    ## [16] "SRR6715676" "SRR6715677" "SRR6715678" "SRR6715679"

``` r
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
```

``` r
plotQualityProfile(fnFs[1:2])
```

![](PRJNA433243-Zinke18_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
out <- filterAndTrim(fnFs, filtFs, trimLeft = 20, maxN = 0, maxEE = 2, truncQ = 2)
head(out)
```

    ##                    reads.in reads.out
    ## SRR6715661_1.fastq   246877    246877
    ## SRR6715662_1.fastq   293932    293932
    ## SRR6715663_1.fastq   359238    359238
    ## SRR6715664_1.fastq   363694    363694
    ## SRR6715665_1.fastq   329859    329859
    ## SRR6715666_1.fastq   161896    161896

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 118977980 total bases in 540809 reads from 2 samples will be used for learning the error rates.
    ## Initializing error rates to maximum possible estimate.
    ## selfConsist step 1 ..
    ##    selfConsist step 2
    ##    selfConsist step 3
    ##    selfConsist step 4
    ##    selfConsist step 5
    ##    selfConsist step 6
    ##    selfConsist step 7
    ## Convergence after  7  rounds.

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](PRJNA433243-Zinke18_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
derepFs <- derepFastq(filtFs, verbose=TRUE)
```

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA433243//filtered/SRR6715661_F_filt.fastq.gz

    ## Encountered 51143 unique sequences from 246877 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA433243//filtered/SRR6715662_F_filt.fastq.gz

    ## Encountered 130964 unique sequences from 293932 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA433243//filtered/SRR6715663_F_filt.fastq.gz

    ## Encountered 159867 unique sequences from 359238 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA433243//filtered/SRR6715664_F_filt.fastq.gz

    ## Encountered 171953 unique sequences from 363694 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA433243//filtered/SRR6715665_F_filt.fastq.gz

    ## Encountered 155229 unique sequences from 329859 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA433243//filtered/SRR6715666_F_filt.fastq.gz

    ## Encountered 70739 unique sequences from 161896 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA433243//filtered/SRR6715667_F_filt.fastq.gz

    ## Encountered 89346 unique sequences from 176938 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA433243//filtered/SRR6715668_F_filt.fastq.gz

    ## Encountered 116426 unique sequences from 237748 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA433243//filtered/SRR6715669_F_filt.fastq.gz

    ## Encountered 93554 unique sequences from 183070 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA433243//filtered/SRR6715670_F_filt.fastq.gz

    ## Encountered 111385 unique sequences from 240469 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA433243//filtered/SRR6715671_F_filt.fastq.gz

    ## Encountered 131764 unique sequences from 265663 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA433243//filtered/SRR6715672_F_filt.fastq.gz

    ## Encountered 122168 unique sequences from 237632 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA433243//filtered/SRR6715673_F_filt.fastq.gz

    ## Encountered 108563 unique sequences from 213783 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA433243//filtered/SRR6715674_F_filt.fastq.gz

    ## Encountered 79581 unique sequences from 162642 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA433243//filtered/SRR6715675_F_filt.fastq.gz

    ## Encountered 91781 unique sequences from 187919 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA433243//filtered/SRR6715676_F_filt.fastq.gz

    ## Encountered 113597 unique sequences from 226551 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA433243//filtered/SRR6715677_F_filt.fastq.gz

    ## Encountered 112800 unique sequences from 243098 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA433243//filtered/SRR6715678_F_filt.fastq.gz

    ## Encountered 73104 unique sequences from 172870 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Ocean-Crust-Synthesis/PRJNA433243//filtered/SRR6715679_F_filt.fastq.gz

    ## Encountered 128410 unique sequences from 260493 total sequences read.

``` r
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
```

``` r
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 246877 reads in 51143 unique sequences.
    ## Sample 2 - 293932 reads in 130964 unique sequences.
    ## Sample 3 - 359238 reads in 159867 unique sequences.
    ## Sample 4 - 363694 reads in 171953 unique sequences.
    ## Sample 5 - 329859 reads in 155229 unique sequences.
    ## Sample 6 - 161896 reads in 70739 unique sequences.
    ## Sample 7 - 176938 reads in 89346 unique sequences.
    ## Sample 8 - 237748 reads in 116426 unique sequences.
    ## Sample 9 - 183070 reads in 93554 unique sequences.
    ## Sample 10 - 240469 reads in 111385 unique sequences.
    ## Sample 11 - 265663 reads in 131764 unique sequences.
    ## Sample 12 - 237632 reads in 122168 unique sequences.
    ## Sample 13 - 213783 reads in 108563 unique sequences.
    ## Sample 14 - 162642 reads in 79581 unique sequences.
    ## Sample 15 - 187919 reads in 91781 unique sequences.
    ## Sample 16 - 226551 reads in 113597 unique sequences.
    ## Sample 17 - 243098 reads in 112800 unique sequences.
    ## Sample 18 - 172870 reads in 73104 unique sequences.
    ## Sample 19 - 260493 reads in 128410 unique sequences.

``` r
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
```

    ## [1]    19 22371

``` r
table(nchar(getSequences(seqtab)))
```

    ## 
    ##   220 
    ## 22371

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 435 bimeras out of 22371 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]    19 21936

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.9797772

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
```

``` r
ps_z18 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))
ps_z18
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 21936 taxa and 19 samples ]
    ## tax_table()   Taxonomy Table:    [ 21936 taxa by 6 taxonomic ranks ]

``` r
colnames(tax_table(ps_z18)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
```

``` r
meta <- read.csv("names.csv", row.names = 1)
nm <- sample_data(meta)
```

``` r
ps_z18d <- merge_phyloseq(ps_z18, nm)
```

Two Samples used in analysis taken from total

``` r
ps_z18_sm <- subset_samples(ps_z18d, sample_names(ps_z18d) == "SRR6715674" | sample_names(ps_z18d)== "SRR6715672")
ps_z18_sm
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 21936 taxa and 2 samples ]
    ## sample_data() Sample Data:       [ 2 samples by 1 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 21936 taxa by 6 taxonomic ranks ]

``` r
final_tab <- otu_table(ps_z18_sm)
final_tax <- tax_table(ps_z18_sm)
zinke_final <- merge_phyloseq(final_tab, final_tax)
zinke_final
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 21936 taxa and 2 samples ]
    ## tax_table()   Taxonomy Table:    [ 21936 taxa by 6 taxonomic ranks ]

``` r
saveRDS(zinke_final, file="zinke-final")
```
