``` r
library(phyloseq)
library(vegan)
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.5-2

``` r
library(ggplot2)
library(plyr)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:plyr':
    ## 
    ##     arrange, count, desc, failwith, id, mutate, rename, summarise,
    ##     summarize

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
setwd("~/Desktop/Ocean-Crust-Synthesis/")
```

Merged Phyloseq object from bar chart analysis

``` r
all_sr <- readRDS(file="~/Desktop/Ocean-Crust-Synthesis/final_RDS/all_sr")
```

``` r
all_sr_fam <- tax_glom(all_sr, taxrank = "Family")
all_fam_rel <- transform_sample_counts(all_sr_fam, function(x) x/sum(x))
#plot_bar(all_fam_rel, fill = "Family")
```

``` r
metasr1 <- read.csv("~/Desktop/Ocean-Crust-Synthesis/final_RDS/meta_allsr1.csv", row.names = 1)
pmsr1 <- sample_data(metasr1)
dim(pmsr1)
```

    ## [1] 116   7

``` r
all_fam_rel <- merge_phyloseq(all_fam_rel, pmsr1)
dim(sample_data(all_fam_rel))
```

    ## [1] 116   7

``` r
all_fam_rel = subset_samples(all_fam_rel, sample_names(all_fam_rel) != "Blank-1")
all_fam_rel = subset_samples(all_fam_rel, sample_names(all_fam_rel) != "Blank-2")
all_fam_rel = subset_samples(all_fam_rel, sample_names(all_fam_rel) != "Blank-3")
all_fam_rel = subset_samples(all_fam_rel, sample_names(all_fam_rel) != "Blank-4")
all_fam_rel = subset_samples(all_fam_rel, sample_names(all_fam_rel) != "JdFBack")
all_fam1 = subset_samples(all_fam_rel, sample_names(all_fam_rel) != "4HCC")
all_fam1
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 679 taxa and 110 samples ]
    ## sample_data() Sample Data:       [ 110 samples by 7 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 679 taxa by 6 taxonomic ranks ]

run this

``` r
#fam_mds <- ordinate(all_fam1, method = "NMDS", distance = "bray")
```

``` r
#fam_DF <- plot_ordination(physeq=all_fam_rel, ordination=fam_mds, justDF = TRUE)
```

``` r
nmdscols <- c("#688fce", "#89903c", "#66b24b", "#7664ca", "#c062bb", "#4aac8b", "#ca5336", "#c85979", "#c98e44")
```

to plot

``` r
#ntest3 <- ggplot(data=fam_DF, aes(y =NMDS2, x = NMDS1))
```

``` r
#ntest3 + geom_point(size = 4, aes(shape=type, color = study, fill=study, alpha=biome)) + scale_color_manual(values=nmdscols) + scale_fill_manual(values=nmdscols)  + scale_shape_manual(values=c(21,22,23,24,25)) + scale_alpha_manual(values=c("subsurface" = 1, "surface" = 0.5 )) + theme_classic()
```
