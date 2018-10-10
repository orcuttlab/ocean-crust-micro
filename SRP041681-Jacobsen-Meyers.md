``` r
library(phyloseq)
```

``` r
setwd("~/Desktop/Ocean-Crust-Synthesis/SRP041681/")
```

See jacobsen-meyers-mothur-logfile.txt for documentation of mothur pipeline.

``` r
shared <- "~/Desktop/Ocean-Crust-Synthesis/SRP041681/mothur-retry/enzymesV6all.subsample.unique.precluster.an.unique_list.shared"
```

``` r
taxa <- "~/Desktop/Ocean-Crust-Synthesis/SRP041681/mothur-retry/enzymesV6all.subsample.unique.precluster.an.unique_list.unique.cons.taxonomy"
```

``` r
loihi_mothur <- import_mothur(mothur_shared_file = shared, mothur_constaxonomy_file = taxa)
```

``` r
loihi_mothur
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 5359 taxa and 6 samples ]
    ## tax_table()   Taxonomy Table:    [ 5359 taxa by 6 taxonomic ranks ]

``` r
saveRDS(loihi_mothur, file="~/Desktop/Ocean-Crust-Synthesis/final_RDS/loihi_mothur")
```
