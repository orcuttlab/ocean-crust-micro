``` r
library(phyloseq)
```

``` r
setwd("~/Desktop/Ocean-Crust-Synthesis/454s/")
```

Mothur logfile: smith-454-mothur-reproduce.txt

``` r
smith_shared <- "~/Desktop/Ocean-Crust-Synthesis/454s/smith-final-shared.txt"
smith_tax <- "~/Desktop/Ocean-Crust-Synthesis/454s/smith-final-taxonomy.txt"
smith454 <- import_mothur(mothur_shared_file = smith_shared, mothur_constaxonomy_file = smith_tax)
smith454
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 16381 taxa and 8 samples ]
    ## tax_table()   Taxonomy Table:    [ 16381 taxa by 6 taxonomic ranks ]

``` r
colnames(tax_table(smith454)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
```

Mothur logfile: jdf-labonte-reclassify-mothur-logfile.txt

``` r
jdf_shared <- "jdf-labonte-alternate-shared.txt"
jdf_tax <- "jdf-alternate-otu-tax.txt"
jdf_labonte <- import_mothur(mothur_shared_file = jdf_shared, mothur_constaxonomy_file = jdf_tax)
```

    ## Warning in readLines(mothur_shared_file): incomplete final line found on
    ## 'jdf-labonte-alternate-shared.txt'

``` r
colnames(tax_table(jdf_labonte)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
jdf_labonte
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 1809 taxa and 12 samples ]
    ## tax_table()   Taxonomy Table:    [ 1809 taxa by 6 taxonomic ranks ]

``` r
saveRDS(smith454, file = "smith454")
```

``` r
saveRDS(jdf_labonte, file = "~/Desktop/Ocean-Crust-Synthesis/final_RDS/jdf_labonte_alt")
```
