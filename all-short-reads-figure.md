Necessary Libraries

``` r
library(phyloseq)
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

Working directory with saved phyloseq objects produced from sequence processing of individual datasets

``` r
setwd("~/Desktop/Ocean-Crust-Synthesis/final_RDS/")
```

Read in RDS files of phyloseq objects created for each data set

``` r
ps_j16 <- readRDS(file ="final_RDS/ps_j16")
ps_m12 <- readRDS(file ="final_RDS/ps_m12")
ps_m14 <- readRDS(file ="final_RDS/loihi_mothur")
ps_r17 <- readRDS(file ="final_RDS/ps_r17")
ps_l16 <- readRDS(file ="final_RDS/ps_l16")
ps_jf16 <- readRDS(file ="final_RDS/ps_jf16")
smith454 <- readRDS(file ="final_RDS/smith454")
jdf_labonte_alt <- readRDS(file ="final_RDS/jdf_labonte_alt")
zinke <- readRDS(file ="final_RDS/zinke-final")
```

merge them into one phyloseq object

``` r
testmerge <- merge_phyloseq(ps_j16, ps_m12)
merge2 <- merge_phyloseq(testmerge, ps_m14)
merge3 <- merge_phyloseq(merge2, ps_r17)
merge4 <- merge_phyloseq(merge3, ps_l16)
merge5 <- merge_phyloseq(merge4, ps_jf16)
merge6 <- merge_phyloseq(merge5, smith454)
merge7 <- merge_phyloseq(merge6, jdf_labonte_alt)
all_sr <- merge_phyloseq(merge7, zinke)
all_sr
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 64366 taxa and 116 samples ]
    ## tax_table()   Taxonomy Table:    [ 64366 taxa by 6 taxonomic ranks ]

Subset taxa without proteobacteria

``` r
sr_noprot <- subset_taxa(all_sr, Phylum != "Proteobacteria")
sr_noprot
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 39824 taxa and 116 samples ]
    ## tax_table()   Taxonomy Table:    [ 39824 taxa by 6 taxonomic ranks ]

Conglomerate at Phylum level

``` r
srglom_noprot <- tax_glom(sr_noprot, taxrank = "Phylum")
tax_table(srglom_noprot) <- tax_table(srglom_noprot)[,c(2)]
colnames(tax_table(srglom_noprot)) <- c("Taxonomy")
#head(tax_table(srglom_noprot))
```

Subset proteobacteria, conglomerate at class level

``` r
prot_sr <- subset_taxa(all_sr, Phylum=="Proteobacteria")
srglom_prot <- tax_glom(prot_sr, taxrank = "Class")
tax_table(srglom_prot) <- tax_table(srglom_prot)[,c(3)]
colnames(tax_table(srglom_prot)) <- c("Taxonomy")
#print(tax_table(srglom_prot))
```

Manually changing proteobacterial names

``` r
write.csv(tax_table(srglom_prot), file="~/Desktop/Ocean-Crust-Synthesis/final_RDS/glom_prot_tax.csv")
```

``` r
srtax_cust <- as.matrix(read.table(file="~/Desktop/Ocean-Crust-Synthesis/final_RDS/prot_tax_cust.txt"))
```

Create new OTU table with customized proteobacterial names

``` r
cust_otu <- otu_table(srglom_prot)
cust_tax = tax_table(srtax_cust)
```

Create new phyloseq object

``` r
srglom_prot_cust = phyloseq(cust_otu, cust_tax)
```

Merge all phyla with customized proteobacterial labels

``` r
srphy_prot <- merge_phyloseq(srglom_noprot, srglom_prot_cust)
srphy_prot
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 92 taxa and 116 samples ]
    ## tax_table()   Taxonomy Table:    [ 92 taxa by 1 taxonomic ranks ]

Convert to relative abundance

``` r
srphy_protrel <- transform_sample_counts(srphy_prot, function(x) x/sum(x))
srphy_protrel
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 92 taxa and 116 samples ]
    ## tax_table()   Taxonomy Table:    [ 92 taxa by 1 taxonomic ranks ]

read in meta data table

``` r
metasr <- read.csv("final_RDS/meta_allsr.csv", row.names = 1)
pmsr <- sample_data(metasr)
```

``` r
class_allsr <- merge_phyloseq(srphy_protrel, pmsr)
dim(sample_data(class_allsr))
```

    ## [1] 116   4

``` r
sample_names(class_allsr) <- pmsr$name
```

Taking out problematic samples

``` r
class_allsr = subset_samples(class_allsr, sample_names(class_allsr) != "Blank-1")
class_allsr = subset_samples(class_allsr, sample_names(class_allsr) != "Blank-2")
class_allsr = subset_samples(class_allsr, sample_names(class_allsr) != "Blank-3")
class_allsr = subset_samples(class_allsr, sample_names(class_allsr) != "Blank-4")
class_allsr = subset_samples(class_allsr, sample_names(class_allsr) != "JdFBack")
class_all1sr = subset_samples(class_allsr, sample_names(class_allsr) != "4HCC")
class_all1sr
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 92 taxa and 110 samples ]
    ## sample_data() Sample Data:       [ 110 samples by 4 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 92 taxa by 1 taxonomic ranks ]

Color Key

``` r
colmatch2 <- c("Acidobacteria" = "#e69878", "Actinobacteria" = "#502ebb", "Aerophobetes" = 
"#60c92c", "Proteobacteria-Alpha" = "#d13bd4", "Proteobacteria-Delta" = "#58ca68", "Proteobacteria-Gamma" = "#a04cde", "Proteobacteria-Zeta" =  "#54a534", "Bacteria" =  "#db60d3", "BRC1" = "#b6bf37", "Calditrichaeota" = "#5764de", "Chlamydiae" =  "#d9ab3e", "Chloroflexi" =  "#742a8c", "Crenarchaeota" =  "#8da945", "Cyanobacteria" = "#da37a2", "Dadabacteria" =  "#64bc87", "Dependentiae" =  "#e92829", "Entotheonellaeota" =  "#53beba", "Epsilonbacteraeota" = "#de4d2b", "Euryarchaeota" = "#D3D3D3", "Firmicutes" = "#6389db", "Fusobacteria" =  "#e48e2e", "Gemmatimonadetes" = "#474691", "Hydrogenedentes" = "#936cc8", "Kiritimatiellaeota" = "#45803d", "Latescibacteria" = "#ca7ad8", "Lentisphaerae" = "#6e7024", "Margulisbacteria" = "#dd316c", "Nanoarchaeaeota" = "#2b5126", "Nitrospinae" =  "#df7fbd", "Nitrospirae" =  "#3e8070", "Omnitrophicaeota" =  "#da4450", "Opisthokonta" =  "#61aad1", "Patescibacteria" =  "#d16d37", "PAUC34f" =  "#4b6891", "Picozoa" = "#a8732c", "Planctomycetes" = "#b199d2", "Rokubacteria" =  "#93341c", "SAR" = "#98ac7c", "Schekmanbacteria" = "#8a2564", "Spirochaetes" = "#ceac79", "Tenericutes" = "#4c315d", "Thaumarchaeota" =  "#9b7952", "unclassified" = "#d55b90", "Verrucomicrobia" =  "#4f4f25", "Acetothermia" =  "#ce6067", "Marinimicrobia_(SAR406_clade)" = "#5d2c24", "Hydrothermarchaeota" = "#d999ac", "Atribacteria" = "#744821", "< 2% Abundance" = "#9b5f84", "Bacteroidetes" = "#7c263a", "Unclassified Proteobacteria" = "#4b6f60", "AEGEAN-245 Proteobacteria" = "#d69d7e", "ARKICE-90 Proteobacteria" = "#59b2d2", "Elev-16S-509 Proteobacteria" = "#9e5f53", "SC3-20 Proteobacteria" = "#7bb6b4", "Candidate division OP3" = "#9b6787", "Cloacimonetes" = "#86b385", "Elusimicrobia" = "#636690", "Gracilibacteria" = "#5b7141", "Parcubacteria" = "#a9a5d0", "SHA-109" = "#856b45", "TM6" = "#4e7a8c", "Woesearchaeota" = "#d4a1ac", "Taxonomic Groups < 0.5%" = "#49917f", "Microgenomates" = "#826e68", "Proteobacteria-Beta" = "#aeab84", "TA 18 Proteobacteria" = "#6c4122", "Aquificae" = "#97ad62", "FCPU426" = "#a567e3", "Hadesarchaeaeota" = "#57a927", "Marinimicrobia" = "#8f41aa", "Thermotogae" = "#4fc758", "WS2" = "#cc47a9", "Asgardaeota" =  "#b5b540", "Bacteria_unclassified" = "#bf6130", "BHI80-139" = "#67af44", "Halanaerobiaeota" = "#787b35", "Aegiribacteria" = "#de3985", "Altiarchaeota"= "#5dae82", "AncK6" = "#ab2e7d", "Archaea_unclassified" = "#6c841f", "Armatimonadetes" = "#aa83e3", "Arthropoda" = "#be9e24", "Ascomycota"= "#598eec", "Basidiomycota" = "#d19421", "Chlorophyta_ph" = "#7353ac", "CK-2C2-2" = "#91a94e", "Deferribacteres" = "#cc79cd", "Deinococcus-Thermus" = "#8b7c1a", "Desantisbacteria" = "#3d66b2", "Diapherotrites" = "#d6852a", "Euglenozoa" = "#3e9dcd", "Eukaryota_unclassified" = "#ba9645", "FBP" = "#964894", "Fibrobacteres" = "#556c27", "GN01" =  "#eb75bc", "Hydrothermae" = "#358760", "LCP-89" = "#ce559b", "MAST-3" = "#286841","MAT-CR-M4-B07" = "#cc7dbb", "Ochrophyta" = "#645f20", "Parabasalia" = "#9a8ed5", "Poribacteria" = "#946427", "Protalveolata" = "#7293d2", "Prymnesiophyceae"  = "#877a3a", "Retaria" = "#61609e", "Synergistetes" = "#959e5e", "TA06" = "#7e4c82", "unknown_unclassified" = "#30aaa7", "WOR-1" = "#a03f6c", "WPS-2" = "#c4915b", "WS1" = "#a372af", "Zixibacteria" = "#dc84ac", "Proteobacteria-Unclassified" = "#dc84ad")
```

To simplify the graph all phyla whose minimum abundance in any sample is greater than 5% is subsetted. This was determined by abundance filtering taxa that have an abundance of 5% in at least one sample. Taxa for Figure 3

``` r
high_5 <- subset_taxa(class_all1sr, Taxonomy=="Acetothermia" | Taxonomy=="Acidobacteria" | Taxonomy=="Actinobacteria" | Taxonomy=="Aerophobetes" | Taxonomy == "Asgardaeota" | Taxonomy=="Atribacteria"  | Taxonomy=="Bacteroidetes" | Taxonomy == "Bacteria_unclassified" |  Taxonomy=="Chloroflexi" | Taxonomy=="Crenarchaeota" | Taxonomy=="Cyanobacteria"| Taxonomy == "Elusimicrobia" | Taxonomy=="Epsilonbacteraeota" | Taxonomy=="Euryarchaeota" | Taxonomy=="Firmicutes" | Taxonomy=="Gemmatimonadetes" | Taxonomy=="Hadesarchaeaeota" | Taxonomy=="Hydrothermarchaeota" | Taxonomy=="Marinimicrobia_(SAR406_clade)" |  Taxonomy=="Nitrospirae" | Taxonomy=="Planctomycetes" | Taxonomy=="Proteobacteria-Alpha" | Taxonomy=="Proteobacteria-Delta" | Taxonomy=="Proteobacteria-Gamma" | Taxonomy == "Proteobacteria-Unclassified" | Taxonomy=="Schekmanbacteria" | Taxonomy=="Tenericutes" | Taxonomy=="Thaumarchaeota" | Taxonomy =="Thermotogae" | Taxonomy=="unknown_unclassified")
high_5
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 30 taxa and 110 samples ]
    ## sample_data() Sample Data:       [ 110 samples by 4 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 30 taxa by 1 taxonomic ranks ]

``` r
himeltf5 <- psmelt(high_5)
himeltf5$Taxonomy <- as.character(himeltf5$Taxonomy) 
himeltf5$fac = factor(himeltf5$type, levels=c("Basalt", "fluids", "sediment", "seawater"))
```

Code to create high abundance bar chart

``` r
#fig3 = ggplot(himeltf5[order(himeltf5$Taxonomy),], aes_string(x = "Sample", y = "Abundance", fill = "Taxonomy"), ordered = TRUE)

#fig3 = fig3 + geom_bar(stat = "identity", 
 #                 position = "stack", color = "black") 

#fig3 = fig3 + scale_fill_manual(values = colmatch2)

#fig3 = fig3 + theme(axis.text.x = element_text(angle = -90, hjust = 0)) #theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())

#fig3 = fig3 + guides(fill = guide_legend(override.aes = list(colour = NULL), reverse=FALSE))  + 
#    theme(legend.key = element_rect(colour = "black")) + #scale_y_continuous(expand=c(0,0), limits = c(0,1)) + theme(panel.background = element_rect(fill=NA), axis.line = element_line(size=0.3), strip.background = element_rect(fill=NA)) + facet_grid(.~type, scales = "free") +
#   ylab("Relative Abundance") + xlab("Sample") 

#fig3 = fig3 + theme(panel.background = element_rect(fill=NA),
#                  axis.text = element_text(size=10),
#                  axis.title = element_text(size=11), axis.line = #element_line(colour = 'black', size=.2)) + 
#   scale_y_continuous(expand=c(0.0002,.0002), limits = c(0,1)) +
#  ylab("Relative Abundance")
  
#fig3 = fig3 + facet_grid(. ~fac, scales = "free_x", space = "free_x", as.ordered(fac))

#fig3
```

All taxa whose maximum abundance in any sample is less than 5% is subsetted to make another graph of low abundance taxa Taxa for Figure 4

``` r
low_5 <- subset_taxa(class_all1sr, Taxonomy=="Aegiribacteria" |  Taxonomy=="Altiarchaeota" | Taxonomy=="AncK6" | Taxonomy=="Aquificae" | Taxonomy=="Archaea_unclassified" | Taxonomy=="Armatimonadetes" | Taxonomy=="Arthropoda" | Taxonomy=="Ascomycota" | Taxonomy=="BRC1"| Taxonomy=="BHI80-139" | Taxonomy=="Basidiomycota" | Taxonomy=="Calditrichaeota" | Taxonomy=="CK-2C2-2" | Taxonomy=="Chlamydiae" | Taxonomy=="Chlorophyta_ph" | Taxonomy=="Cloacimonetes" | Taxonomy=="Dadabacteria" | Taxonomy=="Deferribacteres" | Taxonomy=="Deinococcus-Thermus" | Taxonomy=="Dependentiae" | Taxonomy=="Desantisbacteria" | Taxonomy=="Diapherotrites" | Taxonomy=="Entotheonellaeota" | Taxonomy=="Euglenozoa" | Taxonomy=="Eukaryota_unclassified" | Taxonomy=="FBP" | Taxonomy=="FCPU426" | Taxonomy=="Fibrobacteres" | Taxonomy=="Fusobacteria" | Taxonomy=="GN01" | Taxonomy =="Halanaerobiaeota" | Taxonomy=="Hydrogenedentes" | Taxonomy=="Hydrothermae" | Taxonomy=="Kiritimatiellaeota" | Taxonomy=="LCP-89" | Taxonomy=="Latescibacteria" | Taxonomy=="Lentisphaerae" | Taxonomy=="MAST-3" | Taxonomy=="MAT-CR-M4-B07" | Taxonomy=="Margulisbacteria" | Taxonomy=="Nanoarchaeaeota" | Taxonomy=="Nitrospinae" | Taxonomy=="Ochrophyta" | Taxonomy =="Omnitrophicaeota" | Taxonomy=="PAUC34f" | Taxonomy =="Parabasalia" | Taxonomy =="Poribacteria" | Taxonomy =="Protalveolata" | Taxonomy=="Proteobacteria-Zeta" |Taxonomy == "Patescibacteria" | Taxonomy=="Prymnesiophyceae" | Taxonomy=="Retaria" | Taxonomy=="Rokubacteria" | Taxonomy=="Spirochaetes" | Taxonomy=="Synergistetes" | Taxonomy=="TA06" | Taxonomy=="WOR-1" | Taxonomy=="WPS-2" | Taxonomy=="WS1" | Taxonomy=="WS2" | Taxonomy=="Verrucomicrobia" | Taxonomy=="Zixibacteria")
low_5
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 62 taxa and 110 samples ]
    ## sample_data() Sample Data:       [ 110 samples by 4 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 62 taxa by 1 taxonomic ranks ]

``` r
lomeltf5 <- psmelt(low_5)
lomeltf5$Taxonomy <- as.character(lomeltf5$Taxonomy) 
lomeltf5$fac = factor(lomeltf5$type, levels=c("Basalt", "fluids", "sediment", "seawater"))
```

Code to create low abundance bar chart

``` r
#fig4 = ggplot(lomeltf5[order(lomeltf5$Taxonomy),], aes_string(x = "Sample", y = "Abundance", fill = "Taxonomy"), environment = .e, ordered = TRUE)

#fig4 = fig4 + geom_bar(stat = "identity", 
#                  position = "stack", color = "black") 

#fig4 = fig4 + scale_fill_manual(values = colmatch2)

#fig4 = fig4 + theme(axis.text.x = element_text(angle = -90, hjust = 0)) + theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())

#fig4 = fig4 + guides(fill = guide_legend(override.aes = list(colour = NULL), reverse=FALSE))  + 
#    theme(legend.key = element_rect(colour = "black")) + scale_y_continuous(expand=c(0,0), limits = c(0,1)) + theme(panel.background = element_rect(fill=NA), axis.line = element_line(size=0.3), strip.background = element_rect(fill=NA)) + facet_grid(.~type, scales = "free") +
#   ylab("Relative Abundance") + xlab("Sample") 

#fig4 = fig4 + theme(panel.background = element_rect(fill=NA),
 #                 axis.text = element_text(size=10),
#                  axis.title = element_text(size=11), axis.line = element_line(colour = 'black', size=.2)) + 
#    scale_y_continuous(expand=c(0.0002,.0002), limits = c(0,.1)) +
#  ylab("Relative Abundance")

#fig4 = fig4 + facet_grid(. ~fac, scales = "free_x", space = "free_x", as.ordered(fac))

#fig4
```

code for gammaproteobacterial families

``` r
allsr_rel <- transform_sample_counts(all_sr, function(x) x/sum(x))
all_gamma_rel <- subset_taxa(allsr_rel, Class=="Gammaproteobacteria")
gamma_fam <- tax_glom(all_gamma_rel, taxrank = "Family")
```

``` r
gf_meta <- merge_phyloseq(gamma_fam, pmsr)
```

``` r
sample_names(gf_meta) <- pmsr$name
```

Taking out problematic samples

``` r
gf_meta = subset_samples(gf_meta, sample_names(gf_meta) != "Blank-1")
gf_meta = subset_samples(gf_meta, sample_names(gf_meta) != "Blank-2")
gf_meta = subset_samples(gf_meta, sample_names(gf_meta) != "Blank-3")
gf_meta = subset_samples(gf_meta, sample_names(gf_meta) != "Blank-4")
gf_meta = subset_samples(gf_meta, sample_names(gf_meta) != "JdFBack")
gf_meta1 = subset_samples(gf_meta, sample_names(gf_meta) != "4HCC")
gf_meta1
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 91 taxa and 110 samples ]
    ## sample_data() Sample Data:       [ 110 samples by 4 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 91 taxa by 6 taxonomic ranks ]

Gammaprotoeobacterial families that have an abundance of 2% in at least one sample

``` r
gf_2 <- subset_taxa(gf_meta1, Family=="Acidiferrobacteraceae" | Family=="Burkholderiaceae" | Family=="Aeromonadaceae" | Family=="Alcanivoracaceae" | Family == "Alteromonadaceae" | Family=="Cellvibrionaceae"  | Family=="Colwelliaceae" | Family == "Enterobacteriaceae" |  Family=="Gammaproteobacteria_unclassified" | Family=="Halieaceae" | Family=="Halomonadaceae"| Family == "Idiomarinaceae" | Family=="Methylomonaceae" | Family=="Moraxellaceae" | Family=="Neisseriaceae" | Family=="Nitrosococcaceae" | Family=="Nitrosomonadaceae" | Family=="Porticoccaceae" | Family=="Pseudoalteromonadaceae" |  Family=="Pseudohongiellaceae" | Family=="Pseudomonadaceae" | Family=="Saccharospirillaceae" | Family=="Sedimenticolaceae" | Family=="Solimonadaceae" | Family == "Spongiibacteraceae" | Family=="Thioalkalispiraceae" | Family=="Thioglobaceae" | Family=="Thiomicrospiraceae" | Family =="Thiotrichaceae" | Family=="Unknown_Family" | Family=="Woeseiaceae" | Family=="Xanthomonadaceae")
gf_2
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 32 taxa and 110 samples ]
    ## sample_data() Sample Data:       [ 110 samples by 4 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 32 taxa by 6 taxonomic ranks ]

``` r
meltdg2 <- psmelt(gf_2)
meltdg2$Family <- as.character(meltdg2$Family) 
meltdg2$fac = factor(meltdg2$type, levels=c("Basalt", "fluids", "sediment", "seawater"))
```

code for gammaproteobacterial familiy plot

``` r
#fig5 = ggplot(meltdg2[order(meltdg2$Family),], aes_string(x = "Sample", y = "Abundance", fill = "Family"), environment = .e, ordered = TRUE)

#fig5 = fig5 + geom_bar(stat = "identity", 
#                  position = "stack", color = "black") 

#fig5 = fig5 + scale_fill_manual(values = pal46)

#fig5 = fig5 + theme(axis.text.x = element_text(angle = -90, hjust = 0)) + theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())# - for no names

#fig5 = fig5 + guides(fill = guide_legend(override.aes = list(colour = NULL), reverse=FALSE))  + 
#    theme(legend.key = element_rect(colour = "black")) + scale_y_continuous(expand=c(0,0), limits = c(0,1)) + theme(panel.background = element_rect(fill=NA), axis.line = element_line(size=0.3), strip.background = element_rect(fill=NA)) + facet_grid(.~type, scales = "free") +
 #  ylab("Relative Abundance") + xlab("Sample") 

#fig5 = fig5 + theme(panel.background = element_rect(fill=NA),
#                  axis.text = element_text(size=10),
#                 axis.title = element_text(size=11), axis.line = element_line(colour = 'black', size=.2)) + 
 #   scale_y_continuous(expand=c(0.0002,.0002), limits = c(0,.8)) +
 # ylab("Relative Abundance")

#fig5 = fig5 + facet_grid(. ~fac, scales = "free_x", space = "free_x", as.ordered(fac))

#fig5
```
