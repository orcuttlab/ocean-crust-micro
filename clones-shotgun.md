Libraries

``` r
library(reshape2)
library(ggplot2)
library(RColorBrewer)
```

``` r
setwd("~/Desktop/Ocean-Crust-Synthesis/clone_taxonomy/")
```

Clone library composition calculated as (\# of occurences of a given phyla) / (total number of clones)

``` r
clones <- as.data.frame(read.csv("all-clones.csv", check.names = FALSE))
cloneso <- clones[, order(names(clones))]
```

``` r
clones_long <- melt(cloneso, id.vars = "Study", variable.name = "Taxonomy", value.name = "value", factorsAsStrings = FALSE)
```

``` r
colmatch2 <- c("Acidobacteria" = "#e69878", "Actinobacteria" = "#502ebb", "Aerophobetes" = 
"#60c92c", "Proteobacteria-Alpha" = "#d13bd4", "Proteobacteria-Delta" = "#58ca68", "Proteobacteria-Gamma" = "#a04cde", "Proteobacteria-Zeta" =  "#54a534", "Bacteria" =  "#db60d3", "BRC1" = "#b6bf37", "Calditrichaeota" = "#5764de", "Chlamydiae" =  "#d9ab3e", "Chloroflexi" =  "#742a8c", "Crenarchaeota" =  "#8da945", "Cyanobacteria" = "#da37a2", "Dadabacteria" =  "#64bc87", "Dependentiae" =  "#e92829", "Entotheonellaeota" =  "#53beba", "Epsilonbacteraeota" = "#de4d2b", "Euryarchaeota" = "#D3D3D3", "Firmicutes" = "#6389db", "Fusobacteria" =  "#e48e2e", "Gemmatimonadetes" = "#474691", "Hydrogenedentes" = "#936cc8", "Kiritimatiellaeota" = "#45803d", "Latescibacteria" = "#ca7ad8", "Lentisphaerae" = "#6e7024", "Margulisbacteria" = "#dd316c", "Nanoarchaeaeota" = "#2b5126", "Nitrospinae" =  "#df7fbd", "Nitrospirae" =  "#3e8070", "Omnitrophicaeota" =  "#da4450", "Opisthokonta" =  "#61aad1", "Patescibacteria" =  "#d16d37", "PAUC34f" =  "#4b6891", "Picozoa" = "#a8732c", "Planctomycetes" = "#b199d2", "Rokubacteria" =  "#93341c", "SAR" = "#98ac7c", "Schekmanbacteria" = "#8a2564", "Spirochaetes" = "#ceac79", "Tenericutes" = "#4c315d", "Thaumarchaeota" =  "#9b7952", "unclassified" = "#d55b90", "Verrucomicrobia" =  "#4f4f25", "Acetothermia" =  "#ce6067", "Marinimicrobia_(SAR406_clade)" = "#5d2c24", "Hydrothermarchaeota" = "#d999ac", "Atribacteria" = "#744821", "< 2% Abundance" = "#9b5f84", "Bacteroidetes" = "#7c263a", "Unclassified Proteobacteria" = "#4b6f60", "AEGEAN-245 Proteobacteria" = "#d69d7e", "ARKICE-90 Proteobacteria" = "#59b2d2", "Elev-16S-509 Proteobacteria" = "#9e5f53", "SC3-20 Proteobacteria" = "#7bb6b4", "Candidate division OP3" = "#9b6787", "Cloacimonetes" = "#86b385", "Elusimicrobia" = "#636690", "Gracilibacteria" = "#5b7141", "Parcubacteria" = "#a9a5d0", "SHA-109" = "#856b45", "TM6" = "#4e7a8c", "Woesearchaeota" = "#d4a1ac", "Taxonomic Groups < 0.5%" = "#49917f", "Microgenomates" = "#826e68", "Proteobacteria-Beta" = "#aeab84", "TA 18 Proteobacteria" = "#6c4122", "Aquificae" = "#97ad62", "FCPU426" = "#a567e3", "Hadesarchaeaeota" = "#57a927", "Marinimicrobia" = "#8f41aa", "Thermotogae" = "#4fc758", "WS2" = "#cc47a9", "Asgardaeota" =  "#b5b540", "Bacteria_unclassified" = "#bf6130", "BHI80-139" = "#67af44", "Halanaerobiaeota" = "#787b35", "Aegiribacteria" = "#de3985", "Altiarchaeota"= "#5dae82", "AncK6" = "#ab2e7d", "Archaea_unclassified" = "#6c841f", "Armatimonadetes" = "#aa83e3", "Arthropoda" = "#be9e24", "Ascomycota"= "#598eec", "Basidiomycota" = "#d19421", "Chlorophyta_ph" = "#7353ac", "CK-2C2-2" = "#91a94e", "Deferribacteres" = "#cc79cd", "Deinococcus-Thermus" = "#8b7c1a", "Desantisbacteria" = "#3d66b2", "Diapherotrites" = "#d6852a", "Euglenozoa" = "#3e9dcd", "Eukaryota_unclassified" = "#ba9645", "FBP" = "#964894", "Fibrobacteres" = "#556c27", "GN01" =  "#eb75bc", "Hydrothermae" = "#358760", "LCP-89" = "#ce559b", "MAST-3" = "#286841","MAT-CR-M4-B07" = "#cc7dbb", "Ochrophyta" = "#645f20", "Parabasalia" = "#9a8ed5", "Poribacteria" = "#946427", "Protalveolata" = "#7293d2", "Prymnesiophyceae"  = "#877a3a", "Retaria" = "#61609e", "Synergistetes" = "#959e5e", "TA06" = "#7e4c82", "unknown_unclassified" = "#30aaa7", "WOR-1" = "#a03f6c", "WPS-2" = "#c4915b", "WS1" = "#a372af", "Zixibacteria" = "#dc84ac")
```

``` r
order <- c("Huber-JDF", "Jungbluth-JDF-13", "Orcutt-JDF", "Smith-JDF", "Baquarin-JDF", "Orcutt-Loihi", "Barco-Surface", "Sudek-Vailulu'uSeamount", "Sylvan-LauBasin", "Thorseth-KnipovichRidge", "Santelli-Loihi", "Santelli-EPR", "Fisk-Hawaii","Leysnes-MohnsRidge", "Lee-Dorado", "Nitahara-Takuyo")
```

Same taxa included as in Figure 3

``` r
clones5 <- as.data.frame(read.csv("all-clones-mod5.csv", check.names = FALSE))
clones5o <- clones5[, order(names(clones5))]
```

``` r
clones_long5 <- melt(clones5o, id.vars = "Study", variable.name = "Taxonomy", value.name = "value", factorsAsStrings = FALSE)
```

Code for Figure 3B

``` r
#fig3a = ggplot(clones_long5, aes(x = Study, y = value, fill = Taxonomy), ordered = FALSE)  + scale_x_discrete(limits = order)
#fig3b = fig3b + theme(axis.text.x = element_text(angle = -90, hjust = 0))
#fig3b = fig3b + theme(panel.background = element_rect(fill=NA))
#fig3b = fig3b + geom_bar(stat = "identity", position = "stack", color = "black")
#fig3b = fig3b + scale_fill_manual(values = colmatch2)
#cplot = cplot + theme(legend.key.size = unit(1, "cm"))
#fig3b = fig3b + scale_y_continuous(labels = scales::percent)
#fig3b
```

``` r
orderT <- c("1383C Shallow filter 2012", "1383C Shallow sled TP4", "1383C Shallow sled TP5", "1383C Shallow filter 2014", "1383C Middle sled TP4", "1383C Middle filter 2014", "1383C Deep filter 2012", "1383C Deep sled TP4", "1383C Deep filter 2014", "1382A filter 2012", "1382A sled TP1", "1382A sled TP2", "1382A sled TP3", "1382A sled TP4", "1382A sled TP5", "1382A sled TP6", "1382A sled TP7", "1382A sled TP8", "1382A filter 2014", "Bottomwater filter 2012", "Bottomwater filter 2014", "Jungbluth-JDF-16")
```

``` r
tully5 <- as.data.frame(read.csv("tulley-relab-mod5.csv", check.names = FALSE))
tully5o <- tully5[, order(names(tully5))]
```

``` r
tully_long5 <- melt(tully5o, id.vars = "Sample", variable.name = "Taxonomy")
```

Code for Figure 3C

``` r
#fig3c = ggplot(tully_long5, aes(x = Sample, y = value, fill = Taxonomy, width = 1), ordered = FALSE) + 
#    geom_bar(stat = "identity")
#fig3c = fig3c + theme(axis.text.x = element_text(angle = -90, hjust = 0)) + scale_x_discrete(limits = orderT)
#fig3c = fig3c + theme(panel.background = element_rect(fill=NA))
#fig3c = fig3c + geom_bar(stat = "identity", 
#                  position = "stack", color = "black")
#fig3c = fig3c + scale_fill_manual(values = colmatch2)
#fig3c
```
