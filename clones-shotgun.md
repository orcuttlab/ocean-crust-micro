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
tully5
```

    ##                       Sample Proteobacteria-Alpha Proteobacteria-Beta
    ## 1  1383C Shallow filter 2012             15.06060              7.2188
    ## 2     1383C Shallow sled TP4             49.82580              0.0000
    ## 3     1383C Shallow sled TP5              8.30540              0.0000
    ## 4  1383C Shallow filter 2014             16.64360              0.0000
    ## 5      1383C Middle sled TP4             11.64420              8.5268
    ## 6   1383C Middle filter 2014             11.68620              0.0000
    ## 7     1383C Deep filter 2012              2.56750              0.0000
    ## 8        1383C Deep sled TP4             12.11250              0.5420
    ## 9     1383C Deep filter 2014             15.72080              2.6183
    ## 10         1382A filter 2012             10.61140              0.0000
    ## 11            1382A sled TP1             23.17970              0.0000
    ## 12            1382A sled TP2             26.09700              1.6656
    ## 13            1382A sled TP3             20.26490              0.0000
    ## 14            1382A sled TP4             15.59620              0.0000
    ## 15            1382A sled TP5             45.19030              0.0000
    ## 16            1382A sled TP6             43.98420              0.6082
    ## 17            1382A sled TP7             22.67720              1.1838
    ## 18            1382A sled TP8             32.47070              1.3046
    ## 19         1382A filter 2014             18.59360              0.0000
    ## 20   Bottomwater filter 2012             45.27070              0.0000
    ## 21   Bottomwater filter 2014             12.42200              0.0000
    ## 22          Jungbluth-JDF-16              9.71867              0.0000
    ##    Proteobacteria-Delta Epsilonbacteraeota Proteobacteria-Gamma
    ## 1              0.691600           5.769000             55.85630
    ## 2              0.000000           3.273900             38.69590
    ## 3             12.087200          46.513400             23.41720
    ## 4              9.263900           3.255300             25.44390
    ## 5              0.000000          29.412400             35.63230
    ## 6              0.000000           1.975600             74.53670
    ## 7              0.579000          77.810100              8.03240
    ## 8              1.166400          43.136700             22.01130
    ## 9              2.701900           7.456300             20.24090
    ## 10             3.308600          28.210800             45.38350
    ## 11             0.000000           0.000000             21.26860
    ## 12             2.593000          14.451900             35.43540
    ## 13             0.000000           1.266300             47.72630
    ## 14            49.259100           1.834000             17.19710
    ## 15             0.000000           4.509500             44.51960
    ## 16             0.000000           1.628200             23.02970
    ## 17             5.105000          15.768700             29.14820
    ## 18             2.527900           9.667200             34.21490
    ## 19             0.000000           0.682700             71.76990
    ## 20             7.370900           0.000000             15.01720
    ## 21             5.682400           0.000000             55.69930
    ## 22             3.324808           1.790281             15.08951
    ##    Proteobacteria-Zeta Acetothermia Acidobacteria Actinobacteria
    ## 1             0.000000     0.000000      0.000000       0.000000
    ## 2             0.705600     0.000000      0.000000       0.000000
    ## 3             1.033300     0.000000      0.000000       0.000000
    ## 4             1.487000     0.000000      0.000000       0.000000
    ## 5             0.000000     0.000000      0.000000       0.000000
    ## 6             0.000000     0.000000      0.000000       0.572900
    ## 7             0.787600     0.000000      0.000000       0.000000
    ## 8             4.712300     0.000000      0.000000       0.000000
    ## 9            13.747400     0.000000      0.000000       0.000000
    ## 10            1.007500     0.000000      0.000000       0.000000
    ## 11            0.000000     0.000000      0.000000       0.000000
    ## 12            1.490800     0.000000      0.000000       0.000000
    ## 13            0.000000     0.000000      0.000000       0.000000
    ## 14            0.000000     0.000000      0.000000       0.000000
    ## 15            0.000000     0.000000      0.000000       0.000000
    ## 16            0.000000     0.000000      0.000000       0.000000
    ## 17            5.161300     0.000000      0.000000       0.000000
    ## 18            4.969100     0.000000      0.000000       0.000000
    ## 19            0.507000     0.000000      0.000000       0.880600
    ## 20            0.000000     0.000000      0.512700       2.635400
    ## 21            0.000000     0.000000      0.000000       0.000000
    ## 22            3.324808     5.626598      1.790281       1.790281
    ##    Bacteroidetes Chloroflexi Cyanobacteria Crenarchaeota Euryarchaeota
    ## 1      10.812000    0.000000      0.000000      0.000000       0.00000
    ## 2       0.000000    0.000000      0.000000      0.000000       0.00000
    ## 3       3.201200    0.000000      0.000000      0.000000       0.00000
    ## 4       7.133000    2.237400      0.000000      0.000000       0.00000
    ## 5       5.748900    0.000000      0.669800      0.000000       0.00000
    ## 6       4.389300    0.000000      0.000000      0.000000       0.00000
    ## 7       2.268600    0.000000      0.000000      0.000000       0.00000
    ## 8       5.950300    0.000000      0.000000      0.000000       0.00000
    ## 9      12.810500    0.644400      0.000000      0.000000       0.00000
    ## 10      7.475700    0.000000      0.000000      0.000000       0.00000
    ## 11     47.935800    0.000000      0.000000      0.000000       0.00000
    ## 12      5.594000    0.000000      0.000000      0.000000       0.00000
    ## 13     28.828100    0.000000      0.000000      0.000000       0.00000
    ## 14     13.715900    0.000000      0.000000      0.000000       0.00000
    ## 15      3.576500    0.000000      0.000000      0.000000       0.00000
    ## 16     14.762100    0.000000      0.000000      0.000000       0.00000
    ## 17      3.550800    0.000000      0.000000      0.000000       0.00000
    ## 18      5.357500    0.000000      0.000000      0.000000       0.00000
    ## 19      1.531700    0.000000      0.000000      0.000000       0.00000
    ## 20      4.958700    6.831600      0.000000      0.000000       0.00000
    ## 21      1.413600    4.080100      0.000000      0.000000       1.36880
    ## 22      0.511509    4.347826      0.511509      5.370844      12.78772
    ##    Hydrothermarchaeota Firmicutes Marinimicrobia_(SAR406_clade)
    ## 1             0.000000   0.000000                        0.0000
    ## 2             0.000000   0.000000                        0.0000
    ## 3             0.000000   0.000000                        0.0000
    ## 4             0.000000  17.401900                        0.0000
    ## 5             0.000000   1.455500                        0.0000
    ## 6             0.000000   0.000000                        0.0000
    ## 7             0.000000   0.000000                        0.0000
    ## 8             0.000000   0.000000                        0.0000
    ## 9             0.000000   0.994800                        0.0000
    ## 10            0.000000   0.000000                        0.0000
    ## 11            0.000000   0.000000                        0.0000
    ## 12            0.000000   0.000000                        0.0000
    ## 13            0.000000   0.000000                        0.0000
    ## 14            0.000000   0.000000                        0.0000
    ## 15            0.000000   0.000000                        0.0000
    ## 16            0.000000   0.000000                        0.0000
    ## 17            0.000000   0.958300                        0.0000
    ## 18            0.000000   0.843500                        0.0000
    ## 19            0.000000   0.000000                        0.0000
    ## 20            0.000000   0.000000                        2.7394
    ## 21            0.000000   0.000000                        2.7529
    ## 22            2.557545   1.790281                        0.0000
    ##    Nitrospirae Planctomycetes Tenericutes Thaumarchaeota
    ## 1     0.000000         1.3400   0.0000000         0.0000
    ## 2     0.000000         0.0000   0.0000000         0.0000
    ## 3     0.000000         1.0928   0.0000000         0.9765
    ## 4     0.000000         4.7227   0.0000000         4.5936
    ## 5     0.000000         2.8063   0.0000000         0.0000
    ## 6     0.000000         3.7823   0.0000000         1.4398
    ## 7     0.000000         0.0000   0.0000000         0.0000
    ## 8     0.000000         3.0823   0.0000000         0.5906
    ## 9     0.000000         4.9188   0.0000000         8.5079
    ## 10    0.000000         0.0000   0.0000000         0.0000
    ## 11    0.000000         0.0000   0.0000000         0.0000
    ## 12    0.000000         2.6197   0.0000000         0.0000
    ## 13    0.000000         0.0000   0.0000000         0.0000
    ## 14    0.000000         0.0000   0.0000000         0.0000
    ## 15    0.000000         0.0000   0.0000000         0.0000
    ## 16    0.000000         4.7833   0.0000000         0.0000
    ## 17    0.000000         2.6347   0.0000000         1.3407
    ## 18    0.000000         2.6276   0.0000000         1.3537
    ## 19    0.000000         0.0000   0.0000000         0.5793
    ## 20    0.000000         3.9515   0.0000000         4.7726
    ## 21    0.000000         4.9223   0.0000000         5.2083
    ## 22    3.324808         0.0000   0.2557545         0.0000
    ##    unknown_unclassified         
    ## 1                  0.00       NA
    ## 2                  0.00       NA
    ## 3                  0.00       NA
    ## 4                  0.00       NA
    ## 5                  0.00       NA
    ## 6                  0.00       NA
    ## 7                  0.00       NA
    ## 8                  0.00       NA
    ## 9                  0.00       NA
    ## 10                 0.00       NA
    ## 11                 0.00       NA
    ## 12                 0.00       NA
    ## 13                 0.00       NA
    ## 14                 0.00       NA
    ## 15                 0.00       NA
    ## 16                 0.00       NA
    ## 17                 0.00       NA
    ## 18                 0.00       NA
    ## 19                 0.00       NA
    ## 20                 0.00       NA
    ## 21                 0.00       NA
    ## 22                26.06 99.97304

``` r
tully5o <- tully5[, order(names(tully5))]
tully5o
```

    ##             Acetothermia Acidobacteria Actinobacteria Bacteroidetes
    ## 1        NA     0.000000      0.000000       0.000000     10.812000
    ## 2        NA     0.000000      0.000000       0.000000      0.000000
    ## 3        NA     0.000000      0.000000       0.000000      3.201200
    ## 4        NA     0.000000      0.000000       0.000000      7.133000
    ## 5        NA     0.000000      0.000000       0.000000      5.748900
    ## 6        NA     0.000000      0.000000       0.572900      4.389300
    ## 7        NA     0.000000      0.000000       0.000000      2.268600
    ## 8        NA     0.000000      0.000000       0.000000      5.950300
    ## 9        NA     0.000000      0.000000       0.000000     12.810500
    ## 10       NA     0.000000      0.000000       0.000000      7.475700
    ## 11       NA     0.000000      0.000000       0.000000     47.935800
    ## 12       NA     0.000000      0.000000       0.000000      5.594000
    ## 13       NA     0.000000      0.000000       0.000000     28.828100
    ## 14       NA     0.000000      0.000000       0.000000     13.715900
    ## 15       NA     0.000000      0.000000       0.000000      3.576500
    ## 16       NA     0.000000      0.000000       0.000000     14.762100
    ## 17       NA     0.000000      0.000000       0.000000      3.550800
    ## 18       NA     0.000000      0.000000       0.000000      5.357500
    ## 19       NA     0.000000      0.000000       0.880600      1.531700
    ## 20       NA     0.000000      0.512700       2.635400      4.958700
    ## 21       NA     0.000000      0.000000       0.000000      1.413600
    ## 22 99.97304     5.626598      1.790281       1.790281      0.511509
    ##    Chloroflexi Crenarchaeota Cyanobacteria Epsilonbacteraeota
    ## 1     0.000000      0.000000      0.000000           5.769000
    ## 2     0.000000      0.000000      0.000000           3.273900
    ## 3     0.000000      0.000000      0.000000          46.513400
    ## 4     2.237400      0.000000      0.000000           3.255300
    ## 5     0.000000      0.000000      0.669800          29.412400
    ## 6     0.000000      0.000000      0.000000           1.975600
    ## 7     0.000000      0.000000      0.000000          77.810100
    ## 8     0.000000      0.000000      0.000000          43.136700
    ## 9     0.644400      0.000000      0.000000           7.456300
    ## 10    0.000000      0.000000      0.000000          28.210800
    ## 11    0.000000      0.000000      0.000000           0.000000
    ## 12    0.000000      0.000000      0.000000          14.451900
    ## 13    0.000000      0.000000      0.000000           1.266300
    ## 14    0.000000      0.000000      0.000000           1.834000
    ## 15    0.000000      0.000000      0.000000           4.509500
    ## 16    0.000000      0.000000      0.000000           1.628200
    ## 17    0.000000      0.000000      0.000000          15.768700
    ## 18    0.000000      0.000000      0.000000           9.667200
    ## 19    0.000000      0.000000      0.000000           0.682700
    ## 20    6.831600      0.000000      0.000000           0.000000
    ## 21    4.080100      0.000000      0.000000           0.000000
    ## 22    4.347826      5.370844      0.511509           1.790281
    ##    Euryarchaeota Firmicutes Hydrothermarchaeota
    ## 1        0.00000   0.000000            0.000000
    ## 2        0.00000   0.000000            0.000000
    ## 3        0.00000   0.000000            0.000000
    ## 4        0.00000  17.401900            0.000000
    ## 5        0.00000   1.455500            0.000000
    ## 6        0.00000   0.000000            0.000000
    ## 7        0.00000   0.000000            0.000000
    ## 8        0.00000   0.000000            0.000000
    ## 9        0.00000   0.994800            0.000000
    ## 10       0.00000   0.000000            0.000000
    ## 11       0.00000   0.000000            0.000000
    ## 12       0.00000   0.000000            0.000000
    ## 13       0.00000   0.000000            0.000000
    ## 14       0.00000   0.000000            0.000000
    ## 15       0.00000   0.000000            0.000000
    ## 16       0.00000   0.000000            0.000000
    ## 17       0.00000   0.958300            0.000000
    ## 18       0.00000   0.843500            0.000000
    ## 19       0.00000   0.000000            0.000000
    ## 20       0.00000   0.000000            0.000000
    ## 21       1.36880   0.000000            0.000000
    ## 22      12.78772   1.790281            2.557545
    ##    Marinimicrobia_(SAR406_clade) Nitrospirae Planctomycetes
    ## 1                         0.0000    0.000000         1.3400
    ## 2                         0.0000    0.000000         0.0000
    ## 3                         0.0000    0.000000         1.0928
    ## 4                         0.0000    0.000000         4.7227
    ## 5                         0.0000    0.000000         2.8063
    ## 6                         0.0000    0.000000         3.7823
    ## 7                         0.0000    0.000000         0.0000
    ## 8                         0.0000    0.000000         3.0823
    ## 9                         0.0000    0.000000         4.9188
    ## 10                        0.0000    0.000000         0.0000
    ## 11                        0.0000    0.000000         0.0000
    ## 12                        0.0000    0.000000         2.6197
    ## 13                        0.0000    0.000000         0.0000
    ## 14                        0.0000    0.000000         0.0000
    ## 15                        0.0000    0.000000         0.0000
    ## 16                        0.0000    0.000000         4.7833
    ## 17                        0.0000    0.000000         2.6347
    ## 18                        0.0000    0.000000         2.6276
    ## 19                        0.0000    0.000000         0.0000
    ## 20                        2.7394    0.000000         3.9515
    ## 21                        2.7529    0.000000         4.9223
    ## 22                        0.0000    3.324808         0.0000
    ##    Proteobacteria-Alpha Proteobacteria-Beta Proteobacteria-Delta
    ## 1              15.06060              7.2188             0.691600
    ## 2              49.82580              0.0000             0.000000
    ## 3               8.30540              0.0000            12.087200
    ## 4              16.64360              0.0000             9.263900
    ## 5              11.64420              8.5268             0.000000
    ## 6              11.68620              0.0000             0.000000
    ## 7               2.56750              0.0000             0.579000
    ## 8              12.11250              0.5420             1.166400
    ## 9              15.72080              2.6183             2.701900
    ## 10             10.61140              0.0000             3.308600
    ## 11             23.17970              0.0000             0.000000
    ## 12             26.09700              1.6656             2.593000
    ## 13             20.26490              0.0000             0.000000
    ## 14             15.59620              0.0000            49.259100
    ## 15             45.19030              0.0000             0.000000
    ## 16             43.98420              0.6082             0.000000
    ## 17             22.67720              1.1838             5.105000
    ## 18             32.47070              1.3046             2.527900
    ## 19             18.59360              0.0000             0.000000
    ## 20             45.27070              0.0000             7.370900
    ## 21             12.42200              0.0000             5.682400
    ## 22              9.71867              0.0000             3.324808
    ##    Proteobacteria-Gamma Proteobacteria-Zeta                    Sample
    ## 1              55.85630            0.000000 1383C Shallow filter 2012
    ## 2              38.69590            0.705600    1383C Shallow sled TP4
    ## 3              23.41720            1.033300    1383C Shallow sled TP5
    ## 4              25.44390            1.487000 1383C Shallow filter 2014
    ## 5              35.63230            0.000000     1383C Middle sled TP4
    ## 6              74.53670            0.000000  1383C Middle filter 2014
    ## 7               8.03240            0.787600    1383C Deep filter 2012
    ## 8              22.01130            4.712300       1383C Deep sled TP4
    ## 9              20.24090           13.747400    1383C Deep filter 2014
    ## 10             45.38350            1.007500         1382A filter 2012
    ## 11             21.26860            0.000000            1382A sled TP1
    ## 12             35.43540            1.490800            1382A sled TP2
    ## 13             47.72630            0.000000            1382A sled TP3
    ## 14             17.19710            0.000000            1382A sled TP4
    ## 15             44.51960            0.000000            1382A sled TP5
    ## 16             23.02970            0.000000            1382A sled TP6
    ## 17             29.14820            5.161300            1382A sled TP7
    ## 18             34.21490            4.969100            1382A sled TP8
    ## 19             71.76990            0.507000         1382A filter 2014
    ## 20             15.01720            0.000000   Bottomwater filter 2012
    ## 21             55.69930            0.000000   Bottomwater filter 2014
    ## 22             15.08951            3.324808          Jungbluth-JDF-16
    ##    Tenericutes Thaumarchaeota unknown_unclassified
    ## 1    0.0000000         0.0000                 0.00
    ## 2    0.0000000         0.0000                 0.00
    ## 3    0.0000000         0.9765                 0.00
    ## 4    0.0000000         4.5936                 0.00
    ## 5    0.0000000         0.0000                 0.00
    ## 6    0.0000000         1.4398                 0.00
    ## 7    0.0000000         0.0000                 0.00
    ## 8    0.0000000         0.5906                 0.00
    ## 9    0.0000000         8.5079                 0.00
    ## 10   0.0000000         0.0000                 0.00
    ## 11   0.0000000         0.0000                 0.00
    ## 12   0.0000000         0.0000                 0.00
    ## 13   0.0000000         0.0000                 0.00
    ## 14   0.0000000         0.0000                 0.00
    ## 15   0.0000000         0.0000                 0.00
    ## 16   0.0000000         0.0000                 0.00
    ## 17   0.0000000         1.3407                 0.00
    ## 18   0.0000000         1.3537                 0.00
    ## 19   0.0000000         0.5793                 0.00
    ## 20   0.0000000         4.7726                 0.00
    ## 21   0.0000000         5.2083                 0.00
    ## 22   0.2557545         0.0000                26.06

``` r
tully_long5 <- melt(tully5o, id.vars = "Sample", variable.name = "Taxonomy")
tully_long5
```

    ##                        Sample                      Taxonomy      value
    ## 1   1383C Shallow filter 2012                                       NA
    ## 2      1383C Shallow sled TP4                                       NA
    ## 3      1383C Shallow sled TP5                                       NA
    ## 4   1383C Shallow filter 2014                                       NA
    ## 5       1383C Middle sled TP4                                       NA
    ## 6    1383C Middle filter 2014                                       NA
    ## 7      1383C Deep filter 2012                                       NA
    ## 8         1383C Deep sled TP4                                       NA
    ## 9      1383C Deep filter 2014                                       NA
    ## 10          1382A filter 2012                                       NA
    ## 11             1382A sled TP1                                       NA
    ## 12             1382A sled TP2                                       NA
    ## 13             1382A sled TP3                                       NA
    ## 14             1382A sled TP4                                       NA
    ## 15             1382A sled TP5                                       NA
    ## 16             1382A sled TP6                                       NA
    ## 17             1382A sled TP7                                       NA
    ## 18             1382A sled TP8                                       NA
    ## 19          1382A filter 2014                                       NA
    ## 20    Bottomwater filter 2012                                       NA
    ## 21    Bottomwater filter 2014                                       NA
    ## 22           Jungbluth-JDF-16                               99.9730432
    ## 23  1383C Shallow filter 2012                  Acetothermia  0.0000000
    ## 24     1383C Shallow sled TP4                  Acetothermia  0.0000000
    ## 25     1383C Shallow sled TP5                  Acetothermia  0.0000000
    ## 26  1383C Shallow filter 2014                  Acetothermia  0.0000000
    ## 27      1383C Middle sled TP4                  Acetothermia  0.0000000
    ## 28   1383C Middle filter 2014                  Acetothermia  0.0000000
    ## 29     1383C Deep filter 2012                  Acetothermia  0.0000000
    ## 30        1383C Deep sled TP4                  Acetothermia  0.0000000
    ## 31     1383C Deep filter 2014                  Acetothermia  0.0000000
    ## 32          1382A filter 2012                  Acetothermia  0.0000000
    ## 33             1382A sled TP1                  Acetothermia  0.0000000
    ## 34             1382A sled TP2                  Acetothermia  0.0000000
    ## 35             1382A sled TP3                  Acetothermia  0.0000000
    ## 36             1382A sled TP4                  Acetothermia  0.0000000
    ## 37             1382A sled TP5                  Acetothermia  0.0000000
    ## 38             1382A sled TP6                  Acetothermia  0.0000000
    ## 39             1382A sled TP7                  Acetothermia  0.0000000
    ## 40             1382A sled TP8                  Acetothermia  0.0000000
    ## 41          1382A filter 2014                  Acetothermia  0.0000000
    ## 42    Bottomwater filter 2012                  Acetothermia  0.0000000
    ## 43    Bottomwater filter 2014                  Acetothermia  0.0000000
    ## 44           Jungbluth-JDF-16                  Acetothermia  5.6265980
    ## 45  1383C Shallow filter 2012                 Acidobacteria  0.0000000
    ## 46     1383C Shallow sled TP4                 Acidobacteria  0.0000000
    ## 47     1383C Shallow sled TP5                 Acidobacteria  0.0000000
    ## 48  1383C Shallow filter 2014                 Acidobacteria  0.0000000
    ## 49      1383C Middle sled TP4                 Acidobacteria  0.0000000
    ## 50   1383C Middle filter 2014                 Acidobacteria  0.0000000
    ## 51     1383C Deep filter 2012                 Acidobacteria  0.0000000
    ## 52        1383C Deep sled TP4                 Acidobacteria  0.0000000
    ## 53     1383C Deep filter 2014                 Acidobacteria  0.0000000
    ## 54          1382A filter 2012                 Acidobacteria  0.0000000
    ## 55             1382A sled TP1                 Acidobacteria  0.0000000
    ## 56             1382A sled TP2                 Acidobacteria  0.0000000
    ## 57             1382A sled TP3                 Acidobacteria  0.0000000
    ## 58             1382A sled TP4                 Acidobacteria  0.0000000
    ## 59             1382A sled TP5                 Acidobacteria  0.0000000
    ## 60             1382A sled TP6                 Acidobacteria  0.0000000
    ## 61             1382A sled TP7                 Acidobacteria  0.0000000
    ## 62             1382A sled TP8                 Acidobacteria  0.0000000
    ## 63          1382A filter 2014                 Acidobacteria  0.0000000
    ## 64    Bottomwater filter 2012                 Acidobacteria  0.5127000
    ## 65    Bottomwater filter 2014                 Acidobacteria  0.0000000
    ## 66           Jungbluth-JDF-16                 Acidobacteria  1.7902813
    ## 67  1383C Shallow filter 2012                Actinobacteria  0.0000000
    ## 68     1383C Shallow sled TP4                Actinobacteria  0.0000000
    ## 69     1383C Shallow sled TP5                Actinobacteria  0.0000000
    ## 70  1383C Shallow filter 2014                Actinobacteria  0.0000000
    ## 71      1383C Middle sled TP4                Actinobacteria  0.0000000
    ## 72   1383C Middle filter 2014                Actinobacteria  0.5729000
    ## 73     1383C Deep filter 2012                Actinobacteria  0.0000000
    ## 74        1383C Deep sled TP4                Actinobacteria  0.0000000
    ## 75     1383C Deep filter 2014                Actinobacteria  0.0000000
    ## 76          1382A filter 2012                Actinobacteria  0.0000000
    ## 77             1382A sled TP1                Actinobacteria  0.0000000
    ## 78             1382A sled TP2                Actinobacteria  0.0000000
    ## 79             1382A sled TP3                Actinobacteria  0.0000000
    ## 80             1382A sled TP4                Actinobacteria  0.0000000
    ## 81             1382A sled TP5                Actinobacteria  0.0000000
    ## 82             1382A sled TP6                Actinobacteria  0.0000000
    ## 83             1382A sled TP7                Actinobacteria  0.0000000
    ## 84             1382A sled TP8                Actinobacteria  0.0000000
    ## 85          1382A filter 2014                Actinobacteria  0.8806000
    ## 86    Bottomwater filter 2012                Actinobacteria  2.6354000
    ## 87    Bottomwater filter 2014                Actinobacteria  0.0000000
    ## 88           Jungbluth-JDF-16                Actinobacteria  1.7902813
    ## 89  1383C Shallow filter 2012                 Bacteroidetes 10.8120000
    ## 90     1383C Shallow sled TP4                 Bacteroidetes  0.0000000
    ## 91     1383C Shallow sled TP5                 Bacteroidetes  3.2012000
    ## 92  1383C Shallow filter 2014                 Bacteroidetes  7.1330000
    ## 93      1383C Middle sled TP4                 Bacteroidetes  5.7489000
    ## 94   1383C Middle filter 2014                 Bacteroidetes  4.3893000
    ## 95     1383C Deep filter 2012                 Bacteroidetes  2.2686000
    ## 96        1383C Deep sled TP4                 Bacteroidetes  5.9503000
    ## 97     1383C Deep filter 2014                 Bacteroidetes 12.8105000
    ## 98          1382A filter 2012                 Bacteroidetes  7.4757000
    ## 99             1382A sled TP1                 Bacteroidetes 47.9358000
    ## 100            1382A sled TP2                 Bacteroidetes  5.5940000
    ## 101            1382A sled TP3                 Bacteroidetes 28.8281000
    ## 102            1382A sled TP4                 Bacteroidetes 13.7159000
    ## 103            1382A sled TP5                 Bacteroidetes  3.5765000
    ## 104            1382A sled TP6                 Bacteroidetes 14.7621000
    ## 105            1382A sled TP7                 Bacteroidetes  3.5508000
    ## 106            1382A sled TP8                 Bacteroidetes  5.3575000
    ## 107         1382A filter 2014                 Bacteroidetes  1.5317000
    ## 108   Bottomwater filter 2012                 Bacteroidetes  4.9587000
    ## 109   Bottomwater filter 2014                 Bacteroidetes  1.4136000
    ## 110          Jungbluth-JDF-16                 Bacteroidetes  0.5115090
    ## 111 1383C Shallow filter 2012                   Chloroflexi  0.0000000
    ## 112    1383C Shallow sled TP4                   Chloroflexi  0.0000000
    ## 113    1383C Shallow sled TP5                   Chloroflexi  0.0000000
    ## 114 1383C Shallow filter 2014                   Chloroflexi  2.2374000
    ## 115     1383C Middle sled TP4                   Chloroflexi  0.0000000
    ## 116  1383C Middle filter 2014                   Chloroflexi  0.0000000
    ## 117    1383C Deep filter 2012                   Chloroflexi  0.0000000
    ## 118       1383C Deep sled TP4                   Chloroflexi  0.0000000
    ## 119    1383C Deep filter 2014                   Chloroflexi  0.6444000
    ## 120         1382A filter 2012                   Chloroflexi  0.0000000
    ## 121            1382A sled TP1                   Chloroflexi  0.0000000
    ## 122            1382A sled TP2                   Chloroflexi  0.0000000
    ## 123            1382A sled TP3                   Chloroflexi  0.0000000
    ## 124            1382A sled TP4                   Chloroflexi  0.0000000
    ## 125            1382A sled TP5                   Chloroflexi  0.0000000
    ## 126            1382A sled TP6                   Chloroflexi  0.0000000
    ## 127            1382A sled TP7                   Chloroflexi  0.0000000
    ## 128            1382A sled TP8                   Chloroflexi  0.0000000
    ## 129         1382A filter 2014                   Chloroflexi  0.0000000
    ## 130   Bottomwater filter 2012                   Chloroflexi  6.8316000
    ## 131   Bottomwater filter 2014                   Chloroflexi  4.0801000
    ## 132          Jungbluth-JDF-16                   Chloroflexi  4.3478261
    ## 133 1383C Shallow filter 2012                 Crenarchaeota  0.0000000
    ## 134    1383C Shallow sled TP4                 Crenarchaeota  0.0000000
    ## 135    1383C Shallow sled TP5                 Crenarchaeota  0.0000000
    ## 136 1383C Shallow filter 2014                 Crenarchaeota  0.0000000
    ## 137     1383C Middle sled TP4                 Crenarchaeota  0.0000000
    ## 138  1383C Middle filter 2014                 Crenarchaeota  0.0000000
    ## 139    1383C Deep filter 2012                 Crenarchaeota  0.0000000
    ## 140       1383C Deep sled TP4                 Crenarchaeota  0.0000000
    ## 141    1383C Deep filter 2014                 Crenarchaeota  0.0000000
    ## 142         1382A filter 2012                 Crenarchaeota  0.0000000
    ## 143            1382A sled TP1                 Crenarchaeota  0.0000000
    ## 144            1382A sled TP2                 Crenarchaeota  0.0000000
    ## 145            1382A sled TP3                 Crenarchaeota  0.0000000
    ## 146            1382A sled TP4                 Crenarchaeota  0.0000000
    ## 147            1382A sled TP5                 Crenarchaeota  0.0000000
    ## 148            1382A sled TP6                 Crenarchaeota  0.0000000
    ## 149            1382A sled TP7                 Crenarchaeota  0.0000000
    ## 150            1382A sled TP8                 Crenarchaeota  0.0000000
    ## 151         1382A filter 2014                 Crenarchaeota  0.0000000
    ## 152   Bottomwater filter 2012                 Crenarchaeota  0.0000000
    ## 153   Bottomwater filter 2014                 Crenarchaeota  0.0000000
    ## 154          Jungbluth-JDF-16                 Crenarchaeota  5.3708440
    ## 155 1383C Shallow filter 2012                 Cyanobacteria  0.0000000
    ## 156    1383C Shallow sled TP4                 Cyanobacteria  0.0000000
    ## 157    1383C Shallow sled TP5                 Cyanobacteria  0.0000000
    ## 158 1383C Shallow filter 2014                 Cyanobacteria  0.0000000
    ## 159     1383C Middle sled TP4                 Cyanobacteria  0.6698000
    ## 160  1383C Middle filter 2014                 Cyanobacteria  0.0000000
    ## 161    1383C Deep filter 2012                 Cyanobacteria  0.0000000
    ## 162       1383C Deep sled TP4                 Cyanobacteria  0.0000000
    ## 163    1383C Deep filter 2014                 Cyanobacteria  0.0000000
    ## 164         1382A filter 2012                 Cyanobacteria  0.0000000
    ## 165            1382A sled TP1                 Cyanobacteria  0.0000000
    ## 166            1382A sled TP2                 Cyanobacteria  0.0000000
    ## 167            1382A sled TP3                 Cyanobacteria  0.0000000
    ## 168            1382A sled TP4                 Cyanobacteria  0.0000000
    ## 169            1382A sled TP5                 Cyanobacteria  0.0000000
    ## 170            1382A sled TP6                 Cyanobacteria  0.0000000
    ## 171            1382A sled TP7                 Cyanobacteria  0.0000000
    ## 172            1382A sled TP8                 Cyanobacteria  0.0000000
    ## 173         1382A filter 2014                 Cyanobacteria  0.0000000
    ## 174   Bottomwater filter 2012                 Cyanobacteria  0.0000000
    ## 175   Bottomwater filter 2014                 Cyanobacteria  0.0000000
    ## 176          Jungbluth-JDF-16                 Cyanobacteria  0.5115090
    ## 177 1383C Shallow filter 2012            Epsilonbacteraeota  5.7690000
    ## 178    1383C Shallow sled TP4            Epsilonbacteraeota  3.2739000
    ## 179    1383C Shallow sled TP5            Epsilonbacteraeota 46.5134000
    ## 180 1383C Shallow filter 2014            Epsilonbacteraeota  3.2553000
    ## 181     1383C Middle sled TP4            Epsilonbacteraeota 29.4124000
    ## 182  1383C Middle filter 2014            Epsilonbacteraeota  1.9756000
    ## 183    1383C Deep filter 2012            Epsilonbacteraeota 77.8101000
    ## 184       1383C Deep sled TP4            Epsilonbacteraeota 43.1367000
    ## 185    1383C Deep filter 2014            Epsilonbacteraeota  7.4563000
    ## 186         1382A filter 2012            Epsilonbacteraeota 28.2108000
    ## 187            1382A sled TP1            Epsilonbacteraeota  0.0000000
    ## 188            1382A sled TP2            Epsilonbacteraeota 14.4519000
    ## 189            1382A sled TP3            Epsilonbacteraeota  1.2663000
    ## 190            1382A sled TP4            Epsilonbacteraeota  1.8340000
    ## 191            1382A sled TP5            Epsilonbacteraeota  4.5095000
    ## 192            1382A sled TP6            Epsilonbacteraeota  1.6282000
    ## 193            1382A sled TP7            Epsilonbacteraeota 15.7687000
    ## 194            1382A sled TP8            Epsilonbacteraeota  9.6672000
    ## 195         1382A filter 2014            Epsilonbacteraeota  0.6827000
    ## 196   Bottomwater filter 2012            Epsilonbacteraeota  0.0000000
    ## 197   Bottomwater filter 2014            Epsilonbacteraeota  0.0000000
    ## 198          Jungbluth-JDF-16            Epsilonbacteraeota  1.7902813
    ## 199 1383C Shallow filter 2012                 Euryarchaeota  0.0000000
    ## 200    1383C Shallow sled TP4                 Euryarchaeota  0.0000000
    ## 201    1383C Shallow sled TP5                 Euryarchaeota  0.0000000
    ## 202 1383C Shallow filter 2014                 Euryarchaeota  0.0000000
    ## 203     1383C Middle sled TP4                 Euryarchaeota  0.0000000
    ## 204  1383C Middle filter 2014                 Euryarchaeota  0.0000000
    ## 205    1383C Deep filter 2012                 Euryarchaeota  0.0000000
    ## 206       1383C Deep sled TP4                 Euryarchaeota  0.0000000
    ## 207    1383C Deep filter 2014                 Euryarchaeota  0.0000000
    ## 208         1382A filter 2012                 Euryarchaeota  0.0000000
    ## 209            1382A sled TP1                 Euryarchaeota  0.0000000
    ## 210            1382A sled TP2                 Euryarchaeota  0.0000000
    ## 211            1382A sled TP3                 Euryarchaeota  0.0000000
    ## 212            1382A sled TP4                 Euryarchaeota  0.0000000
    ## 213            1382A sled TP5                 Euryarchaeota  0.0000000
    ## 214            1382A sled TP6                 Euryarchaeota  0.0000000
    ## 215            1382A sled TP7                 Euryarchaeota  0.0000000
    ## 216            1382A sled TP8                 Euryarchaeota  0.0000000
    ## 217         1382A filter 2014                 Euryarchaeota  0.0000000
    ## 218   Bottomwater filter 2012                 Euryarchaeota  0.0000000
    ## 219   Bottomwater filter 2014                 Euryarchaeota  1.3688000
    ## 220          Jungbluth-JDF-16                 Euryarchaeota 12.7877238
    ## 221 1383C Shallow filter 2012                    Firmicutes  0.0000000
    ## 222    1383C Shallow sled TP4                    Firmicutes  0.0000000
    ## 223    1383C Shallow sled TP5                    Firmicutes  0.0000000
    ## 224 1383C Shallow filter 2014                    Firmicutes 17.4019000
    ## 225     1383C Middle sled TP4                    Firmicutes  1.4555000
    ## 226  1383C Middle filter 2014                    Firmicutes  0.0000000
    ## 227    1383C Deep filter 2012                    Firmicutes  0.0000000
    ## 228       1383C Deep sled TP4                    Firmicutes  0.0000000
    ## 229    1383C Deep filter 2014                    Firmicutes  0.9948000
    ## 230         1382A filter 2012                    Firmicutes  0.0000000
    ## 231            1382A sled TP1                    Firmicutes  0.0000000
    ## 232            1382A sled TP2                    Firmicutes  0.0000000
    ## 233            1382A sled TP3                    Firmicutes  0.0000000
    ## 234            1382A sled TP4                    Firmicutes  0.0000000
    ## 235            1382A sled TP5                    Firmicutes  0.0000000
    ## 236            1382A sled TP6                    Firmicutes  0.0000000
    ## 237            1382A sled TP7                    Firmicutes  0.9583000
    ## 238            1382A sled TP8                    Firmicutes  0.8435000
    ## 239         1382A filter 2014                    Firmicutes  0.0000000
    ## 240   Bottomwater filter 2012                    Firmicutes  0.0000000
    ## 241   Bottomwater filter 2014                    Firmicutes  0.0000000
    ## 242          Jungbluth-JDF-16                    Firmicutes  1.7902813
    ## 243 1383C Shallow filter 2012           Hydrothermarchaeota  0.0000000
    ## 244    1383C Shallow sled TP4           Hydrothermarchaeota  0.0000000
    ## 245    1383C Shallow sled TP5           Hydrothermarchaeota  0.0000000
    ## 246 1383C Shallow filter 2014           Hydrothermarchaeota  0.0000000
    ## 247     1383C Middle sled TP4           Hydrothermarchaeota  0.0000000
    ## 248  1383C Middle filter 2014           Hydrothermarchaeota  0.0000000
    ## 249    1383C Deep filter 2012           Hydrothermarchaeota  0.0000000
    ## 250       1383C Deep sled TP4           Hydrothermarchaeota  0.0000000
    ## 251    1383C Deep filter 2014           Hydrothermarchaeota  0.0000000
    ## 252         1382A filter 2012           Hydrothermarchaeota  0.0000000
    ## 253            1382A sled TP1           Hydrothermarchaeota  0.0000000
    ## 254            1382A sled TP2           Hydrothermarchaeota  0.0000000
    ## 255            1382A sled TP3           Hydrothermarchaeota  0.0000000
    ## 256            1382A sled TP4           Hydrothermarchaeota  0.0000000
    ## 257            1382A sled TP5           Hydrothermarchaeota  0.0000000
    ## 258            1382A sled TP6           Hydrothermarchaeota  0.0000000
    ## 259            1382A sled TP7           Hydrothermarchaeota  0.0000000
    ## 260            1382A sled TP8           Hydrothermarchaeota  0.0000000
    ## 261         1382A filter 2014           Hydrothermarchaeota  0.0000000
    ## 262   Bottomwater filter 2012           Hydrothermarchaeota  0.0000000
    ## 263   Bottomwater filter 2014           Hydrothermarchaeota  0.0000000
    ## 264          Jungbluth-JDF-16           Hydrothermarchaeota  2.5575448
    ## 265 1383C Shallow filter 2012 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 266    1383C Shallow sled TP4 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 267    1383C Shallow sled TP5 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 268 1383C Shallow filter 2014 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 269     1383C Middle sled TP4 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 270  1383C Middle filter 2014 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 271    1383C Deep filter 2012 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 272       1383C Deep sled TP4 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 273    1383C Deep filter 2014 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 274         1382A filter 2012 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 275            1382A sled TP1 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 276            1382A sled TP2 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 277            1382A sled TP3 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 278            1382A sled TP4 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 279            1382A sled TP5 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 280            1382A sled TP6 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 281            1382A sled TP7 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 282            1382A sled TP8 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 283         1382A filter 2014 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 284   Bottomwater filter 2012 Marinimicrobia_(SAR406_clade)  2.7394000
    ## 285   Bottomwater filter 2014 Marinimicrobia_(SAR406_clade)  2.7529000
    ## 286          Jungbluth-JDF-16 Marinimicrobia_(SAR406_clade)  0.0000000
    ## 287 1383C Shallow filter 2012                   Nitrospirae  0.0000000
    ## 288    1383C Shallow sled TP4                   Nitrospirae  0.0000000
    ## 289    1383C Shallow sled TP5                   Nitrospirae  0.0000000
    ## 290 1383C Shallow filter 2014                   Nitrospirae  0.0000000
    ## 291     1383C Middle sled TP4                   Nitrospirae  0.0000000
    ## 292  1383C Middle filter 2014                   Nitrospirae  0.0000000
    ## 293    1383C Deep filter 2012                   Nitrospirae  0.0000000
    ## 294       1383C Deep sled TP4                   Nitrospirae  0.0000000
    ## 295    1383C Deep filter 2014                   Nitrospirae  0.0000000
    ## 296         1382A filter 2012                   Nitrospirae  0.0000000
    ## 297            1382A sled TP1                   Nitrospirae  0.0000000
    ## 298            1382A sled TP2                   Nitrospirae  0.0000000
    ## 299            1382A sled TP3                   Nitrospirae  0.0000000
    ## 300            1382A sled TP4                   Nitrospirae  0.0000000
    ## 301            1382A sled TP5                   Nitrospirae  0.0000000
    ## 302            1382A sled TP6                   Nitrospirae  0.0000000
    ## 303            1382A sled TP7                   Nitrospirae  0.0000000
    ## 304            1382A sled TP8                   Nitrospirae  0.0000000
    ## 305         1382A filter 2014                   Nitrospirae  0.0000000
    ## 306   Bottomwater filter 2012                   Nitrospirae  0.0000000
    ## 307   Bottomwater filter 2014                   Nitrospirae  0.0000000
    ## 308          Jungbluth-JDF-16                   Nitrospirae  3.3248082
    ## 309 1383C Shallow filter 2012                Planctomycetes  1.3400000
    ## 310    1383C Shallow sled TP4                Planctomycetes  0.0000000
    ## 311    1383C Shallow sled TP5                Planctomycetes  1.0928000
    ## 312 1383C Shallow filter 2014                Planctomycetes  4.7227000
    ## 313     1383C Middle sled TP4                Planctomycetes  2.8063000
    ## 314  1383C Middle filter 2014                Planctomycetes  3.7823000
    ## 315    1383C Deep filter 2012                Planctomycetes  0.0000000
    ## 316       1383C Deep sled TP4                Planctomycetes  3.0823000
    ## 317    1383C Deep filter 2014                Planctomycetes  4.9188000
    ## 318         1382A filter 2012                Planctomycetes  0.0000000
    ## 319            1382A sled TP1                Planctomycetes  0.0000000
    ## 320            1382A sled TP2                Planctomycetes  2.6197000
    ## 321            1382A sled TP3                Planctomycetes  0.0000000
    ## 322            1382A sled TP4                Planctomycetes  0.0000000
    ## 323            1382A sled TP5                Planctomycetes  0.0000000
    ## 324            1382A sled TP6                Planctomycetes  4.7833000
    ## 325            1382A sled TP7                Planctomycetes  2.6347000
    ## 326            1382A sled TP8                Planctomycetes  2.6276000
    ## 327         1382A filter 2014                Planctomycetes  0.0000000
    ## 328   Bottomwater filter 2012                Planctomycetes  3.9515000
    ## 329   Bottomwater filter 2014                Planctomycetes  4.9223000
    ## 330          Jungbluth-JDF-16                Planctomycetes  0.0000000
    ## 331 1383C Shallow filter 2012          Proteobacteria-Alpha 15.0606000
    ## 332    1383C Shallow sled TP4          Proteobacteria-Alpha 49.8258000
    ## 333    1383C Shallow sled TP5          Proteobacteria-Alpha  8.3054000
    ## 334 1383C Shallow filter 2014          Proteobacteria-Alpha 16.6436000
    ## 335     1383C Middle sled TP4          Proteobacteria-Alpha 11.6442000
    ## 336  1383C Middle filter 2014          Proteobacteria-Alpha 11.6862000
    ## 337    1383C Deep filter 2012          Proteobacteria-Alpha  2.5675000
    ## 338       1383C Deep sled TP4          Proteobacteria-Alpha 12.1125000
    ## 339    1383C Deep filter 2014          Proteobacteria-Alpha 15.7208000
    ## 340         1382A filter 2012          Proteobacteria-Alpha 10.6114000
    ## 341            1382A sled TP1          Proteobacteria-Alpha 23.1797000
    ## 342            1382A sled TP2          Proteobacteria-Alpha 26.0970000
    ## 343            1382A sled TP3          Proteobacteria-Alpha 20.2649000
    ## 344            1382A sled TP4          Proteobacteria-Alpha 15.5962000
    ## 345            1382A sled TP5          Proteobacteria-Alpha 45.1903000
    ## 346            1382A sled TP6          Proteobacteria-Alpha 43.9842000
    ## 347            1382A sled TP7          Proteobacteria-Alpha 22.6772000
    ## 348            1382A sled TP8          Proteobacteria-Alpha 32.4707000
    ## 349         1382A filter 2014          Proteobacteria-Alpha 18.5936000
    ## 350   Bottomwater filter 2012          Proteobacteria-Alpha 45.2707000
    ## 351   Bottomwater filter 2014          Proteobacteria-Alpha 12.4220000
    ## 352          Jungbluth-JDF-16          Proteobacteria-Alpha  9.7186701
    ## 353 1383C Shallow filter 2012           Proteobacteria-Beta  7.2188000
    ## 354    1383C Shallow sled TP4           Proteobacteria-Beta  0.0000000
    ## 355    1383C Shallow sled TP5           Proteobacteria-Beta  0.0000000
    ## 356 1383C Shallow filter 2014           Proteobacteria-Beta  0.0000000
    ## 357     1383C Middle sled TP4           Proteobacteria-Beta  8.5268000
    ## 358  1383C Middle filter 2014           Proteobacteria-Beta  0.0000000
    ## 359    1383C Deep filter 2012           Proteobacteria-Beta  0.0000000
    ## 360       1383C Deep sled TP4           Proteobacteria-Beta  0.5420000
    ## 361    1383C Deep filter 2014           Proteobacteria-Beta  2.6183000
    ## 362         1382A filter 2012           Proteobacteria-Beta  0.0000000
    ## 363            1382A sled TP1           Proteobacteria-Beta  0.0000000
    ## 364            1382A sled TP2           Proteobacteria-Beta  1.6656000
    ## 365            1382A sled TP3           Proteobacteria-Beta  0.0000000
    ## 366            1382A sled TP4           Proteobacteria-Beta  0.0000000
    ## 367            1382A sled TP5           Proteobacteria-Beta  0.0000000
    ## 368            1382A sled TP6           Proteobacteria-Beta  0.6082000
    ## 369            1382A sled TP7           Proteobacteria-Beta  1.1838000
    ## 370            1382A sled TP8           Proteobacteria-Beta  1.3046000
    ## 371         1382A filter 2014           Proteobacteria-Beta  0.0000000
    ## 372   Bottomwater filter 2012           Proteobacteria-Beta  0.0000000
    ## 373   Bottomwater filter 2014           Proteobacteria-Beta  0.0000000
    ## 374          Jungbluth-JDF-16           Proteobacteria-Beta  0.0000000
    ## 375 1383C Shallow filter 2012          Proteobacteria-Delta  0.6916000
    ## 376    1383C Shallow sled TP4          Proteobacteria-Delta  0.0000000
    ## 377    1383C Shallow sled TP5          Proteobacteria-Delta 12.0872000
    ## 378 1383C Shallow filter 2014          Proteobacteria-Delta  9.2639000
    ## 379     1383C Middle sled TP4          Proteobacteria-Delta  0.0000000
    ## 380  1383C Middle filter 2014          Proteobacteria-Delta  0.0000000
    ## 381    1383C Deep filter 2012          Proteobacteria-Delta  0.5790000
    ## 382       1383C Deep sled TP4          Proteobacteria-Delta  1.1664000
    ## 383    1383C Deep filter 2014          Proteobacteria-Delta  2.7019000
    ## 384         1382A filter 2012          Proteobacteria-Delta  3.3086000
    ## 385            1382A sled TP1          Proteobacteria-Delta  0.0000000
    ## 386            1382A sled TP2          Proteobacteria-Delta  2.5930000
    ## 387            1382A sled TP3          Proteobacteria-Delta  0.0000000
    ## 388            1382A sled TP4          Proteobacteria-Delta 49.2591000
    ## 389            1382A sled TP5          Proteobacteria-Delta  0.0000000
    ## 390            1382A sled TP6          Proteobacteria-Delta  0.0000000
    ## 391            1382A sled TP7          Proteobacteria-Delta  5.1050000
    ## 392            1382A sled TP8          Proteobacteria-Delta  2.5279000
    ## 393         1382A filter 2014          Proteobacteria-Delta  0.0000000
    ## 394   Bottomwater filter 2012          Proteobacteria-Delta  7.3709000
    ## 395   Bottomwater filter 2014          Proteobacteria-Delta  5.6824000
    ## 396          Jungbluth-JDF-16          Proteobacteria-Delta  3.3248082
    ## 397 1383C Shallow filter 2012          Proteobacteria-Gamma 55.8563000
    ## 398    1383C Shallow sled TP4          Proteobacteria-Gamma 38.6959000
    ## 399    1383C Shallow sled TP5          Proteobacteria-Gamma 23.4172000
    ## 400 1383C Shallow filter 2014          Proteobacteria-Gamma 25.4439000
    ## 401     1383C Middle sled TP4          Proteobacteria-Gamma 35.6323000
    ## 402  1383C Middle filter 2014          Proteobacteria-Gamma 74.5367000
    ## 403    1383C Deep filter 2012          Proteobacteria-Gamma  8.0324000
    ## 404       1383C Deep sled TP4          Proteobacteria-Gamma 22.0113000
    ## 405    1383C Deep filter 2014          Proteobacteria-Gamma 20.2409000
    ## 406         1382A filter 2012          Proteobacteria-Gamma 45.3835000
    ## 407            1382A sled TP1          Proteobacteria-Gamma 21.2686000
    ## 408            1382A sled TP2          Proteobacteria-Gamma 35.4354000
    ## 409            1382A sled TP3          Proteobacteria-Gamma 47.7263000
    ## 410            1382A sled TP4          Proteobacteria-Gamma 17.1971000
    ## 411            1382A sled TP5          Proteobacteria-Gamma 44.5196000
    ## 412            1382A sled TP6          Proteobacteria-Gamma 23.0297000
    ## 413            1382A sled TP7          Proteobacteria-Gamma 29.1482000
    ## 414            1382A sled TP8          Proteobacteria-Gamma 34.2149000
    ## 415         1382A filter 2014          Proteobacteria-Gamma 71.7699000
    ## 416   Bottomwater filter 2012          Proteobacteria-Gamma 15.0172000
    ## 417   Bottomwater filter 2014          Proteobacteria-Gamma 55.6993000
    ## 418          Jungbluth-JDF-16          Proteobacteria-Gamma 15.0895141
    ## 419 1383C Shallow filter 2012           Proteobacteria-Zeta  0.0000000
    ## 420    1383C Shallow sled TP4           Proteobacteria-Zeta  0.7056000
    ## 421    1383C Shallow sled TP5           Proteobacteria-Zeta  1.0333000
    ## 422 1383C Shallow filter 2014           Proteobacteria-Zeta  1.4870000
    ## 423     1383C Middle sled TP4           Proteobacteria-Zeta  0.0000000
    ## 424  1383C Middle filter 2014           Proteobacteria-Zeta  0.0000000
    ## 425    1383C Deep filter 2012           Proteobacteria-Zeta  0.7876000
    ## 426       1383C Deep sled TP4           Proteobacteria-Zeta  4.7123000
    ## 427    1383C Deep filter 2014           Proteobacteria-Zeta 13.7474000
    ## 428         1382A filter 2012           Proteobacteria-Zeta  1.0075000
    ## 429            1382A sled TP1           Proteobacteria-Zeta  0.0000000
    ## 430            1382A sled TP2           Proteobacteria-Zeta  1.4908000
    ## 431            1382A sled TP3           Proteobacteria-Zeta  0.0000000
    ## 432            1382A sled TP4           Proteobacteria-Zeta  0.0000000
    ## 433            1382A sled TP5           Proteobacteria-Zeta  0.0000000
    ## 434            1382A sled TP6           Proteobacteria-Zeta  0.0000000
    ## 435            1382A sled TP7           Proteobacteria-Zeta  5.1613000
    ## 436            1382A sled TP8           Proteobacteria-Zeta  4.9691000
    ## 437         1382A filter 2014           Proteobacteria-Zeta  0.5070000
    ## 438   Bottomwater filter 2012           Proteobacteria-Zeta  0.0000000
    ## 439   Bottomwater filter 2014           Proteobacteria-Zeta  0.0000000
    ## 440          Jungbluth-JDF-16           Proteobacteria-Zeta  3.3248082
    ## 441 1383C Shallow filter 2012                   Tenericutes  0.0000000
    ## 442    1383C Shallow sled TP4                   Tenericutes  0.0000000
    ## 443    1383C Shallow sled TP5                   Tenericutes  0.0000000
    ## 444 1383C Shallow filter 2014                   Tenericutes  0.0000000
    ## 445     1383C Middle sled TP4                   Tenericutes  0.0000000
    ## 446  1383C Middle filter 2014                   Tenericutes  0.0000000
    ## 447    1383C Deep filter 2012                   Tenericutes  0.0000000
    ## 448       1383C Deep sled TP4                   Tenericutes  0.0000000
    ## 449    1383C Deep filter 2014                   Tenericutes  0.0000000
    ## 450         1382A filter 2012                   Tenericutes  0.0000000
    ## 451            1382A sled TP1                   Tenericutes  0.0000000
    ## 452            1382A sled TP2                   Tenericutes  0.0000000
    ## 453            1382A sled TP3                   Tenericutes  0.0000000
    ## 454            1382A sled TP4                   Tenericutes  0.0000000
    ## 455            1382A sled TP5                   Tenericutes  0.0000000
    ## 456            1382A sled TP6                   Tenericutes  0.0000000
    ## 457            1382A sled TP7                   Tenericutes  0.0000000
    ## 458            1382A sled TP8                   Tenericutes  0.0000000
    ## 459         1382A filter 2014                   Tenericutes  0.0000000
    ## 460   Bottomwater filter 2012                   Tenericutes  0.0000000
    ## 461   Bottomwater filter 2014                   Tenericutes  0.0000000
    ## 462          Jungbluth-JDF-16                   Tenericutes  0.2557545
    ## 463 1383C Shallow filter 2012                Thaumarchaeota  0.0000000
    ## 464    1383C Shallow sled TP4                Thaumarchaeota  0.0000000
    ## 465    1383C Shallow sled TP5                Thaumarchaeota  0.9765000
    ## 466 1383C Shallow filter 2014                Thaumarchaeota  4.5936000
    ## 467     1383C Middle sled TP4                Thaumarchaeota  0.0000000
    ## 468  1383C Middle filter 2014                Thaumarchaeota  1.4398000
    ## 469    1383C Deep filter 2012                Thaumarchaeota  0.0000000
    ## 470       1383C Deep sled TP4                Thaumarchaeota  0.5906000
    ## 471    1383C Deep filter 2014                Thaumarchaeota  8.5079000
    ## 472         1382A filter 2012                Thaumarchaeota  0.0000000
    ## 473            1382A sled TP1                Thaumarchaeota  0.0000000
    ## 474            1382A sled TP2                Thaumarchaeota  0.0000000
    ## 475            1382A sled TP3                Thaumarchaeota  0.0000000
    ## 476            1382A sled TP4                Thaumarchaeota  0.0000000
    ## 477            1382A sled TP5                Thaumarchaeota  0.0000000
    ## 478            1382A sled TP6                Thaumarchaeota  0.0000000
    ## 479            1382A sled TP7                Thaumarchaeota  1.3407000
    ## 480            1382A sled TP8                Thaumarchaeota  1.3537000
    ## 481         1382A filter 2014                Thaumarchaeota  0.5793000
    ## 482   Bottomwater filter 2012                Thaumarchaeota  4.7726000
    ## 483   Bottomwater filter 2014                Thaumarchaeota  5.2083000
    ## 484          Jungbluth-JDF-16                Thaumarchaeota  0.0000000
    ## 485 1383C Shallow filter 2012          unknown_unclassified  0.0000000
    ## 486    1383C Shallow sled TP4          unknown_unclassified  0.0000000
    ## 487    1383C Shallow sled TP5          unknown_unclassified  0.0000000
    ## 488 1383C Shallow filter 2014          unknown_unclassified  0.0000000
    ## 489     1383C Middle sled TP4          unknown_unclassified  0.0000000
    ## 490  1383C Middle filter 2014          unknown_unclassified  0.0000000
    ## 491    1383C Deep filter 2012          unknown_unclassified  0.0000000
    ## 492       1383C Deep sled TP4          unknown_unclassified  0.0000000
    ## 493    1383C Deep filter 2014          unknown_unclassified  0.0000000
    ## 494         1382A filter 2012          unknown_unclassified  0.0000000
    ## 495            1382A sled TP1          unknown_unclassified  0.0000000
    ## 496            1382A sled TP2          unknown_unclassified  0.0000000
    ## 497            1382A sled TP3          unknown_unclassified  0.0000000
    ## 498            1382A sled TP4          unknown_unclassified  0.0000000
    ## 499            1382A sled TP5          unknown_unclassified  0.0000000
    ## 500            1382A sled TP6          unknown_unclassified  0.0000000
    ## 501            1382A sled TP7          unknown_unclassified  0.0000000
    ## 502            1382A sled TP8          unknown_unclassified  0.0000000
    ## 503         1382A filter 2014          unknown_unclassified  0.0000000
    ## 504   Bottomwater filter 2012          unknown_unclassified  0.0000000
    ## 505   Bottomwater filter 2014          unknown_unclassified  0.0000000
    ## 506          Jungbluth-JDF-16          unknown_unclassified 26.0600000

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
