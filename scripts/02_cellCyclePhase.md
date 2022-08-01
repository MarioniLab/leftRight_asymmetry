---
title: "<span style='font-size: 28px'>Left-right asymmetric gene expression in the mouse developing heart</style>"
date: '01 August, 2022'
output:
  html_document:
    keep_md: true
    fig_width: 5
    fig_height: 5
    fig_caption: yes
    code_folding: hide
    toc: true
    toc_depth: 4
    toc_float: 
      collapsed: false
---



### Cell cycle stage analysis

During the characterisation of the mouse embryonic heart atlas data (Tyser et al., Science, 2021) we used `cyclone` to infer the cell cycle stage of each single cell. Let's explore this in the context of the left and right sides.


```r
## normalised counts
sce <- readRDS(paste0(dir, "data/sce_goodQual.NORM.clusters.Rds"))
row.names(sce) <- uniquifyFeatureNames(ID=row.names(rowData(sce)), names = rowData(sce)$gene)

## inferred phase from previous work
phase <- read.table(paste0(dir, "data/cellCyclePhase.tsv"))
# add to coldata
sce$phase <- phase[match(colnames(sce), phase$V1),2]

## retain only cells with side information
sce <- sce[,-which(sce$stage %in% c(-1,"LHT"))]

umap <- reducedDim(sce)
umap$phase <- as.factor(sce$phase)
o <- sample(1:nrow(umap), replace = FALSE)
plot(umap$x[o], umap$y[o], 
     col=ifelse(umap$phase[o] == "G1", "grey60",
                ifelse(umap$phase[o] == "S", "burlywood2","burlywood4")),
     pch=16, cex=0.75, axes=FALSE, xlab="", ylab=""); box(bty="l")
legend("bottomright", legend = c("G1", "S", "G2/M"),
       pch=16, col=c("grey60","burlywood2","burlywood4"))
```

![](02_cellCyclePhase_files/figure-html/data-1.png)<!-- -->

When looking at the proportion of cells in each phase, for each cluster, we observe a reduction of cells in **G1** phase for the cells from the right side in the cardiac mesoderm progenitor clusters (Me4-8). There are no differences in the cardiomyocytes from Me3.


```r
df <- data.frame(cluster=sce$clusterAnn, side=sce$side, phase=sce$phase)

left <- prop.table(table(df[df$side=="Left",]$cluster, df[df$side=="Left",]$phase),1)*100
right <- prop.table(table(df[df$side=="Right",]$cluster, df[df$side=="Right",]$phase),1)*100

df <- as.data.frame(rbind(left, right))
df$side <- rep(c("left", "right"), each = 12)
df$cluster <- substr(row.names(df), 1, 3)
df <- melt(df)

df$variable <- factor(df$variable, levels = rev(c("G1", "S", "G2M")))

ggplot(df, aes(side, value, fill=variable)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = rev(c("grey60","burlywood2","burlywood4"))) +
  facet_wrap(~cluster) +
  ylab("% cells") +
  labs(fill="phase") +
  th
```

![](02_cellCyclePhase_files/figure-html/proportions-1.png)<!-- -->

We confirm these observations by testing directly whether the proportion of cells that are not in G1 phase is significantly greater in the right side. This is indeed the case for clusters `Me4-7`.


```r
## collect numbers
df <- data.frame(cluster=sce$clusterAnn, side=sce$side, phase=sce$phase)

left <- table(df[df$side=="Left",]$cluster, df[df$side=="Left",]$phase)
left <- data.frame(G1 = left[,1], non_G1 = rowSums(left[,2:3]))
right <- table(df[df$side=="Right",]$cluster, df[df$side=="Right",]$phase)
right <- data.frame(G1 = right[,1], non_G1 = rowSums(right[,2:3]))

## take the left as expected proportions
props <- prop.table(as.matrix(left), 1)

## test if right is equal
pvals <- sapply(7:12, function(i) binom.test(x = right[i,'non_G1'], n = sum(right[i,]),
                                             p = props[i,'non_G1'], alternative = "greater")$p.value)
pvals <- data.frame(cluster = row.names(props)[7:12], p_value = pvals)
pvals
```

```
##   cluster     p_value
## 1     Me3 0.917393828
## 2     Me4 0.025108714
## 3     Me5 0.014072063
## 4     Me6 0.056911448
## 5     Me7 0.008996448
## 6     Me8 0.325828802
```

Despite the overall number of G1 cells in progenitor clusters being quite low, observing the effect in all clusters increases our confidence in the effect. Furthermore, since the left and right cells were collected from the same embryos, sampling biases are less a concern as if the groups came from independent samples.

---

Full numbers.

Left:


```r
cbind(
  table(sce[,sce$side == "Left"]$clusterAnn, sce[,sce$side == "Left"]$phase),
  prop.table(table(sce[,sce$side == "Left"]$clusterAnn, sce[,sce$side == "Left"]$phase),1)*100
)
```

```
##      G1 G2M   S        G1      G2M        S
## Ec1  10  46  47  9.708738 44.66019 45.63107
## Ec2   4  30  27  6.557377 49.18033 44.26230
## En1  19  87  64 11.176471 51.17647 37.64706
## En2  49  28  29 46.226415 26.41509 27.35849
## Me1   0   3   6  0.000000 33.33333 66.66667
## Me2   3  17   8 10.714286 60.71429 28.57143
## Me3 145 117  66 44.207317 35.67073 20.12195
## Me4  16  73  27 13.793103 62.93103 23.27586
## Me5  15 101  46  9.259259 62.34568 28.39506
## Me6   9  19  11 23.076923 48.71795 28.20513
## Me7  29 127 100 11.328125 49.60938 39.06250
## Me8   4  37  17  6.896552 63.79310 29.31034
```

Right:


```r
cbind(
  table(sce[,sce$side == "Right"]$clusterAnn, sce[,sce$side == "Right"]$phase),
  prop.table(table(sce[,sce$side == "Right"]$clusterAnn, sce[,sce$side == "Right"]$phase),1)*100
)
```

```
##      G1 G2M  S        G1      G2M        S
## Ec1   8  35 33 10.526316 46.05263 43.42105
## Ec2   0  18 10  0.000000 64.28571 35.71429
## En1  15  89 61  9.090909 53.93939 36.96970
## En2  65  25 32 53.278689 20.49180 26.22951
## Me1   1  21 13  2.857143 60.00000 37.14286
## Me2   3  13  7 13.043478 56.52174 30.43478
## Me3 151 118 46 47.936508 37.46032 14.60317
## Me4   6  51 34  6.593407 56.04396 37.36264
## Me5   8 118 51  4.519774 66.66667 28.81356
## Me6   1  15  2  5.555556 83.33333 11.11111
## Me7  11  97 78  5.913978 52.15054 41.93548
## Me8   1  19 13  3.030303 57.57576 39.39394
```




```r
sessionInfo()
```

```
## R version 4.1.0 (2021-05-18)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] ggpubr_0.4.0                RColorBrewer_1.1-3         
##  [3] reshape2_1.4.4              scater_1.20.1              
##  [5] ggplot2_3.3.5               scran_1.20.1               
##  [7] scuttle_1.2.1               SingleCellExperiment_1.14.1
##  [9] SummarizedExperiment_1.22.0 Biobase_2.52.0             
## [11] GenomicRanges_1.44.0        GenomeInfoDb_1.28.4        
## [13] IRanges_2.26.0              S4Vectors_0.30.2           
## [15] BiocGenerics_0.38.0         MatrixGenerics_1.4.3       
## [17] matrixStats_0.62.0         
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-7              backports_1.4.1          
##  [3] tools_4.1.0               bslib_0.3.1              
##  [5] utf8_1.2.2                R6_2.5.1                 
##  [7] irlba_2.3.5               vipor_0.4.5              
##  [9] DBI_1.1.2                 colorspace_2.0-3         
## [11] withr_2.5.0               gridExtra_2.3            
## [13] tidyselect_1.1.2          compiler_4.1.0           
## [15] cli_3.2.0                 BiocNeighbors_1.10.0     
## [17] DelayedArray_0.18.0       labeling_0.4.2           
## [19] sass_0.4.1                scales_1.2.0             
## [21] stringr_1.4.0             digest_0.6.29            
## [23] rmarkdown_2.13            XVector_0.32.0           
## [25] pkgconfig_2.0.3           htmltools_0.5.2          
## [27] sparseMatrixStats_1.4.2   highr_0.9                
## [29] fastmap_1.1.0             limma_3.48.3             
## [31] rlang_1.0.2               rstudioapi_0.13          
## [33] DelayedMatrixStats_1.14.3 farver_2.1.0             
## [35] jquerylib_0.1.4           generics_0.1.2           
## [37] jsonlite_1.8.0            BiocParallel_1.26.2      
## [39] car_3.0-12                dplyr_1.0.8              
## [41] RCurl_1.98-1.6            magrittr_2.0.3           
## [43] BiocSingular_1.8.1        GenomeInfoDbData_1.2.6   
## [45] Matrix_1.4-1              Rcpp_1.0.8.3             
## [47] ggbeeswarm_0.6.0          munsell_0.5.0            
## [49] fansi_1.0.3               abind_1.4-5              
## [51] viridis_0.6.2             lifecycle_1.0.1          
## [53] stringi_1.7.6             yaml_2.3.5               
## [55] edgeR_3.34.1              carData_3.0-5            
## [57] zlibbioc_1.38.0           plyr_1.8.7               
## [59] grid_4.1.0                dqrng_0.3.0              
## [61] crayon_1.5.1              lattice_0.20-45          
## [63] beachmat_2.8.1            locfit_1.5-9.5           
## [65] metapod_1.0.0             knitr_1.38               
## [67] pillar_1.7.0              igraph_1.3.0             
## [69] ggsignif_0.6.3            ScaledMatrix_1.0.0       
## [71] glue_1.6.2                evaluate_0.15            
## [73] vctrs_0.4.1               tidyr_1.2.0              
## [75] gtable_0.3.0              purrr_0.3.4              
## [77] assertthat_0.2.1          xfun_0.30                
## [79] rsvd_1.0.5                broom_0.8.0              
## [81] rstatix_0.7.0             viridisLite_0.4.0        
## [83] tibble_3.1.6              beeswarm_0.4.0           
## [85] cluster_2.1.3             bluster_1.2.1            
## [87] statmod_1.4.36            ellipsis_0.3.2
```

