---
title: "<span style='font-size: 28px'>Left-right asymmetric gene expression in the mouse developing heart</style>"
date: '30 July, 2022'
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



### Atlas of mouse heart development

In previous work, we collected mouse embryos at different stages of development, and dissected out the region comprising the cardiac crescent. Sampling spanned stages from just before the left and right portions of the prospective cardiac crescent fuse, up to the linear heart tube (LHT) stage. The cells from the isolated regions were used to construct single-cell RNA-seq libraries. 


```r
## normalised counts
sce <- readRDS(paste0(dir, "data/sce_goodQual.NORM.clusters.Rds"))

## UMAP coordinates
umap <- reducedDim(sce)
umap$cluster <- sce$clusterAnn
```

The analysis of these data, comprising 3104 cells, revealed multiple cell populations from all three germ layers. The majority were of mesodermal origin, and comprised two differentiation trajectories from cardiac progenitor cells (`Me5` and `Me7`) to mature cardiomyocytes (`Me3`).


```r
## define median cluster position for labels
text.x <- sapply(unique(umap$cluster), function(x) median(umap[umap$cluster==x,1]))
text.y <- sapply(unique(umap$cluster), function(x) median(umap[umap$cluster==x,2]))

## plot
o <- sample(1:nrow(umap), replace = FALSE)
plot(umap$x[o], umap$y[o], col=sce$clusterCol[o], pch=16, cex=0.85, axes=FALSE, xlab="", ylab=""); box(bty="l")
text(x = text.x, y = text.y, labels = names(text.x), cex=0.9)
```

![](01_left_vs_right_files/figure-html/umap-1.png)<!-- -->

### Left *versus* right gene expression {.tabset}

Upon dissection, for the embryos from stages 0 to 3, the left and right halves of the cardiac crescent were separated, and for each cell it was recorded whether it came from the left or right sides. Thus, with these data we can study asymmetric expression during early heart development.

Roughly equal numbers of cells from each side were collected at every stage, although stages 0 and 3 are biased 60:40 towards the left side.


```r
## retain only cells with side information
sce <- sce[,-which(sce$stage %in% c(-1,"LHT"))]
umap <- reducedDim(sce)

table(stage=sce$stage, side=sce$side)
```

```
##      side
## stage Left Right
##     0  303   179
##     1  422   407
##     2  463   487
##     3  248   196
```

On the UMAP, cells from the left and right sides intermingle well, suggesting they have very similar transcriptomes. Most clusters are composed of roughly equal numbers of cells from each side; the clusters that deviate from this behaviour are generally small clusters with few cells (the width of the bars are proportional to cluster size), and thus the fluctuations are likely the result from sampling effects.


```r
par(mfrow=c(1,2), mar=c(4,2,2,2))
plot(umap$x[o], umap$y[o], col=as.factor(sce$side)[o], pch=16, axes=FALSE, xlab="", ylab=""); box(bty="l")
legend("topright", legend = levels(as.factor(sce$side)), pch=16, col=1:2, cex=0.85)

plot(table(sce$clusterAnn, sce$side), col=1:2, main=""); abline(h=0.5, lty=2)
mtext(side=2, line=0.25, text=c(0,0.5,1), at=c(0,0.5,1), las=2, cex=0.75)
```

![](01_left_vs_right_files/figure-html/props-1.png)<!-- -->

We know that certain genes are expressed in an asymmetric fashion. Some of the best characterised, like *Lefty2* and *Pitx2*, show clear bias to be expressed preferentially on the cells from the left side.


```r
plots <- list()

tmp <- data.frame(side=sce$side, expr=logcounts(sce)[which(rowData(sce)$gene == "Lefty2"),])
plots[[1]] <- ggplot(tmp, aes(side, expr)) + geom_boxplot() + ylab(expression('log'[2]*' normalised counts')) + ggtitle("Lefty2") + th

tmp <- data.frame(side=sce$side, expr=logcounts(sce)[which(rowData(sce)$gene == "Pitx2"),])
plots[[2]] <- ggplot(tmp, aes(side, expr)) + geom_boxplot() + ylab(expression('log'[2]*' normalised counts')) + ggtitle("Pitx2") + th

# tmp <- data.frame(side=sce$side, expr=logcounts(sce)[which(rowData(sce)$gene == "Nodal"),])
# plots[[3]] <- ggplot(tmp, aes(side, expr)) + geom_boxplot() + ylab(expression('log'[2]*' normalised counts')) + ggtitle("Nodal") + th

ggarrange(plotlist = plots)
```

![](01_left_vs_right_files/figure-html/leftGenes-1.png)<!-- -->

We can use differential expression to identify other genes expressed unequally between the left and right side cells. We analyse clusters from the three germ layers separately.


```r
y <- convertTo(sce, "edgeR")

## save DE results in a list
resAdj <- list()
```

#### Mesoderm

For the mesodermal lineage analysis, we exclude the blood and endothelial clusters (`Me1` and `Me2`) because they have very few cells. For all others, we test for differences between the left and right cells, within each cluster. However, we test all comparisons together, so the p-value indicates whether the gene is asymmetrically expressed in *any* of the clusters.


```r
cluster <- "meso"
clusters <- paste0("Me",3:8)
cells.meso <- which(y$samples$clusterAnn %in% clusters)

y.tmp <- y
y.tmp$counts <- y.tmp$counts[,cells.meso]
y.tmp$samples <- y.tmp$samples[cells.meso,]
y.tmp$counts <- y.tmp$counts[rowMeans(y.tmp$counts)>1,] ## filter lowly expressed genes
y.tmp$genes <- y.tmp$genes[row.names(y.tmp$counts),]
stopifnot(identical(row.names(y.tmp$samples), colnames(y.tmp$counts)))

y.tmp$samples$batch <- droplevels(y.tmp$samples$batch)
## use a combination of cluster and side as groups
y.tmp$samples$group <- paste(y.tmp$samples$clusterAnn, y.tmp$samples$side, sep=".")

design <- model.matrix(~0+y.tmp$samples$group+y.tmp$samples$batch)
colnames(design) <- substr(colnames(design),20,40)

## edgeR workflow
y.tmp <- estimateDisp(y.tmp, design)
# plotBCV(y.tmp)
fit <- glmQLFit(y.tmp, design, robust=TRUE)
# plotQLDisp(fit)

my.contrasts <- makeContrasts(Me3 = Me3.Left - Me3.Right,
                              Me4 = Me4.Left - Me4.Right,
                              Me5 = Me5.Left - Me5.Right,
                              Me6 = Me6.Left - Me6.Right,
                              Me7 = Me7.Left - Me7.Right,
                              Me8 = Me8.Left - Me8.Right, levels=design)
test <- glmQLFTest(fit, contrast = my.contrasts)
resAdj[[cluster]] <- as.data.frame(topTags(test, n=nrow(y.tmp$counts)))
```

We consider differentially expressed genes those with an FDR < 0.05. To retain meaningful changes in expression, we also require a minimum fold-change of 1.5 in at least one cluster. And to remove artefacts from very few cells expressing the gene, clusters with fold-changes > 1.5 must have mean log2 expression of at least 1, to ensure the fold-change estimate is robust.


```r
## mesoderm DEGs
de.meso <- as.data.frame(resAdj[[cluster]])
de.meso$sig <- ifelse(de.meso$FDR<0.05, "DEG", "notSignificant")
de.meso <- de.meso[de.meso$sig=="DEG",]

## define clusters with fold-changes > 1.5
perCluster <- ifelse(abs(de.meso[,3:8]) > log2(1.5), 1, 0)

## one problem is that some clusters have expression in very few cells, and this results in unstable and non-informative fold-changes. So these should be ignored
## to filter out clusters with very low expression, we use mean expression per cluster as an estimate
means <- matrix <- matrix(nrow = nrow(de.meso), ncol=6)
colnames(means) <- paste0("Me",3:8)
row.names(means) <- row.names(de.meso)
for(c in paste0("Me",3:8)){
  cells <- which(sce$clusterAnn == c)
  means[,c] <- rowMeans(logcounts(sce)[row.names(de.meso),cells])
}

## remove 'large' fold-changes for cluster with very low expression
perCluster[means[row.names(perCluster),] < 1] <- NA

## filter genes that don't reach the minimum fold-change
perCluster <- as.data.frame(perCluster)
perCluster$count <- rowSums(perCluster, na.rm=TRUE)
perCluster <- perCluster[perCluster$count>0,]

## restrict DEGs
de.meso <- de.meso[row.names(perCluster),]

tmp <- as.data.frame(resAdj[[cluster]])
tmp$sig <- ifelse(row.names(tmp) %in% row.names(de.meso), "DEG", "notSignificant")
## to get a summary fold-change for plotting, take the mean of all clusters as an approximation
tmp$logFC <- rowMeans(tmp[,3:8])
## for DEGs, get the max FC from clusters that pass the expression filter
tmp[tmp$sig=="DEG",]$logFC <- unlist(sapply(1:nrow(tmp[tmp$sig=="DEG",]), function(x)
  tmp[tmp$sig=="DEG",][x,names(which.max(abs(tmp[tmp$sig=="DEG",][x,3:8][means[row.names(tmp[tmp$sig=="DEG",])[x],]>1])))]))
tmp <- tmp[order(tmp$sig, abs(tmp$logFC), decreasing=TRUE),]

ggplot(tmp, aes(logCPM, logFC, label=gene)) + 
  geom_point(aes(col=sig)) + 
  scale_color_manual(name = "", values = c("red", "black")) + 
  geom_hline(yintercept = c(-log2(1.5), log2(1.5)), lty=2) + 
  ggtitle("Cardiac mesoderm") + 
  xlab(expression('log'[2]*' mean expression')) + 
  ylab(expression('log'[2]*' fold-change (left / right)')) + 
  geom_text_repel(data=tmp[tmp$sig=="DEG",][1:20,]) + 
  th 
```

```
## Warning: ggrepel: 4 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
```

![](01_left_vs_right_files/figure-html/mesoderm_degs-1.png)<!-- -->

There are 229 DE genes, with many more (two thirds) expressed higher in the cells from the left side (positive fold-change). Among these are the expected genes (*Lefty1*, *Lefty2*, *Pitx2*).

Generally, when a gene is significantly different in more than one cluster, all clusters tend to be biased towards the same side. That is, for two thirds of the genes significant in two or more clusters, all clusters are biased to the left or the right. The remaining third of the genes have conflicting patterns, with one cluster to the left and one cluster to the right, and often involve the small clusters (particularly `Me6` and `Me8`).


```r
## now that we have filtered out the DEGs, we can use a less stringent filter to determine which fold-changes we consider. We want to only remove clusters with little to no expression.
## however, mean expression depends on the size of the cluster; same mean in a small and big clusters mean the big cluster probably expresses the gene in many more cells.
sf <- 100/as.numeric(table(y.tmp$samples$clusterAnn))
means.scaled <- t(t(means)/sf)

perCluster <- ifelse(abs(de.meso[,3:8]) > log2(1.5), 1, 0)
perCluster[means.scaled[row.names(perCluster),]<0.3] <- NA

## filter genes that don't reach the minimum fold-change
perCluster <- as.data.frame(perCluster)
perCluster$count <- rowSums(perCluster, na.rm=TRUE)
perCluster <- perCluster[perCluster$count>0,]
perCluster$gene <- de.meso[row.names(perCluster),]$gene

## remove fold-change from results df
de.meso[,3:8] <- t(sapply(1:nrow(de.meso), function(x) 
  ifelse(is.na(perCluster[x,1:6]), NA, de.meso[x,3:8][!is.na(perCluster[x,1:6])])))
# use largest fold-change as the summary value
de.meso$logFC <- apply(de.meso[,3:8], 1, function(x) x[which.max(abs(x))])

## define the direction of bias
bias.meso <- ifelse(de.meso[,3:8]>0, "left", "right")
bia.mesos <- bias.meso[match(row.names(perCluster), row.names(bias.meso)),]
bias.meso[is.na(perCluster[,1:6])] <- "notExpr" # gene not expressed in cluster
bias.meso[perCluster[,1:6]==0] <- "notSig" ## fold-change not large enough
bias.meso <- as.data.frame(bias.meso)
bias.meso$gene <- de.meso[row.names(bias.meso),]$gene

bias.meso$left <- apply(bias.meso[,1:6], 1, function(x) sum(x=="left", na.rm = TRUE))
bias.meso$right <- apply(bias.meso[,1:6], 1, function(x) sum(x=="right", na.rm = TRUE))
bias.meso$count <- bias.meso$left+bias.meso$right
# table(left=bias.meso$left, right=bias.meso$right)
#     right
# left  0  1  2  3  4
#    0  0 26 14  6  1
#    1 33 26  7  3  0
#    2 36 19  6  0  0
#    3 23  5  1  0  0
#    4 10  4  0  0  0
#    5  7  2  0  0  0

## flag the genes with bias to both sides as not-congruent
bias.meso$congruent <- ifelse(bias.meso$left>0 & bias.meso$right>0, 0, 1)

# i=1
# gene <- bias.meso[bias.meso$left>0 & bias.meso$right > 0,][i,7]
# plotGeneSide(cells=cells.meso, gene=gene, bias = bias.meso)
# i<-i+1
## sometimes, the incongruent genes show opposing pattern between first and second heart field clusters.
```

The 5 most significant genes are plotted below, plus *Pitx2* for reference.


```r
plots <- list()
for(i in c(1:5,7)){
  id <- row.names(de.meso)[i]
  gene <- de.meso[i,1]
  df <- data.frame(x=umap$x, y=umap$y,
                   side=sce$side, cluster=sce$clusterAnn, 
                   expr=logcounts(sce)[id,], gene = gene)
  df <- df[df$cluster %in% paste0("Me", 3:8),]
  df <- df[order(df$expr),]

  plots[[i]] <- ggplot(df, aes(x, y, colour=expr)) +
    geom_point(size=1, alpha=0.8) +
    scale_color_gradientn(colours = c("grey", brewer.pal(n=9, "Blues"))) +
    ylab(gene) + xlab("") +
    facet_wrap(~side) +
    th + theme(axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
}
ggarrange(plotlist = plots[c(7,1:5)], ncol=2, nrow=3)
```

![](01_left_vs_right_files/figure-html/plot_meso-1.png)<!-- -->

#### Endoderm

Now we do the same for the two endoderm clusters.


```r
cluster <- "endo"
clusters <- paste0("En",1:2)
cells.endo <- which(y$samples$clusterAnn %in% clusters)

y.tmp <- y
y.tmp$counts <- y.tmp$counts[,cells.endo]
y.tmp$samples <- y.tmp$samples[cells.endo,]
y.tmp$counts <- y.tmp$counts[rowMeans(y.tmp$counts)>1,] ## filter lowly expressed genes
y.tmp$genes <- y.tmp$genes[row.names(y.tmp$counts),]
stopifnot(identical(row.names(y.tmp$samples), colnames(y.tmp$counts)))

y.tmp$samples$batch <- droplevels(y.tmp$samples$batch)
## use a combination of cluster and side as groups
y.tmp$samples$group <- paste(y.tmp$samples$clusterAnn, y.tmp$samples$side, sep=".")

design <- model.matrix(~0+y.tmp$samples$group+y.tmp$samples$batch)
colnames(design) <- substr(colnames(design),20,40)

## edgeR workflow
y.tmp <- estimateDisp(y.tmp, design)
# plotBCV(y.tmp)
fit <- glmQLFit(y.tmp, design, robust=TRUE)
# plotQLDisp(fit)

my.contrasts <- makeContrasts(En1 = En1.Left - En1.Right,
                              En2 = En2.Left - En2.Right, levels=design)
test <- glmQLFTest(fit, contrast = my.contrasts)
resAdj[[cluster]] <- as.data.frame(topTags(test, n=nrow(y.tmp$counts)))
```

We relax the fold-change threshold to 1.3 because we have less power due to fewer cells.


```r
## endoderm DEGs
de.endo <- as.data.frame(resAdj[[cluster]])
de.endo$sig <- ifelse(de.endo$FDR<0.05, "DEG", "notSignificant")
de.endo <- de.endo[de.endo$sig=="DEG",]

## define mean expression per-cluster for DEGs
means <- matrix <- matrix(nrow=nrow(resAdj[[cluster]]), ncol=2)
colnames(means) <- paste0("En",1:2)
row.names(means) <- row.names(resAdj[[cluster]])
for(c in paste0("En",1:2)){
  cells <- which(sce$clusterAnn == c)
  means[,c] <- rowMeans(logcounts(sce)[row.names(resAdj[[cluster]]),cells])
}

## define clusters with large fold-changes
perCluster <- ifelse(abs(de.endo[,3:4]) > log2(1.3), 1, 0)

## remove 'large' fold-changes for cluster with very low expression
perCluster[means[row.names(perCluster),]<1] <- NA

## filter genes that don't reach the minimum fold-change
perCluster <- as.data.frame(perCluster)
perCluster$count <- rowSums(perCluster, na.rm=TRUE)
perCluster <- perCluster[perCluster$count>0,]
perCluster$gene <- de.endo[row.names(perCluster),]$gene

## restrict DEGs
de.endo <- de.endo[row.names(perCluster),]

tmp <- as.data.frame(resAdj[[cluster]])
## to get a summary fold-change for plotting, take the largest fold-change as an approximation
tmp$logFC <- rowMeans(tmp[,3:4])
tmp$sig <- ifelse(row.names(tmp) %in% row.names(de.endo), "DEG", "notSignificant")
## for DEGs, we only consider FCs when mean expression in that cluster is >1. Thus, get the max FC from only pertinent clusters
tmp[tmp$sig=="DEG",]$logFC <- unlist(sapply(1:nrow(tmp[tmp$sig=="DEG",]), function(x) tmp[tmp$sig=="DEG",][x,names(which.max(abs(tmp[tmp$sig=="DEG",][x,3:4][means[row.names(tmp[tmp$sig=="DEG",])[x],]>1])))]))
tmp <- tmp[order(tmp$sig, abs(tmp$logFC), decreasing=TRUE),]

ggplot(tmp, aes(logCPM, logFC, label=gene)) + geom_point(aes(col=sig)) + scale_color_manual(name = "", values = c("red", "black")) + geom_hline(yintercept = c(-log2(1.3), log2(1.3)), lty=2) + xlab(expression('log'[2]*' mean expression')) + ylab(expression('log'[2]*' fold-change (left / right)')) + geom_text_repel(data=tmp[tmp$sig=="DEG",][1:20,]) + th + ggtitle("Endoderm")
```

```
## Warning: ggrepel: 1 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
```

![](01_left_vs_right_files/figure-html/endoderm_degs-1.png)<!-- -->

There are 23 genes differentially expressed, with nearly three quarters expressed higher on the right-side cells.


```r
sf <- 100/as.numeric(table(y.tmp$samples$clusterAnn))
means.scaled <- t(t(means)/sf)

perCluster <- ifelse(abs(de.endo[,3:4]) > log2(1.3), 1, 0)
perCluster[means.scaled[row.names(perCluster),]<0.3] <- NA

## filter genes that don't reach the minimum fold-change
perCluster <- as.data.frame(perCluster)
perCluster$count <- rowSums(perCluster, na.rm=TRUE)
perCluster <- perCluster[perCluster$count>0,]
perCluster$gene <- de.endo[row.names(perCluster),]$gene

## remove fold-change from results df
de.endo[,3:4] <- t(sapply(1:nrow(de.endo), function(x) 
  ifelse(is.na(perCluster[x,1:2]), NA, de.endo[x,3:4][!is.na(perCluster[x,1:2])])))
# use largest fold-change as the summary value
de.endo$logFC <- apply(de.endo[,3:4], 1, function(x) x[which.max(abs(x))])

## define the direction of bias
bias.endo <- ifelse(de.endo[,3:4]>0, "left", "right")
bias.endo <- bias.endo[match(row.names(perCluster), row.names(bias.endo)),]
bias.endo[is.na(perCluster[,1:2])] <- "notExpr" # gene not expressed in cluster
bias.endo[perCluster[,1:2]==0] <- "notSig" ## fold-change not large enough
bias.endo <- as.data.frame(bias.endo)
bias.endo$gene <- de.endo[row.names(bias.endo),]$gene

bias.endo$left <- apply(bias.endo[,1:2], 1, function(x) sum(x=="left", na.rm = TRUE))
bias.endo$right <- apply(bias.endo[,1:2], 1, function(x) sum(x=="right", na.rm = TRUE))
bias.endo$count <- bias.endo$left+bias.endo$right
# table(left=bias.endo$left, right=bias.endo$right)
#     right
# left 0 1 2
#    0 0 4 2
#    1 4 1 0
#    2 1 0 0

## flag the genes with bias to both sides as not-congruent
bias.endo$congruent <- ifelse(bias.endo$left>0 & bias.endo$right>0, 0, 1)
```

Seven genes are significant in both clusters; five are right-biased, one is left-biased and one shows conflicting direction between the two clusters.

The 6 most significant genes are plotted below.


```r
plots <- list()
for(i in 1:6){
  id <- row.names(de.endo)[i]
  gene <- de.endo[i,1]
  df <- data.frame(x=umap$x, y=umap$y,
                   side=sce$side, cluster=sce$clusterAnn, 
                   expr=logcounts(sce)[id,], gene = gene)
  df <- df[df$cluster %in% paste0("En", 1:2),]
  df <- df[order(df$expr),]

  plots[[i]] <- ggplot(df, aes(x, y, colour=expr)) +
    geom_point(size=1, alpha=0.8) +
    scale_color_gradientn(colours = c("grey", brewer.pal(n=9, "Blues"))) +
    ylab(gene) + xlab("") +
    facet_wrap(~side) +
    th + theme(axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
}
ggarrange(plotlist = plots, ncol=2, nrow=3)
```

![](01_left_vs_right_files/figure-html/plot_endo-1.png)<!-- -->


#### Ectoderm

Finally, we repeat for the two ectoderm clusters, with the same 1.3 fold-change threshold.


```r
cluster <- "ecto"
clusters <- paste0("Ec",1:2)
cells.ecto <- which(y$samples$clusterAnn %in% clusters)

y.tmp <- y
y.tmp$counts <- y.tmp$counts[,cells.ecto]
y.tmp$samples <- y.tmp$samples[cells.ecto,]
y.tmp$counts <- y.tmp$counts[rowMeans(y.tmp$counts)>1,] ## filter lowly expressed genes
y.tmp$genes <- y.tmp$genes[row.names(y.tmp$counts),]
stopifnot(identical(row.names(y.tmp$samples), colnames(y.tmp$counts)))

y.tmp$samples$batch <- droplevels(y.tmp$samples$batch)
## use a combination of cluster and side as groups
y.tmp$samples$group <- paste(y.tmp$samples$clusterAnn, y.tmp$samples$side, sep=".")

design <- model.matrix(~0+y.tmp$samples$group+y.tmp$samples$batch)
colnames(design) <- substr(colnames(design),20,40)

## edgeR workflow
y.tmp <- estimateDisp(y.tmp, design)
# plotBCV(y.tmp)
fit <- glmQLFit(y.tmp, design, robust=TRUE)
# plotQLDisp(fit)

my.contrasts <- makeContrasts(Ec1 = Ec1.Left - Ec1.Right,
                              Ec2 = Ec2.Left - Ec2.Right, levels=design)
test <- glmQLFTest(fit, contrast = my.contrasts)
resAdj[[cluster]] <- as.data.frame(topTags(test, n=nrow(y.tmp$counts)))
```



```r
## ectoderm DEGs
de.ecto <- as.data.frame(resAdj[[cluster]])
de.ecto$sig <- ifelse(de.ecto$FDR<0.05, "DEG", "notSignificant")
de.ecto <- de.ecto[de.ecto$sig=="DEG",]

## define mean expression per-cluster for DEGs
means <- matrix <- matrix(nrow=nrow(resAdj[[cluster]]), ncol=2)
colnames(means) <- paste0("Ec",1:2)
row.names(means) <- row.names(resAdj[[cluster]])
for(c in paste0("Ec",1:2)){
  cells <- which(sce$clusterAnn == c)
  means[,c] <- rowMeans(logcounts(sce)[row.names(resAdj[[cluster]]),cells])
}

## define clusters with large fold-changes
perCluster <- ifelse(abs(de.ecto[,3:4]) > log2(1.3), 1, 0)

## remove 'large' fold-changes for cluster with very low expression
perCluster[means[row.names(perCluster),]<1] <- NA

## filter genes that don't reach the minimum fold-change
perCluster <- as.data.frame(perCluster)
perCluster$count <- rowSums(perCluster, na.rm=TRUE)
perCluster <- perCluster[perCluster$count>0,]
perCluster$gene <- de.ecto[row.names(perCluster),]$gene

## restrict DEGs
de.ecto <- de.ecto[row.names(perCluster),]

tmp <- as.data.frame(resAdj[[cluster]])
## to get a summary fold-change for plotting, take the largest fold-change as an approximation
tmp$logFC <- rowMeans(tmp[,3:4])
tmp$sig <- ifelse(row.names(tmp) %in% row.names(de.ecto), "DEG", "notSignificant")
## for DEGs, we only consider FCs when mean expression in that cluster is >1. Thus, get the max FC from only pertinent clusters
tmp[tmp$sig=="DEG",]$logFC <- unlist(sapply(1:nrow(tmp[tmp$sig=="DEG",]), function(x) tmp[tmp$sig=="DEG",][x,names(which.max(abs(tmp[tmp$sig=="DEG",][x,3:4][means[row.names(tmp[tmp$sig=="DEG",])[x],]>1])))]))
tmp <- tmp[order(tmp$sig, abs(tmp$logFC), decreasing=TRUE),]

ggplot(tmp, aes(logCPM, logFC, label=gene)) + geom_point(aes(col=sig)) + scale_color_manual(name = "", values = c("red", "black")) + geom_hline(yintercept = c(-log2(1.3), log2(1.3)), lty=2) + xlab(expression('log'[2]*' mean expression')) + ylab(expression('log'[2]*' fold-change (left / right)')) + geom_text_repel(data=tmp[tmp$sig=="DEG",][1:20,]) + th + ggtitle("Ectoderm")
```

![](01_left_vs_right_files/figure-html/ectoderm_degs-1.png)<!-- -->

There are 60 differentially expressed genes, with roughly equal numbers of genes expressed higher on the left or right sides.


```r
sf <- 100/as.numeric(table(y.tmp$samples$clusterAnn))
means.scaled <- t(t(means)/sf)

perCluster <- ifelse(abs(de.ecto[,3:4]) > log2(1.3), 1, 0)
perCluster[means.scaled[row.names(perCluster),]<0.3] <- NA

## filter genes that don't reach the minimum fold-change
perCluster <- as.data.frame(perCluster)
perCluster$count <- rowSums(perCluster, na.rm=TRUE)
perCluster <- perCluster[perCluster$count>0,]
perCluster$gene <- de.ecto[row.names(perCluster),]$gene

## remove fold-change from results df
de.ecto[,3:4] <- t(sapply(1:nrow(de.ecto), function(x) 
  ifelse(is.na(perCluster[x,1:2]), NA, de.ecto[x,3:4][!is.na(perCluster[x,1:2])])))
# use largest fold-change as the summary value
de.ecto$logFC <- apply(de.ecto[,3:4], 1, function(x) x[which.max(abs(x))])

## define the direction of bias
bias.ecto <- ifelse(de.ecto[,3:4]>0, "left", "right")
bias.ecto <- bias.ecto[match(row.names(perCluster), row.names(bias.ecto)),]
bias.ecto[is.na(perCluster[,1:2])] <- "notExpr" # gene not expressed in cluster
bias.ecto[perCluster[,1:2]==0] <- "notSig" ## fold-change not large enough
bias.ecto <- as.data.frame(bias.ecto)
bias.ecto$gene <- de.ecto[row.names(bias.ecto),]$gene

bias.ecto$left <- apply(bias.ecto[,1:2], 1, function(x) sum(x=="left", na.rm = TRUE))
bias.ecto$right <- apply(bias.ecto[,1:2], 1, function(x) sum(x=="right", na.rm = TRUE))
bias.ecto$count <- bias.ecto$left+bias.ecto$right
# table(left=bias.ecto$left, right=bias.ecto$right)
#     right
# left  0  1  2
#    0  0 14 16
#    1 16  8  0
#    2  6  0  0

## flag the genes with bias to both sides as not-congruent
bias.ecto$congruent <- ifelse(bias.ecto$left>0 & bias.ecto$right>0, 0, 1)
```

Half of all genes are significant in both clusters; half of these are biased to the right and the other half is split between left and incongruent biases.

The 6 most significant are below.


```r
plots <- list()
for(i in 1:6){
  id <- row.names(de.ecto)[i]
  gene <- de.ecto[i,1]
  df <- data.frame(x=umap$x, y=umap$y,
                   side=sce$side, cluster=sce$clusterAnn, 
                   expr=logcounts(sce)[id,], gene = gene)
  df <- df[df$cluster %in% paste0("Ec", 1:2),]
  df <- df[order(df$expr),]

  plots[[i]] <- ggplot(df, aes(x, y, colour=expr)) +
    geom_point(size=1, alpha=0.8) +
    scale_color_gradientn(colours = c("grey", brewer.pal(n=9, "Blues"))) +
    ylab(gene) + xlab("") +
    xlim(3.9,6) +
    facet_wrap(~side) +
    th + theme(axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
}
ggarrange(plotlist = plots, ncol=2, nrow=3)
```

![](01_left_vs_right_files/figure-html/plot_ecto-1.png)<!-- -->

### Recurrence

The sets of genes differentially expressed in each germ layer are mostly independent. There are only six genes DE in the mesoderm clusters that are also DE in the endoderm (Eif3b, Polr2a, Gm13050) or in the ectoderm (Gm18821, Efhd2, Klf6).


```r
# intersect(de.meso$gene, de.endo$gene) # "Eif3b"   "Polr2a"  "Gm13050"
# intersect(de.meso$gene, de.ecto$gene) # "Gm18821" "Efhd2"   "Klf6"   
```

In the case of the cardiac mesoderm, the asymmetrically expressed genes have a diversity of expression patterns, with most being expressed selectively in some clusters. We can cluster these mesoderm DE genes into six clusters, based on their expression patterns. Some genes are expressed preferentially in the cells from clusters with a first heart field (FHF) signature (clusters 1 and 6 in the heatmap below) or second heart field (SHF) signature (clusters 2 and 3), while others are more prevalent in the mature cardiomyocytes from Me3 (cluster 4), or have a mixed expression pattern (cluster 5).


```r
## check how the asymmetric genes are expressed in the different cardiac clusters
means <- matrix <- matrix(nrow=nrow(de.meso), ncol=6)
colnames(means) <- paste0("Me",3:8)
row.names(means) <- row.names(de.meso)
for(c in paste0("Me",3:8)){
  cells <- which(sce$clusterAnn == c)
  means[,c] <- rowMeans(logcounts(sce)[row.names(de.meso),cells])
}
tmp <- t(apply(means, 1, function(x) (x-mean(x))/sd(x)))
set.seed(1827)
km = kmeans(tmp[,c(3:1,4:6)], centers = 6)
c <- factor(paste0("cluster", km$cluster), levels = paste0("cluster", c(6,1,4,3,2,5)))
names(c) <- row.names(tmp)

hm <- Heatmap(tmp[,c(3:1,4:6)], cluster_columns = FALSE, col=colorRamp2(breaks = c(-2,0,2), colors = c("steelblue","white", "indianred3")), name = "z-score", show_row_names = FALSE, split=c, cluster_row_slices = FALSE)
hm <- draw(hm)
```

![](01_left_vs_right_files/figure-html/cluster_meso-1.png)<!-- -->

```r
## FHF: 1,6; SHF: 3,2; cardiomyo: 4; progenitors: 5
```

The mean behaviour for each cluster clearly exemplifies the patterns of expression.


```r
plots <- list()
for(cluster in levels(c)){
  df <- data.frame(x=reducedDim(sce)[,1], y=reducedDim(sce)[,2], expr=colMeans(logcounts(sce)[names(c[c==cluster]),]))
  df <- df[order(df$expr),]
  plots[[cluster]] <- ggplot(df, aes(x,y,col=expr)) + geom_point() + scale_color_gradientn(colors=brewer.pal(n=9, "Blues")) + ggtitle(cluster) + th + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) + xlab("") + ylab("")
}

ggarrange(plotlist = plots, ncol=3, nrow=2, common.legend = TRUE, legend = "bottom")
```

![](01_left_vs_right_files/figure-html/meanExpr-1.png)<!-- -->

Interestingly, the genes with most significant p-values in the differential expression test are enriched for genes from cluster 3, which are primarily expressed in SHF cells (a third of all genes in the top 20 or 30 DE set). And the p-values for all genes in this cluster are much lower compared to the distributions in the rest of the clusters.


```r
stats <- data.frame(gene = de.meso$gene, FDR = de.meso$FDR, cluster = c[row.names(de.meso)], congruent = bias.meso[row.names(de.meso),]$congruent)

# round(table(stats$cluster, stats$congruent)/rowSums(table(stats$cluster, stats$congruent))*100,2)
## cluster2 has lower congruency rate than all others, followed by 4
## cluster5 and 1 have the highest

plot(density(stats$FDR), ylim=c(0,75), main="FDR per cluster", bty="l", lwd=2)
for(i in 1:6){
  lines(density(stats[stats$cluster == paste0("cluster",i),]$FDR), col=i, lwd=2)
}
legend("topright", legend = c(paste0("cluster",1:6), "all_DEGs"), lwd=2, col=c(1:6,"black"), cex=0.85)
```

![](01_left_vs_right_files/figure-html/stats-1.png)<!-- -->

```r
# ## FHF
# bias.meso[names(c[c=="cluster1"]),]
# bias.meso[names(c[c=="cluster6"]),]
# ## cardiomyocytes
# bias.meso[names(c[c=="cluster4"]),]
# ## SHF
# bias.meso[names(c[c=="cluster3"]),]
# bias.meso[names(c[c=="cluster2"]),]
# 
# bias.meso[names(c[c=="cluster5"]),]
```

This suggests that some of the strongest differences between the left and right cells arise in the SHF, whereas the asymmetries observed in the FHF cells are more subtle. 

This is supported by looking at the fold-changes of the genes in each cluster; the fold-change used is that of the cluster(s) with strongest expression, to avoid taking large fold-changes that arise from clusters with very low expression. That is:

- Cluster 6: fold-change of Me5 cells.
- Cluster 1: largest fold-change of Me3/4/5 cells.
- Cluster 4: fold-change of Me3 cells.
- Cluster 3: fold-change of Me6 cells.
- Cluster 2: largest fold-change of Me6/7/8 cells.
- Cluster 5: largest fold-change of all cells. This will retrieve meaningless fold-changes, and thus we ignore this cluster.

Genes from clusters 2 and 3, which are expressed in SHF cells (Me6/7/8) have much larger fold-changes, compared to the rest.


```r
stats$FC <- ifelse(stats$cluster == "cluster6", de.meso$logFC.Me5, ifelse(stats$cluster == "cluster1", apply(abs(de.meso[,3:5]), 1, max), ifelse(stats$cluster == "cluster4", de.meso$logFC.Me3, ifelse(stats$cluster == "cluster3", de.meso$logFC.Me6, ifelse(stats$cluster == "cluster2", apply(abs(de.meso[,6:8]), 1, max), apply(abs(de.meso[,3:8]), 1, max))))))

tmp <- stats[stats$cluster != "cluster5",]
ggplot(tmp, aes(cluster, abs(FC))) + geom_violin() + geom_boxplot(width=0.1) + th + xlab("") + ylab("|fold-change|")
```

```
## Warning: Removed 5 rows containing non-finite values (stat_ydensity).
```

```
## Warning: Removed 5 rows containing non-finite values (stat_boxplot).
```

![](01_left_vs_right_files/figure-html/fold-change-1.png)<!-- -->



```r
write.table(resAdj[['meso']], file=paste0(dir, "results/01_leftRight_DEGs_mesoderm.tsv"), 
            quote = FALSE, sep="\t")
write.table(resAdj[['endo']], file=paste0(dir, "results/01_leftRight_DEGs_endoderm.tsv"), 
            quote = FALSE, sep="\t")
write.table(resAdj[['ecto']], file=paste0(dir, "results/01_leftRight_DEGs_ectoderm.tsv"), 
            quote = FALSE, sep="\t")

write.table(de.meso, file=paste0(dir, "results/01_leftRight_DEGs_signif_mesoderm.tsv"), 
            quote = FALSE, sep="\t")
write.table(de.endo, file=paste0(dir, "results/01_leftRight_DEGs_signif_endoderm.tsv"), 
            quote = FALSE, sep="\t")
write.table(de.ecto, file=paste0(dir, "results/01_leftRight_DEGs_signif_ectoderm.tsv"), 
            quote = FALSE, sep="\t")

write.table(bias.meso, file=paste0(dir, "results/01_leftRight_bias_mesoderm.tsv"), 
            quote = FALSE, sep="\t")
write.table(bias.endo, file=paste0(dir, "results/01_leftRight_bias_endoderm.tsv"), 
            quote = FALSE, sep="\t")
write.table(bias.ecto, file=paste0(dir, "results/01_leftRight_bias_ectoderm.tsv"), 
            quote = FALSE, sep="\t")

write.table(stats, file=paste0(dir, "results/01_leftRight_DEGs_mesoderm_clusters.tsv"), 
            quote = FALSE, sep="\t")
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
##  [1] grid      parallel  stats4    stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] circlize_0.4.14             ComplexHeatmap_2.8.0       
##  [3] ggrepel_0.9.1               ggpubr_0.4.0               
##  [5] ggplot2_3.3.5               RColorBrewer_1.1-3         
##  [7] edgeR_3.34.1                limma_3.48.3               
##  [9] scran_1.20.1                scuttle_1.2.1              
## [11] SingleCellExperiment_1.14.1 SummarizedExperiment_1.22.0
## [13] Biobase_2.52.0              GenomicRanges_1.44.0       
## [15] GenomeInfoDb_1.28.4         IRanges_2.26.0             
## [17] S4Vectors_0.30.2            BiocGenerics_0.38.0        
## [19] MatrixGenerics_1.4.3        matrixStats_0.62.0         
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_2.0-3          ggsignif_0.6.3           
##  [3] rjson_0.2.21              ellipsis_0.3.2           
##  [5] bluster_1.2.1             XVector_0.32.0           
##  [7] GlobalOptions_0.1.2       BiocNeighbors_1.10.0     
##  [9] clue_0.3-60               rstudioapi_0.13          
## [11] farver_2.1.0              fansi_1.0.3              
## [13] splines_4.1.0             codetools_0.2-18         
## [15] sparseMatrixStats_1.4.2   doParallel_1.0.17        
## [17] knitr_1.38                jsonlite_1.8.0           
## [19] Cairo_1.5-15              broom_0.8.0              
## [21] cluster_2.1.3             png_0.1-7                
## [23] compiler_4.1.0            dqrng_0.3.0              
## [25] backports_1.4.1           assertthat_0.2.1         
## [27] Matrix_1.4-1              fastmap_1.1.0            
## [29] cli_3.2.0                 BiocSingular_1.8.1       
## [31] htmltools_0.5.2           tools_4.1.0              
## [33] rsvd_1.0.5                igraph_1.3.0             
## [35] gtable_0.3.0              glue_1.6.2               
## [37] GenomeInfoDbData_1.2.6    dplyr_1.0.8              
## [39] Rcpp_1.0.8.3              carData_3.0-5            
## [41] jquerylib_0.1.4           vctrs_0.4.1              
## [43] iterators_1.0.14          DelayedMatrixStats_1.14.3
## [45] xfun_0.30                 stringr_1.4.0            
## [47] beachmat_2.8.1            lifecycle_1.0.1          
## [49] irlba_2.3.5               statmod_1.4.36           
## [51] rstatix_0.7.0             zlibbioc_1.38.0          
## [53] scales_1.2.0              yaml_2.3.5               
## [55] gridExtra_2.3             sass_0.4.1               
## [57] stringi_1.7.6             highr_0.9                
## [59] foreach_1.5.2             ScaledMatrix_1.0.0       
## [61] BiocParallel_1.26.2       shape_1.4.6              
## [63] rlang_1.0.2               pkgconfig_2.0.3          
## [65] bitops_1.0-7              evaluate_0.15            
## [67] lattice_0.20-45           purrr_0.3.4              
## [69] labeling_0.4.2            cowplot_1.1.1            
## [71] tidyselect_1.1.2          magrittr_2.0.3           
## [73] R6_2.5.1                  magick_2.7.3             
## [75] generics_0.1.2            metapod_1.0.0            
## [77] DelayedArray_0.18.0       DBI_1.1.2                
## [79] pillar_1.7.0              withr_2.5.0              
## [81] abind_1.4-5               RCurl_1.98-1.6           
## [83] tibble_3.1.6              crayon_1.5.1             
## [85] car_3.0-12                utf8_1.2.2               
## [87] rmarkdown_2.13            GetoptLong_1.0.5         
## [89] locfit_1.5-9.5            digest_0.6.29            
## [91] tidyr_1.2.0               munsell_0.5.0            
## [93] bslib_0.3.1
```

