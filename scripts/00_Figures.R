## Code to reproduce the figures
library(scran)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

# dir <- "/Users/ibarra01/OneDrive - CRUK Cambridge Institute/github/leftRight_asymmetry/"
dir <- "/Users/xi629080/Documents/P_Dev/Tyser/leftRight_asymmetry/"
out <- "/Users/xi629080/Documents/P_Dev/Tyser/leftRight_asymmetry/figs/"
# out <- "/Users/ibarra01/OneDrive - CRUK Cambridge Institute/WRITING/HEARTleft_right/Figures/figureElements/"

palette(brewer.pal(n=12, "Set3")[-c(1:2)])

th <- theme_bw() + theme(axis.text.x = element_text(size=12), 
                         axis.title.x = element_text(size=12), 
                         axis.text.y = element_text(size=12), 
                         axis.title.y = element_text(size=12), 
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), 
                         axis.line = element_line(colour = "black"), 
                         panel.border = element_blank(), 
                         plot.title = element_text(face="bold", hjust = 0.5))

## cluster colours
cols <- c(Ec1 = "#ec6646", Ec2 = "#af4424", En1 = "#3c537c", En2 = "#768ba5",
          Me1 = "#bf9a77", Me2 = "#debc95", Me3 = "#556dad", Me4 = "#f28a31", 
          Me5 = "#729f3c", Me6 = "#fbba14", Me7 = "#5fa398", Me8 = "#9FD3C5")
cols2 <- c(Ec1 = "#F48F7A", Ec2 = "#C1634D", En1 = "#5378AD", En2 = "#96B0CC",
          Me1 = "#E5B893", Me2 = "#F7D9BC", Me3 = "#7490D1", Me4 = "#F7A466", 
          Me5 = "#A6D366", Me6 = "#F4C969", Me7 = "#7DD3C2", Me8 = "#BDEFE1")


### Figure 1 ###############
sce <- readRDS(paste0(dir, "data/sce_goodQual.NORM.clusters.Rds"))
sce <- sce[,-which(sce$stage %in% c(-1,"LHT"))]

## C: UMAP with clusters ======
umap <- reducedDim(sce)
umap$cluster <- sce$clusterAnn
umap$side <- factor(sce$side, levels=c("Right", "Left"))
# umap$col <- cols[umap$cluster]
# umap[umap$side=="Right",]$col <- cols2[umap[umap$side == "Right",]$cluster]

bg <- umap
bg$side <- NULL
o <- sample(1:nrow(umap), replace = FALSE)
pdf(paste0(out, "Fig1c_UMAP.pdf"), width = 7, height = 3.5, useDingbats = FALSE)
ggplot(umap[o,], aes(x, y)) +
  geom_point(data=bg, colour="grey80", size=0.5) +
  geom_point(aes(colour=cluster), size=0.35) +
  scale_color_manual(values = cols) +
  facet_wrap(~side) +
  xlab("") +
  ylab("") +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  th + theme(axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks = element_blank(),
             strip.background = element_rect(fill=NA, size = 1),
             strip.text = element_text(face="bold", size=11))
dev.off()

# pdf(paste0(out, "Fig1c_UMAP.pdf"), width = 7, height = 7, useDingbats = FALSE)
# plot(umap$x[o], umap$y[o], 
#      col="black", lwd=0.5,
#      bg=umap$col[o],
#      pch=umap$pch, cex=0.85, 
#      axes=FALSE, 
#      xlab="", ylab="")
# box(bty="l")
# dev.off()

## D: left-side genes ======
plots <- list()
for(gene in c("Lefty2", "Pitx2", "Nodal")){
  id <- row.names(rowData(sce)[rowData(sce)$gene == gene,])
  umap$expr <- logcounts(sce)[id, row.names(umap)]
  o <- order(umap$expr)
  plots[[gene]] <- ggplot(umap[o,], aes(x, y, colour=expr)) +
  geom_point(size=0.75, alpha=0.8) +
  scale_color_gradientn(colours = c("grey", brewer.pal(n=9, "Blues"))) +
  ylab(gene) + xlab("") +
  facet_wrap(~side) +
  th + theme(axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks = element_blank(),
             panel.border = element_rect(fill=NA, size=1),
             strip.background = element_rect(fill=NA, size = 1),
             strip.text = element_text(face="bold", size=11))
}

pdf(paste0(out, "Fig1d_leftGenes.pdf"), width = 7, height = 10.5, useDingbats = FALSE)
ggarrange(plotlist = plots, nrow=3, align="hv")
dev.off()


### Figure 2 ###############

## A: mesoderm DE analysis ======
## mesoderm DEGs
de.meso <- read.table(paste0(dir, "results/01_leftRight_DEGs_mesoderm.tsv"))
# plot  average fold-change for all clusters
de.meso$logFC <- rowMeans(de.meso[,3:8]) 
# for DEGs, plot the largest fold-change from all clusters
de.meso.sig <- read.table(paste0(dir, "results/01_leftRight_DEGs_signif_mesoderm.tsv"))
de.meso$sig <- ifelse(row.names(de.meso) %in% row.names(de.meso.sig), "DEG", "notSignificant")
de.meso[row.names(de.meso.sig),]$logFC <- de.meso.sig$logFC
de.meso$label <- ifelse(de.meso$gene %in% c("Nodal", "Pitx2", "Lefty2", "Gal",
                                        "Acss1", "Mmp9", "Mecom", "Sema3a", "Srgap1",
                                        "Eif3b", "Ybx3"), de.meso$gene, "")

pdf(paste0(out, "Fig1f_volcano_meso.pdf"), width = 7, height = 6, useDingbats = FALSE)
ggplot(de.meso, aes(logFC, -log10(FDR), label=label)) +
  geom_point(aes(col=sig)) + 
  scale_color_manual(name = "", values = c("red", "grey")) + 
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), lty=2) + 
  geom_hline(yintercept = -log10(0.05), lty=2) +
  xlab(expression('log'[2]*' fold-change (left / right)')) + 
  ylab(expression('-log'[10]*' p-value')) + 
  ylim(c(0,40)) +
  geom_point(data=de.meso[de.meso$label != "",], aes(logFC, -log10(FDR)), colour="black", pch=21) +
  geom_text_repel(min.segment.length = 0, direction="both", max.overlaps = 50) + 
  th 
# + scale_x_continuous(breaks = c(-6,-4,-2,0,2),labels = c(6,4,2,0,-2))
dev.off()


## C: DE genes ======
plots <- list()
for(gene in c("Acss1","Sema3a","Srgap1","Ybx3")){
  p <- list()
  id <- row.names(rowData(sce)[rowData(sce)$gene == gene,])
  umap$expr <- logcounts(sce)[id, row.names(umap)]
  
  ## only Me3-8
  tmp <- umap[row.names(umap[umap$cluster %in% paste0("Me", 3:8),]),]
  
  ## violin
  p[[paste0("1_",gene)]] <- ggplot(tmp, aes(side, expr, fill=expr)) + 
    geom_violin() + 
    geom_jitter(width=0.05, pch=21, colour="grey20") +
    stat_summary(fun = median,
                 fun.min = function(z) { quantile(z,0.25) },
                 fun.max = function(z) { quantile(z,0.75) },
                 geom="pointrange", color="red") + 
    ggtitle(gene) +
    xlab("") +
    ylab(expression('log'[2]*' expression')) + 
    scale_fill_gradientn(colours = c("grey", brewer.pal(n=9, "Blues"))) +
    th + theme(plot.title = element_text(face="italic", hjust = 0.5),
               legend.position = "none")
  
  ## proportion cells > 0
  df <- as.data.frame(prop.table(table(tmp$side, tmp$expr>0), 1)*100)
  colnames(df) <- c("side", "expr", "prop")
  
  p[[paste0("2_", gene)]] <- ggplot(df, aes(side, prop, fill=expr)) + 
    geom_bar(stat="identity") +
    ggtitle(gene) +
    xlab("") +
    ylab("% cells with counts > 0") + 
    scale_fill_manual(values = c("grey", "steelblue4")) +
    th + theme(plot.title = element_text(face="italic", hjust = 0.5),
               legend.position = "none")
  
  plots[[gene]] <- ggarrange(plotlist = p, widths = c(0.65,0.35), align="hv")
}


pdf(paste0(out, "Fig2_DEgenes_violin.pdf"), width = 18, height = 3, useDingbats = FALSE)
ggarrange(plotlist = plots, nrow=1, ncol=4, align="hv")
dev.off()


### Figure 3 ###############

## A,D: DE genes ======
plots <- list()
for(gene in c("Mecom","Mmp9")){
  p <- list()
  id <- row.names(rowData(sce)[rowData(sce)$gene == gene,])
  umap$expr <- logcounts(sce)[id, row.names(umap)]
  
  ## only Me3-8
  tmp <- umap[row.names(umap[umap$cluster %in% paste0("Me", 3:8),]),]
  
  ## violin
  p[[paste0("1_", gene)]] <- ggplot(tmp, aes(side, expr, fill=expr)) + 
    geom_violin() + 
    geom_jitter(width=0.1, pch=21, colour="grey20") +
    stat_summary(fun = median,
                 fun.min = function(z) { quantile(z,0.25) },
                 fun.max = function(z) { quantile(z,0.75) },
                 geom="pointrange", color="red") +
    ylab(gene) + xlab("") +
    scale_fill_gradientn(colours = c("grey", brewer.pal(n=9, "Blues"))) + 
    facet_wrap(~cluster, ncol=6) +
    th + theme(panel.border = element_rect(fill=NA, size=1),
               strip.background = element_rect(fill=NA, size = 1),
               strip.text = element_text(face="bold", size=11),
               legend.position = "none")
  
  ## proportion cells > 0
  df <- rbind(
    cbind(aggregate(expr ~ side + cluster, data = tmp, FUN = \(x) mean(x > 0)), expr="TRUE"),
    cbind(aggregate(expr ~ side + cluster, data = tmp, FUN = \(x) 1-mean(x > 0)), expr="FALSE")
  )
  colnames(df) <- c("side", "cluster", "prop", "expr")
  
  p[[paste0("2_", gene)]] <- ggplot(df, aes(side, prop*100, fill=expr)) + 
    geom_bar(stat="identity") +
    xlab("") +
    ylab("% cells with counts > 0") + 
    scale_fill_manual(values = c("grey", "steelblue4")) +
    facet_wrap(~cluster, ncol=6) +
    th + theme(panel.border = element_rect(fill=NA, size=1),
               strip.background = element_blank(),
               strip.text.x = element_blank(),
               legend.position = "none")
  plots[[gene]] <- ggarrange(plotlist = p, nrow=2, ncol=1, align="hv")
}

pdf(paste0(out, "Fig3_DEgenes_violin.pdf"), width = 12, height = 8, useDingbats = FALSE)
ggarrange(plotlist = plots, ncol=1, nrow=2, align="hv")
dev.off()



### Figure S1 ###############

## C: Cell-cycle proportions ======
phase <- read.table(paste0(dir, "data/cellCyclePhase.tsv"))
sce$phase <- phase[match(colnames(sce), phase$V1),2]
df <- data.frame(cluster=sce$clusterAnn, side=sce$side, phase=sce$phase)

left <- prop.table(table(df[df$side=="Left",]$cluster, df[df$side=="Left",]$phase),1)*100
right <- prop.table(table(df[df$side=="Right",]$cluster, df[df$side=="Right",]$phase),1)*100

df <- as.data.frame(rbind(left, right))
df$side <- rep(c("left", "right"), each = 12)
df$cluster <- substr(row.names(df), 1, 3)
df <- reshape2::melt(df)

df$variable <- factor(df$variable, levels = rev(c("G1", "S", "G2M")))

pdf(paste0(out, "FigX_cellCycle.pdf"), width = 7, height = 4, useDingbats = FALSE)
ggplot(df, aes(side, value, fill=variable)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = rev(c("grey60","burlywood2","burlywood4"))) +
  facet_wrap(~cluster, ncol=6) +
  ylab("% cells") +
  labs(fill="phase") +
  th + theme(strip.background = element_rect(fill=NA))
dev.off()


### Figure S2 ###############

## A: endoderm DE analysis ======
de.endo <- read.table(paste0(dir, "results/01_leftRight_DEGs_endoderm.tsv"))
# plot  average fold-change for all clusters
de.endo$logFC <- rowMeans(de.endo[,3:4]) 
# for DEGs, plot the largest fold-change from all clusters
de.endo.sig <- read.table(paste0(dir, "results/01_leftRight_DEGs_signif_endoderm.tsv"))
de.endo$sig <- ifelse(row.names(de.endo) %in% row.names(de.endo.sig), "DEG", "notSignificant")
de.endo[row.names(de.endo.sig),]$logFC <- de.endo.sig$logFC
de.endo$label <- ifelse(de.endo$gene %in% c("Eif3b", "Polr2a", "Gm13050"), de.endo$gene, "")

pdf(paste0(out, "FigS2a_volcano_endo.pdf"), width = 7, height = 6, useDingbats = FALSE)
ggplot(de.endo, aes(logFC, -log10(FDR), label=label)) +
  geom_point(aes(col=sig)) + 
  scale_color_manual(name = "", values = c("red", "grey")) + 
  geom_vline(xintercept = c(-log2(1.2), log2(1.2)), lty=2) + 
  geom_hline(yintercept = -log10(0.05), lty=2) +
  xlab(expression('log'[2]*' fold-change (left / right)')) + 
  ylab(expression('-log'[10]*' p-value')) + 
  ggtitle("Endoderm") +
  ylim(c(0,8)) +
  geom_point(data=de.endo[de.endo$label != "",], aes(logFC, -log10(FDR)), colour="black", pch=21) +
  geom_text_repel(min.segment.length = 0, direction="both", max.overlaps = 50) +
  th 
dev.off()

## B: ectoderm DE analysis ======
de.ecto <- read.table(paste0(dir, "results/01_leftRight_DEGs_ectoderm.tsv"))
# plot  average fold-change for all clusters
de.ecto$logFC <- rowMeans(de.ecto[,3:4]) 
# for DEGs, plot the largest fold-change from all clusters
de.ecto.sig <- read.table(paste0(dir, "results/01_leftRight_DEGs_signif_ectoderm.tsv"))
de.ecto$sig <- ifelse(row.names(de.ecto) %in% row.names(de.ecto.sig), "DEG", "notSignificant")
de.ecto[row.names(de.ecto.sig),]$logFC <- de.ecto.sig$logFC
de.ecto$label <- ifelse(de.ecto$gene %in% c("Gm18821", "Efhd2", "Klf6"), de.ecto$gene, "")

pdf(paste0(out, "FigS2b_volcano_ecto.pdf"), width = 7, height = 6, useDingbats = FALSE)
ggplot(de.ecto, aes(logFC, -log10(FDR), label=label)) +
  geom_point(aes(col=sig)) + 
  scale_color_manual(name = "", values = c("red", "grey")) + 
  geom_vline(xintercept = c(-log2(1.2), log2(1.2)), lty=2) + 
  geom_hline(yintercept = -log10(0.05), lty=2) +
  xlab(expression('log'[2]*' fold-change (left / right)')) + 
  ylab(expression('-log'[10]*' p-value')) + 
  ggtitle("Ectoderm") +
  ylim(c(0,8)) +
  geom_point(data=de.ecto[de.ecto$label != "",], aes(logFC, -log10(FDR)), colour="black", pch=21) +
  geom_text_repel(min.segment.length = 0, direction="both", max.overlaps = 50) +
  th 
dev.off()

## C: expression pattern of mesoderm DEGs ======
means <- matrix <- matrix(nrow=nrow(de.meso.sig), ncol=6)
colnames(means) <- paste0("Me",3:8)
row.names(means) <- row.names(de.meso.sig)
for(c in paste0("Me",3:8)){
  cells <- which(sce$clusterAnn == c)
  means[,c] <- rowMeans(logcounts(sce)[row.names(de.meso.sig),cells])
}
tmp <- t(apply(means, 1, function(x) (x-mean(x))/sd(x)))
set.seed(1827)
km = kmeans(tmp[,c(3:1,4:6)], centers = 6)
c <- factor(paste0("cluster", km$cluster), levels = paste0("cluster", c(6,1,4,3,2,5)))
names(c) <- row.names(tmp)

hm <- Heatmap(tmp[,c(3:1,4:6)], cluster_columns = FALSE, 
              col=colorRamp2(breaks = c(-2,0,2), colors = c("steelblue","white", "indianred3")), 
              name = "z-score", show_row_names = FALSE, split=c, cluster_row_slices = FALSE)

pdf(paste0(out, "FigS2c_exprPatterns_meso.pdf"), width = 5, height = 5, useDingbats = FALSE)
hm <- draw(hm)
dev.off()

## plot average expression
plots <- list()
for(cluster in levels(c)){
  df <- data.frame(x=reducedDim(sce)[,1], y=reducedDim(sce)[,2], expr=colMeans(logcounts(sce)[names(c[c==cluster]),]))
  df <- df[order(df$expr),]
  plots[[cluster]] <- ggplot(df, aes(x,y,col=expr)) + 
    geom_point(size=0.5) + 
    scale_color_gradientn(colors=brewer.pal(n=9, "Blues")) + 
    ggtitle(cluster) + 
    xlab("") + ylab("") +
    th + theme(axis.ticks = element_blank(), 
               axis.text.x = element_blank(), 
               axis.text.y = element_blank())
}

pdf(paste0(out, "FigS2c_meanExpr_meso.pdf"), width = 6, height = 4, useDingbats = FALSE)
ggarrange(plotlist = plots, ncol=3, nrow=2, common.legend = TRUE, legend = "bottom")
dev.off()


### Figure S6 ###############

## A: DE genes ======
plots <- list()
for(gene in c("Acss1","Sema3a","Srgap1","Ybx3","Mmp9","Mecom")){
    id <- row.names(rowData(sce)[rowData(sce)$gene == gene,])
    umap$expr <- logcounts(sce)[id, row.names(umap)]
    o <- order(umap$expr)
    plots[[gene]] <- ggplot(umap[o,], aes(x, y, colour=expr)) +
      geom_point(size=0.75, alpha=0.8) +
      scale_color_gradientn(colours = c("grey", brewer.pal(n=9, "Blues"))) +
      ylab(gene) + xlab("") +
      facet_wrap(~side) +
      th + theme(axis.text.x = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks = element_blank(),
                 panel.border = element_rect(fill=NA, size=1),
                 strip.background = element_rect(fill=NA, size = 1),
                 strip.text = element_text(face="bold", size=11))
  }

pdf(paste0(out, "FigX_leftRight_UMAPs.pdf"), width = 7, height = 21, useDingbats = FALSE)
ggarrange(plotlist = plots, nrow=6, align="hv")
dev.off()


