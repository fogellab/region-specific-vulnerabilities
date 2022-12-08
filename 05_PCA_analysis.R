#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(factoextra)
library(pca3d)
library(limma)
library(umap)

# Read in data
data=read.csv("oligoTRAP_rlog_matrix.csv", row.names = 1)
head(data)
data=t(data)
res.pca <- prcomp(data, scale = T)
pdf(paste0(gsub(".csv","", args[2]),'.pdf'))

# Box and whisker plot
par(mar=c(7,4,2,1))
title=args[2]
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)

# expression value distribution plot
par(mar=c(4,4,2,1))
plotDensities(ex, main=title, legend=F)

# hierarchical clustering
sampleTree = hclust(dist(scale(data)))#, method = "average");
par(mar = c(0,4,2,0))
par(oma=c(1,1,1,1))
plot(sampleTree, main = title, sub="", xlab="", cex.lab = 1.5,cex.main=2)

# Principal component analysis
fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             col.ind = "black", # Color by the quality of representation
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
                )
groups <- as.factor(c(rep('Cerebellum',6),rep('Frontal cortex',7),rep('Spinal cord',5),rep('Corpus callosum',6)))

fviz_pca_ind(res.pca,
             #col.ind = groups, # color by groups
             #palette = c("#00AFBB",  "#FC4E07"),
             habillage=groups,
	     addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             ellipse.level = 0.95,
	     legend.title = "Region",
             label="none"
             )

fviz_pca_ind(res.pca, #use this one
	     axes=c(1,2),
	     axes.linetype=NA,
             palette = c("cyan3","chartreuse2", "red3", "blue"),#c("green4","gold2", "red3", "deepskyblue3"),#c("cyan3","green3", "gold2", "purple"),
             habillage=groups,invisible="quali",
             legend.title = "Region",
             label="none"
             )+
        xlim(-120, 105)+
        ylim(-130, 70)+
        scale_shape_manual(values=c(15,16,17,18))+scale_size_manual(values=rep(100,4))+
	       theme(
          panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  	     labs(title=args[2], x = "PC1 (19.6%)", y = "PC2 (12.2%)")
dev.off()

ggsave("PCA_plot.pdf", height = 3, width = 6, device='pdf')
