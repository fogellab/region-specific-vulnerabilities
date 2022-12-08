#!/usr/bin/env Rscript

path="/Users/Dario/Desktop/Research/Fogel_Lab/Data/4.32_oligoTRAP_settings4.25t4.26e(html_done_final)/Input"
files=grep("*Aligned.sortedByCoord.out.bam.LastLinesRem.txt", list.files(path), value=T)
Countfiles=lapply(files, function(x) read.table(paste0(path, "/",  x), sep="\t", header=F))
Count <- Reduce(function(x, y) {
  merge(x, y, all=TRUE, by="V1")
}, Countfiles)
rownames(Count)=Count$V1
Count$V1=NULL
colnames(Count)=gsub("Aligned.sortedByCoord.out.bam.LastLinesRem.txt", "", files)
write.csv(Count, "oligoTRAP_count_matrix.csv", row.names=F)

# Remove genes with zero variance
data=t(Count)
data=data[, which(apply(data, 2, var) != 0)] 

# Log transform
library(DESeq2)
rlog_data=rlog(t(data))
write.csv(rlog_data, "oligoTRAP_rlog_matrix.csv", row.names=T)
matrix=cor(rlog_data)
head(matrix)

colors=rev(c(rep("blue",5), rep("chartreuse2", 3), rep("red3", 2), 
             rep("chartreuse2", 1), rep("red3", 5), rep("chartreuse2", 2), 
             rep("cyan3", 6)))

library(ggdendro)
library(ggplot2)
library(dplyr)
library(ape)
library(dendextend)

dend1 <- matrix %>% scale %>% dist %>% 
   hclust %>% as.dendrogram %>%
   set("branches_k_color", colors) #%>% #set("branches_lwd", 1.2) %>%
   #set("labels_colors") %>% 
   #set("labels_cex", 0.9) %>% 
   #set("leaves_pch", 19) %>% 
  #set("leaves_col", colors)

ggd1=as.ggdend(dend1)

dend2 <- hclust(dist(scale(matrix)), method="complete") %>% as.dendrogram #%>%
  #set("branches_k_color",rev(c(rep("blue",5),rep("cyan3",6),rep("chartreuse2",6),rep("red3",7)))) %>% set("branches_lwd", 1.2) %>%
  #set("labels_colors") %>% set("labels_cex", c(1)) %>% 
  #set("leaves_pch", 19) #%>% set("leaves_col", c("blue", "cyan3","chartreuse2","red3"))

ggd2=as.ggdend(dend2)

pdf("dendrogram.pdf", height=5, width=8)
ggplot(ggd1, horiz = TRUE, theme = NULL)+ 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	      panel.background = element_blank(), axis.text=element_text(color="black"), 
	      axis.line = element_line(colour = "black"))
dev.off()

pdf("dendrogram_uncolored.pdf", height=5, width=8)
ggplot(ggd2, horiz = TRUE, theme = NULL)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
