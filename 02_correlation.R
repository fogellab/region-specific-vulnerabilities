#!/usr/bin/env Rscript

# Load rlog tranformed counts
rlog_data = read.csv("oligoTRAP_rlog_matrix.csv", row.names = 1)

# Load DESeq count matrix
path='/Users/Dario/Desktop/Research/Fogel_Lab/Data/4.32_oligoTRAP_settings4.25t4.26e(html_done_final)/Input'
files=grep("Aligned.sortedByCoord.out.bam.LastLinesRem.txt", 
           list.files(path), value=T)
counts=lapply(files, function(x) read.table(paste0(path, "/", x), header=F, sep="\t"))
count_data=Reduce(function(x, y) {
  merge(x, y, all=TRUE, by="V1")
}, Countfiles)
rnames=count_data$V1
count_data=count_data[,-1]
colnames(count_data)=gsub("Aligned.sortedByCoord.out.bam.LastLinesRem.txt", "", files)
rownames(count_data)=rnames
#write.csv(count_data, "oligoTRAP_count_matrix.csv")
count_data=read.csv("oligoTRAP_count_matrix.csv", row.names=1)

# Dependencies
library(ggpubr)
library(stringr)
library(scales)
library(DESeq2)

# given two groups, plot correlation and DE genes
flamePlot = function(rlog_data, count_data, group1, group2, color1, color2){
  
  # If group1 and group2 are the same, then randomly assign replicates
  if(group1==group2){
    conditions=colnames(count_data)
    selected=conditions[startsWith(conditions, group1)]
    n.replicates=length(selected)
    g1=sample(selected, n.replicates/2)
    g2=selected[!selected %in% g1]
    conditions[conditions %in% g1]=gsub(group1, paste0(group1, "r1"), conditions[conditions %in% g1])
    conditions[conditions %in% g2]=gsub(group1, paste0(group2, "r2"), conditions[conditions %in% g2])
    colnames(count_data)=conditions
    colnames(rlog_data)=conditions
    group1=paste0(group1, "r1")
    group2=paste0(group2, "r2")
  }
  
  # Differential gene expression
  conditions=colnames(count_data)
  subset_counts=cbind(count_data[,startsWith(conditions, group1)], 
                      count_data[,startsWith(conditions, group2)])
  group = str_split_fixed(colnames(subset_counts), "_", 2)[,1]
  group = factor(group, levels=unique(group))
  dds1 = DESeqDataSetFromMatrix(subset_counts, data.frame(condition=group), ~ condition)
  dds = DESeq(dds1)
  res = results(dds)
  resOrdered <- as.data.frame(res[order(res$pvalue),])
  
  # Plot correlation
  mean_data=cbind(rowMeans(rlog_data[,startsWith(conditions, group1)]), rowMeans(rlog_data[,startsWith(conditions, group2)]))
  colnames(mean_data)=c(group1, group2)
  rownames(resOrdered)=toupper(rownames(resOrdered))
  diff.1b=resOrdered[rownames(resOrdered) %in% rownames(mean_data),]
  layer1=merge(mean_data, diff.1b, by=0, all=TRUE)
  #layer1=cbind(mean_data, resOrdered[match(rownames(mean_data), toupper(rownames(resOrdered))),])
  
  plot=ggplot(layer1, aes(x=eval(parse(text = group1)), y=eval(parse(text = group2))))+
    xlab(group1)+
    ylab(group2)+
    scale_x_continuous(limits=c(-1, 18), labels = math_format(2^.x))+
    scale_y_continuous(limits=c(-1, 18), labels = math_format(2^.x))+
    #scale_color_identity(name = 'Legend', guide = 'legend', labels = c(args[4], args[5], args[6])) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "white", colour = "black", size=0.75), 
          axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))+
    #theme(plot.title = element_text(hjust = 0.5, size=35), axis.text=element_text(size=28), axis.title=element_text(size=32))+
    #theme_classic()+
    geom_point(data = subset(layer1, abs(log2FoldChange)<1 & !padj<0.01), alpha = 0.04, size=0.5)+
    geom_point(data = subset(layer1, log2FoldChange>1 & padj<0.01), color = color2, alpha=1, size=0.5)+
    geom_point(data = subset(layer1, log2FoldChange<(-1) & padj<0.01), color = color1, alpha=1, size=0.5)+
    geom_abline(intercept=0, slope = 1, size=0.75)+
    #geom_label_repel(data=layer1 %>% filter(abs(log2FoldChange)>10 & padj<0.01), aes(label=Row.names), seed=3, fill="white", size=7, max.overlaps = 100)+
    coord_cartesian(clip = "off") +
    annotate(geom="text", x=15, y=0, label=paste0("r = ", sprintf("%.3f", round(cor(layer1[group1], layer1[group2]), 3))), size=4, color="black")+
    annotate(geom="text", x=4, y=16, label=paste0(nrow(layer1[abs(layer1$log2FoldChange)>1 & layer1$padj<0.01 & !is.na(layer1$padj),]), " DEGs"), size=4, color="black")
  
  return(plot)
}

# Experimental
WM_cortex=flamePlot(rlog_data, count_data, "WM", "cortex", "chartreuse2", "red3")
CBL_cortex=flamePlot(rlog_data, count_data, "CBL", "cortex", "cyan3", "red3")
SC_cortex=flamePlot(rlog_data, count_data, "SC", "cortex", "blue", "red3")
CBL_WM=flamePlot(rlog_data, count_data, "CBL", "WM", "cyan3", "chartreuse2")
SC_WM=flamePlot(rlog_data, count_data, "SC", "WM", "blue", "chartreuse2")
SC_CBL=flamePlot(rlog_data, count_data, "SC", "CBL", "blue", "cyan3")

pdf("experimental_flame.pdf", height=5, width=7.5)
ggarrange(WM_cortex, CBL_cortex, SC_cortex, CBL_WM, SC_WM, SC_CBL, nrow=2, ncol=3)
dev.off()

# Controls
set.seed(1)
SC_SC=flamePlot(rlog_data, count_data, "SC", "SC", "blue", "blue")
CBL_CBL=flamePlot(rlog_data, count_data, "CBL", "CBL", "cyan3", "cyan3")
WM_WM=flamePlot(rlog_data, count_data, "WM", "WM", "chartreuse2", "chartreuse2")
cortex_cortex=flamePlot(rlog_data, count_data, "cortex", "cortex", "red3", "red3")

pdf("control_flame.pdf", height=5, width=5)
ggarrange(cortex_cortex, WM_WM, CBL_CBL, SC_SC, nrow=2, ncol=2)
dev.off()
