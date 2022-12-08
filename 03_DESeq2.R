#!/usr/bin/env Rscript

# Example call: DESeq2.R cKO WT ../Expression

args = commandArgs(trailingOnly=TRUE)

directory <- getwd()

treatment <- grep(args[1],list.files(directory),value=TRUE)
control <- grep(args[2],list.files(directory),value=TRUE)

sampleFiles <- c(treatment,control)
sampleFiles

condition <- c(rep('Treatment',length(treatment)), rep('Control',length(control)))

# If correcting for technical covariates, use: 
# batch<- c(1,1,1,1,1,1,1,1,2,2,2,2)
# batch <- c(1,1,2,2,2,1,1,1,2,2)

sampleTable <- data.frame(sampleName = sampleFiles,
                      fileName = sampleFiles,
                      condition = condition)
		     # batch = batch)
sampleTable

library("DESeq2")

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                   directory = directory,
                                   design= ~condition) # ~batch + condition)

dds <- DESeq(ddsHTSeq)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
head(resOrdered)

write.table(resOrdered, paste0(args[3],"/",args[1],"_vs_",args[2],".txt"), sep="\t")

#save.image("deseq2.RData")

