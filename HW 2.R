library(ggplot2)
library(sp)
library(DESeq2)
install.packages("gplots")
library(gplots)
library(RColorBrewer)

#**Obtain Data

getwd()
setwd("C:/Users/sharm/OneDrive/Desktop/CS10/HW2 R")
raw_counts = read.table("GSE190524_Acipenser_ruthenus_non-normalized_counts.csv", sep = ",", header = TRUE, row.names = 1)
# raw_counts = raw_counts[2:46028, 2:17]

meta_subset = read.table("Series Matrix_accession-condition-time-replicate.csv", sep = ",", header = TRUE, row.names = 1)

meta_subset = as.data.frame(meta_subset)

#**Separate Conditions

control = raw_counts[,1:4]
CPT = raw_counts[,5:10]
CPT48 = raw_counts[,5:7]
CPT72 = raw_counts[,8:10]
OLA = raw_counts[,11:16]
OLA24 = raw_counts[,11:13]
OLA48 = raw_counts[,14:16]

# Normalize data set
summed = as.numeric(unlist(sum(raw_counts[,3])))

Normalized = raw_counts

colsnum = as.numeric(unlist(ncol(Normalized)))

summed = as.numeric(unlist(sum(Normalized[,1])))

i = 0
for (i in 1:colsnum){
  i = as.numeric(unlist(i))
  Normalized[,i] = as.numeric(Normalized[,i])/as.numeric(unlist(sum(Normalized[,i])))
  
}

#reads per million

Normalized_per_mil = Normalized*1000000

#**Statistical Tests

#t-test() is for only normal datasets, and Wilcox is for  datasets with more than 30 samples

wilcox.test(unlist(Normalized_per_mil[11,5:10]),unlist(Normalized_per_mil[11,11:16]))

#**Violin Plot

raw_counts_labeled = cbind(meta_subset$treatment.condition, as.data.frame(t(raw_counts)))

colnames(raw_counts_labeled)[1] = "Condition"

ggplot(raw_counts_labeled, aes(x=Condition,y=ABA_contig_00009)) + geom_violin(fill = "pink") + geom_boxplot(width=.2, fill = "light blue")

#**Spearman Correlation

cor.test(CPT[,1], OLA[,1])

cor(CPT, OLA)

#**Scatter plot

plot(CPT[,1], OLA[,1])
 
abline(lm(CPT[,1]~OLA[,1]))


#**Histogram

#Please use a column of a data frame instead of a column of a matrix
hist(log(as.matrix(raw_counts)+1), breaks = 100)

#**Shapiro Test

shapiro.test(sample(as.matrix(raw_counts),5000))
t.test(unlist(CPT),unlist(OLA))
wilcox.test(unlist(Rads),unlist(Controls))

#**Heatmap

CorMat = cor(CPT, OLA)

heatmap.2(CorMat, margins = c(12,10), cexRow = .6, trace = "none", col = brewer.pal(9,"RdBu"), cexCol = .7, Rowv = FALSE, Colv = FALSE, colsep = 8, sepcolor = "black")

#**PCA Plot or volcano plot

?DESeqDataSetFromMatrix

MyDESeq = DESeqDataSetFromMatrix(countData = raw_counts, colData = meta_subset, design = ~treatment.condition)

MyDESeq = DESeq(MyDESeq)

MyDESeqResult = results(MyDESeq)

GeneStates = as.data.frame(MyDESeqResult@listData)

MyDESeqTime = DESeqDataSetFromMatrix(countData = raw_counts, colData = meta_subset, design = ~treatment.time)

MyDESeqTime = DESeq(MyDESeqTime)

MyDESeqTimeResult = results(MyDESeqTime)

GeneStatesTime = as.data.frame(MyDESeqTimeResult@listData)

