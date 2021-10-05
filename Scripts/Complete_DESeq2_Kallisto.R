#Load Packages
library(DESeq2)
library(RColorBrewer)
library(pheatmap)

#Import Kallisto Count Matrix
Count_Matrix <- read.table("featureCounts_count_matrix.txt", header = T, comment.char = "")

#Remove length and eff length columns. N denotes the final count #column
Count_Matrix <- Count_Matrix[,c(1,4:N)]

#Format matrix, make IDs into rownames for DESeq2
rownames(Count_Matrix) <- Count_Matrix$Gene

#Keep only Transcripts with non critically low coverage 
my.good <- which(apply(Count_Matrix[,-1]>0, 1, sum) >= 6)
my.filtered.matrix <- Count_Matrix[my.good,-1]

#Prerequisite: round the counts to nearest integer. Unnecessary #when using featureCounts with unique counts only.
rounded_matrix <- round(my.filtered.matrix)

#Set variables to be investigated, ex. age
my.Var1  <- c(rep(Condition1, Number of Occurrences), rep(Condition1, Number of Occurrences))

#Design matrix
dataDesign = data.frame( row.names = colnames( rounded_matrix ), 
                         condition = my.Var1)

#Get matrix using the condition(s) as a modeling covariate
dds <- DESeqDataSetFromMatrix(countData = rounded_matrix,
                              colData = dataDesign,
                              design = ~ condition)

#Run DESeq2 normalizations and export results
dds.deseq <- DESeq(dds)

#Model with my.Var1 as a variable
res.condition <- results(dds.deseq, name= "condition")

#Write results file with both genes and TEs included
res.condition.df <- as.data.frame(res.condition)

res.condition.df$ID <- rownames(res.condition.df)

write.table(res.condition.df, file = "TEs_Genes_Diff_Exp.txt", quote = F, row.names = F, sep = '\t')

#Write results file with only TEs included
res.condition.df.TE <- subset(res.condition.df, grepl("NotFur", res.condition.df$ID))

write.table(res.condition.df.TE, file = "TEs_Only_Diff_Exp.txt", quote = F, row.names = F, sep = '\t')

#Normalized expression value
norm.cts <- getVarianceStabilizedData(dds.deseq)

#Prepare TE Heatmap

#Get the heatmap of conditional changes at FDR5; exclude NA
res.condition <- res.condition[!is.na(res.condition$padj),]

#Keep significantly differentially expressed transcripts
transcripts.condition <- rownames(res.condition)[res.condition$padj < 0.05]

transcript_list <- as.data.frame(transcripts.condition)

#Exclude any genes; retain only TEs
TE_list <- subset(transcript_list, grepl("NotFur", transcript_list$transcripts.condition))

TE_chars <- as.character(TE_list$transcripts.condition)

#Create heatmap
pheatmap(norm.cts[TE_chars,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = "Title", cellwidth = 15, border = NA)