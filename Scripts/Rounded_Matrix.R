#Format matrix, make IDs into rownames for DESeq2
rownames(Count_Matrix) <- Count_Matrix$Gene

#Keep only Transcripts with non critically low coverage 
my.good <- which(apply(Count_Matrix[,-1]>0, 1, sum) >= 6)
my.filtered.matrix <- Count_Matrix[my.good,-1]

#Prerequisite: round the counts to nearest integer. Unnecessary when using featureCounts with unique counts only.
rounded_matrix <- round(my.filtered.matrix)