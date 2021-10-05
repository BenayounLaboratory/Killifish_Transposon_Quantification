#Set variables to be investigated, ex. age
my.Var1  <- c(rep(Condition1, Number of Occurrences), rep(Condition1, Number of Occurrences))

#Design matrix
dataDesign = data.frame( row.names = colnames( rounded_matrix ), 
                         condition = my.Var1)

20. Run DESeq2 test and normalization.

#Get matrix using the condition(s) as a modeling covariate
dds <- DESeqDataSetFromMatrix(countData = rounded_matrix,
                              colData = dataDesign,
                              design = ~ condition)
     
#Run DESeq2 normalizations and export results
dds.deseq <- DESeq(dds)

#Model with my.Var1 as a variable
res.condition <- results(dds.deseq, name= "condition")