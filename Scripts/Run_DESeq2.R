#Get matrix using the condition(s) as a modeling covariate
dds <- DESeqDataSetFromMatrix(countData = rounded_matrix,
                              colData = dataDesign,
                              design = ~ condition)
     
#Run DESeq2 normalizations and export results
dds.deseq <- DESeq(dds)

#Model with my.Var1 as a variable
res.condition <- results(dds.deseq, name= "condition")
