#Write results file with both genes and TEs included
res.condition.df <- as.data.frame(res.condition)

res.condition.df$ID <- rownames(res.condition.df)

write.table(res.condition.df, file = "TEs_Genes_Diff_Exp.txt", quote = F, row.names = F, sep = '\t')

#Write results file with only TEs included
res.condition.df.TE <- subset(res.condition.df, grepl("NotFur", res.condition.df$ID))

write.table(res.condition.df.TE file = "TEs_Only_Diff_Exp.txt", quote = F, row.names = F, sep = '\t')