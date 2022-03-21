#Prepare TE Heatmap
     
#Get the heatmap of conditional changes at FDR5; exclude NA
res.condition <- res.condition[!is.na(res.condition$padj),]

#Keep significantly differentially expressed transcripts
transcripts.condition <- rownames(res.condition)[res.condition$padj < 0.05]

transcript_list <- as.data.frame(transcripts.condition)

#Exclude any genes; retain only TEs
TE_list <- subset(transcript_list, grepl("NotFur", transcript_list$transcripts.condition))

TE_chars <- as.character(TE_list$transcripts.condition)
