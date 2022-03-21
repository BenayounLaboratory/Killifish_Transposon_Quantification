#Constructing the 2020 Genome GTF file with the embedded TEs from fishTEDB.org
#load packages

library(gdata)

#set working directory to the direcotry in which 1) the Repeatmasker .out file is stored and 2) the genome genic .gtf file is stored

setwd("directory_with_Repeatmasker.out_file_and_Genic_GTF_file")

# before importing the repeatmasker outfile, within a text editor, 
# 1) replace spaces with tabs (find "\t" ; replace " ")
# 2) remove * column (find "\t*" ; replace "")

outfile <- read.table("path_to_modified_(above)_Repeatmasker_.out_file", skip = 3, sep = "\t", comment.char = "")

#change "complementary" to "-"
outfile$V10 <- gsub("C", "-",outfile$V10)

#keep only positinal info and names; strand
outfile <- outfile[,c(6:8,10:12)]
colnames(outfile) <- c("chr", "start", "stop", "strand", "name", "family")

#generate all the formatted gene, transcript, exon names for each TE
outfile$gene <- paste0('gene_id "',outfile$name,'"; transcript_id ""; description ""; gbkey "Gene"; gene_biotype "protein_coding"; locus_tag "',outfile$name,'"; note "";')

outfile$transcript <- paste0('gene_id "',outfile$name,'"; transcript_id "',outfile$name,'"; gbkey "mRNA"; locus_tag "',outfile$name,'"; orig_protein_id "',outfile$name,'"; orig_transcript_id "',outfile$name,'"; product ""; transcript_biotype "mRNA";')

outfile$exon <- paste0('gene_id "',outfile$name,'"; transcript_id "',outfile$name,'"; locus_tag "',outfile$name,'"; orig_protein_id "',outfile$name,'"; orig_transcript_id "',outfile$name,'"; product ""; transcript_biotype "mRNA"; exon_number "1";')

#now add other needed metrics (position, type)

outfile$pos1 <- "."
outfile$pos2 <- "."
outfile$TE <- "TE"

outfile <- outfile[,c("chr","TE","start","stop","pos1","strand","pos2", "gene", "transcript", "exon")]

#now interleave to make a gene, transcript, and exon line for each TE location

gene <- outfile[,c(1:8)]
gene$type <- "gene"
gene <- gene[,c("chr","TE", "type","start","stop","pos1","strand","pos2","gene")]
colnames(gene) <- c("chr","TE", "type","start","stop","pos1","strand","pos2","descriptor")

transcript <- outfile[,c(1:7,9)]
transcript$type <- "transcript"
transcript <- transcript[,c("chr","TE", "type","start","stop","pos1","strand","pos2","transcript")]
colnames(transcript) <- c("chr","TE", "type","start","stop","pos1","strand","pos2","descriptor")

exon <- outfile[,c(1:7,10)]
exon$type <- "exon"
exon <- exon[,c("chr","TE", "type","start","stop","pos1","strand","pos2","exon")]
colnames(exon) <- c("chr","TE", "type","start","stop","pos1","strand","pos2","descriptor")

#remove outfile to reduce memory burden
rm(outfile)

#create a TE only GTF file
TE_gtf <- gdata::interleave(gene, transcript, exon)

#remove temporary files to reduce memory burden
rm(gene)
rm(transcript)
rm(exon)

#import genic GTF file and create the same headers
genic <- read.table("genic_gtf_file.gtf", sep = "\t", quote = "")
colnames(genic) <- c("chr","TE", "type","start","stop","pos1","strand","pos2","descriptor")

#combine gtfs
complete_gtf <- rbind(genic, gtf)

#write tables for a TE-specific GTF and a whole genome GTF
write.table(TE_gtf, file = "TE_only.gtf", sep = "\t", col.names = F, row.names = F, quote = F)

write.table(complete_gtf, file = "Gene_and_TE.gtf", sep = "\t", col.names = F, row.names = F, quote = F)