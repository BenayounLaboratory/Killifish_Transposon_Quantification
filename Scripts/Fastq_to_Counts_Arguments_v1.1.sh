# This file is a list of arguments to be passed to Fastq_to_Counts_v1.1.sh
# Save as Fastq_to_Counts_Filled_Arguments.sh
#Enter within quotes if quotes are present
####################################################

#name the directory where the fastq.gz files are;
FastqDir=""

#is your data is paired-end, please enter "Paired". Otherwise assumed to be single-end data
End=""

#name the path to fastxtrimmer, specifying the program i.e. xxx/bin/fastx_trimmer
FastXDir=""

#name the path to trimgalore i.e. xxx/bin/trimgalore
TrimGalDir=""

#number of nucleotides to remove from the 5' end of the read
Front=

#length of the reads
Length=

#path to transcriptome i.e. xxx/xxx/transcriptome.fa
trans=""

#path to kallisto
kallisto=""

#for single end only- what is the expected trimmed length
trim_len=

#for single end only- what is the expected st. dev in trimmed length
sd=

#do you want to run kallisto? Y for yes, N for no
Kal_run=""

#do you want to run STAR mapping? Y for yes, N for no
STAR_run=""

#STAR path
Path_to_STAR=""

#If using STAR, specify the path to the directory where you want the indexed genome reference stored
Indexed_Genome_Dir=""

#If using STAR, specify the path to the genome_TE genomic reference
genome_fasta=""

#If using STAR, specify the path to the genome_TE GTF
genomic_gtf=""

#genomic output directory
Genomic_output_directory=""

#Run featurecounts? "Y" for yes
FtCts=""

#path to featurecounts
FtCts_Path=""

#name of featurecounts output
FtCts_output=""

#For featurecounts, count fractionally? "Y" for yes
Fraction=""