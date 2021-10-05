README - Transposable Element Differential Expression Analysis Protocol from RNA-seq Data
###########################

**************************
*** R Package versions ***

DESeq2 v1.30.1 
Rcolorbrewer v1.1-2
pheatmap v1.0.12

**************************
***   Other software   ***

fastx_toolkit v0.0.13 (fastx_trimmer)
trim_galore v0.6.7
kallisto v0.46.1
cutadapt v3.4
FastQC
R v4.0.3
RepeatMasker v4.1.2
SubRead v2.0.3  (featurecounts)
STAR v2.7.0e

**************************

* Bash_Scripts

 - FastXTrimmer.sh: Hardtrim reads with fastx_trimmer
 - TrimGalore_SE.sh: Remove adapters using trim_galore for single-end reads
 - TrimGalore_PE.sh: Remove adapters using trim_galore for paired-end reads
 - Kallisto_Index.sh: Create a kallisto index reference from a transcriptome input
 - Kallisto_Quant_SE.sh: Run kallisto quant to pseudoalign reads and create counts for single-end data
 - Kallisto_Quant_PE.sh: Run kallisto quant to pseudoalign reads and create counts for paired-end data
 - RepeatMasker_FishTEDB.sh: Run repeatmasker using RMBlast to create a FishTEDB masked genomic reference
 - STAR_Index.sh: Index genome for STAR mapping
 - STAR_Map_SE.sh: Map reads to STAR genomic reference allowing for up to 200 multimappers using single-end data
 - STAR_Map_PE.sh: Map reads to STAR genomic reference allowing for up to 200 multimappers using paired-end data
 - featureCounts_SE_Fractional.sh: Generate count matrix using featureCounts using fractional counting for multimappers on single-end data
 - featureCounts_SE_Unique.sh: Generate count matrix using featureCounts without allowing multimappers on single-end data
 - featureCounts_PE_Fractional.sh: Generate count matrix using featureCounts using fractional counting for multimappers on paired-end data
 - featureCounts_PE_Unique.sh: Generate count matrix using featureCounts without allowing multimappers on paired-end data
 - Fastq_to_Counts.sh: Takes zipped fastq files as input and 1) hardtrims 2) removes adapters and either i) uses Kallisto to generate counts or ii) uses STAR and featureCounts to generate counts. Takes arguments from "Fastq_to_Counts_Filled_Arguments.sh", which is the completed result of "Fastq_to_Counts_Arguments.sh"
 - Fastq_to_Counts_Arguments.sh: List of prompts to be filled in by user as instructed. To be saved as "Fastq_to_Counts_Filled_Arguments.sh" and used as input for "Fastq_to_Counts.sh"

* Perl_Scripts

 - collapse_perl.pl: Concatenate count results from Kallisto into a single count matrix
 
* R_Scripts

 - Import_Libraries.R: Import necessary libraries into R for differential expression analysis and heatmap creation
 - Format_Kallisto_Matrix.R: Format count matrix from Kallisto for analysis
 - Format_featureCounts_Matrix.R: Format count matrix from featureCounts for analysis
 - Rounded_Matrix.R: Round count matrix for use in DESeq2 pipeline
 - Data_Design.R: Create design matrix for DESeq2
 - Run_DESeq2.R: Run DESeq2 differential expression analysis
 - Write_Outfiles.R: Write results files to text files
 - VST_Transform_Data.R: Transform data for heatmap plotting
 - Isolate_TE_Data.R: Keep only TE data for heatmap visualization
 - Plot_Heatmap.R: Generate TE heatmap
 - Complete_DESeq2_Kallisto.R: Takes Kalliso count matrix input after "collapse_perl.pl" concantenation and run full DESeq2 analysis with heatmap creation
 - Complete_DESeq2_featureCounts.R: Takes featureCounts count matrix and runs full DESeq2 analysis with heatmap creation