######################################
###########     README     ###########
######################################


Bryan Teefy, Matthew Malone and Bérénice A. Benayoun


Scripts and reference files to analyze transcription from transposons 
in the African turquoise killifish genome

Also see companion figshare repository for reference files: https://doi.org/10.6084/m9.figshare.19394723


**************************
*** R Package versions ***

DESeq2 v1.30.1 
Rcolorbrewer v1.1-2
pheatmap v1.0.12

**************************
***   Other software   ***

fastx_toolkit (0.0.13, fastx_trimmer)
kallisto (0.46.2)
cutadapt (3.3, dependency for TrimGalore)
R (3.5.1)
RepeatMasker (4.1.2-p1)
Subread (2.0.2; for the featureCounts function)
STAR (2.7.0e)
trim_galore (0.6.7, TrimGalore)

**************************

* Bash_Scripts

 - FastXTrimmer.sh: Hardtrim reads with fastx_trimmer
 - Fastq_to_Counts_v1.1.sh: Takes zipped fastq files as input and 1) hardtrims 2) removes adapters and either i) uses Kallisto to generate counts or ii) uses STAR and featureCounts to generate counts. Takes arguments from "Fastq_to_Counts_Filled_Arguments.sh", which is the completed result of "Fastq_to_Counts_Arguments.sh"
 - Fastq_to_Counts_Arguments_v1.1.sh: List of prompts to be filled in by user as instructed. To be saved as "Fastq_to_Counts_Filled_Arguments.sh" and used as input for "Fastq_to_Counts_v1.1.sh"
 - featureCounts_SE_Fractional_v1.1.sh: Generate count matrix using featureCounts using fractional counting for multimappers on single-end data
 - featureCounts_SE_Unique_v1.1.sh: Generate count matrix using featureCounts without allowing multimappers on single-end data
 - featureCounts_PE_Fractional_v1.1.sh: Generate count matrix using featureCounts using fractional counting for multimappers on paired-end data
 - featureCounts_PE_Unique_v1.1.sh: Generate count matrix using featureCounts without allowing multimappers on paired-end data
 - Kallisto_Index_v1.1.sh: Create a kallisto index reference from a transcriptome input
 - Kallisto_Quant_SE.sh: Run kallisto quant to pseudoalign reads and create counts for single-end data
 - Kallisto_Quant_PE.sh: Run kallisto quant to pseudoalign reads and create counts for paired-end data
 - perl_collapse_v1.1.sh: Run collapse_perl.pl
 - RepeatMasker_FishTEDB_v1.1.sh: Run repeatmasker using RMBlast to create a FishTEDB softmasked genomic reference
 - STAR_Index.sh: Index genome for STAR mapping
 - STAR_Map_SE_v1.1.sh: Map reads to STAR genomic reference allowing for up to 200 multimappers using single-end data
 - STAR_Map_PE_v1.1.sh: Map reads to STAR genomic reference allowing for up to 200 multimappers using paired-end data
 - TrimGalore_SE.sh: Remove adapters using trim_galore for single-end reads
 - TrimGalore_PE.sh: Remove adapters using trim_galore for paired-end reads
 
* Perl_Scripts

 - collapse_perl.pl: Concatenate count results from Kallisto into a single count matrix
 
* R_Scripts

 - Complete_DESeq2_Kallisto_v1.1.R: Takes Kalliso count matrix input after "collapse_perl.pl" concantenation and run full DESeq2 analysis with heatmap creation
 - Complete_DESeq2_featureCounts_v1.1.R: Takes featureCounts count matrix and runs full DESeq2 analysis with heatmap creation
 - Data_Design.R: Create design matrix for DESeq2
 - Format_Kallisto_Matrix.R: Format count matrix from Kallisto for analysis
 - Format_featureCounts_Matrix.R: Format count matrix from featureCounts for analysis
 - GTF_Generation_v1.1.R: Create GTF file of genes and TEs from the Repeatmasker .out file and the genome-associated genic GTF file
 - Import_Libraries.R: Import necessary libraries into R for differential expression analysis and heatmap creation
 - Isolate_TE_Data.R: Keep only TE data for heatmap visualization
 - Plot_Heatmap.R: Generate TE heatmap
 - Rounded_Matrix.R: Round count matrix for use in DESeq2 pipeline
 - Run_DESeq2.R: Run DESeq2 differential expression analysis
 - VST_Transform_Data.R: Transform data for heatmap plotting
 - Write_Outfiles.R: Write results files to text files
