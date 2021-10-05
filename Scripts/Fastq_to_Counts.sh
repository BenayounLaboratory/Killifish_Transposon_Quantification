#!/usr/bin/bash

## trimmer

. Fastq_to_Counts_Filled_Arguments.sh

fastxdir=$FastqDir/"fastxtrimmed"

mkdir $fastxdir

cd $FastqDir

for f in $(find "." -name '*.fastq.gz')
do
        f2=$(basename "${f}");
        of=$(basename "${f}" | sed 's/\.fastq\.gz/_hardtrim\.fq\.gz/g');

        gzcat $f2 | $FastXDir -f $Front -l $Length -z -i - -Q33 -o $fastxdir/$of
        
        echo "Done hard-trimming $f2"
        
done

echo "Done with hard trimming!"

#trim_galore

cd $fastxdir

trimgalout=$fastxdir/"quality"

mkdir $trimgalout


if [ $End = "P" ] 
then
    for f in $(find "." -name '*1_hardtrim.fq.gz')
	do
        f1=$(basename "${f}");
        f2=$(basename "${f}" | sed 's/1_hardtrim\.fq\.gz/2_hardtrim\.fq\.gz/g');
		$TrimGalDir --paired $fastxdir/$f1 $fastxdir/$f2 -o $trimgalout  
	done
else
echo "I'm in single end mode"
        for f in $(find "." -name '*1_hardtrim.fq.gz')
	do
        f1=$(basename "${f}");
		$TrimGalDir $fastxdir/$f1 -o $trimgalout 
	done		
fi

echo "Done trimming adapters!"


#run kallisto

if [ $Kal_run = "Y" ] 
then

	#indexing
	
	echo "Indexing the transcriptome"
	
	cd $trimgalout
	
	$kallisto index -i $trimgalout/transcriptome_index.idx $trans
	
	echo "Done Indexing!"
	
	#counting
	
	if [ $End = "P" ] 
	then
	    for f in $(find "." -name '*1_hardtrim_val_1.fq.gz')
		do
	        f1=$(basename "${f}");
	        of=$(basename "${f}" | sed 's/1_hardtrim_val_1\.fq\.gz/counts/g');
	        f2=$(basename "${f}" | sed 's/1_hardtrim_val_1\.fq\.gz/2_hardtrim_val_2\.fq\.gz/g');
	
			kallisto quant -i transcriptome_index.idx $f1 $f2 -o $of --plaintext
	        
		done
		
		echo "Done Generating Kallisto Counts! Please concatenate counts using the perl script: collapse_perl.pl"
	else
	    for f in $(find "." -name '*_1_hardtrim_trimmed.fq.gz')
	do
	        f1=$(basename "${f}");
	        of=$(basename "${f}" | sed 's/_1_hardtrim_trimmed\.fq\.gz/_counts/g');
	
			kallisto quant -i transcriptome_index.idx --single -l $trim_len -s $sd $f1 -o $of --plaintext
	done
	
		echo "Done Generating Kallisto Counts!"
	
	fi
	
	else
	
	echo "Goodbye!"
fi	
	
if [ $STAR_run = "Y" ] 
then
	
	#indexing the genome
	
	echo "Indexing the Genome"
	
	cd $trimgalout
	
	$Path_to_STAR --runThreadN 6 --limitGenomeGenerateRAM 60000000000 --runMode genomeGenerate --genomeDir $Indexed_Genome_Dir --genomeFastaFiles $genome_fasta --sjdbGTFfile $genomic_gtf
	
	#counting
	
		if [ $End = "P" ] 
	then
	    for f in $(find "." -name '*_1_hardtrim_val_1.fq.gz')
			do
    			f2=$(basename "${f}" | sed 's/_1_hardtrim_val_1\.fq\.gz/_2_hardtrim_val_2\.fq\.gz/g');
    			inf2="${trimgalout}/${f2}"
    			of=$(basename "${f}" | sed 's/_1_hardtrim_val_1\.fq\.gz//g');
    			oFname="${Genomic_output_directory}/${of}"
    			STAR --genomeDir $Indexed_Genome_Dir --readFilesIn $f $inf2 --readFilesCommand gzcat --runThreadN 2 --outFilterMultimapNmax 200 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --alignEndsProtrude 10 ConcordantPair --limitGenomeGenerateRAM 60000000000 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $oFname
			done
			
			echo "Done Mapping Reads to Genome with STAR"
	else
		for f in $(find "." -name '*_1_hardtrim_trimmed.fq.gz')
			do
		    	of=$(basename "${f}" | sed 's/_1_hardtrim_trimmed\.fq\.gz//g');
		    	oFname="${Genomic_output_directory}/${of}"
		    	STAR --genomeDir $Indexed_Genome_Dir --readFilesIn $f --readFilesCommand gzcat --runThreadN 2 --outFilterMultimapNmax 200 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --alignEndsProtrude 10 ConcordantPair --limitGenomeGenerateRAM 60000000000 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $oFname
			done
			
			echo "Done Mapping Reads to Genome with STAR"
	fi
	
	else
	
	echo "Goodbye!"
	
fi

if [ $FtCts = "Y" ] 
then

cd $Genomic_output_directory

bam_names=$(ls -p | grep "Aligned.sortedByCoord.out.bam" | tr '\n' ' ')

	if [ $End = "P" ] 
	then
		if [ $Fraction = "Y" ]
		then
		$FtCts_Path -O -p -M --fraction -a $genomic_gtf -o $FtCts_output $bam_names
		else
		$FtCts_Path -O -p -a $genomic_gtf -o $FtCts_output $bam_names
		fi
	else
		if [ $Fraction = "Y" ]
		then
		$FtCts_Path -O -M --fraction -a $genomic_gtf -o $FtCts_output $bam_names
		else
		$FtCts_Path -O -a $genomic_gtf -o $FtCts_output $bam_names
		fi
	fi
else
echo "Done Without Count Matrix Generation"
fi