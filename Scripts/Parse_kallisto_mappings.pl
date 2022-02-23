#!/usr/bin/perl

use warnings;
use strict;

# 2021-09-17
# parse kallisto maps into one file

unless (@ARGV >= 2) {
	die "parse_kallisto_mappings.pl <outfile name> <kallisto 1> <kallisto 2> ....";
}

my $outname = shift @ARGV;

# to get gene length on the first file
my $start = 1;

my %outresults = ();

# prepare header
my $header = "Gene\tlength\teff_length\t";
$header .= join("\t",@ARGV);
$header .= "\n";

# loop on files
foreach my $file (@ARGV) {
	
	open(FILE,$file) or die "Could not open $file: $!\n";
	
	# skip first line
	my $headerline = <FILE>;
	
	if ($start == 1) {
		
		while (my $line = <FILE>) {
			### target_id	length	eff_length	est_counts	tpm
			### lcl|CM025008.1_mrna_1	5234	5235	0	0
			### NotFur1_Unknown_966#DNA/DNA	907	908	0	0

			my @linedata = get_line_data($line);
			#print "@linedata\n:";

			push(@{$outresults{$linedata[0]}}, ($linedata[1],$linedata[2],$linedata[3]) );
		}
		
		$start = 0;
	
		
	} else {
	
		while (my $line = <FILE>) {
			my @linedata = get_line_data($line);
			push(@{$outresults{$linedata[0]}}, $linedata[3] );
		}
		
	}
	
	close FILE;
}



#### output results
open(OUT,'>',$outname) or die "Could not open $outname: $!\n";

print OUT $header;

foreach my $gene (sort keys %outresults) {
	
	#print $gene."\n";
	
	# ENSMUSG00000000094|ENSMUST00000000096|Tbx4|protein_coding|11
	# "GeneSymbol\tGene_Type\tlength\teff_length\t"

	my $outline = $gene."\t".join("\t",@{$outresults{$gene}})."\n";
	print OUT $outline;
}


close OUT;

exit;

###########################################################
# SUBROUTINES
###########################################################

###########################################################
# a subroutine that separates fields from a data line and
# returns them in an array

sub get_line_data {

    my $line = $_[0];
    
    chomp $line;  

    my @linedata = split(/\t/, $line);
        
    return @linedata;
}