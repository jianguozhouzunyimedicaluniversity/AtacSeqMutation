#!/usr/bin/perl
# Small script that takes in a GTF file of TSS (e.g. ENSEMBL) and will output a list of TSS in BED format
# Starting position is zero based, end position is one based
# Created 06/29/12 by Ewy Mathe
# Edited 11/02/15 by Elizabeth Baskin
# Edited 3/15/18 by Tara Eicher 
use strict;
use warnings;
use Pod::Usage;

pod2usage("Usage: $0 ensembl.gtf tss.bed") if ((@ARGV != 2));

#Open files.
open(my $ref, "<", $ARGV[0]) or die "Can't open $ARGV[0] $!!";
open(my $out, ">", $ARGV[1]) or die "Can't open $ARGV[1] $!";

#This is a hash that keeps track of TSS.
#We will be making a hash of arrays below.
#The hash will consist of the gene ID followed by the transcript ID.
my %tss; 

#Remove the header.
<$ref>;
<$ref>;
<$ref>;
<$ref>;
<$ref>;

#Loop through each line in the GTF file. If the line is an exon with exon number 1, add it to the hash.
#Once a line is in the hash, print it.
while(<$ref>) {
	#Remove new lines and other end of line characters. 
	chomp($_);
	
	#Seperate the fields by a comma.
	my @fields = split(/,/,$_);
	
	#Place fields in variables for easier readability.
	my $category = $fields[2];
    my $exonnumber = $fields[14];
    my $strand = $fields[6];
	my $chrom = $fields[0];
    my $start = $fields[3];
    my $stop = $fields[4];
    my $gene = $fields[8];
    my $transc = $fields[9];
    my $genename = $fields[11];
    my $transcname = $fields[13];
	#If the strand is equal to a plus sign (positive strand), has category exon, and has exon number 1, add it to the hash and print.
	if ($strand eq "-" && $category eq "exon" && $exonnumber eq "exon_number1") {
		my @tssarray = ($chrom, $start, $strand, $gene, $transc, $genename, $transcname, $category, $exonnumber);
	
		#If this transcript already exists in the hash, skip. 
		if (exists $tss{($gene."_".$transc)}) {
			next;
		}
		#Otherwise, this transcript does not exist in the hash. We add it and print the info we want to the output.
		else {
			$tss{($gene."_".$transc)}=\@tssarray;
			print $out $chrom."\t".$start."\t".($start+1)."\t".$category."\t".$exonnumber."\t".$gene."\t".$transc."\t".$genename."\t".$transcname."\t".$strand."\n";
		}
	}
	##If the strand is equal to a negative sign (negative strand), has category exon, and exon number 1, add it to the hash and print.
	elsif ($strand eq "+" && $category eq "exon" && $exonnumber eq "exon_number1") {
		my @tssarray = ($chrom, $start, $strand, $gene, $transc, $genename, $transcname, $category, $exonnumber);

		#If this transcript already exists in the hash, skip. 
		if (exists $tss{($gene."_".$transc)}) {
			next;
		}
		#Otherwise, this transcript does not exist in the hash. We add it and print the info we want to the output.
		else {
			$tss{($gene."_".$transc)}=\@tssarray;
			print $out $chrom."\t".($stop-1)."\t".($stop)."\t".$category."\t".$exonnumber."\t".$gene."\t".$transc."\t".$genename."\t".$transcname."\t".$strand."\n";
		}
	}

	#We are not interested in non-exonic regions or exon numbers apart from 1.
	#Skip these lines.
	else { next;}

}
    
#Close the files and exit the program.
close($ref) or die "Can't close $ref\n";
close($out) or die "Can't close $out\n";
exit;
