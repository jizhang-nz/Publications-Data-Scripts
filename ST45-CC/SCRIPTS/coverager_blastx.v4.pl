#!/usr/bin/perl
#coverager: caculate the percentage of an aa sequence can be covered by blastx alignment of a nucleotide sequence

#Copy (C) 2014-2015  Helsinki University.  All rights reserved
#Written by Ji Zhang, MD, PhD

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY. See the GNU General Public License for 
#more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#Revision notes
#Created 13.1.2015

#Options: 
#-r list of the reference files (*.faa).
#-g list of the sequence files (*.fna).
#-h: Show help
#-v: Show progress
#-i: identity threshold (default 0)
#-e: evalue (default 0.0001)

use strict;
use Getopt::Std;
#use SVG;
use List::Util qw/sum/;

getopts('g:r:d:i:e:h v');
our ($opt_g, $opt_r, $opt_d, $opt_i, $opt_e, $opt_h, $opt_v);
my $distance_d = $opt_d || "20000";
my $identidy_i = $opt_i || "0";
my $evalue_e = $opt_e || "0.0001";
###################
if($opt_h==1){
print ("
coverager: caculate the percentage of an aa sequence can be covered by blastx alignment of a nucleotide sequence.

Options: 
-r list of the reference files (*.faa).
-g list of the sequence files (*.fna).
-h: Show help
-v: Show progress
-i: identity threshold (pident >= 70)

");
exit;
}
print ("Making preparations for the analysis. It is going to take a while ...") if ($opt_v == 1);
###################
my($counter);
open(GENOME, "<$opt_g") or die "Cannot open genome list file!";
open(OUT, ">genome_list.tmp");
my $counter = 0;
while(<GENOME>){
	if($_ =~ /^\s/){
		next;
	}
	else{
		$counter++;
		chomp;
		print OUT "$_\n";
	}
}
close GENOME;
close OUT;

my($counter);
open(REF, "<$opt_r") or die "Cannot open the list of the references";
open(OUT, ">ref_list.tmp");
my $counter = 0;
while(<REF>){
	if($_ =~ /^\s/){
		next;
	}
	else{
		$counter++;
		chomp;
		print OUT "$_\n";
	}
}
close REF;
close OUT;
###################combine the contigs
my($genome, $seq, $counter, $spacer);
open(GENOME, "<genome_list.tmp");
while(<GENOME>){
	chomp;
	$genome = $_;
	open(SEQIN, "<$genome") or die "Cannot open genome sequence $genome!";
	open(OUT, ">>$genome.combined.fas");
	$/=">";
	print OUT ">", "$genome\n";
	while(<SEQIN>){
		$spacer = ();
		$counter = 0;
		while($counter<$distance_d){
			$spacer = "$spacer"."N";
			$counter++;
		}
		
		if($_ =~ /^>/){
			next;
		}
		else{
			my $seq = $_;
			$seq =~ s/.*\n//;
			$seq =~ s/>//;
			$seq =~ s/\n//g;
			print OUT "$seq", "\n$spacer\n";
		}
	}
	close OUT;
	close SEQIN;
	$/="\n";
}
close GENOME;
print ("... done!\n") if ($opt_v == 1);
###################
my($reference, $counter, $gene_total, $aaseq, $aalength, $sn);
my(@percent_all);

open(REFERENCE, "<ref_list.tmp");

while(<REFERENCE>){
	chomp;
	$reference = $_;
	$counter = 0;
	@percent_all = ();
	
	open(RESULTS, ">>output.$reference.txt");
	
	open(REFIN, "<$reference");
	open(TMP, ">$reference.tmp");
	while(<REFIN>){
		chomp;
		if($_ =~ />/){
			print TMP "$_   \n";
		}
		else{
			print TMP "$_\n";
		}
	}
	close TMP;
	close REFIN;
	
	open(REFIN, "<$reference.tmp");
	open(TMP, ">$reference.tmp.tmp");
	$/= undef;
	while(<REFIN>){
		$_ =~ s/\n//gs;
		$_ =~ s/   /\n/gs;
		$_ =~ s/>/\n>/gs;
		$_ =~ s/^\n//gs;
		print TMP "$_";
	}
	close TMP;
	close REFIN;
	$/="\n";
	
	open(REFIN, "<$reference.tmp.tmp");
	open(REFOUT, ">$reference.aa.tmp");
	while(<REFIN>){
		chomp;
		if($_ =~ /^>/){
			$percent_all["$counter"] = 0;
			print REFOUT ">gene$counter";
			$counter++;
			$aaseq = ();
		}
		else{
			$aaseq = $_;
			$aalength = length($aaseq);
			print REFOUT "_$aalength\n$aaseq\n";
		}
	}
	close REFIN;
	close REFOUT;
	$gene_total = $counter + 1;
	
	print RESULTS "g";
	$sn = 1;
	while($sn < $gene_total){
		print RESULTS "\t$sn";
		$sn++;
	}
	print RESULTS "\n";
	
	system("makeblastdb -in $reference.aa.tmp -dbtype prot -logfile makeblastdb.log");
	system("rm -f makeblastdb.log");
	###################
	my($genome, $aa_len, $aa_len2, $gene_number, $gene_number2, $s_start, $s_end, $counter, $sum, 
			$percentage, $line, $gene_counter, $pident, $sseqid, $sseqid2, $sum_percent, $coverage);
	my ($svg, $x);
	my (@blastx, @inline, @aa, @percent, @percent_all2);
	open(GENOME, "<genome_list.tmp");
	while(<GENOME>){
		chomp;
		$genome = $_;
		@aa=(); @percent=();
		@percent_all2 = @percent_all;
		$sseqid2 = "NA";
		open(SEQIN, "<$genome");
		print ("blastx: ", "$genome", " vs. ", "$reference", " ...\n") if ($opt_v == 1);
		@blastx = readpipe("blastx -query $genome.combined.fas -db $reference.aa.tmp -max_target_seqs 99999 -evalue $evalue_e -seg no -num_threads 4 -query_gencode 11 -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qframe sframe\"");
		if(@blastx){
			$line=0;
			$gene_counter = 0;
			
			foreach(@blastx){
				chomp;
				$line++;
				@inline = split (/\t/,$_);
				$sseqid = $inline[1];
				$aa_len = $inline[1];
				$aa_len =~ m/_[0-9]+/;
				$aa_len = $&;
				$aa_len =~ s/_//;
				$gene_number = $inline[1];
				$gene_number =~ s/_.*//g;
				$gene_number =~ s/gene//;
				$pident = $inline[2];
				$s_start = $inline[8];
				$s_end = $inline[9];
				
				if($line == 1){
					$sseqid2 = $sseqid;
				}
				else{};
					
				if($pident >= $identidy_i){
					if($sseqid2 eq $sseqid){
						$gene_number2 = $gene_number;
						$counter = $s_start - 1;
						$aa_len2 = $aa_len;
						while($counter<$s_end){
							$aa["$counter"] = 1;
							$counter++;
						}
					}
					else{
						$sum = sum @aa;
						if($sum==0){
							$percentage = 0;
						}
						else{
							$percentage = $sum/$aa_len2;
						}
						push @percent, $percentage;
						$percent_all2["$gene_number2"] = $percentage;
						@aa = ();
						
						$gene_number2 = $gene_number;
						$counter = $s_start - 1;
						$aa_len2 = $aa_len;
						while($counter<$s_end){
							$aa["$counter"] = 1;
							$counter++;
						}
						$sseqid2 = $sseqid;
					}
				}
				else{}
			}
			$sum = sum @aa;
			if($sum==0){
				$percentage = 0;
			}
			else{
				$percentage = $sum/$aa_len2;
			}
			push @percent, $percentage;
			$percent_all2["$gene_number2"] = $percentage;
			
			$sum_percent = sum @percent;
			if($sum_percent==0){
				$coverage = 0;
			}
			else{
				$coverage = $sum_percent/$gene_total;
			}
#			print RESULTS "$genome\t$reference\t$coverage\t@percent_all2\n";
			print RESULTS "$genome";
			foreach(@percent_all2){
				print RESULTS "\t$_";
			}
			print RESULTS "\n";
		}
		else{
			print RESULTS "$genome\t$reference\tM\n";
		}
	}
}
close REFERENCE;
###################
system ("rm -f *.combined.fas");
system ("rm -f *.tmp");
system ("rm -f *.*hr");
system ("rm -f *.*in");
system ("rm -f *.*sq");

#system ("mkdir tmp");
#system ("mkdir combined_contigs");
#system ("mv *.combined.fas ./combined_contigs/");
#system ("mv *.tmp ./tmp/");
