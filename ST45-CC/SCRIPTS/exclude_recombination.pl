#!/usr/bin/perl
use strict;

#author: Zhang Ji
#email: ji.zhang@helsinki.fi
#last update: 30.7.2014
#function: get read of recombination sites predicted using BNG and print new alignment in *_screened.fas
#usage: perl exclude_recombination.pl BNG_tabular_format.txt multi-FASTA.fas


my $pozfile = $ARGV[0];
my $seqfile = $ARGV[1];

my ($line_counter, $inline, $start, $end, $isolate1, $isolate2, $iso_counter);
my (@line, @isolates);
open(IN, "<$pozfile");

$line_counter = 0;
$iso_counter = 0;
while(<IN>){
	chomp;
	$line_counter++;
	if($line_counter<3){
		next;
	}
	else{
		$inline = $_;
		$inline =~ s/  / /g;
		@line = split(" ",$inline);
		$start = $line[0];
		$end = $line[1];
		$isolate1 = $line[-1];
		if("$isolate1" eq "$isolate2"){
			print POZ "$start\n", "$end\n";
		}
		else{
			$iso_counter++;
			open (POZ, ">>$iso_counter.poz.tmp");
			print POZ "$start\n", "$end\n";
			push @isolates, $isolate1;
		}
		$isolate2 = $isolate1;
	}
}
close IN;

my($counter, $isolate, $seq);
my(@seqs);
my(%isolates_hash);
open(IN, "<$seqfile");
$/=">";
$counter = 0;
@seqs = <IN>;
foreach(@seqs){
	if("$counter"==0){
		$counter++;
	}
	else{
		$isolate = $_;
		$isolate =~ s/\n.*//g;
		$isolate =~ s/\n//g;
		$seq = $_;
		$seq =~ s/>//;
		$seq =~ s/.*\n//;
		$seq =~ s/\n//g;
		open (OUT, ">$counter.seq.tmp");
		print OUT "\n>$isolate\n", "$seq\n";
		$isolates_hash{$isolate} = $counter;
		close OUT;
		$counter++;
	}
}
$/="\n";
close IN;

my ($poz1, $poz2, $seqname, $seq, $full_len, $poz_counter, $start, $end, $len, $extract, $roller);
my (@poz);
$counter = 0;
foreach(@isolates){
	$counter++;
	$isolate = $_;
	$seqname = $isolates_hash{$isolate};
	open (SEQ, "<$seqname.seq.tmp");
	open (POZ, "<$counter.poz.tmp");
	open (OUT, ">>$counter.extrction.seq.tmp");
	
	while(<SEQ>){
		chomp;
		if($_ =~ />/){
			print OUT "\n$_\n";
		}
		else{
			$seq = $_;
			$full_len = length($seq) + 1;
		}
	}
	@poz = <POZ>;
	push @poz, $full_len;
	unshift @poz, "0";
	
	$poz_counter = 0;
	$roller = 1;
	foreach(@poz){
		chomp;
		$poz_counter++;
		if($roller == 0){
			$start = $_ - 1;
			$end = $poz["$poz_counter"];
			$end =~ s/\n//;
			$len = $end - $start;
			$extract = substr $seq, $start, $len;
			$extract =~ s/./-/g;
			$roller = 1;
		}
		else{
			$start = $_;
			$end = $poz["$poz_counter"];
			$end =~ s/\n//;
			$len = $end - $start - 1;
			if($len<0){
				$extract = ();
				print "OK\n";
			}
			else{
			$extract = substr $seq, $start, $len;
			}
			$roller = 0;
		}
		print OUT "$extract";
	}
	system ("rm -f $seqname.seq.tmp");
}
close SEQ;
close POZ;
close OUT;

system ("cat *.seq.tmp >$seqfile.screened.fas");
system ("rm -f *.tmp");
