#!/usr/bin/env perl

############################################
#Authored by: Sarah E. Schmedes
#This script is used to quality trim and filter reads for hidskinplex targeted sequencing data
#Release Date: 9/1/17
#############################################
use strict;
use warnings;
use Getopt::Long;

my $samples;
my $datapath;
my $outpath;
my $fastqc; #path to fastqc output
#Set flags
GetOptions('samplelist=s' => \$samples,
	   'datapath=s' => \$datapath,
	   'outpath=s' => \$outpath,
	   'fastqc=s'=> \$fastqc)
    or die "Bodysite flag is required for run\n";
my $outpath;

open IN, "<$samples" or die "Could not open sample list for reading $!\n";

while (<IN>) {
    my $file1 = $_;
    chomp $file1;
    my $file2 = <IN>;
    chomp $file2;
    my @samplename = split /\./, $file1;
    my $sample = $samplename[0];
    #cutadapt parameters
    my $qual = 20; #can change quality score threshold
    my $length = 50; #can change length filter cutoff

    if (! -f "$outpath/$sample\_cutadapt_report.txt") {
	print "Cutadapt running on $file1\n";
	if (system ("cutadapt -n 2 --trim-n -q $qual -m $length -o $outpath/$file1 -p $outpath/$file2 $datapath/$file1 $datapath/$file2 > $outpath/$sample\_cutadapt_report.txt"))
	{
	    die "Cutadapt failed\nError is: $!\n";
	}
	if (system("fastqc -o $fastqc -t 20 -f fastq $outpath/$file1 $outpath/$file2")) {
	    die "Error running fastqc $!\n";
	}
    }
}
close(IN);
