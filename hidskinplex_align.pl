#!/usr/bin/env perl

########################################
#Created by Sarah Schmedes
#Release date: 9/1/2017
#######################################

use strict;
use warnings;
use Genomix qw(:constants :meta);
use Getopt::Long;

#select samples to run through metaphlan2
my $QCdatapath;
my $outdir;
my $samples;
#Set flags
GetOptions('samplelist=s' => \$samples,
	   'QCdatapath=s' => \$QCdatapath,
	   'outdir=s' => \$outdir)
    or die "Must declare samplie flag to run script\n";
open IN, "<$samples" or die "Cannot open samplelist for reading $!\n";

while (<IN>) {
    my $file1 = $_;
    chomp $file1;
    my $file2 = <IN>;
    chomp $file2;
    my @filename = split /\./, $file1;
    my @samplename = split /_/, $filename[0];
    my $samplename = $samplename[0];
    
    my $filesin = "$QCdatapath/$file1 $QCdatapath/$file2"; 
    
#run metaphlan2 on all samples
#all metaphlan2/strainphlan scripts must be in the environmental path
    if (! -f "$outdir/$samplename\_profile.txt") {
	print "Metaphlan running on $samplename\n";
	if (system("bash", "-c", "metaphlan2.py --input_type fastq <(zcat $filesin) --mpa_pkl /home/sarah/strainphlan/db_v20/mpa_v20_m200.pkl --output_file $outdir/$samplename\_profile.txt --bowtie2out $outdir/$samplename\_bowtie2.txt --samout $outdir/$samplename\.sam.bz2 --nproc 20"))
	{
	    die "Metaphlan2 ERROR: $!";
	}
    }    
}









