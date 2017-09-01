#!/usr/bin/env perl

############################################################
#Created by Sarah Schmedes
#Script will pull out reads aligning to markers in targeted panel across samples in a bodysite
#and create mpileup files, identify variants, and calculate theta pi across markers
#(parse through an .mpileup file (run with -Os flags))
############################################################

use strict;
use warnings;
use Getopt::Long;
use Genomix qw(:constants :meta);

my ($bodysite, $help);
my $samfilepath = "/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results"; #file path to .sam.bz2 files
my $pileuppath = "/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/mpileup_results";
my (%markerND, %markerNDlength, %finalmarkers);
#setflags
GetOptions('bodysite=s' => \$bodysite, #group of SRSs to run for bodysite
	   'help' => \$help)
    or die "Error with flags. Designatre -bodysite or see -help for more information\n";
if ($help) {
    die "Must designate bodysite for targeted panel analysis.\n";
} 

#Get desired samples to analyze
open SAMPLES, "</mnt/blsm/sarah/targetedstudy/miseqdata/fastqs/fastqlist.txt" or die "Could not open sample list for reading! $!\n";
my %samplenames;
while (<SAMPLES>) {
    my $line = $_;
    chomp $line;
    my $line2 = <SAMPLES>;
    chomp $line2;
    my @partialname = split /\./, $line;
    my @samplename = split /_/, $partialname[0];
    if ($samplename[0] =~ /^SwabBlank/) {
	next;
    }
    $samplenames{$samplename[0]} = 1;
}
close(SAMPLES);
my @sortednames = sort {$a cmp $b} keys %samplenames;
my %bodysamplename;
my @bodysortednames;
foreach my $sample (@sortednames) {
    my ($sampleid, $body, $time, $rep) = split /-/, $sample;
    if ($body eq $bodysite) {
	$bodysamplename{$sample} = 1;
    }
}
if ($bodysite eq "all") {
    @bodysortednames = @sortednames;
} else {
    @bodysortednames = sort {$a cmp $b} keys %bodysamplename;
}
#Get taxon path and length for each marker
my ($taxamarker_ref, $taxalevel_ref, $markerlength_ref) = gettaxa();
my %taxamarker = %$taxamarker_ref; #hash of each marker and associated full taxon path
my %taxalevel = %$taxalevel_ref; #hash of each full taxon path and associated
                                 #leaf clade tab separated corresponding taxon level
my %markerlength = %$markerlength_ref; #hash of each marker and associated length
my @totalmarkers = sort keys %taxamarker;

#Run samtools-1.3.1 mpileup on designated samples
my @files;
foreach my $file (@bodysortednames) {
    push(@files, "$samfilepath/$file\.sorted.bam");
}
my $filestring = join("\t", @files);
if (! -f "$pileuppath/$bodysite\.mpileup") {
    if (system("samtools-1.3.1 mpileup -Os -d 8000 -f /mnt/blsm/sarah/targetedstudy/skinplex_panel_db/MarkersForTargetPanel.fasta $filestring > $pileuppath/$bodysite\.mpileup")) {
	die "Samtools-1.3.1 mpileup error: $!";
    }
}
open MPILE, "<$pileuppath/$bodysite\.mpileup" or die "Could not open input file for reading\n";
open OUT, ">$pileuppath/nonuniversal/$bodysite\_ND.txt\n" or die "Could not open output file for writing\n";
my $linenum = 1;
my $samplelabel;
my $nsamples = scalar(@bodysortednames);

while (my $record = <MPILE>) {
    my ($site, $total);
    my (@sites, @bases, @coverage);
    my @samplelabels;
    my (%basestring, %bases, %stotal, %totalcount, %samplecount);
    my $lineoutput;
#    warn "$.\n";

    chomp $record;
    #Store reference fields
    my ($marker, $pos, $refbase, @sampledata) = split/\t/, $record;
    $site = join('_', $marker, $pos);
    push(@sites, $site);

    #Determine number of samples in input file/correct any empty string or * for samples
    my $data = scalar(@sampledata);
    if ($data % 5 != 0) {
	if ($sampledata[$#sampledata] == 0) { #correct empty strings
	    for (1..4) {
		push(@sampledata, '');
	    }
	    $data = scalar(@sampledata);
	} else {
	    die "First:There are not 5 data columns per sample!\nLine number $linenum\n$record\n";
	}
    }
    if ($data % 5 != 0) { #double check the corrected empty strings
	die "Second:There are not 5 data columns per sample!\nLine number $linenum\n$record\n";
    }

    if ($linenum == 1) { #Print header line for output
	my $sampleheader;
	foreach my $samplename (@bodysortednames) {
	    $sampleheader .= "$samplename\t$samplename"."ND\t$samplename"."TotalCov\t";
	}
	print OUT "Marker\tPosition\tTaxonPath\tTaxon\tLevel\t$sampleheader\n";
    }
    my $coverage;
    foreach my $sample (@bodysortednames) { #hash of all bases for each sample
	my @fields = splice (@sampledata, 0, 5);
	if ($fields[0] < 1) {
	    $coverage++;
	}
	$basestring{$sample} = $fields[1];
    }

    if (defined $coverage && ($coverage == $nsamples)) { #skip over any sites that did not align for all samples with at least a coverage of > 0
	$linenum++;
	next;
    }
    
    #count bases for each sample and across all samples
    foreach my $sample (keys %basestring) {#hash of base strings for each sample
	my $samplebases = $basestring{$sample};
	$samplebases =~ s/[+-](\d+)(??{"[ACGTN]{$1}"})//gi;
	my @samplebases = split //, $samplebases;
	foreach my $base (@samplebases) {
	    if ($base eq '.' || $base eq ',') {
		$base = $refbase;
		$samplecount{$base}{$sample}++;
		$totalcount{$base}++;
		$stotal{$sample}++;
		$total++;
	    }
	    elsif ($base eq 'A' || $base eq 'a') {
		$base = 'A';
		$samplecount{$base}{$sample}++;
		$totalcount{$base}++;
		$stotal{$sample}++;
		$total++;
	    }
	    elsif ($base eq 'T' || $base eq 't') {
		$base = 'T';
		$samplecount{$base}{$sample}++;
		$totalcount{$base}++;
		$stotal{$sample}++;
		$total++;
	    }
	    elsif ($base eq 'C' || $base eq 'c') {
		$base = 'C';
		$samplecount{$base}{$sample}++;
		$totalcount{$base}++;
		$stotal{$sample}++;
		$total++;
	    }
	    elsif ($base eq 'G' || $base eq 'g') {
		$base = 'G';
		$samplecount{$base}{$sample}++;
		$totalcount{$base}++;
		$stotal{$sample}++;
		$total++;
	    }
	    else { #not interested in other characters
		next;
	    }
	}
    }
    my @basecalls = keys %totalcount; #Bases observed across samples
    if (!@basecalls) {
	$linenum++;
	next;
    }

#    my $samplecheck = 0;
#    foreach my $sample (@bodysortednames) {
#	if (! exists $stotal{$sample}) {
#	    $samplecheck++;
#	} else {
#	    next;
#	}
#    }
#    if ($samplecheck > 0) {
#	$linenum++;
#	next;
#    }
    
    #Calculate nucleotide diversity for each sample compared to reference allele
    my %nucdiv;
    my $check;
    foreach my $basecall (@basecalls) {
	if ($basecall eq $refbase) {
	    foreach my $sample (@bodysortednames) {
		my $p;
		if (exists $samplecount{$basecall}{$sample}) {
		    $p = $samplecount{$basecall}{$sample}/$stotal{$sample};
		} else {
		    $p = 0;
		}
		my $q = 1-$p;
		$nucdiv{$sample} = 2*$p*$q;
		$markerND{$marker}{$sample} +=$nucdiv{$sample};
		$markerNDlength{$marker}{$sample}++;
		$finalmarkers{$marker}++;
	    } 
	} else {
	    $check++;
	    next;
	}
    }
    #Check if no basecalls matched reference call
    if (defined $check && $check == scalar(@basecalls)) {
	foreach my $sample (@bodysortednames) {
	    $nucdiv{$sample} = 0;
	    $markerND{$marker}{$sample} +=$nucdiv{$sample};
	    $markerNDlength{$marker}{$sample}++;
	    $finalmarkers{$marker}++;
	}
    }
    
    #Output ThetaPi values for each sample at each loci/site
    #Get taxon for each loci
    my $taxonpath = $taxamarker{$marker};
    my $cladelevel = $taxalevel{$taxonpath};
    my $sampleout;

    foreach my $sample (@bodysortednames) {
	if (! exists $stotal{$sample}) {
	    $stotal{$sample} = 0;
	}    
	$sampleout .= "\t$sample\t$nucdiv{$sample}\t$stotal{$sample}";
    }
    
    $lineoutput = "$marker\t$pos\t$taxonpath\t$cladelevel$sampleout";
    print OUT "$lineoutput\n";
    $linenum++;
#    die "$lineoutput\n";
}
close(OUT);
close(MPILE);
