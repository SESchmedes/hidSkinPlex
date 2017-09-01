#!/usr/bin/env perl

#############################################
#Created by Sarah E. Schmedes
#Release Date: 9/1/2017
#############################################

use strict;
use warnings;
use Genomix qw(:constants :meta);
use Getopt::Long;

my ($samplelist, $datapath, $outpath, $bedpath, $qcdatapath, $amppath);
my (%rawtotalseq, %ampsize, %postqctotalseq, %mappedreads, %avglength, %maxlength, %avgqual,
    %bedcov, %avgmarkercov);
#Set flags
GetOptions('samplelist=s' => \$samplelist,
    'datapath=s' => \$datapath,
    'outpath=s' => \$outpath,
    'bedpath=s' => \$bedpath,
    'qcdatapath=s' => \$qcdatapath,
    'amppath=s' => \$amppath)
    or die "Error with flags. Designate path to samplelist\n";

#Get marker info
my ($taxa_ref, $taxalevel_ref, $marker_ref) = gettaxa();
my %taxamarker = %$taxa_ref;
my %taxalevel = %$taxalevel_ref;
my %markerlength = %$marker_ref;
my @metaphlanmarkers = sort {$a cmp $b} keys %taxamarker;

#Get amplicon info
open AMP, "<$amppath" or die "Cannot open amplicon size file for reading $!\n";
while (<AMP>) {
    my $line = $_;
    chomp $line;
    if ($line =~ /^Marker/) {
	next;
    }
    my ($marker, $taxonpath, $clade, $level, $markerlength, $ampliconsize) = split /\t/, $line;
    $ampsize{$marker} = $ampliconsize;
}
close(AMP);
my @panelmarkers = sort {$a cmp $b} keys %ampsize;

#Get sample names and file names
open IN, "<$samplelist" or die "Could not open sample list for reading $!\n";
my @samplenames;
while (<IN>) {
    my $file1 = $_;
    chomp $file1;
    my $file2 = <IN>;
    chomp $file2;
    my @name = split /\./, $file1;
    my @name2 = split /_/, $name[0];
    my $samplename = $name2[0];
    push(@samplenames, $samplename);
    my $lines1 = `zcat $datapath/$file1 | wc -l`;
    my $wc1 = $lines1/4;
    my $lines2 = `zcat $datapath/$file2 | wc -l`;
    my $wc2 = $lines2/4;
    if ($lines1 % 4 != 0 || $lines2 % 4 !=0) {
	die "corrupted fastq file $file1 or $file2\n";
    }
    my $qclines1 = `zcat $qcdatapath/$file1 | wc -l`;
    my $qcwc1 = $qclines1/4;
    my $qclines2 = `zcat $qcdatapath/$file2 | wc -l`;
    my $qcwc2 = $qclines2/4;
    if ($qclines1 % 4 != 0 || $qclines2 % 4 != 0) {
	die "corrupted fastq file $file1 or $file2\n";
    }
    $rawtotalseq{$samplename} = $wc1 + $wc2;
    $postqctotalseq{$samplename} = $qcwc1 + $qcwc2;
}
close(IN);

#Create .sorted.bam files containing only markers in HIDSkinPlex bed file, align stats, and coverage at each marker
foreach my $sample (@samplenames) {
    if (! -f "$outpath/$sample\.sorted.bam") {
	if (system("bzcat $outpath/$sample\.sam.bz2 | samtools-1.3.1 view -h -u -L $bedpath - | samtools-1.3.1 sort -m 10G -o $outpath/$sample\.sorted.bam -")) {
	    die "Error creating .sorted.bam file $!\n";
	}
    }
    if (! -f "$outpath/alignstats/$sample\.alignstats") {
	if (system("samtools-1.3.1 stats $outpath/$sample\.sorted.bam > $outpath/alignstats/$sample\.alignstats")) {
	    die "Error running samtools stats $!\n";
	}
    }
    if (! -f "$outpath/$sample.sorted.bam.bai") {
	if (system("samtools-1.3.1 index $outpath/$sample\.sorted.bam")) {
	    die "Error creating .sorted.bam.bai $!\n";
	}
    }
    if (! -f "$outpath/alignstats/$sample\.marker.bedcov") {
	if (system("samtools-1.3.1 bedcov $bedpath $outpath/$sample\.sorted.bam > $outpath/alignstats/$sample\.marker.bedcov")) {
	    die "Error running samtools bedcov $!\n";
	}
    }
    if ($sample =~ /^Undetermined/) {
	next;
    }
#Store read/alignment stats for each sample
#parse alignstats files and bedcov (use amplicon size)
    open STATS, "<$outpath/alignstats/$sample\.alignstats" or die "Could not open alignstats file for reading for $sample $!\n";
    my $linenum = 1;
    while (<STATS>) {
	my $line = $_;
	chomp $line;
	if ($linenum == 14) {
	    my @mapped = split /\s+/, $line;
	    $mappedreads{$sample} = $mapped[3];
	    $linenum++;
	    next;
	} elsif ($linenum == 30) {
	    my @averagelength = split /\s+/, $line;
	    $avglength{$sample} = $averagelength[3];
	    $linenum++;
	    next;
	} elsif ($linenum == 31) {
	    my @maxlength = split /\s+/, $line;
	    $maxlength{$sample} = $maxlength[3];
	    $linenum++;
	    next;
	} elsif ($linenum == 32) {
	    my @qual = split /\s+/, $line;
	    $avgqual{$sample} = $qual[3];
	    $linenum++;
	    next;
	} else {
	    $linenum++;
	    next;
	}
    }
    close(STATS);
    
    open COV, "<$outpath/alignstats/$sample\.marker.bedcov" or die "Cannot open bedcov file for $sample for reading $!\n";
    while (<COV>) {
	my $line = $_;
	chomp $line;
	my ($marker, $start, $end, $depth) = split /\t/, $line;
	$bedcov{$marker}{$sample} = $depth;
    }
    close(COV);
    #Get average amplicon coverage
    foreach my $amp (@panelmarkers) {
	foreach my $refmarker (@metaphlanmarkers) {
	    if ($amp eq $refmarker) {
		my $taxonpath = $taxamarker{$amp};
		my $cladelevel = $taxalevel{$taxonpath};
		my ($clade, $level) = split /\t/, $cladelevel;
		my $avgampcov = $bedcov{$amp}{$sample}/$ampsize{$amp};
		$avgmarkercov{$amp}{$sample} = $avgampcov;
	    }
	}
    }
}
#Generate output with all info
open OUT, ">$outpath/alignstats/hidskinplex_output.txt" or die "Could not open output file for writing $!\n";
print OUT "Marker\tClade\tLevel\tMarkerSizeDatabase\tAmpliconSize\tPercentMarkerCovbyAmp\tSamplename\tSampleID\tBodySite\tTimePoint\tReplicate\tAvgMarkerCov\tTotalRawSeq\tPostQCtotalReads\tTotalMappedReads\tAvgQuality\tMaxLength\tAvgLength\n";
foreach my $amp (@panelmarkers) {
    my $taxonpath = $taxamarker{$amp};
    my $cladelevel = $taxalevel{$taxonpath};
    my ($clade, $level) = split /\t/, $cladelevel;
    my $markerdatalength = $markerlength{$amp};
    my $ampsize = $ampsize{$amp};
    my $percentmarkercov = ($ampsize/$markerdatalength)*100;
    foreach my $sample (@samplenames) {
	my ($sampleid, $bodysite, $time, $replicate);
	if ($sample =~ /^Undetermined/) {
	    next;
	}
	if ($sample =~ /^SwabBlank/) {
	    $sampleid = $sample;
	    $bodysite = "swabblank";
	    $time = "NA";
	    $replicate = "R1";
	} else {
	    ($sampleid, $bodysite, $time, $replicate) = split /-/, $sample;
	}
	my $out = "$amp\t$clade\t$level\t$markerdatalength\t$ampsize\t$percentmarkercov\t$sample\t$sampleid\t$bodysite\t$time\t$replicate\t$avgmarkercov{$amp}{$sample}\t$rawtotalseq{$sample}\t$postqctotalseq{$sample}\t$mappedreads{$sample}\t$avgqual{$sample}\t$maxlength{$sample}\t$avglength{$sample}";
	print OUT "$out\n";
    }
}
close(OUT);
