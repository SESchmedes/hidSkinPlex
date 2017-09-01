#!/usr/bin/env perl

use strict;
use warnings;
use Genomix qw(:constants :meta);

my @bodysites = ("bodyclass");
my @classifiers = ("logistic", "knn", "attselectLogistic", "attselectknn");
my @thresholds = (2, 10, 25, 50, 100, 150, 200);
open OUT, ">/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/mpileup_results/nonuniversal/for_accuracy_table_nonuniversal_bodyclass.txt" or die "Could not open file for writing\n";
print OUT "Bodysite\tThreshold\tNumSamples\tNumIndividuals\tMarkrersAligned\tNumMarkers\tClassifier\tNumCorrectClass\tNumIncorrectClass\tAccuracy\tRandomAccuracy\tKappastatistic\n";

open WIN, ">/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/mpileup_results/nonuniversal/for_accuracy_table_nonuniversal_bodyclass_wide.txt" or die "Could not open table file for writing\n";
print WIN "Bodysite\tThreshold\tNumSamples\tNumIndividuals\tMarkersAligned\tNumMarkers\tNumFSmarkers\tLogistic\tLogKappa\tKNN\tKNNkappa\tFSLogistic\tFSlogkappa\tFSknn\tFSknnkappa\tRandomAccuracy\n";

foreach my $bodysite (@bodysites) {
    foreach my $threshold (@thresholds) {
	my (%attributes, %correctclass, %accuracy, %incorrectclass, %bodymarkers, %kappastat);
	my $natt;
	my $wc = `wc -l /mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/mpileup_results/nonuniversal/$bodysite\_markerND_fv_R_gteq$threshold\.txt`;
	my @wc = split /\s+/, $wc;
	my $nsamples = $wc[0] - 1;
	my $nind = $nsamples/3;
	if ($bodysite eq "all") {
	    $nind = $nsamples/9;
	}
	open FV, "</mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/mpileup_results/nonuniversal/$bodysite\_markerND_fv_R_gteq$threshold\.txt" or die "Could not open FV for reading\n";
	my $linenum = 1;
	my $markersnum;
	while (<FV>) {
	    my $line = $_;
	    chomp $line;
	    if ($linenum == 1) {
		my @markers = split /\t/, $line;
		pop @markers;
		$markersnum = scalar(@markers);
		last;
	    }
	}
	close(FV);
	open NDFV, "</mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/mpileup_results/nonuniversal/$bodysite\_markerND_fv_gteq$threshold\.txt" or die "Could not open FV for reading\n";
	my $nline = 1;
	my $initialmarkersnum;
	while (<NDFV>) {
	    my $head = $_;
	    chomp $head;
	    if ($nline == 1) {
		my @ndmarkers = split /\t/, $head;
		pop @ndmarkers;
		$initialmarkersnum = scalar(@ndmarkers);
		last;
	    }
	}
	close(NDFV);
	foreach my $class (@classifiers) {
	    open CLASS, "</mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/mpileup_results/nonuniversal/$bodysite\_weka_$class\_markerND_gteq$threshold\.txt" or die "Cannot open weka results for reading. $bodysite $!\n";
	    if ($class eq "logistic" || $class eq "knn") {
		while (<CLASS>) {
		    my $line = $_;
		    chomp $line;
		    if ($line =~ /=== Stratified cross-validation ===/) {
			my $blank = <CLASS>;
			chomp $blank;
			my $correctline = <CLASS>;
			chomp $correctline;
			my @correctinfo = split /\s+/, $correctline;
			$correctclass{$class} = $correctinfo[3];
			$accuracy{$class} = $correctinfo[4];
			my $incorrectline = <CLASS>;
			chomp $incorrectline;
			my @incorrectinfo = split /\s+/, $incorrectline;
			$incorrectclass{$class} = $incorrectinfo[3];
			my $kappaline = <CLASS>;
			chomp $kappaline;
			my @kappastat = split /\s+/, $kappaline;
			$kappastat{$class} = $kappastat[2];
		    }
		}
	    } else {
		while (<CLASS>) {
		    my $line = $_;
		    chomp $line;
		    if ($line =~ /Selected attributes:/) {
			my @att = split /\s+/, $line;
			$natt = $att[$#att];
		    }
		    if ($line =~ /=== Stratified cross-validation ===/) {
			my $blank = <CLASS>;
			chomp $blank;
			my $correctline = <CLASS>;
			chomp $correctline;
			my @correctinfo = split /\s+/, $correctline;
			$correctclass{$class} = $correctinfo[3];
			$accuracy{$class} = $correctinfo[4];
			my $incorrectline = <CLASS>;
			chomp $incorrectline;
			my @incorrectinfo = split /\s+/, $incorrectline;
			$incorrectclass{$class} = $incorrectinfo[3];
			my $kappaline = <CLASS>;
			chomp $kappaline;
			my @kappastat = split /\s+/, $kappaline;
			$kappastat{$class} = $kappastat[2];
		    }
		}
	    }
	    my $nmarkers;
	    close(CLASS);
	    if ($class eq "knn" || $class eq "logistic") {
		$nmarkers = $markersnum;
	    } else {
		$nmarkers = $natt;
	    }
	    my $random = (23/($nsamples-1))*100;
	    print OUT "$bodysite\t$threshold\t$nsamples\t$nind\t$initialmarkersnum\t$nmarkers\t$class\t$correctclass{$class}\t$incorrectclass{$class}\t$accuracy{$class}\t$random\t$kappastat{$class}\n";
	}
	my $classout;
	foreach my $class (@classifiers) {
	    $classout .= "$accuracy{$class}\t$kappastat{$class}\t";
	}
	my $random = (23/($nsamples-1))*100;
	print WIN "$bodysite\t$threshold\t$nsamples\t$nind\t$initialmarkersnum\t$markersnum\t$natt\t$classout$random\n";   
    }
}
close(OUT);
close(WIN);
