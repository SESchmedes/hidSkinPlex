#!/usr/bin/env perl

use strict;
use warnings;
use Genomix qw(:constants :meta);

my @bodysites = ("all", "Fb", "Hp", "Mb");
my $pilepath = "/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/mpileup_results";
my (%ampsize, %amppaths);
my $threshold = 10;
#Get amplicon info
open AMP, "</mnt/blsm/sarah/targetedstudy/skinplex_panel_db/targetedpanel_ampliconsize.txt" or die "Cannot open amplicon size file for reading $!\n";
while (<AMP>) {
    my $line = $_;
    chomp $line;
    if ($line =~ /^Marker/) {
	next;
    }
    my ($marker, $taxonpath, $clade, $level, $markerlength, $ampliconsize) = split /\t/, $line;
    $amppaths{$taxonpath}++;
    $ampsize{$marker} = $ampliconsize;
}
close(AMP);
my @panelmarkers = sort {$a cmp $b} keys %ampsize;

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
my ($taxa_ref, $taxalevel_ref, $markerlength_ref) = gettaxa();
my %taxamarker = %$taxa_ref;
my %taxalevel = %$taxalevel_ref;
my %markerlength = %$markerlength_ref;
my @totalDBmarkers = sort {$a cmp $b} keys %taxamarker;
#Run pileparseMarkerND.pl to produce new mpileup, ND calculations, and master output
foreach my $bodysite (@bodysites) {
    my %bodysortednames;
    my @bodysortednames;
    if (! -f "$pilepath/$bodysite\_ND.txt") {
	if (system("perl /home/sarah/src/hidskinplex_pileparseMarkerND.pl -bodysite $bodysite")) {
	    die "problem executing hidskinplex_pileparseMarkerND.pl. $!\n";
	}
    }
    foreach my $sample (@sortednames) {
	my ($sampleid, $body, $time, $rep) = split /-/, $sample;
	if ($body eq $bodysite) {
	    $bodysortednames{$sample} = 1;
	}
    } 
    if ($bodysite eq "all") {
	@bodysortednames = @sortednames;
    }
    if ($bodysite ne "all") {
	@bodysortednames = sort {$a cmp $b} keys %bodysortednames;
    }
    my $nsamples = scalar(@bodysortednames);
    
    open IN, "<$pilepath/$bodysite\_ND.txt" or die "Could not open ND.txt file for reading\n";
    my %markercov;
    my %samplemarkercov;
    my %sampleND;
    my %samplemarkerlengthcov;
    while (<IN>) {
	my $line = $_;
	chomp $line;
	if ($line =~ /^Marker/) {
	    next;
	}
	my ($marker, $pos, $taxonpath, $taxon, $level, @sampledata) = split /\t/, $line;
	my $covcheck = 0;
	#save if you want to try for only common markers
	for (my $i=2; $i<=$#sampledata; $i+=3) {
	    if ($sampledata[$i] < $threshold) {
		$covcheck++;
	    }
	}
	if ($covcheck > 0) {
	    next;
	}
	$markercov{$marker}++;
	foreach my $sample (@bodysortednames) {
	    my ($sampleID, $sampleND, $samplecov) = splice (@sampledata, 0, 3);
	    if ($sample eq $sampleID) {
#		if ($samplecov >= 100) {
		    $samplemarkercov{$marker}{$sample} += $samplecov;
		    $sampleND{$marker}{$sample} += $sampleND;
#		    $samplemarkerlengthcov{$sample}++;
#		}
	    } else {
		die "Sample do not align\n";
	    }
	}
    }
    close(IN);
    
    my @markers = sort {$a cmp $b} keys %markercov;
    
    open FV, ">$pilepath/$bodysite\_markerND_fv_gteq$threshold\.txt" or die "Cannot open fv file for writing\n";
    my $markerstring = join("\t", @markers);
    print FV "$markerstring\tIndividual\n";
    
    foreach my $sample (@bodysortednames) {
	my $out;
	foreach my $marker (@markers) {
	    my $ND = $sampleND{$marker}{$sample}/$markercov{$marker};	    
	    $out .= "$ND\t";
	}
	my ($sampleid, $body, $time, $rep) = split /-/, $sample;
	$out .= $sampleid;
	print FV "$out\n";
    }
    close(FV);
    
    
    
#Run R script to create .arff
    if(system("Rscript /home/sarah/src/hidskinplex_markerND_fv_to_arff.R --bodysite $bodysite --threshold $threshold")) {
	die "Error executing rscript to create .arff file\n";
    }
    
#weka commands here
#Logistic
    if(system("java -Xmx8g -classpath /home/sarah/weka-3-8-1/weka.jar weka.classifiers.functions.Logistic -R 1.0E-8 -M -1 -num-decimal-places 4 -x $nsamples -t $pilepath/$bodysite\_markerND_fv_gteq$threshold\.arff > $pilepath/$bodysite\_weka_logistic_markerND_gteq$threshold\.txt")) {
	die "Error executing java and weka logistic. $!\n";
    }
    
#KNN
    if(system("java -Xmx8g -classpath /home/sarah/weka-3-8-1/weka.jar weka.classifiers.lazy.IBk -K 1 -x $nsamples -t $pilepath/$bodysite\_markerND_fv_gteq$threshold\.arff > $pilepath/$bodysite\_weka_knn_markerND_gteq$threshold\.txt")) {
	die "Error executing java and weka knn. $!\n";
    }
    
#AttSelectLog
    if(system("java -Xmx8g -classpath /home/sarah/weka-3-8-1/weka.jar weka.classifiers.meta.AttributeSelectedClassifier -x $nsamples -t $pilepath/$bodysite\_markerND_fv_gteq$threshold\.arff -W weka.classifiers.functions.Logistic -- -R 1.0E-8 -M -1 -num-decimal-places 4 > $pilepath/$bodysite\_weka_attselectLogistic_markerND_gteq$threshold\.txt")) {
	die "Error executing java and weka attselect logistic. $!\n";
    }
    
#AttSelectKNN
    if(system("java -Xmx8g -classpath /home/sarah/weka-3-8-1/weka.jar weka.classifiers.meta.AttributeSelectedClassifier -x $nsamples -t $pilepath/$bodysite\_markerND_fv_gteq$threshold\.arff -W weka.classifiers.lazy.IBk -- -K 1 > $pilepath/$bodysite\_weka_attselectknn_markerND_gteq$threshold\.txt")) {
	die "Error executing java and weka attselect knn euclidean. $!\n";
    }
    
    my @classifiers =("logistic", "knn", "attselectLogistic", "attselectknn");
    my (%attributes, %correctclass, %accuracy, %incorrectclass);
    foreach my $class (@classifiers) {
	open CLASS, "</$pilepath/$bodysite\_weka_$class\_markerND_gteq$threshold\.txt" or die "Cannot open weka results for reading. $!\n";
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
		}
	    }
	} else {
	    while (<CLASS>) {
		my $line = $_;
		chomp $line;
		if ($line =~ /Selected attributes:/) {
		    my @att = split /\s+/, $line;
		    my $natt = $att[$#att];
		    for my $n (1..$natt) {
			my $m = <CLASS>;
			chomp $m;
			my @marker = split /\s+/, $m;
			$attributes{$marker[1]}{$class}++;
		    }
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
		}
	    }
	}	
	close(CLASS);
    }
#Get all taxon paths for each feature vector marker
#and number of markers per clade
    my %markersperpath;
    foreach my $marker (@markers) {
	my $taxonpath = $taxamarker{$marker};
	$markersperpath{$taxonpath}++;
    }
    my @paths = sort {$a cmp $b} keys %markersperpath;

#Get number of FS markers per clade
    my %featurecheck;
    foreach my $marker (@markers) {
	my $revise = $marker;
	$revise =~ s/[|:-]/./g;
	foreach my $class (@classifiers) {
	    if (exists $attributes{$revise}{$class}) {
		$featurecheck{$marker}++;
	    } else {
		next;
	    }
	}
    }
    my @FSmarkers = sort {$a cmp $b} keys %featurecheck;
    foreach my $marker (@FSmarkers) {
	if ($featurecheck{$marker} != 2) {
	    die "There are differences in the marker attributes selected for each AttSelect classifier\n";
	}
    }
    my %fsmarkcount;
    foreach my $fsmark (@FSmarkers) {
	my $path = $taxamarker{$fsmark};
	$fsmarkcount{$path}++;
    }
    foreach my $path (@paths) {
	if (! exists $fsmarkcount{$path}) {
	    $fsmarkcount{$path} = 0;
	}
    }
    
#Prepare line out
    open RESULTS, ">/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/mpileup_results/$bodysite\_markerND_weka_compare_MASTER_gteq$threshold\.txt" or die "Cannot open output for writing. $!\n";
    print RESULTS "Marker\tAmpSizeDB\tMarkerCoverage\tTaxonPath\tTaxon\tLevel\tTargetMarkerCountPanel\tSampleID\tSampleAvgMarkerCov\tBodysite\tTotalCladeMarkersAligned\tSelectedFeatures\tLogistic\tKNN\tAttSelectLog\tAttSelectKNN\n";
    open SHORT, ">/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/mpileup_results/$bodysite\_markerND_weka_compare_short_gteq$threshold\.txt" or die "Cannot open short output for writing. $!\n";
    print SHORT "Bodysite\tTaxonPath\tTaxon\tLevel\tNumMarkersinPanel\tNumofMarkersAligned\tNumofMarkersFSClade\tLogistic\tKNN\tAttSelectLog\tAttSelectKNN\n";
    my $accuracyout;
    foreach my $class (@classifiers) {
	$accuracyout .= "\t$accuracy{$class}";
    }
    foreach my $marker (@markers) {
	my $FS;
	my $length = $ampsize{$marker};
	my $cov = $markercov{$marker};
	my $taxonpath = $taxamarker{$marker};
	my $cladelevel = $taxalevel{$taxonpath};
	my ($clade, $level) = split /\t/, $cladelevel;
	if (exists $featurecheck{$marker}) {
	    $FS = 1;
	} else {
	    $FS = 0;
	}
	foreach my $sample (@bodysortednames) {
	    my $avgsitecov = $samplemarkercov{$marker}{$sample}/$cov;
	    my $clademarkers = $markersperpath{$taxonpath};
	    my $lineout = "$marker\t$length\t$cov\t$taxonpath\t$clade\t$level\t$amppaths{$taxonpath}\t$sample\t$avgsitecov\t$bodysite\t$clademarkers\t$FS$accuracyout";
	    print RESULTS "$lineout\n";
	}
    }
#Condensed output
    foreach my $path (@paths) {
	my $cladelevel = $taxalevel{$path};
	my ($clade, $level) = split /\t/, $cladelevel;
	my $markersforclade = $markersperpath{$path};
	my $fsmarkerforclade = $fsmarkcount{$path};
	my $out = "$bodysite\t$path\t$clade\t$level\t$amppaths{$path}\t$markersforclade\t$fsmarkerforclade$accuracyout";
	print SHORT "$out\n";
    }
    close(RESULTS);
    close(SHORT);
}
