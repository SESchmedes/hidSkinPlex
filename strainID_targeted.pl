#!/usr/bin/env perl
########################################################
#Created by Sarah Schmedes
#Will run entire strainphlan workflow, including images
#Designate -clade flag
#Date Released: 9/1/2017 
#########################################################

use strict;
use warnings;
use Genomix qw(:constants :meta);
use Getopt::Long;

#select samples to run through strainphlan
my $metafile = $OH_METAFILE;
my $sampath; #path to sam files from metaphlan2 run
my $straindir; #output directory for all strainphlan results
my $samplelist;
my ($clade, $refgen);
#my @bodysites = ("all", "Fb", "Hp", "Mb");
my @bodysites = ("Fb");
#set flags
GetOptions('clade=s' => \$clade,
#	   'refgen=s' => \$refgen)
	   'sampath=s' => \$sampath,
	   'straindir=s' => \$straindir,
	   'samplelist=s' => \$samplelist) 
    or die "Must use analysis flag and provide 'clade' argument\n"; 
#Get samplelist
open SAMPLES, "<$samplelist" or die "Could not open sample list for reading! $!\n";
my %bodysamples;
while (<SAMPLES>) {
    my $file1 = $_;
    chomp $file1;
    my $file2 = <SAMPLES>;
    chomp $file2;
    my @name = split /\./, $file1;
    my @name2 = split /_/, $name[0];
    my $samplename = $name2[0];
    if ($samplename =~ /^SwabBlank/) {
	next;
    }
    my ($ind, $bod, $time, $rep) = split /-/, $samplename;
    $bodysamples{$samplename} = $bod;
}
close(SAMPLES);
my @sortednames = sort {$a cmp $b} keys %bodysamples;

#run strainphlan on all samples
#all metaphlan2/strainphlan scripts must be in the environmental path
#Also make metadata sheet at the same time
foreach my $bodysite (@bodysites) {
    open META, ">$straindir/$bodysite/$bodysite\_metadata.txt" or die "Could not open metadata file for writing $!\n";
    print META "Sample\tSubjectID\tTimePoint\tBodysite\n";
    foreach my $sample (@sortednames) {
	my ($ind, $bod, $time, $rep) = split /-/, $sample;
	if (($bodysite ne "all") && ($bod ne $bodysite)) {
	    next;
	}
	print META "$sample\t$ind\t$time\t$bod\n";
	#Generate a marker file for each sample.
	#The marker file contains the consensus of unique marker genes for each species found in the sample
	#This marker file will be used for SNP profiling
	if (! -f "$straindir/$bodysite/consensus_markers/$sample\.markers") {
	    if (system("sample2markers.py --ifn_samples $sampath/$sample\.sam.bz2 --input_type sam --output_dir $straindir/$bodysite/consensus_markers --nprocs 20 1> $straindir/$bodysite/consensus_markers/$sample\_log.txt 2> $straindir/$bodysite/consensus_markers/$sample\_error.txt")) {
		die "Strainphlan sample2markers.py ERROR: $!\n";
	    }
	}
    }
    print META "Pacnes_GCF_000217615.1\tReference Genome\tT1\tReference\n";
    close(META);
    #Run strainphlan to identify clades that were detected in all samples
    #providing the marker files generated in the prior step
    #to see which clades can be SNP-profiled
    if (! -f "$straindir/$bodysite/clades/$bodysite\_clades.txt") {
	if (system("strainphlan.py --ifn_samples $straindir/$bodysite/consensus_markers/*.markers --output_dir $straindir/$bodysite/clades --print_clades_only --nprocs_main 20 1> $straindir/$bodysite/clades/$bodysite\_clades.txt 2> $straindir/$bodysite/clades/$bodysite\_errorlog.txt")) {
	    die "strainphlan.py ERROR: $!\n";
	} 
    }
    
    if (! -f "$straindir/db_markers/$clade\.markers.fasta") {
	#Build reference database for the designated clade
	#This step only needs to be done once for each species for all projects
	if (system("extract_markers.py --mpa_pkl /home/sarah/strainphlan/db_v20/mpa_v20_m200.pkl --ifn_markers $straindir/db_markers/all_markers.fasta --clade $clade --ofn_markers $straindir/db_markers/$clade\.markers.fasta")) {
	    die "Strainphlan extract_markers.py ERROR: $!\n";
	}
    }
    else  {
	print "$clade\.markers.fasta already exits\n";
    }
    
#Build the multiple sequence alignment and phylogenetic tree
#Will align and clean sample-reconstructed strains (stored in .markers)
#and reference-genome-reconstructed strains (from clade.markers.fasta)
#Builds tree using RAxML
#If a reference genome is not specified or if no clade is specified then
#strainphlan.py will build the tree for all species it can detect
#    if (system("strainphlan.py --mpa_pkl /home/sarah/strainphlan/db_v20/mpa_v20_m200.pkl --ifn_samples $straindir/$bodysite/consensus_markers/*.markers --ifn_markers $straindir/db_markers/$clade\.markers.fasta --ifn_ref_genomes $straindir/reference_genomes/$refgen\.fna.bz2 --output_dir $straindir/$bodysite/output --relaxed_parameters2 --nprocs_main 5 --clades $clade 1> $straindir/$bodysite/output/log_full.txt 2> $straindir/$bodysite/output/error_full.txt"))
#    {
#	die "strainphlan.py ERROR: $!\n";
 #   }
    if (system("strainphlan.py --mpa_pkl /home/sarah/strainphlan/db_v20/mpa_v20_m200.pkl --ifn_samples $straindir/$bodysite/consensus_markers/*.markers --ifn_markers $straindir/db_markers/$clade\.markers.fasta --output_dir $straindir/$bodysite/output --marker_in_clade 0.1 --relaxed_parameters2 --nprocs_main 20 --clades $clade 1> $straindir/$bodysite/output/log_full.txt 2> $straindir/$bodysite/output/error_full.txt"))
    {
	die "strainphlan.py ERROR: $!\n";
    }
    
#Add metadata to the tree
#must of metadata file in strainphlan group directory
#multiple trees and multiple metadata files can be used (space separated, and wild card can be used)
#metadata file (tab separated, can have multiple columns)
    my $metadata = "SubjectID"; #change based on what metadata you want listed on the tree
    if (system("add_metadata_tree.py --ifn_trees $straindir/$bodysite/output/RAxML_bestTree.$clade\.tree --ifn_metadatas $straindir/$bodysite/$bodysite\_metadata.txt --metadatas $metadata"))
    {
	die "Strainphlan add_metadata_tree.py ERROR: $!\n";
    }
    
#Plot tree using Graphlan (graphlan scripts must be in path)
    if (system("plot_tree_graphlan.py --ifn_tree $straindir/$bodysite/output/RAxML_bestTree.$clade\.tree.metadata --colorized_metadata $metadata --leaf_marker_size 60 --legend_marker_size 60"))
    {
	die "Graphlan plot_tree_graphlan.py ERROR: $!\n";
    }
    
#Create dendrogram using ggtree script
#breadcrumbs directory must be in path
    if (system("strainphlan_ggtree_Mod.R $straindir/$bodysite/output/RAxML_bestTree.$clade\.tree $straindir/$bodysite/$bodysite\_metadata.txt $straindir/$bodysite/output/$clade\.fasta $straindir/$bodysite/output/$bodysite\_$clade\_ggtree_1.png $straindir/$bodysite/output/$bodysite\_$clade\_ggtree_2.png"))
    {
	die "strainphlan_ggtree_Mod.R ERROR: $!\n";
    }
    
#Create a distance matrix
    if (system("distmat -sequence $straindir/$bodysite/output/$clade\.fasta -nucmethod 2 -outfile $straindir/$bodysite/output/$clade\.distmat"))
    {
	die "distmat ERROR: $!\n"; 
    }
    if (system("strainphlan_ordination_Mod.R $straindir/$bodysite/output/$clade\.distmat $straindir/$bodysite/$bodysite\_metadata.txt $straindir/$bodysite/output/$bodysite\_$clade\_strainord.png"))
    {
	die "strainphlan_ordination_Mod.R ERROR: $!\n";
    }
}
