#!/usr/bin/env Rscript

##############################################
#Created by Sarah E. Schmedes
#Release Date: 9/1/2017
##############################################
#This script is used to make figures for align stats for hidskinplex samples
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape)

input.file <- "/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/alignstats/hidskinplex_output.txt"
tbl <- read.table(input.file, header = T, sep = "\t")

#amp.tbl <- unique(cbind.data.frame(tbl$Marker, tbl$AmpliconSize))
tbl <- arrange(tbl, Clade, AmpliconSize)
tbl$Marker <- factor(tbl$Marker, levels=unique(tbl$Marker))
tbl.sample <- aggregate(AvgMarkerCov ~ BodySite + Clade + Marker + SampleID, FUN = "mean", data=tbl)
tbl.body <- aggregate(AvgMarkerCov ~ BodySite + Clade + Marker, FUN ="mean", data = tbl)
#Make coverage per marker per bodysite and subject
ggplot(tbl.sample[tbl.sample$BodySite=="Mb",], aes(x=Marker, y=AvgMarkerCov, fill=Clade)) + geom_col() + facet_wrap(~SampleID, nrow=2, ncol=4) + theme(axis.text.x=element_blank(), legend.position = "top") + scale_y_log10() + ylab("Average Read Depth\nPer Marker")
ggplot(tbl.sample[tbl.sample$BodySite=="Hp",], aes(x=Marker, y=AvgMarkerCov, fill=Clade)) + geom_col() + facet_wrap(~SampleID, nrow=2, ncol=4) + theme(axis.text.x=element_blank(), legend.position= "top") + scale_y_log10() + ylab("Average Read Depth\nPer Marker")
ggplot(tbl.sample[tbl.sample$BodySite=="Fb",], aes(x=Marker, y=AvgMarkerCov, fill=Clade)) + geom_col() + facet_wrap(~SampleID, nrow=2, ncol=4) + theme(axis.text.x=element_blank(), legend.position = "top") + scale_y_log10() + ylab("Average Read Depth\nPer Marker")
ggplot(tbl.body[tbl.body$BodySite!="swabblank",], aes(x=Marker, y=AvgMarkerCov, fill=Clade)) + geom_col() + facet_wrap(~BodySite) + theme(axis.text.x=element_blank(), legend.position = "top", axis.title.x = element_text(size=18), axis.title.y = element_text(size=18)) + scale_y_log10(breaks=c(0.01, 1, 10, 100, 1000, 10000), labels=c("0.01", "1", "10", "100", "1000", "10000")) + ylab("Average Read Depth\nPer Marker")
#ggplot(tbl[tbl$SampleID=="S001",], aes(x=Marker, y=AvgMarkerCov, fill=Clade)) + geom_col() + facet_wrap(~BodySite) + theme(axis.title.x = element_blank(), legend.position = "top")
#ggplot(tbl[tbl$SampleID=="S002",], aes(x=Marker, y=AvgMarkerCov, fill=Clade)) + geom_col() + facet_wrap(~BodySite) + theme(axis.title.x = element_blank(), legend.position = "top")
ggsave("/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/alignstats/figures/BodysiteCov.png", scale=1, width = 55, height = 25, units = "cm", dpi = 300, limitsize = T)
#Make read count figure
read.tbl <- select(tbl, Samplename, BodySite, TotalRawSeq, PostQCtotalReads)
read.tbl.melt <- melt(read.tbl, c("Samplename", "BodySite"))
colnames(read.tbl.melt)[3] = "SequenceReads"
colnames(read.tbl.melt)[4] = "TotalReads"
ggplot(read.tbl.melt, aes(x=Samplename, y=TotalReads, fill=SequenceReads)) + geom_col(position = "dodge")

#Make quality figure
ggplot(tbl, aes(x=Samplename, y=AvgQuality, fill=Primer)) + geom_col(position = "dodge") + facet_wrap(~Bodysite)

#Save images for coverage per marker plots
ggsave("/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/alignstats/figures/coverage.png", scale=1, width = 55, height = 25, units = "cm", dpi = 300, limitsize = T)
#Save images for accuracy per sample figure
ggsave("/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/alignstats/figures/Quality.png", scale=1, width=42, height =20, units = "cm", dpi = 300, limitsize = T)
