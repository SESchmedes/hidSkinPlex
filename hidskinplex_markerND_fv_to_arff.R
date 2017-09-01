#!/usr/bin/Rscript

#This script normalizes marker nucleotide diversity data
#and converts the .txt file to an .arff file
#Arguments input file, bodysite

library(foreign)
library(GetoptLong)

#set flags
GetoptLong(
    "bodysite=s", "bodysite symbol, character, mandatory option",
    "threshold=s", "threshold, character, mandatory option")
#Read in data
input.file <- sprintf("/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/mpileup_results/%s_markerND_fv_gteq%s.txt", bodysite, threshold)
fv <- read.table(input.file, header=T, sep="\t")

#Remove any columns with sd < 1e-6
alldevs <- apply(fv,2,sd)
alldevs['Individual'] <- 1
fv2 <- subset(fv, T, select=(alldevs > 1e-6))
fv2$Individual <- as.factor(fv2$Individual)

#Output .arff file
write.table(fv2, file=sprintf("/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/mpileup_results/%s_markerND_fv_R_gteq%s.txt", bodysite, threshold), sep="\t", quote=F, eol = "\n", row.names = F)
write.arff(fv2, file=sprintf("/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/mpileup_results/%s_markerND_fv_gteq%s.arff", bodysite, threshold))
