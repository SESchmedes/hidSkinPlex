#!/usr/bin/env Rscript

###############################################
#Created by Sarah E. Schmedes
#Release date: 9/1/2017
##############################################

library(ggplot2)

input.file.common <- "/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/for_accuracy_table_common.txt"
tbl.common <- read.table(input.file.common, header=T, sep="\t")

input.file.nonuniversal <- "/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/mpileup_results/nonuniversal/for_accuracy_table_nonuniversal.txt"
tbl.nonuniversal <- read.table(input.file.nonuniversal, header=T, sep="\t")

#tbl.common$MarkerType <- rep("common", 112)
#tbl.nonuniversal$MarkerType <- rep("nonuniversal", 112)
colnames(tbl.common)[10] = "CommonAccuracy"
colnames(tbl.nonuniversal)[10] = "NonuniversalAccuracy"
tbl.merge <- merge(tbl.common, tbl.nonuniversal, by=c("Bodysite", "Threshold", "Classifier"))
tbl.merge <- tbl.merge[, c("Bodysite", "Threshold", "Classifier", "CommonAccuracy", "NonuniversalAccuracy")]
tbl.merge$Threshold <- as.factor(tbl.merge$Threshold)

#tbl.common$CommonAccuracy <- tbl.common$Accuracy
#tbl.common$NonuniversalAccuracy <- as.numeric(rep(0, 112))
#tbl.nonuniversal$CommonAccuracy <- as.numeric(rep(0, 112))
#tbl.nonuniversal$NonuniversalAccuracy <- tbl.nonuniversal$Accuracy

#tbl.both <- rbind.data.frame(tbl.common, tbl.nonuniversal)

ggplot(tbl.merge, aes(x=CommonAccuracy, y=NonuniversalAccuracy, color=Classifier, shape=Threshold)) + geom_point(size=3, position="jitter") + geom_abline() + facet_wrap(~Bodysite) + xlim(0, 101) + ylim(0,101) +
scale_shape_manual(values=c(15,16,17,18,7,8,9)) + xlab("Classification Accuracy using Universal Markers") + ylab("Classification Accuracy using Non-universal Markers") +
scale_color_discrete(name="Classifier", breaks=c("attselectknn", "attselectLogistic", "knn", "logistic"), labels=c("1NN w/AttSelect", "RMLR w/AttSelect", "1NN", "RMLR"))
  
ggsave("/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/alignstats/figures/AccuracyComparison.png", scale=1, width=30, height =20, units = "cm", dpi = 300, limitsize = T)
