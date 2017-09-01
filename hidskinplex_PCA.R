#!/usr/bin/Rscript

library(ggplot2)
library(dplyr)
library(glmnet)
library(RColorBrewer)
library(cowplot)

dat.all <- read.table("/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/mpileup_results/all_markerND_fv_R_gteq10.txt", header=T, sep="\t")
dat.all.body <- dat.all
dat.all.body$body <- rep(c(rep("Fb", 3), rep("Hp", 3), rep("Mb", 3)), 4) 
dat.all.matrix <- as.matrix(select(dat.all, -Individual))
pcafit.all <- prcomp(dat.all.matrix, center=T, scale.=T)
names(pcafit.all)
var.all <- pcafit.all$sdev^2
var.all.percent <- var.all/sum(var.all)
df.all <- data.frame(dat.all$Individual, dat.all.body$body, pcafit.all$x[,1], pcafit.all$x[,2])
colnames(df.all) <- c("Individual", "Bodysite", "PC1", "PC2")
plot.all <- ggplot(df.all, aes(x=PC1, y=PC2, color=Individual)) + geom_point(size=3) + xlab("PC1 (46% variance explained)") + ylab("PC2 (14% variance explained)") + theme(legend.position = "none")
#+ scale_color_discrete(name="Body Site", breaks=c("Fb", "Hp","Mb"), labels=c("Fb", "Hp", "Mb"))
ggsave("/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/alignstats/figures/all_bodysite_gteq10_PCA.png", plot=plot.all, scale = 1, width=25, height=20, units=c("cm"),dpi=300, limitsize=T)

dat.fb <- read.table("/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/mpileup_results/Fb_markerND_fv_R_gteq10.txt", header=T, sep="\t")
dat.fb.matrix <- as.matrix(select(dat.fb, -Individual))
pcafit.fb <- prcomp(dat.fb.matrix, center=T, scale.=T)
names(pcafit.fb)
var.fb <- pcafit.fb$sdev^2
var.fb.percent <- var.fb/sum(var.fb)
df.fb <- data.frame(dat.fb$Individual, pcafit.fb$x[,1], pcafit.fb$x[,2])
colnames(df.fb) <- c("Individual", "PC1", "PC2")
plot.fb <- ggplot(df.fb, aes(x=PC1, y=PC2, color=Individual)) + geom_point(size=3) + xlab("PC1 (56% variance explained)") + ylab("PC2 (10% variance explained)") 

dat.hp <- read.table("/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/mpileup_results/Hp_markerND_fv_R_gteq10.txt", header=T, sep="\t")
dat.hp.matrix <- as.matrix(select(dat.hp, -Individual))
pcafit.hp <- prcomp(dat.hp.matrix, center=T, scale.=T)
names(pcafit.hp)
var.hp <- pcafit.hp$sdev^2
var.hp.percent <- var.hp/sum(var.hp)
df.hp <- data.frame(dat.hp$Individual, pcafit.hp$x[,1], pcafit.hp$x[,2])
colnames(df.hp) <- c("Individual", "PC1", "PC2")
plot.hp <- ggplot(df.hp, aes(x=PC1, y=PC2, color=Individual)) + geom_point(size=3) + xlab("PC1 (38% variance explained)") + ylab("PC2 (23% variance explained)") + theme(legend.position = "none")

dat.mb <- read.table("/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/mpileup_results/Mb_markerND_fv_R_gteq10.txt", header=T, sep="\t")
dat.mb.matrix <- as.matrix(select(dat.mb, -Individual))
pcafit.mb <- prcomp(dat.mb.matrix, center=T, scale.=T)
names(pcafit.mb)
var.mb <- pcafit.mb$sdev^2
var.mb.percent <- var.mb/sum(var.mb)
df.mb <- data.frame(dat.mb$Individual, pcafit.mb$x[,1], pcafit.mb$x[,2])
colnames(df.mb) <- c("Individual", "PC1", "PC2")
plot.mb <- ggplot(df.mb, aes(x=PC1, y=PC2, color=Individual)) + geom_point(size=3) + xlab("PC1 (36% variance explained)") + ylab("PC2 (26% variance explained)") 
#dat.ch.fs <- read.table("/mnt/blsm/sarah/OhData/metaphlan_results/mpileup_results/Ch/Ch_markerND_fv_shared50-10_gteq5_featureselection.txt", header=T, sep="\t")
#dat.ch.fs.matrix <- as.matrix(select(dat.ch.fs, -Individual))
#pcafit.ch.fs <- prcomp(dat.ch.fs.matrix, center=T, scale.=T)
#names(pcafit.ch.fs)
#var.ch.fs <- pcafit.ch.fs$sdev^2
#var.ch.fs.percent <- var.ch.fs/sum(var.ch.fs)
#df.ch.fs <- data.frame(dat.ch.fs$Individual, pcafit.ch.fs$x[,1], pcafit.ch.fs$x[,2])
#colnames(df.ch.fs) <- c("Individual", "PC1", "PC2")
#df.ch.fs$PC1 <- df.ch.fs$PC1*-1
#plot.ch.fs <- ggplot(df.ch.fs, aes(x=PC1, y=PC2, color=Individual)) + geom_point(size=3) + xlab("PC1 (40.19% variance explained)") + ylab("PC2 (19.96% variance explained)")

plot <- plot_grid(plot.all, plot.fb, plot.hp, plot.mb, nrow=2, ncol=2, labels=c("A", "B", "C", "D"), scale=1, rel_widths = c(1,1.17,1,1.17))

ggsave("/mnt/blsm/sarah/targetedstudy/analysis/metaphlan_results/alignstats/figures/PCA_cowplot.png", plot=plot, scale = 1, width=30, height=20, units=c("cm"),dpi=300, limitsize=T)
