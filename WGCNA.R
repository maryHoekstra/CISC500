# WGCNA
# build gene network, extract MEs, and test MEs for DE 
# plot results using boxplots

# last modified: April 3rd, 2018

library(WGCNA)
library(genefilter)
library(cluster)
library(ggplot2)
library(reshape2)
options(stringsAsFactors = FALSE)

# get collapsed expression matrix from collapseProbes.R
datExprs <- t(collapsedExprs)
samples <- rownames(datExprs)

### determine soft-thresholding power to build network 
# choose set of soft-thresholding powers
candidatePowers = c(c(1:10), seq(from = 12, to=32, by=2))

# call network topology analysis function
softThreshold = pickSoftThreshold(datExprs,networkType="signed",corFnc="bicor", powerVector = candidatePowers, verbose = 5)

# plot results
sizeGrWindow(10,7)
par(mfrow = c(1,2))
cex1 = 0.9

# scale-free topology fit index as a function of the soft-thresholding power
plot(softThreshold$fitIndices[,1], -sign(softThreshold$fitIndices[,3])*softThreshold$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(softThreshold$fitIndices[,1], -sign(softThreshold$fitIndices[,3])*softThreshold$fitIndices[,2],
     labels=candidatePowers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")

# mean connectivity as a function of the soft-thresholding power
plot(softThreshold$fitIndices[,1], softThreshold$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(softThreshold$fitIndices[,1], softThreshold$fitIndices[,5], labels=candidatePowers, cex=cex1,col="red")

### find modules 
# we have 16GB laptop and can run all genes at once

# getModules takes an expression matrix and soft-thresholding power and returns a network
# note: reduce maxBlockSize if running on a computer with <16GB memory
getModules <- function(datExprs,sfPower) {
  bwnet = blockwiseModules(datExprs, maxBlockSize = 15000,
                           power = sfPower, networkType= "signed", TOMType = "signed", minModuleSize = 30,
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = TRUE,
                           corType = "bicor", maxPOutliers = 0.05,
                           verbose = 3)
  return (bwnet)
}

### comparing expression of module eigengenes
# use default soft-thresholding power of 12
net <- getModules(datExprs,12) 
# view number of modules and size of modules 
table(net$colors)
moduleLabelsAutomatic=net$colors
# convert labels to colors for plotting
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
# get data frame with module eigengenes 
MEsAutomatic=net$MEs

allDays <- filtered.eset$Day
allCodes <- filtered.eset$code

# create expression set with only early and late samples
early <- c("day 1")
late <- c("day 14")
earlyLate.eset <- getESet(filtered.eset,early,late)
groupA.eset = earlyLate.eset[,earlyLate.eset$code=="A"]
groupB.eset = earlyLate.eset[,earlyLate.eset$code=="B"]

# group A
title = "Expression of MEs in control group"
A <- which(allCodes=="A")
early <- which(allDays=="day 1")
late <- which(allDays=="day 14")
inds <- intersect(A,c(early,late))
timepointLabels <- sub("day 14","Late",groupA.eset$Day)
timepointLabels <- sub("day 1","Early",timepointLabels)

# group B
title = "Expression of MEs in lactoferrin group"
B <- which(allCodes=="B")
early <- which(allDays=="day 1")
late <- which(allDays=="day 14")
inds <- intersect(B,c(early,late))
timepointLabels <- sub("day 14","Late",groupB.eset$Day)
timepointLabels <- sub("day 1","Early",timepointLabels)

# select the MEs from group A or group B 
MEsAutomatic=net$MEs
MEsAutomatic <- MEsAutomatic[inds,] # use inds from groupA OR groupB
# optionally, only look at MEs of interest
MEsAutomatic <- MEsAutomatic[,c("ME7","ME10","ME12")]
MEsAutomatic$Timepoint <- factor(timepointLabels) # use timepoint labels from groupA OR groupB
MEsAutomatic$Timepoint <- factor(groupA.eset$Day) # use groupA.eset OR groupB.eset

# create boxplots
df <- melt(MEsAutomatic)
ggplot(data=df) + geom_boxplot(aes(x=Timepoint,y=value)) + 
  facet_wrap(~variable,scales = "free") + 
  ggtitle(title) + scale_fill_brewer(palette = "Accent") + 
  ylab("Value") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))

# identify DE in MEs using t-tests
numModules <- length(MEsAutomatic)-1
ctt <- colttests(data.matrix(MEsAutomatic),MEsAutomatic$Timepoint,tstatOnly = FALSE)

# plot p-values
title = "DE of MEs between Day 1 and Day 14 - Control"
barplot(ctt$p.value[1:numModules],names.arg = rownames(ctt)[1:numModules],ylim = c(0,1),main = title,xlab = "Module Eigengenes",ylab = "P-values")
abline(h=0.05,col="red")
legend(x = "topright",legend = "P = 0.05",lty = 1,col = "red")

