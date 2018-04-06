# HEATMAPS
# create heatmap comparing expression from early to late timepoints in all samples 

# last modified: April 5th, 2018

library(gplots)
library(RColorBrewer)

# define early and late timepoints 
earlyDays <- c("day 1")
lateDays <- c("day 14")
earlyVsLateESet <- getESet(filtered.eset,earlyDays,lateDays)

# find DE probes
earlyLateVec <- character(ncol(earlyVsLateESet))
earlyIndices <- logical(ncol(earlyVsLateESet))
for (i in 1:length(earlyDays)) {
  earlyIndices <- earlyIndices | earlyVsLateESet$Day==earlyDays[i]
}
earlyLateVec[earlyIndices]="early"
lateIndices <- logical(ncol(earlyVsLateESet))
for (i in 1:length(lateDays)) {
  lateIndices <- lateIndices | earlyVsLateESet$Day==lateDays[i]
}
earlyLateVec[lateIndices]="late"
designVec <- factor(earlyLateVec)
designMatrix <- model.matrix(~0+designVec)
colnames(designMatrix) <- levels(designVec)
diffs <- runLimma(earlyVsLateESet,designMatrix,"early-late")
diffs$ID <- as.factor(rownames(diffs))
diffs <- merge(diffs, dplyr::select(annoData, ID, Entrez.Gene))

# take only the DE probes
topProbes <- as.character(diffs$ID)
exprsData <- exprs(earlyVsLateESet[topProbes,])
dayLabels <- sub("day 14","khaki1",earlyVsLateESet$Day)
dayLabels <- sub("day 1","mediumpurple",dayLabels)

# use superheat package
library(superheat)
heatmapTitle = "Differentially expressed genes between Day 1 and 14"
superheat(exprsData,
          clustering.method = "hierarchical",
          pretty.order.rows = TRUE,
          pretty.order.cols = TRUE,
          left.label.size = 0.05,
          bottom.label.size = 0.1,
          col.dendrogram = TRUE,
          row.dendrogram = FALSE,
          row.title = "Probes",
          column.title = "Samples",
          left.label.text.size = 2,
          bottom.label.text.size = 2,
          bottom.label.text.angle = 90,
          legend.breaks = c(5,7.5, 10,12.5, 15),
          title = heatmapTitle,
          bottom.label.col = dayLabels,
          grid.vline = FALSE)
