# EARLY VS LATE
# find DE between early vs. late samples in group A and group B

# last modified: January 15th, 2018

# getESet subsets an expression set based on samples from specific days
# earlyDays and lateDays are arrays containing the desired timepoints, ex. c("day 1","day 3")
getESet <- function(filteredESet,earlyDays,lateDays) {
  
  # create expression set with only early and late samples 
  earlyLateDays = c(earlyDays,lateDays)
  earlyLateIndices <- logical(ncol(filteredESet))
  for (i in 1:length(earlyLateDays)) {
    earlyLateIndices <- earlyLateIndices | filteredESet$Day==earlyLateDays[i]
  }
  earlyLate.eset <- filteredESet[,earlyLateIndices]
  return (earlyLate.eset)
}

# createDesign forms a design matrix using the smaller expression set
# earlyDays and lateDays are arrays containing the desired timepoints, ex. c("day 1","day 3")
createDesign <- function(earlyLateESet,earlyDays,lateDays) {  
  earlyLateVec <- character(ncol(earlyLateESet))
  earlyIndices <- logical(ncol(earlyLateESet))
  for (i in 1:length(earlyDays)) {
    earlyIndices <- earlyIndices | earlyLateESet$Day==earlyDays[i]
  }
  earlyLateVec[earlyIndices]="early"
  lateIndices <- logical(ncol(earlyLateESet))
  for (i in 1:length(lateDays)) {
    lateIndices <- lateIndices | earlyLateESet$Day==lateDays[i]
  }
  earlyLateVec[lateIndices]="late"
  # create design matrix with the four groups we are interested in (A.early, A.late, etc) 
  designVec <- factor(paste(earlyLateESet$code,earlyLateVec,sep="."))
  designMatrix <- model.matrix(~0+designVec)
  colnames(designMatrix) <- levels(designVec)
  return (designMatrix)
}

# getDiffs gets a table of top DE probes for a given group from limma analysis
# groupCode is a character and design is an nx4 matrix 
# note: "runLimma" function is found in limmaTest.R
getDiffs <- function(earlyLateESet,groupCode,design) {
  # get diffs
  contrastsString <- paste(groupCode,".early-",groupCode,".late",sep="")
  diffs <- runLimma(earlyLateESet,design,contrastsString)
  # add entrez gene IDs to diffs table 
  diffs$ID <- as.factor(rownames(diffs))
  diffsWithAnno <- merge(diffs, dplyr::select(annoData, ID, Entrez.Gene))
  return(diffsWithAnno)
}

# getEntrezIDs returns a list of entrez IDs from the annotated diffs table 
getEntrezIDs <- function(diffsWithAnno) {
  geneEntrez <- dplyr::select(diffsWithAnno, Entrez.Gene)
  geneEntrez <- unique(as.character(geneEntrez[,1]))
  return (geneEntrez)
}

# find DE genes between Day 1 and Day 14 
earlyDays <- c("day 1")
lateDays <- c("day 14")
earlyVsLateESet <- getESet(filtered.eset,earlyDays,lateDays)
design <- createDesign(earlyVsLateESet,earlyDays,lateDays)

# specify where to write gene lists
path = "/Users/maryhoekstra/Desktop/"

# get gene list from group A
diffsA <- getDiffs(earlyVsLateESet,groupCode="A",design)
geneListA <- getEntrezIDs(diffsA)
write(geneListA,file=paste(path,"A_day1vs14_entrez.txt",sep=""))

# get gene list from group B
diffsB <- getDiffs(earlyVsLateESet,groupCode="B",design)
geneListB <- getEntrezIDs(diffsB)
write(geneListB,file=paste(path,"B_day1vs14_entrez.txt",sep=""))

# determine which genes are common to both groups
commonGenes <- intersect(geneListA,geneListB)
write(commonGenes,file=paste(path,"common.txt",sep=""))

# determine which genes are unique to A or B
uniqueToA <- setdiff(geneListA,geneListB)
write(uniqueToA,file=paste(path,"uniqueA.txt",sep=""))
uniqueToB <- setdiff(geneListB,geneListA)
write(uniqueToB,file=paste(path,"uniqueB.txt",sep=""))
