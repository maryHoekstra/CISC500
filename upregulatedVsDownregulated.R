# UPREGULATED VS DOWNREGULATED
# determine which genes are upregulated or downregulated for a given diffs table 

# last modified: June 2017 

# subsetByRegulation returns a list of genes that are upregulated or downregulated
# diffsWithAnno is the annotated diffs table from earlyVsLate.R
# directionString is either "up" or "down"
subsetByRegulation <- function(diffsWithAnno,directionString) {
  if (directionString=="up") {
    diffs <- diffsWithAnno[diffsWithAnno$logFC>0,]
  }
  else {
    diffs <- diffsWithAnno[diffsWithAnno$logFC<0,]
  }
  
  geneList <- getEntrezIDs(diffs)
  return(geneList)
}

# these commands assume that diffsA and diffsB have already been created in earlyVsLate.R
upregulatedGeneListA <- subsetByRegulation(diffsA,"up")
upregulatedGeneListB <- subsetByRegulation(diffsB,"up")
downregulatedGeneListA <- subsetByRegulation(diffsA,"down")
downregulatedGeneListB <- subsetByRegulation(diffsB,"down")

# find interest
upregulatedIntersect <- intersect(upregulatedGeneListA,upregulatedGeneListB)
downregulatedIntersect <- intersect(downregulatedGeneListA,downregulatedGeneListB)

# find unique genes 
downregulatedGeneListB_unique <- setdiff(downregulatedGeneListB,downregulatedGeneListA)
upregulatedGeneListB_unique <- setdiff(upregulatedGeneListB,upregulatedGeneListA)
downregulatedGeneListA_unique <- setdiff(downregulatedGeneListA,downregulatedGeneListB)
upregulatedGeneListA_unique <- setdiff(upregulatedGeneListA,upregulatedGeneListB)

# specify where to write gene lists
path = "/Users/maryhoekstra/Desktop/"

write(downregulatedGeneListA_unique,file=paste(path,"downregulatedGeneListA_unique.txt",sep=""))
write(downregulatedGeneListB_unique,file=paste(path,"downregulatedGeneListB_unique.txt",sep=""))
write(downregulatedIntersect,file=paste(path,"downregulatedIntersect.txt",sep=""))
write(upregulatedGeneListA_unique,file=paste(path,"upregulatedGeneListA_unique.txt",sep=""))
write(upregulatedGeneListB_unique,file=paste(path,"upregulatedGeneListB_unique.txt",sep=""))
write(upregulatedIntersect,file=paste(path,"upregulatedIntersect.txt",sep=""))

