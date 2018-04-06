# GET MODULE GENES
# extract genes from modules and write them to files for enrichment analysis

# last modified: April, 2018.

# get list of all gene IDs
allGenes <- colnames(datExprs)

# getModuleGenes returns a list of gene IDs from a given module 
# geneList contains all genes in the expression matrix
# geneColours contains an array of module colour assignments for each gene in the expression matrix
# moduleColour is the module of interest's colour 
# note: geneColours may contain numbers instead of colours, in which case moduleColour must also be a number
getModuleGenes <- function(geneList,geneColours,moduleColour) {
  moduleGenes <- (geneColours==moduleColour)
  moduleIDs <- geneList[moduleGenes]
  return (moduleIDs)
}

# modules of interest
moduleNumbers <- c(12,7,10)
# for each module of interest, get IDs and write them to a file
# moduleLabelsAutomatic comes from building the network in WGCNA.R
path <- "/Users/maryhoekstra/Desktop/"
for (i in 1:length(moduleNumbers)) {
  moduleIDs <- getModuleGenes(geneList=allGenes,geneColours = moduleLabelsAutomatic,moduleColour = moduleNumbers[i])
  write(moduleIDs,file=paste(path,moduleNumbers[i],".txt",sep=""))
}
