# COLLAPSE PROBES
# reduce expression matrix of probes to expression matrix of genes

# last modified: August, 2017

library(WGCNA)

# collapseProbes returns a collapsed expression matrix with expression values for genes
# filteredExprs is the expression matrix that's been filtered by variance 
# annoData is the dataframe containing the mapping of probe IDs to gene IDs 
collapseProbes <- function(filteredExprs,annoData) {
  
  # find indices of our chosen probes in the original annotation matrix
  probeIndices <- match(rownames(filteredExprs),annoData$ID)
  
  #find indices of our chosen probes in the expanded annotation matrix 
  allIDS <- as.character(annoData$Entrez.Gene)
  # get the gene symbols which correspond to the chosen probes
  IDsubset <- allIDS[probeIndices]
  
  # choose probe with highest mean expression when multiple probes map to a gene
  # use only the symbols from the filtered expression set
  collapsedData <- collapseRows(filteredExprs,IDsubset,rownames(filteredExprs),method="MaxMean")
  collapsedExprs <- collapsedData$datETcollapsed
  
  # remove row with "---" as gene symbol
  collapsedExprs <- collapsedExprs[-1,]
  
  return(collapsedExprs) 
}

# get filteredExprs and annoData from running preprocessData.R
collapsedExprs <- collapseProbes(filteredExprs,annoData)
