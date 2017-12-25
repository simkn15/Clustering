library(transclustr)
library(ggplot2)
library(plyr)
library(cluster)
library(gtools)
library(ClusterR)
library(stringi)
# 
# tclustAndSetClasses <- function(simMatrix, threshold, proteins) {
#     tclustResult <- tclust(simmatrix = simMatrix, convert_dissimilarity_to_similarity = FALSE, threshold = threshold)
#     tclustResultDataFrame <- data.frame(protein = proteins, cluster = tclustResult$clusters[[1]])
#     ## Add 1 to all classes to match first class = 1 in gold standard
#     clusteringResultDataFrame = tclustResultDataFrame$cluster + 1
#     return(list(clusteringResultDataFrame, tclustResult$cost[[1]]))
# }

clusteringWithTclust <- function(simMatrix, proteinLabels, threshold) {
    tclustResult <- tclust(simmatrix = simMatrix, convert_dissimilarity_to_similarity = FALSE, threshold = threshold)
    tclustResultDataFrame <- data.frame(protein = proteinLabels, cluster = tclustResult$clusters[[1]])
    tclustResultDataFrame$cluster = tclustResultDataFrame$cluster + 1
    
    return(tclustResultDataFrame)
}

getProteinLabelsFromClustering <- function(tclustRes) {
    amountOfClasses <- length(table(tclustRes$cluster))
    
    # Initialize each cluster object
    proteinLabels <- list()
    for (i in 1:amountOfClasses) {
        proteinLabels[[i]] <- list(cid = "", proteins = c(), startIndex = 0, endIndex = 0, height = 0)
    }
    
    # Collect labels for all clusters
    for (i in 1:length(tclustRes$cluster)) {
        labelCurrentProtein <- tclustRes$protein[i] # change 'tclustRes$protein[i]' to 'i', to get protein number instead of label
        classCurrentProtein <- tclustRes$cluster[i]
        proteinLabels[[classCurrentProtein]]$proteins <- c(proteinLabels[[classCurrentProtein]]$proteins, c(as.character(labelCurrentProtein)))
    }
    
    return(proteinLabels)
}

buildSimilarityMatrix <- function(proteins, df_table) {
    if (length(proteins) < 2) { stop("Cannot build similarity matrix with less than 2 proteins")}
    # Create an empty matrix
    simMatrix = matrix(0, nrow = length(proteins), ncol = length(proteins))

    # Set row and column names for readability
    row.names(simMatrix) = proteins
    colnames(simMatrix) = proteins
    
    # Extract the pairwise similarities
    for(x in 1:(length(proteins)-1)){
        # Select all rows containing p1
        p1 = proteins[x]
        presec = df_table[df_table[,1] == p1 | df_table[,2] == p1, ]
        
        for(y in (x+1):length(proteins)){
            p2 = proteins[y]
            # We have to check for both directions p1 -> p2 and p2 <- p1. Rule is, to be more conservative, we take the minimum of both values
            if(length(presec[presec[,1] == p2, 3])==0 | length(presec[presec[,2] == p2, 3])==0){
                #Allright, one value was missing, we take 0 as fallback (since there is no meaningful minimum)
                simMatrix[x,y] = simMatrix[y,x] = 0
            }else{
                #Take the minimum of both directions
                simMatrix[x,y] = simMatrix[y,x] = min(presec[presec[,1] == p2, 3], presec[presec[,2] == p2, 3])
            }
        }
    }
    # Set the diagonal to the highest observed value (not actually necessary, since self-similarity is never queried)
    # diag(simMatrix) = max(simMatrix)
    
    return(simMatrix)
}