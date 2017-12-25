library(transclustr)
library(ggplot2)
library(plyr)
library(cluster)
library(gtools)
library(ClusterR)
library(stringi)

#library(pryr) # get size of objects with object_size(l)

source("fmeasure.R")
source("utilities.R")
####################################################################################################
# Read in data and build similarity matrix
####################################################################################################
table <- read.table("sfld_brown_et_al_amidohydrolases_protein_similarities_for_beh.txt", sep = "", skip = 5)
df_table <- as.data.frame(table)

# Get the involved proteins
proteins = levels(df_table[,1])

# Build similarity matrix
simMatrix <- buildSimilarityMatrix(proteins, df_table)

####################################################################################################
# Cluster
####################################################################################################
######
# Hierarchical clustering
######
# hcObject <- list(merge = matrix(0, nrow = 0, ncol = 2), height = c(), order = c(), labels = proteins)

maxThreshold <- 324
threshold <- -1

proteinLabelsOfEachCluster <- list(c1 = list(cid = 1, proteins = c(proteins), startIndex = 1, endIndex = length(proteins), height = -1))

amountOfTotalClusters <- 1
totalClusters <- proteinLabelsOfEachCluster

countSplits <- 0
countSingletons <- 0

charC <- "c"

# hc object
merge <- list()
height <- c() # maxThreshold - threshold
order <- proteins

# while (length(countSingletons) < length(proteins)) {
while (threshold < 2) {
    cat("\n")
    print("##############################")
    print(paste("Threshold =", threshold))
    print(paste("Clusters to tclust =", length(proteinLabelsOfEachCluster)))
    # print(names(proteinLabelsOfEachCluster)) # NULL, since we are losing names. Need unique names for this to work
    
    clustersLargerThanOne <- list()
    
    # For every level: Do tclust on all clusters with current threshold
    for (i in 1:length(proteinLabelsOfEachCluster)) {
        currentTclustCluster <- proteinLabelsOfEachCluster[[i]]# Cluster to tclust
        
        if (length(currentTclustCluster$proteins) > 1) { # Not a singleton
            
            print(paste0("currentTclustCluster: cid = ", currentTclustCluster$cid, " | size = ", length(currentTclustCluster$proteins), " | startIndex = ", currentTclustCluster$startIndex, " | endIndex = ", currentTclustCluster$endIndex, " | height = ", currentTclustCluster$height))
            print(currentTclustCluster$proteins)
            
            simMatrix <- buildSimilarityMatrix(currentTclustCluster$proteins, df_table)
            
            tclustResultDataFrame <- clusteringWithTclust(simMatrix, currentTclustCluster$proteins, threshold)
            amountOfNewClusters <- length(table(tclustResultDataFrame$cluster))

            # Get labels for each new cluster
            # TODO check if clusters == 1, no need to getLabels, since currentTclustCluster holds all labels of cluster. Do the check inside the function
            if (amountOfNewClusters > 1) {
                tempProteinLabelsOfEachCluster <- getProteinLabelsFromClustering(tclustResultDataFrame)
                for (j in 1:amountOfNewClusters) { # Update height on cluster object
                    # Update height
                    tempProteinLabelsOfEachCluster[[j]]$height <- maxThreshold - threshold + 1
                    # Update cid
                    amountOfTotalClusters <- amountOfTotalClusters + 1
                    cid <- paste0(charC, amountOfTotalClusters)
                    tempProteinLabelsOfEachCluster[[j]]$cid <- cid
                }
            }
            else {
                tempProteinLabelsOfEachCluster <- getProteinLabelsFromClustering(tclustResultDataFrame)
                # Update height
                tempProteinLabelsOfEachCluster[[1]]$height <- currentTclustCluster$height
                # Update cid
                tempProteinLabelsOfEachCluster[[1]]$cid <- currentTclustCluster$cid
            }
            
            # Update order
            currentStartIndex <- currentTclustCluster$startIndex # Set start of first new cluster = start of parent cluster
            for (i in 1:amountOfNewClusters) { # For every new cluster in parent cluster
                labelsI <- tempProteinLabelsOfEachCluster[[i]]$proteins
                lengthLabelsI <- length(labelsI)
                
                startLabelsI <- currentStartIndex
                endLabelsI <- currentStartIndex + lengthLabelsI - 1
                tempProteinLabelsOfEachCluster[[i]]$startIndex <- startLabelsI
                tempProteinLabelsOfEachCluster[[i]]$endIndex <- endLabelsI
                
                for (j in 1:lengthLabelsI) { # Update positions in order
                    order[startLabelsI + j - 1] <- tempProteinLabelsOfEachCluster[[i]]$proteins[j]
                }
                # Done with labelsI, update currentStartIndex so start at labelsI.endIndex + 1
                currentStartIndex <- endLabelsI + 1
            }
            
            # Add clusters for next iteration
            if (amountOfNewClusters == 1) {
                print(paste0("There was no split, threshold = ", threshold))
                cat("\n")
                tempCurrentCluster <- tempProteinLabelsOfEachCluster[[1]]
                clustersLargerThanOne <- c(clustersLargerThanOne, list(tempCurrentCluster)) # We lose name at this line
            }
            else {
                print(paste0("Was split into ", length(table(tclustResultDataFrame$cluster)), " clusters, threshold = ", threshold))
                
                height <- c(height, c(maxThreshold - threshold + 1))
                countSplits <- countSplits + 1
                
                currentMerge <- c()
                for (j in 1:amountOfNewClusters) { # For each new cluster
                    tempCurrentCluster <- tempProteinLabelsOfEachCluster[[j]]
                    
                    # Update totalClusters
                    tempCurrentClusterList <- list(tempCurrentCluster)
                    names(tempCurrentClusterList) <- tempCurrentCluster$cid
                    totalClusters <- c(totalClusters, tempCurrentClusterList)
                    
                    currentMerge <- c(currentMerge, c(tempCurrentCluster$cid))
                    
                    if (length(tempCurrentCluster$proteins) == 1) { # singleton found
                        countSingletons <- countSingletons + 1
                        # clustersLargerThanOne <- c(clustersLargerThanOne, list(tempCurrentCluster))
                    }
                    if (length(tempCurrentCluster$proteins) > 1) { # Add cluster for next threshold
                        clustersLargerThanOne <- c(clustersLargerThanOne, list(tempCurrentCluster)) # We lose name at this line
                    }
                    # cat("\n")
                    # print(paste0("Cluster: size = ", length(tempCurrentCluster$proteins)))
                    # print(tempCurrentCluster)
                }
                merge <- c(merge, list(currentMerge))
            }
        }
    }
    
    proteinLabelsOfEachCluster <- clustersLargerThanOne
    
    threshold <- threshold + 1
}
print(paste0("Total splits: ", countSplits, " | length(height): ", length(height)))
as.numeric(levels(as.factor(height)))
print(paste0("Total clusters = ", amountOfTotalClusters))
