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
# table <- read.table("sfld_brown_et_al_amidohydrolases_protein_similarities_for_beh.txt", sep = "", skip = 5)
# df_table <- as.data.frame(table)
# 
# # Get the involved proteins
# proteins = levels(df_table[,1])
# 
# # Build similarity matrix
# simMatrix <- buildSimilarityMatrix(proteins, df_table)

table <- read.table("simBig.txt", sep = "")
df_table <- as.data.frame(table)

# Get the involved proteins
proteins = levels(df_table[,1])

# Build similarity matrix
# simMatrix <- buildSimilarityMatrix(proteins, df_table)
simMatrix <- buildSimilarityMatrixFromBlast(proteins, df_table)

####################################################################################################
# Hierarchical clustering
####################################################################################################
# hcObject <- list(merge = matrix(0, nrow = 0, ncol = 2), height = c(), order = c(), labels = proteins)

maxThreshold <- round(max(simMatrix)) + 1
threshold <- -1

proteinLabelsOfEachCluster <- list(c1 = list(cid = "c1", proteins = c(proteins), startIndex = 1, endIndex = length(proteins), height = -1, parent = ""))

amountOfTotalClusters <- 1
totalClusters <- proteinLabelsOfEachCluster

countSplits <- 0
countSingletons <- 0

charC <- "c"

# hc object
mergeMatrix <- matrix(0, nrow = 0, ncol = 2)
merge <- list()
height <- c() # maxThreshold - threshold
order <- proteins
labels <- proteins
step <- 30
while (countSingletons < length(proteins)) {
# while (threshold < 2) {
    cat("\n")
    print("##############################")
    print(paste("Threshold =", threshold))
    print(paste("Clusters to tclust =", length(proteinLabelsOfEachCluster)))
    # print(names(proteinLabelsOfEachCluster)) # NULL, since we are losing names. Need unique names for this to work
    
    clustersLargerThanOne <- list()
    
    # For every level: Do tclust on all clusters with current threshold
    if (length(proteinLabelsOfEachCluster) > 0) {
        for (i in 1:length(proteinLabelsOfEachCluster)) {
            currentTclustCluster <- proteinLabelsOfEachCluster[[i]]# Cluster to tclust
            
            if (length(currentTclustCluster$proteins) > 1) { # Not a singleton
                
                print(paste0("currentTclustCluster: cid = ", currentTclustCluster$cid, " | size = ", length(currentTclustCluster$proteins), " | startIndex = ", currentTclustCluster$startIndex, " | endIndex = ", currentTclustCluster$endIndex, " | height = ", currentTclustCluster$height, " | parent = ", currentTclustCluster$parent))
                print(currentTclustCluster$proteins)
                
                simMatrix <- buildSimilarityMatrix(currentTclustCluster$proteins, df_table)

                tclustResultDataFrame <- clusteringWithTclust(simMatrix, currentTclustCluster$proteins, threshold)
                amountOfNewClusters <- length(table(tclustResultDataFrame$cluster))
                
                # Get labels for each new cluster
                # TODO check if clusters == 1, no need to getLabels, since currentTclustCluster holds all labels of cluster. Do the check inside the function
                if (amountOfNewClusters > 1) {
                    tempProteinLabelsOfEachCluster <- getProteinLabelsFromClustering(tclustResultDataFrame)
                    for (j in 1:amountOfNewClusters) { # Update height on cluster object
                        tempProteinLabelsOfEachCluster[[j]]$height <- maxThreshold - threshold + 1
                        amountOfTotalClusters <- amountOfTotalClusters + 1
                        cid <- paste0(charC, amountOfTotalClusters)
                        tempProteinLabelsOfEachCluster[[j]]$cid <- cid
                        tempProteinLabelsOfEachCluster[[j]]$parent <- currentTclustCluster$cid
                    }
                }
                else { # Same cluster as currentTclustCluster
                    tempProteinLabelsOfEachCluster <- getProteinLabelsFromClustering(tclustResultDataFrame)
                    # Update height
                    tempProteinLabelsOfEachCluster[[1]]$height <- currentTclustCluster$height
                    # Update cid
                    tempProteinLabelsOfEachCluster[[1]]$cid <- currentTclustCluster$cid
                    tempProteinLabelsOfEachCluster[[1]]$parent <- currentTclustCluster$parent
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
    }
    
    proteinLabelsOfEachCluster <- clustersLargerThanOne
    
    threshold <- threshold + step
}
print(paste0("Total splits: ", countSplits, " | length(height): ", length(height)))
as.numeric(levels(as.factor(height)))
print(paste0("Total clusters = ", amountOfTotalClusters))

singletons <- c()
clusters <- c()
for (i in 1:length(totalClusters)) {
    # clust <- paste0("c", i)
    if (length(totalClusters[[i]]$proteins) == 1) {singletons <- c(singletons, totalClusters[[i]]$cid)}
    else {clusters <- c(clusters, totalClusters[[i]]$cid)}
}

####################################################################################################
# hc$merge
####################################################################################################
mergeMatrix <- matrix(0, nrow = 0, ncol = 2)
mergeLookUpList <- list()
mergeHeights <- c()

which(proteins == totalClusters$c324$proteins[[1]])

# Move in descending order!
for (i in length(merge):1) { # For all parents
    parent <- ""
    for (j in 2:length(merge[[i]])) { # For all children of parent i
        # TODO Simplify!
        # Think about if (j == 1) && (j > 1). We are doing pretty much the same thing
        
        mergeHeights <- c(mergeHeights, height[i])
        
        if (j == 2) { # initialize run -> a cluster cannot be on previous line
            c1 <- merge[[i]][j-1]
            c2 <- merge[[i]][j]
            c1.isSingleton <- length(totalClusters[[c1]]$proteins) == 1
            c2.isSingleton <- length(totalClusters[[c2]]$proteins) == 1
            c1.isCluster <- !c1.isSingleton
            c2.isCluster <- !c2.isSingleton
            parent <- totalClusters[[c1]]$parent
            
            if (c1.isSingleton && c2.isSingleton) { # (s,s)
                rowLabelsC1 <- - which(labels == totalClusters[[c1]]$proteins[1])
                rowLabelsC2 <- - which(labels == totalClusters[[c2]]$proteins[1])
                mergeMatrix <- rbind(mergeMatrix, c(rowLabelsC1, rowLabelsC2))
            }
            
            if (c1.isCluster && c2.isSingleton) { # (c,s)
                mergeMatrixRowOfC1 <- mergeLookUpList[[c1]]
                rowLabelsC2 <- - which(labels == totalClusters[[c2]]$proteins[1])
                mergeMatrix <- rbind(mergeMatrix, c(mergeMatrixRowOfC1, rowLabelsC2))
            }
            
            if (c1.isSingleton && c2.isCluster) { # (s,c)
                mergeMatrixRowOfC2 <- mergeLookUpList[[c2]]
                rowLabelsC1 <- - which(labels == totalClusters[[c1]]$proteins[1])
                mergeMatrix <- rbind(mergeMatrix, c(rowLabelsC1, mergeMatrixRowOfC2))
            }
            
            if (c1.isCluster && c2.isCluster) { # (c,c)
                mergeMatrixRowOfC1 <- mergeLookUpList[[c1]]
                mergeMatrixRowOfC2 <- mergeLookUpList[[c2]]
                mergeMatrix <- rbind(mergeMatrix, c(mergeMatrixRowOfC1, mergeMatrixRowOfC2))
            }
        }
        if (j > 2) { # Merge next element with cluster from previous line in mergeMatrix
            c1 <- nrow(mergeMatrix)
            c2 <- merge[[i]][j]
            c1.isCluster <- TRUE
            c2.isSingleton <- length(totalClusters[[c2]]$proteins) == 1
            c2.isCluster <- !c2.isSingleton
            
            if (c1.isCluster && c2.isSingleton) {
                mergeMatrixRowOfC1 <- c1
                rowLabelsC2 <- - which(labels == totalClusters[[c2]]$proteins[1])
                mergeMatrix <- rbind(mergeMatrix, c(mergeMatrixRowOfC1, rowLabelsC2))
            }
            
            if (c1.isCluster && c2.isCluster) {
                mergeMatrixRowOfC1 <- c1
                mergeMatrixRowOfC2 <- mergeLookUpList[[c2]]
                if (is.null(mergeMatrixRowOfC2)) {
                    stop("j > 2, c1 = c2 = isCluster. c2 was not found in mergeLookUpList. All clusters must be singletons at this state, otherwise this error can occur.")
                }
                mergeMatrix <- rbind(mergeMatrix, c(mergeMatrixRowOfC1, mergeMatrixRowOfC2))
            }
        }
    }
    # Update mergeLookUpList
    newCluster <- list(nrow(mergeMatrix)) # Make a counter, this is too many ops
    names(newCluster) <- parent # The children have been marge into its parent.
    mergeLookUpList <- c(mergeLookUpList, newCluster)
}

# Can't plot without making a hc object and reassign
require(graphics)
hc <- hclust(dist(USArrests), "ave")
hc$merge <- mergeMatrix
hc$height <- mergeHeights
hc$order <- order
hc$labels <- labels
plot(hc, xlab = "protein", ylab = "threshold")
plot(hc, hang = -1, xlab = "Protein", ylab = "Threshold")

# Fails
# hct <- list(merge = mergeMatrix, height = height, order = order, labels = labels)
# plot(hct, xlab = "Protein", ylab = "Threshold")
# plot(hct, hang = -1)

# hclust original dataset
# hc <- hclust((max(as.dist(simMatrix)) + 1) - as.dist(simMatrix))

# Approx optimal threshold
# > mean(as.vector(as.dist(simMatrix)))
#[1] 55.91576


# multidimensional scaling R -> Google
# mydata <- (max(as.dist(simMatrix)) + 1) - as.dist(simMatrix)
# d <- dist((max(as.dist(simMatrix)) + 1) - as.dist(simMatrix))



