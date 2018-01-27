# Still has DEBUG
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(transclustr)
library(ggplot2)
library(plyr)
library(cluster)
library(gtools)
library(ClusterR)
library(stringi)
require(graphics)

source("utilities.R")

# NOTE!
# The Hierarchical Clustering approach is divisive.
# The information for the dendrogram is converted into agglomerative, thus being able to use the plotting defined for a hclust() object
# Required to make a hclust() object, doing a clustering with the method, and reassigning the variables on the object.

####################################################################################################
# Hierarchical Method
####################################################################################################
# Returns a hclust() object: ?hclust
# simMatrix: Similarity matrix
# proteins: All protein names
# step: The steps of thresholds used for clustering
# minSplit: Minimum number of splits when doing a binary search
# maxSplit: Maximum number of splits when doing a binary search
# minTreshold: The starting threshold, gets incremented by 'step'
# maxThreshold: The maximum threshold where all clusters are singletons, max(simMatrix) + 1
# binarySearch: TRUE will enable the binary search
# DEBUG and fileName: Used for debugging
hclustDivisive <- function(simMatrix, proteins, step = 1, minSplit = 2, maxSplit = 10, minThreshold, maxThreshold, binarySearch = FALSE, DEBUG = FALSE, fileName = "") {
    if (DEBUG)  {
        if (fileName == "") { stop("DEBUG mode requires a filename for the output.") } 
        fileName <- paste0(fileName, minSplit, "Max", maxSplit, "Step", step, ".txt")
        file <- file(fileName, "w")
    }
    if (binarySearch && step == 1) { stop("Binary search requires step >= 2.") }
    
    # List of all clusters that needs to be clustered with tclust
    # cid: Cluster id
    # proteins: Proteins in the given cluster
    # startIndex, endIndex: The range in 'order' for the proteins belonging to this cluster. 
    # height: The height where the children clusters would be merged into this cluster. Meaning that this cluster is split into its children at height-1.
    # nextThreshold: Next threshold to run tclust on this cluster.
    proteinLabelsOfEachCluster <- list(c1 = list(cid = "c1", proteins = c(proteins), startIndex = 1, endIndex = length(proteins), height = -1, parent = "", nextThreshold = minThreshold))
    
    # Total number of clusters which has occured during the run. (Includes singleton clusters). Used to get proper cid on all clusters.
    amountOfTotalClusters <- 1
    
    # Holds all information for every cluster during the run. totalClusters$c1 is information about the first cluster. Required to make the merge ordering for the dendrogram.
    totalClusters <- proteinLabelsOfEachCluster
    
    countSplits <- 0 # Total number of splits for the entire run
    countSingletons <- 0 # Equals the amount of proteins when run is done
    
    charC <- "c" # used to make var cid
    
    # Variables used to create the dendrogram with a hclust() object
    mergeMatrix <- matrix(0, nrow = 0, ncol = 2)
    merge <- list() # Hold the order of all splits. merge[[1]] returns which clusters was made from the first split.
    height <- c() # maxThreshold - threshold
    order <- proteins # Holds the ordered proteins for the dendrogram
    
    labels <- proteins # The labels of all proteins
    while (countSingletons < length(proteins)) {
        if(DEBUG) {
            cat("\n")
            print("##############################")
            print(paste("Clusters to tclust =", length(proteinLabelsOfEachCluster)))   
        }
        
        clustersLargerThanOne <- list() # All clusters returned from this iteration of tclust(Not singletons). These will be added to 'proteinLabelsOfEachCluster' to be clustered in the next iteration.
        
        # For every level: Do tclust on all clusters with current threshold
        if (length(proteinLabelsOfEachCluster) > 0) {
            for (i in 1:length(proteinLabelsOfEachCluster)) {
                currentStep <- step
                currentTclustCluster <- proteinLabelsOfEachCluster[[i]]# Cluster to tclust
                
                if (DEBUG) {
                    print(paste("currentTclustCluster$nextThreshold =", currentTclustCluster$nextThreshold))
                    write("", file, append = TRUE)
                    string <- paste0("currentTclustCluster = ", currentTclustCluster$cid, " | nextThreshold = ", currentTclustCluster$nextThreshold)
                    write(string, file, append = TRUE) 
                }
                
                if (length(currentTclustCluster$proteins) > 1) { # Not a singleton. This check might not be necessary, as we dont add singletons for the next iteration
                    if (DEBUG) {
                        print(paste0("currentTclustCluster: cid = ", currentTclustCluster$cid, " | size = ", length(currentTclustCluster$proteins), " | startIndex = ", currentTclustCluster$startIndex, " | endIndex = ", currentTclustCluster$endIndex, " | height = ", currentTclustCluster$height, " | parent = ", currentTclustCluster$parent))   
                    }
                    
                    simMatrix <- buildSimilarityMatrix(currentTclustCluster$proteins, df_table)
                    
                    tclustResultDataFrame <- clusteringWithTclust(simMatrix, currentTclustCluster$proteins, currentTclustCluster$nextThreshold)
                    amountOfClustersInTclustResult <- length(table(tclustResultDataFrame$cluster))
                    
                    # Get labels for each new cluster
                    # TODO check if clusters == 1, no need to getLabels, since currentTclustCluster holds all labels of cluster. Do the check inside the function
                    fallBack <- FALSE # If binary search did not minimize the number of splits, reset to the initial threshold/result
                    if (amountOfClustersInTclustResult > 1) {
                        
                        ### BEGIN binary search
                        if (binarySearch && amountOfClustersInTclustResult >= maxSplit) { # Too many new clusters. Can't do binary search if base step is 1
                            tempTclustResultDataFrame <- tclustResultDataFrame
                            tempAmountOfClustersInTclustResult <- amountOfClustersInTclustResult
                            upperBound <- step
                            lowerBound <- 0
                            
                            if (DEBUG) { 
                                string <- paste0("Beginning binary search with cluster: ", currentTclustCluster$cid, " | amountOfClustersInTclustResult = ", amountOfClustersInTclustResult)
                                print(string)
                                write(string, file, append = TRUE) 
                            }
                            
                            while(amountOfClustersInTclustResult >= maxSplit) {
                                currentStep <- lowerBound + ceiling(abs(upperBound - lowerBound) / 2)
                                # print(paste0(currentTclustCluster$cid, " | Decreasing step -> ", currentStep))
                                
                                if (DEBUG) { 
                                    string <- paste0("Decreasing to step = ", currentStep, " | lowerBound = ", lowerBound, " | upperBound = ", upperBound)
                                    print(string)
                                    write(string, file, append = TRUE) 
                                }
                                
                                tclustResultDataFrame <- clusteringWithTclust(simMatrix, currentTclustCluster$proteins, currentTclustCluster$nextThreshold - step + currentStep)
                                amountOfClustersInTclustResult <- length(table(tclustResultDataFrame$cluster))
                                
                                if (DEBUG) { 
                                    string <- paste0("amountOfClustersInTclustResult = ", amountOfClustersInTclustResult, " | currentStep = ", currentStep)
                                    print(string)
                                    write(string, file, append = TRUE) 
                                }
                                
                                while (amountOfClustersInTclustResult < minSplit && currentStep > 1 && currentStep != step) {
                                    if (DEBUG) { 
                                        string <- paste0("Inner loop | Updating lowerBound from ", lowerBound, " to ", currentStep)
                                        print(string)
                                        write(string, file, append = TRUE) 
                                    }
                                    
                                    lowerBound <- currentStep
                                    currentStep <- lowerBound + ceiling(abs(upperBound - lowerBound) / 2)
                                    
                                    if (DEBUG) { 
                                        string <- paste0("Increasing to step = ", currentStep, " | lowerBound = ", lowerBound, " | upperBound = ", upperBound)
                                        print(string)
                                        write(string, file, append = TRUE) 
                                    }
                                    
                                    if (currentStep == upperBound) { # Did not find a better threshold
                                        if (DEBUG) {
                                            string <- paste0("fallBack mode at inner, currentStep ", currentStep, " == ", upperBound, " upperBound")
                                            print(string)
                                            write(string, file, append = TRUE)
                                        }
                                        
                                        fallBack = TRUE
                                        break
                                    }
                                    
                                    tclustResultDataFrame <- clusteringWithTclust(simMatrix, currentTclustCluster$proteins, currentTclustCluster$nextThreshold - step + currentStep)
                                    amountOfClustersInTclustResult <- length(table(tclustResultDataFrame$cluster))
                                    
                                    if (DEBUG) { 
                                        string <- paste0("amountOfClustersInTclustResult = ", amountOfClustersInTclustResult, " | currentStep = ", currentStep)
                                        print(string)
                                        write(string, file, append = TRUE)
                                    }
                                }
                                
                                if (amountOfClustersInTclustResult >= maxSplit && (upperBound - lowerBound == 1)) {
                                    if (DEBUG) {
                                        string <- paste0("fallBack mode at outer, currentStep ", currentStep, " == 1", " | amountOfClustersInTclustResult >= ", maxSplit)
                                        print(string)
                                        write(string, file, append = TRUE)
                                    }
                                    
                                    fallBack = TRUE
                                }
                                if (fallBack) {
                                    break                            
                                }
                                if (amountOfClustersInTclustResult >= maxSplit) { # currentStep did not descrease the amount of clusters enough
                                    if (DEBUG) {
                                        string <- paste0("Outer loop | amountOfClustersInTclustResult >= ", maxSplit, " | Updating upperBound from ", upperBound, " to ", currentStep)
                                        print(string)
                                        write(string, file, append = TRUE)
                                    }
                                    
                                    upperBound <- currentStep
                                }
                            }
                        }
                        
                        if (fallBack) {
                            if (DEBUG) {
                                string <- paste0("fallBack mode caught")
                                print(string)
                                write(string, file, append = TRUE)
                                string <- paste0("amountOfClustersInTclustResult ", amountOfClustersInTclustResult, " > ", tempAmountOfClustersInTclustResult, " tempAmountOfClustersInTclustResult")
                                write(string, file, append = TRUE)
                            }
                            
                            
                            if ((amountOfClustersInTclustResult > tempAmountOfClustersInTclustResult) || 
                                (amountOfClustersInTclustResult == 1 && (amountOfClustersInTclustResult < tempAmountOfClustersInTclustResult))) {
                                # Either the first tclustResult was better, or the first tclustResult gave the only split
                                if (DEBUG) {
                                    string <- paste0("######## Resetting into pre-binary-search")
                                    print(string)
                                    write(string, file, append = TRUE)
                                }
                                
                                tclustResultDataFrame <- tempTclustResultDataFrame
                                amountOfClustersInTclustResult <- tempAmountOfClustersInTclustResult
                            }
                            
                            if (DEBUG) {
                                string <- paste0("amountOfClustersInTclustResult = ", amountOfClustersInTclustResult)
                                print(string)
                                write(string, file, append = TRUE) 
                            }
                        }
                        ### END binary search
                        
                        tempProteinLabelsOfEachCluster <- getProteinLabelsFromClustering(tclustResultDataFrame)
                        for (j in 1:amountOfClustersInTclustResult) { # Update height on cluster object
                            tempProteinLabelsOfEachCluster[[j]]$height <- maxThreshold - (currentTclustCluster$nextThreshold - step + currentStep - 1)
                            amountOfTotalClusters <- amountOfTotalClusters + 1
                            cid <- paste0(charC, amountOfTotalClusters)
                            tempProteinLabelsOfEachCluster[[j]]$cid <- cid
                            tempProteinLabelsOfEachCluster[[j]]$parent <- currentTclustCluster$cid
                            
                            nextThreshold <- currentTclustCluster$nextThreshold + currentStep
                            if (nextThreshold > maxThreshold) {
                                nextThreshold <- maxThreshold
                            }
                            tempProteinLabelsOfEachCluster[[j]]$nextThreshold <- nextThreshold
                        }
                    }
                    else { # Same cluster as currentTclustCluster, no split occured
                        tempProteinLabelsOfEachCluster <- getProteinLabelsFromClustering(tclustResultDataFrame)
                        tempProteinLabelsOfEachCluster[[1]]$height <- currentTclustCluster$height
                        tempProteinLabelsOfEachCluster[[1]]$cid <- currentTclustCluster$cid
                        tempProteinLabelsOfEachCluster[[1]]$parent <- currentTclustCluster$parent
                        
                        nextThreshold <- currentTclustCluster$nextThreshold + currentStep
                        if (nextThreshold > maxThreshold) {
                            nextThreshold <- maxThreshold
                        }
                        tempProteinLabelsOfEachCluster[[1]]$nextThreshold <- nextThreshold
                    }
                    
                    # Update order
                    currentStartIndex <- currentTclustCluster$startIndex # Set start of first new cluster = start of parent cluster
                    for (k in 1:amountOfClustersInTclustResult) { # For every new cluster in parent cluster
                        labelsK <- tempProteinLabelsOfEachCluster[[k]]$proteins
                        lengthLabelsK <- length(labelsK)
                        
                        startLabelsK <- currentStartIndex
                        endLabelsK <- currentStartIndex + lengthLabelsK - 1
                        tempProteinLabelsOfEachCluster[[k]]$startIndex <- startLabelsK
                        tempProteinLabelsOfEachCluster[[k]]$endIndex <- endLabelsK
                        
                        for (n in 1:lengthLabelsK) { # Update positions in order
                            order[startLabelsK + n - 1] <- tempProteinLabelsOfEachCluster[[k]]$proteins[n]
                        }
                        # Done with labelsK, update currentStartIndex so start at labelsK.endIndex + 1
                        currentStartIndex <- endLabelsK + 1
                    }
                    
                    # Add clusters for next iteration
                    if (amountOfClustersInTclustResult == 1) {
                        if (DEBUG) {
                            print(paste0("There was no split, threshold = ", currentTclustCluster$nextThreshold))
                            cat("\n")                            
                        }
                        tempCurrentCluster <- tempProteinLabelsOfEachCluster[[1]]
                        clustersLargerThanOne <- c(clustersLargerThanOne, list(tempCurrentCluster)) # We lose name at this line
                    }
                    else {
                        if (DEBUG)  {
                            print(paste0("Was split into ", length(table(tclustResultDataFrame$cluster)), " clusters, threshold = ", currentTclustCluster$nextThreshold - step + currentStep))                            
                        }
                        
                        height <- c(height, c(maxThreshold - (currentTclustCluster$nextThreshold - step + currentStep - 1)))
                        countSplits <- countSplits + 1
                        
                        currentMerge <- c()
                        for (j in 1:amountOfClustersInTclustResult) { # For each new cluster
                            tempCurrentCluster <- tempProteinLabelsOfEachCluster[[j]]
                            
                            # Update totalClusters
                            tempCurrentClusterList <- list(tempCurrentCluster)
                            names(tempCurrentClusterList) <- tempCurrentCluster$cid
                            totalClusters <- c(totalClusters, tempCurrentClusterList)
                            
                            currentMerge <- c(currentMerge, c(tempCurrentCluster$cid))
                            
                            if (length(tempCurrentCluster$proteins) == 1) { # singleton found
                                countSingletons <- countSingletons + 1
                            }
                            if (length(tempCurrentCluster$proteins) > 1) { # Add cluster for next iteration
                                clustersLargerThanOne <- c(clustersLargerThanOne, list(tempCurrentCluster)) # We lose name at this line
                            }
                        }
                        merge <- c(merge, list(currentMerge))
                    }
                }
            }
        }
        
        # All clusters to tclust on next iteration
        proteinLabelsOfEachCluster <- clustersLargerThanOne
    }
    if (DEBUG) {
        print(paste0("Total splits: ", countSplits, " | length(height): ", length(height)))
        as.numeric(levels(as.factor(height)))
        print(paste0("Total clusters = ", amountOfTotalClusters))        
    }
    
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
    mergeMatrix <- matrix(0, nrow = 0, ncol = 2) # Negative values are singletons. Positive values are clusters, where the value corresponds to the row at which this cluster came from.
    mergeLookUpList <- list() # Holds the row numbers for 'merge' for a given merge. mergeLookUpList$c1 returns the row in 'merge' where cluster c1 was made.
    mergeHeights <- c() # Holds the height for all merges
    
    # 'merge' in ascending order corresponds to the ordering of the splits (Divisive)
    # 'merge' in descending order corresponds to the ordering of the merges (Agglomerative)
    # The hierarchical clustering is done divisive, but in order to use the plotting function for a 'hclust() object' we need to see it as agglomerative.
    for (i in length(merge):1) { # For all parents
        parent <- ""
        for (j in 2:length(merge[[i]])) { # For all children of parent i
            # length(merge[[i]]) == 2, A two-split occured -> merging 2 clusters into one
            # length(merge[[i]]) > 2, split resulted in more than 2 clusters -> merging multiple cluster into one
            
            mergeHeights <- c(mergeHeights, height[i])
            
            if (j == 2) {
                c1 <- merge[[i]][j-1]
                c2 <- merge[[i]][j]
                c1.isSingleton <- length(totalClusters[[c1]]$proteins) == 1
                c2.isSingleton <- length(totalClusters[[c2]]$proteins) == 1
                c1.isCluster <- !c1.isSingleton
                c2.isCluster <- !c2.isSingleton
                parent <- totalClusters[[c1]]$parent
                
                if (c1.isSingleton && c2.isSingleton) { # (s,s), on first iteration there are no clusters, only singletons
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
                c1.isCluster <- TRUE # cluster from previous line in mergeMatrix. This will always be a cluster, since we handle merging of 2 singletons above(j == 2, initalize run)
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
        names(newCluster) <- parent # cid for parent. The children have been merged into its parent.
        mergeLookUpList <- c(mergeLookUpList, newCluster)
    }
    
    if (DEBUG) { close(file) }
    
    # Can't plot without making a hclust() object and reassign
    hc <- hclust(dist(USArrests), "ave")
    hc$merge <- mergeMatrix
    hc$height <- mergeHeights
    hc$order <- order
    hc$labels <- labels
    
    return(hc)
}