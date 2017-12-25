library(transclustr)
library(ggplot2)
library(plyr)
library(cluster)
library(gtools)
library(ClusterR)
library(stringi)
require(graphics)

#library(pryr) # get size of objects with object_size(l)

source("fmeasure.R")
source("utilities.R")
####################################################################################################
# Read in data and build similarity matrix
####################################################################################################
step <- 30
minSplit <- 2
maxSplit <- 10

readSmallDataSet <- FALSE
readBigDataSet <- !readSmallDataSet

writeToFile <- TRUE

if (readSmallDataSet) {
    table <- read.table("../data/brown/sfld_brown_et_al_amidohydrolases_protein_similarities_for_beh.txt", sep = "", skip = 5)
    df_table <- as.data.frame(table)
    proteins = levels(df_table[,1])
    simMatrix <- buildSimilarityMatrix(proteins, df_table)
    
    minThreshold <- min(simMatrix) - 1
    maxThreshold <- round(max(simMatrix)) + 1
    
    if (writeToFile) {
        fileName <- paste0("SavingWorkSpaceHcOutputSmallDataSetMin", minSplit, "Max", maxSplit, "Step", step, ".txt")
        file <- file(fileName, "w") 
    }
}
if (readBigDataSet) {
    table <- read.table("simBig.txt", sep = "") 
    df_table <- as.data.frame(table)
    proteins = levels(df_table[,1])
    simMatrix <- buildSimilarityMatrixFromBlast(proteins, df_table)
    
    minThreshold <- min(simMatrix) - 1
    maxThreshold <- round(max(simMatrix)) + 1

    if (writeToFile) {
        fileName <- paste0("DEBUGhcTestOutputBigDataSetMin", minSplit, "Max", maxSplit, "Step", step, ".txt")
        file <- file(fileName, "w") 
    }
}

####################################################################################################
# Hierarchical clustering
####################################################################################################
# hcObject <- list(merge = matrix(0, nrow = 0, ncol = 2), height = c(), order = c(), labels = proteins)

proteinLabelsOfEachCluster <- list(c1 = list(cid = "c1", proteins = c(proteins), startIndex = 1, endIndex = length(proteins), height = -1, parent = "", nextThreshold = minThreshold))

amountOfTotalClusters <- 1
totalClusters <- proteinLabelsOfEachCluster

countSplits <- 0
countSingletons <- 0

charC <- "c" # Is this used anywhere ?

# hc object
mergeMatrix <- matrix(0, nrow = 0, ncol = 2)
merge <- list()
height <- c() # maxThreshold - threshold
order <- proteins
labels <- proteins
while (countSingletons < length(proteins)) {
    cat("\n")
    print("##############################")
    print(paste("Clusters to tclust =", length(proteinLabelsOfEachCluster)))
    
    clustersLargerThanOne <- list()
    
    # For every level: Do tclust on all clusters with current threshold
    if (length(proteinLabelsOfEachCluster) > 0) {
        for (i in 1:length(proteinLabelsOfEachCluster)) {
            currentStep <- step
            currentTclustCluster <- proteinLabelsOfEachCluster[[i]]# Cluster to tclust
            print(paste("currentTclustCluster$nextThreshold =", currentTclustCluster$nextThreshold))
            
            if (writeToFile) {
                write("", file, append = TRUE)
                string <- paste0("currentTclustCluster = ", currentTclustCluster$cid, " | nextThreshold = ", currentTclustCluster$nextThreshold)
                write(string, file, append = TRUE) 
            }
            
            if (length(currentTclustCluster$proteins) > 1) { # Not a singleton
                
                print(paste0("currentTclustCluster: cid = ", currentTclustCluster$cid, " | size = ", length(currentTclustCluster$proteins), " | startIndex = ", currentTclustCluster$startIndex, " | endIndex = ", currentTclustCluster$endIndex, " | height = ", currentTclustCluster$height, " | parent = ", currentTclustCluster$parent))
                # print(currentTclustCluster$proteins)
                
                simMatrix <- buildSimilarityMatrix(currentTclustCluster$proteins, df_table)

                tclustResultDataFrame <- clusteringWithTclust(simMatrix, currentTclustCluster$proteins, currentTclustCluster$nextThreshold)
                amountOfClustersInTclustResult <- length(table(tclustResultDataFrame$cluster))
                
                # Get labels for each new cluster
                # TODO check if clusters == 1, no need to getLabels, since currentTclustCluster holds all labels of cluster. Do the check inside the function
                fallBack <- FALSE
                if (amountOfClustersInTclustResult > 1) {
                    
                    ### BEGIN New stuff
                    if (amountOfClustersInTclustResult >= maxSplit && step > 1) { # Too many new clusters. Can't do binary search if base step is 1
                        tempTclustResultDataFrame <- tclustResultDataFrame
                        tempAmountOfClustersInTclustResult <- amountOfClustersInTclustResult
                        upperBound <- step
                        lowerBound <- 0
                        string <- paste0("Beginning binary search with cluster: ", currentTclustCluster$cid, " | amountOfClustersInTclustResult = ", amountOfClustersInTclustResult)
                        print(string)
                        if (writeToFile) { write(string, file, append = TRUE) }
                        
                        while(amountOfClustersInTclustResult >= maxSplit) {
                            currentStep <- lowerBound + ceiling(abs(upperBound - lowerBound) / 2)
                            # print(paste0(currentTclustCluster$cid, " | Decreasing step -> ", currentStep))
                            
                            string <- paste0("Decreasing to step = ", currentStep, " | lowerBound = ", lowerBound, " | upperBound = ", upperBound)
                            print(string)
                            if (writeToFile) { write(string, file, append = TRUE) }
                            
                            tclustResultDataFrame <- clusteringWithTclust(simMatrix, currentTclustCluster$proteins, currentTclustCluster$nextThreshold - step + currentStep)
                            amountOfClustersInTclustResult <- length(table(tclustResultDataFrame$cluster))
                            
                            string <- paste0("amountOfClustersInTclustResult = ", amountOfClustersInTclustResult, " | currentStep = ", currentStep)
                            print(string)
                            if (writeToFile) { write(string, file, append = TRUE) }
                            
                            while (amountOfClustersInTclustResult < minSplit && currentStep > 1 && currentStep != step) {
                                string <- paste0("Inner loop | Updating lowerBound from ", lowerBound, " to ", currentStep)
                                print(string)
                                if (writeToFile) { write(string, file, append = TRUE) }
                                
                                lowerBound <- currentStep
                                currentStep <- lowerBound + ceiling(abs(upperBound - lowerBound) / 2)
                                
                                string <- paste0("Increasing to step = ", currentStep, " | lowerBound = ", lowerBound, " | upperBound = ", upperBound)
                                print(string)
                                if (writeToFile) { write(string, file, append = TRUE) }
                                
                                if (currentStep == upperBound) { # Can't find better threshold
                                    string <- paste0("fallBack mode at inner, currentStep ", currentStep, " == ", upperBound, " upperBound")
                                    print(string)
                                    if (writeToFile) { write(string, file, append = TRUE) }
                                    
                                    fallBack = TRUE
                                    break
                                }
                                
                                tclustResultDataFrame <- clusteringWithTclust(simMatrix, currentTclustCluster$proteins, currentTclustCluster$nextThreshold - step + currentStep)
                                amountOfClustersInTclustResult <- length(table(tclustResultDataFrame$cluster))
                                
                                string <- paste0("amountOfClustersInTclustResult = ", amountOfClustersInTclustResult, " | currentStep = ", currentStep)
                                print(string)
                                if (writeToFile) { write(string, file, append = TRUE) }
                            }
                            
                            if (amountOfClustersInTclustResult >= maxSplit && (upperBound - lowerBound == 1)) {
                                string <- paste0("fallBack mode at outer, currentStep ", currentStep, " == 1", " | amountOfClustersInTclustResult >= ", maxSplit)
                                print(string)
                                if (writeToFile) { write(string, file, append = TRUE) }
                                
                                fallBack = TRUE
                            }
                            if (fallBack) {
                                break                            
                            }
                            if (amountOfClustersInTclustResult >= maxSplit) { # currentStep did not descrease the amount of clusters enough
                                string <- paste0("Outer loop | amountOfClustersInTclustResult >= ", maxSplit, " | Updating upperBound from ", upperBound, " to ", currentStep)
                                print(string)
                                if (writeToFile) { write(string, file, append = TRUE) }
                                
                                upperBound <- currentStep
                            }
                        }
                    }
                    
                    if (fallBack) {
                        string <- paste0("fallBack mode caught")
                        print(string)
                        if (writeToFile) { write(string, file, append = TRUE) }
                        
                        print(paste0("amountOfClustersInTclustResult ", amountOfClustersInTclustResult, " > ", tempAmountOfClustersInTclustResult, " tempAmountOfClustersInTclustResult"))
                        
                        if ((amountOfClustersInTclustResult > tempAmountOfClustersInTclustResult) || 
                            (amountOfClustersInTclustResult == 1 && (amountOfClustersInTclustResult < tempAmountOfClustersInTclustResult))) {
                            # Either the first tclustResult was better, or the first tclustResult gave the only split
                            string <- paste0("######## Resetting into pre-binary-search")
                            print(string)
                            if (writeToFile) { write(string, file, append = TRUE) }
                            
                            tclustResultDataFrame <- tempTclustResultDataFrame
                            amountOfClustersInTclustResult <- tempAmountOfClustersInTclustResult
                        }
                        
                        string <- paste0("amountOfClustersInTclustResult = ", amountOfClustersInTclustResult)
                        print(string)
                        if (writeToFile) { write(string, file, append = TRUE) }
                    }
                    ### END New stuff

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
                else { # Same cluster as currentTclustCluster
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
                    print(paste0("There was no split, threshold = ", currentTclustCluster$nextThreshold))
                    cat("\n")
                    tempCurrentCluster <- tempProteinLabelsOfEachCluster[[1]]
                    clustersLargerThanOne <- c(clustersLargerThanOne, list(tempCurrentCluster)) # We lose name at this line
                }
                else {
                    print(paste0("Was split into ", length(table(tclustResultDataFrame$cluster)), " clusters, threshold = ", currentTclustCluster$nextThreshold - step + currentStep))
                    
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

# which(proteins == totalClusters$c324$proteins[[1]])

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
hc <- hclust(dist(USArrests), "ave")
hc$merge <- mergeMatrix
hc$height <- mergeHeights
hc$order <- order
hc$labels <- labels
plot(hc, xlab = "protein", ylab = "height")
plot(hc, hang = -1, xlab = "Protein", ylab = "Height")

# Fails
# hct <- list(merge = mergeMatrix, height = height, order = order, labels = labels)
# plot(hct, xlab = "Protein", ylab = "Threshold")
# plot(hct, hang = -1)

# hclust original dataset
# hc <- hclust((max(as.dist(simMatrix)) + 1) - as.dist(simMatrix))
# plot(hc, xlab = "protein", ylab = "threshold")
# plot(hc, hang = -1, xlab = "Protein", ylab = "Threshold")

# Approx optimal threshold
# > mean(as.vector(as.dist(simMatrix)))
#[1] 55.91576

# hc <- hclust(dist(USArrests), "ave")
# plot(hc)
# plot(hc, hang = -1)
