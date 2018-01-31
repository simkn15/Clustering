# Not the proper way
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
source("randomization.R")

# NOTE!
# The Hierarchical Clustering approach is divisive.
# The information for the dendrogram is converted into agglomerative, thus being able to use the plotting defined for a hclust() object
# Required to make a hclust() object, doing a clustering with the method, and reassigning the variables on the object.

############# One vs. One plotting of costs with biggest gap and foldChange
plotCostWithGap <- function(df_original, df_random, gap, dimensions = 5, seed = 42, ranApproach = 4) {
    maxGapThreshold <- gap[1]
    maxGapFoldChange <- gap[2]
    dimensionsAsChar <- as.character(dimensions)
    df_plot <- data.frame(threshold = integer(0), cost = integer(0), dimension = character())
    names <- c("threshold", "cost")
    
    df_original$dimension <- rep("2(Original data)", nrow(df_original))
    df_plot <- rbind(df_plot, df_original)
    
    df_random$dimension <- rep(dimensionsAsChar, nrow(df_random))
    df_plot <- rbind(df_plot, df_random)
    
    title <- paste0("Big Data Set : Gap threshold = ", maxGapThreshold, " | foldChange = ", maxGapFoldChange)
    g <- ggplot(df_plot, aes(x = threshold, y = cost, color = dimension)) + 
        ggtitle(title) + geom_line() + 
        geom_vline(xintercept = maxGapThreshold, linetype="dotted") +
        scale_y_continuous(name = "Cost", labels=function(n){format(n, scientific = FALSE)}) +
        labs(color = "Dimensions")
    return(g)
}

getMaxGap <- function(df_original, df_random) {
    rows <- min(nrow(df_original), nrow(df_random))
    maxGapThreshold <- -1
    maxGapCost <- 0
    fold <- 0
    for (row in 1:rows) {
        costOriginal <- df_original[row,2]
        costRandom <- df_random[row,2]
        gap <- abs(costRandom - costOriginal)
        if (gap > maxGapCost) {
            fold <- foldchange(costRandom, costOriginal)
            maxGapThreshold <- df_original[row,1]
            maxGapCost <- gap
        }
    }
    return(c(maxGapThreshold, round(fold, digits = 4)))
}

plotHistogramSimilarityDistribution <- function(simMatrixOriginal, simMatrixRandom) {
    colNames <- c("Similarity", "Dataset")
    df_plot <- data.frame(similarity = integer(0), dataset = character())
    
    rows <- length(proteins)^2
    df_plot <- as.data.frame(cbind(as.vector(simMatrixOriginal), rep("Original", rows)))
    names(df_plot) <-  colNames
    
    df_random <- as.data.frame(cbind(as.vector(simMatrixRandom), rep("Random", rows)))
    names(df_random) <- colNames
    df_plot <- rbind(df_plot, df_random)
    df_plot[,1] <- as.numeric(as.character(df_plot[,1]))
    
    g <- ggplot(data = df_plot, aes(x = Similarity, fill = Dataset)) + 
        geom_histogram(binwidth = 10,
                       position = "identity", colour = "black", boundary = 0) +
        scale_x_continuous(name = "Similarity", breaks = seq(0, max(df_plot$Similarity)+20, 20)) +
        scale_y_log10(name = "Frequency", labels=function(n){format(n, scientific = FALSE)}) +
        ggtitle("Frequency of similarities") +
        theme_bw() +
        theme(axis.line = element_line(size=1, colour = "black"),
              panel.grid.major = element_line(colour = "#d3d3d3"),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), panel.background = element_blank(),
              axis.text.x=element_text(colour="black", size = 9, angle = 90),
              axis.text.y=element_text(colour="black", size = 9),
              legend.position = "none") +
        facet_grid(. ~Dataset , scales = "free")
    
    return(g)
}

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
hclustDivisiveRandom <- function(simMatrix, proteins, step = 1, minSplit = 2, maxSplit = 10, minThreshold, maxThreshold, binarySearch = FALSE, GAP = FALSE, dimensions = 5) {
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
        clustersLargerThanOne <- list() # All clusters returned from this iteration of tclust(Not singletons). These will be added to 'proteinLabelsOfEachCluster' to be clustered in the next iteration.
        
        # For every level: Do tclust on all clusters with current threshold
        if (length(proteinLabelsOfEachCluster) > 0) { # TODO: Remove this check ?
            for (i in 1:length(proteinLabelsOfEachCluster)) {
                currentStep <- step
                currentTclustCluster <- proteinLabelsOfEachCluster[[i]]# Cluster to tclust
                
                if (length(currentTclustCluster$proteins) > 1) { # Not a singleton. This should not be necessary, as we dont add singletons for the next iteration
                    # simMatrix <- buildSimilarityMatrix(currentTclustCluster$proteins, df_table)
                    # 
                    # tclustResultDataFrame <- clusteringWithTclust(simMatrix, currentTclustCluster$proteins, currentTclustCluster$nextThreshold)
                    # amountOfClustersInTclustResult <- length(table(tclustResultDataFrame$cluster))
                    
                    # TODO: Implement randomization
                    # Is this the right spot ?
                    # Or do we need to rethink the "whole" algorithm?
                    if (!GAP) {
                        simMatrixTemp <- simMatrix[currentTclustCluster$proteins, currentTclustCluster$proteins]

                        tclustResultDataFrame <- clusteringWithTclust(simMatrixTemp, currentTclustCluster$proteins, currentTclustCluster$nextThreshold)
                        amountOfClustersInTclustResult <- length(table(tclustResultDataFrame$cluster))
                    }
                    # bad: dim > proteins - 1
                    # good: dim <= proteins - 1
                    # TODO Split this into two if statements, such that it jumps to the right if, when length-1 is not >= dimensions
                    if (GAP && length(currentTclustCluster$proteins) - 1 >= dimensions) { # Apparently we need to have amount of proteins > dim+1, else errors and warnings from isoMDS()/cmdscale()
                        # TODO: Add seed and dim as parameters
                        simMatrixTemp <- simMatrix[currentTclustCluster$proteins, currentTclustCluster$proteins]
                        
                        t_min <- currentTclustCluster$nextThreshold
                        # t_max <- round(sqrt(mean(0simMatrixTemp[simMatrixTemp > 0])))
                        # t_max <- t_min + round(mean(simMatrixTemp))
                        t_max <- t_min + 20
                        if (t_max > maxThreshold) {
                            t_max <- maxThreshold
                        }
                        
                        df_original <- data.frame(threshold = integer(0), cost = integer(0))
                        df_random <- data.frame(threshold = integer(0), cost = integer(0))
                        
                        print(paste0("Doing tclust on range: ", t_min, " to ", t_max, " with ", length(currentTclustCluster$proteins), " proteins"))
                        for (t in seq(t_min, t_max, 1)) {
                            tclustResult <- tclust(simmatrix = simMatrixTemp, convert_dissimilarity_to_similarity = FALSE, threshold = t)
                            df_original <- rbind(df_original, data.frame(threshold = t, cost = tclustResult$costs))
                        }
                        
                        simMatrixRandom <- buildRandomSimMatrixAp4(currentTclustCluster$proteins, simMatrixTemp, k = dimensions, seed = 42)
                        for (t in seq(t_min, t_max, 1)) {
                            tclustResult <- tclust(simmatrix = simMatrixRandom, convert_dissimilarity_to_similarity = FALSE, threshold = t)
                            df_random <- rbind(df_random, data.frame(threshold = t, cost = tclustResult$costs))
                        }
                        
                        gap <- getMaxGap(df_original, df_random)
                        if (gap[1] != -1) {
                            costPlot <- plotCostWithGap(df_original, df_random, gap, dimensions = dimensions, seed = 42, ranApproach = 4)
                            simPlot <- plotHistogramSimilarityDistribution(simMatrixTemp, simMatrixRandom)
                            
                            fileName <- paste0("Vap3", currentTclustCluster$cid, "Threshold", gap[1], ".pdf")
                            pdf(fileName)
                            print(costPlot)
                            print(simPlot)
                            dev.off() # close output file for plots
                            
                            print(paste0("Biggest gap was at threshold:" , gap[1], " | foldChange: ", gap[2]))
                            currentTclustCluster$nextThreshold <- gap[1]
                            tclustResultDataFrame <- clusteringWithTclust(simMatrixTemp, currentTclustCluster$proteins, gap[1])
                            amountOfClustersInTclustResult <- length(table(tclustResultDataFrame$cluster))                            
                        }
                        else {
                            print(paste0("No gap was found"))
                            currentTclustCluster$nextThreshold <- t_max
                            tclustResultDataFrame <- clusteringWithTclust(simMatrixTemp, currentTclustCluster$proteins, currentTclustCluster$nextThreshold)
                            amountOfClustersInTclustResult <- length(table(tclustResultDataFrame$cluster))                            
                        }
                    }
                    else {
                        print(paste0("Amount of proteins are lower than the bound, no randomization will occur for cluster: ", currentTclustCluster$cid))
                        simMatrixTemp <- simMatrix[currentTclustCluster$proteins, currentTclustCluster$proteins]
                        
                        tclustResultDataFrame <- clusteringWithTclust(simMatrixTemp, currentTclustCluster$proteins, currentTclustCluster$nextThreshold)
                        amountOfClustersInTclustResult <- length(table(tclustResultDataFrame$cluster))
                    }
                    
                    # Get labels for each new cluster
                    # TODO check if clusters == 1, no need to getLabels, since currentTclustCluster holds all labels of cluster. Do the check inside the function
                    fallBack <- FALSE # If binary search did not minimize the number of splits, reset to the initial threshold/result
                    if (amountOfClustersInTclustResult > 1) {
                        
                        ### BEGIN binary search
                        if (binarySearch && amountOfClustersInTclustResult >= maxSplit) {
                            tempTclustResultDataFrame <- tclustResultDataFrame
                            tempAmountOfClustersInTclustResult <- amountOfClustersInTclustResult
                            upperBound <- step
                            lowerBound <- 0
                            
                            while(amountOfClustersInTclustResult >= maxSplit) {
                                currentStep <- lowerBound + ceiling(abs(upperBound - lowerBound) / 2)
                                tclustResultDataFrame <- clusteringWithTclust(simMatrixTemp, currentTclustCluster$proteins, currentTclustCluster$nextThreshold - step + currentStep)
                                amountOfClustersInTclustResult <- length(table(tclustResultDataFrame$cluster))
                                
                                while (amountOfClustersInTclustResult < minSplit && currentStep > 1 && currentStep != step) {
                                    lowerBound <- currentStep
                                    currentStep <- lowerBound + ceiling(abs(upperBound - lowerBound) / 2)

                                    if (currentStep == upperBound) { # Did not find a better threshold
                                        fallBack = TRUE
                                        break
                                    }
                                    
                                    tclustResultDataFrame <- clusteringWithTclust(simMatrixTemp, currentTclustCluster$proteins, currentTclustCluster$nextThreshold - step + currentStep)
                                    amountOfClustersInTclustResult <- length(table(tclustResultDataFrame$cluster))
                                }
                                
                                if (amountOfClustersInTclustResult >= maxSplit && (upperBound - lowerBound == 1)) {
                                    fallBack = TRUE
                                }
                                if (fallBack) {
                                    break                            
                                }
                                if (amountOfClustersInTclustResult >= maxSplit) { # currentStep did not descrease the amount of clusters enough
                                    # Decrease upperBound
                                    upperBound <- currentStep
                                }
                            }
                            if (fallBack) { # Resetting back to first result
                                if ((amountOfClustersInTclustResult > tempAmountOfClustersInTclustResult) || 
                                    (amountOfClustersInTclustResult == 1 && (amountOfClustersInTclustResult < tempAmountOfClustersInTclustResult))) {
                                    # Either the first tclustResult was better, or the first tclustResult gave the only split
                                    tclustResultDataFrame <- tempTclustResultDataFrame
                                    amountOfClustersInTclustResult <- tempAmountOfClustersInTclustResult
                                }
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
                        tempCurrentCluster <- tempProteinLabelsOfEachCluster[[1]]
                        clustersLargerThanOne <- c(clustersLargerThanOne, list(tempCurrentCluster)) # We lose name at this line
                    }
                    else {
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
    
    singletons <- c()
    clusters <- c()
    for (i in 1:length(totalClusters)) {
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

####################################################################################################
# Read in data and build similarity matrix
# hclust
####################################################################################################
# Adjust variables to desired bounds for binary search
# The binary search uses a upperBound and a lowerBound to find the proper threshold in the current iteration.
step <- 1 # Default step for threshold. Binary search is skipped if 'step <- 1'. 30 is used for the binary split
minSplit <- 5 # Minimum number of splits. Minimum bound will increase if splits are < 2
maxSplit <- 10 # Maximum number of splits. Upper bound will decrease if splits are >= maxSplit.
readSmallDataSet <- FALSE
readBigDataSet <- !readSmallDataSet
DEBUG <- FALSE

if (readSmallDataSet) {
    table <- read.table("./data/brown/sfld_brown_et_al_amidohydrolases_protein_similarities_for_beh.txt", sep = "", skip = 5)
    df_table <- as.data.frame(table)
    proteins = levels(df_table[,1])
    simMatrix <- buildSimilarityMatrix(proteins, df_table)

    minThreshold <- min(simMatrix) - 1 # The starting threshold. -1 such that the first run returns 1 cluster.
    maxThreshold <- round(max(simMatrix)) + 1 # Threshold where all clusters will be singletons

    if (DEBUG) { fileName <- "hcSmallData" }
}
if (readBigDataSet) {
    table <- read.table("./data/big/simBig.txt", sep = "")
    df_table <- as.data.frame(table)
    proteins = levels(df_table[,1])
    simMatrix <- buildSimilarityMatrixFromBlast(proteins, df_table)

    minThreshold <- min(simMatrix) - 1
    maxThreshold <- round(max(simMatrix)) + 1

    if (DEBUG) { fileName <- "hcBigData" }
}

hc <- hclustDivisiveRandom(simMatrix, proteins, step, minSplit, maxSplit, minThreshold, maxThreshold, binarySearch = FALSE, GAP = TRUE, dimensions = 5)
plot(hc, xlab = "protein", ylab = "height")
plot(hc, hang = -1, xlab = "Protein", ylab = "Height")


