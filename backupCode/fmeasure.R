########## TODO ##########
# Implement method to count how many proteins has wrong "class"
## 
# Clean Code
##
##########################

# clustering : tclust result
# gsClustring : gold standard in same format as clustering
getQualityOfClustering <- function(clustering, gsClustering) {
    fMeasureClustering <- 0.0
    totalProteins <- length(gsClustering)
    gsClusteringTable <- table(gsClustering)
    amountOfGsClusters <- length(gsClusteringTable)
    amountOfClusters <- length(table(clustering))
    gsMatch <- matrix(nrow = 0, ncol = 4, dimnames = list(c(), c("cluster", "size", "fmeasure", "common")))
    
    clusteringTable <- table(clustering)
    for (j in 1:amountOfGsClusters) {
        measure <- list(cluster = 0, fmeasure = 0.0, common = 0)
        
        gsProteinsOfClassJ <- getProteinsOfClass(gsClustering, j)
        
        measure <- fMeasureClass(gsProteinsOfClassJ, clustering)

        gsMatch <- rbind(gsMatch, c(measure$cluster, clusteringTable[measure$cluster], measure$fmeasure, measure$common))
    }
    
    # Find best candidate, delete all other candidates
    for (j in 1:amountOfClusters) {
        candidates <- which(gsMatch[,1] == j, arr.ind = TRUE) # Find candidates to gsCluster j
        max <- list("index", "fMeasureTotal")
        max$index <- 0
        max$fMeasureTotal <- 0.0
        
        if (length(candidates) > 1) {
            for (i in 1:length(candidates)) { # Find the best match out of candidates, delete others
                row <- candidates[i]
                fMeasureTotal <- gsMatch[row,2] * gsMatch[row,3]

                if (fMeasureTotal > max$fMeasureTotal) { # Update new max
                    gsMatch[max$index, 1] <- 0
                    max$index <- row
                    max$fMeasureTotal <- fMeasureTotal
                } else { # Delete the lower match
                    gsMatch[row, 1] <- 0
                }
            }
        }
    }
    
    sumCommon <- 0
    for (j in 1:amountOfGsClusters) {
        if (gsMatch[j,1] > 0) {
            fMeasureClustering <- fMeasureClustering + (gsMatch[j,2] * gsMatch[j,3])
            sumCommon <- sumCommon + gsMatch[j,4]
        }
    }
    # print(paste("Sum of common:", sum))
    
    # Divide by totalProteins to get the mean F-measure
    fMeasureClustering <- fMeasureClustering / totalProteins
    # return(fMeasureClustering)
    returnObject <- list(fmeasure = fMeasureClustering, common = sumCommon)
    return(returnObject)
}

fMeasureClass <- function(gsProteinsOfClassJ, clustering) {
    fmeasureCluster <- list(cluster = 0, fmeasure = 0.0, common = 0)
    bestMatch <- list(cluster = 0, common = 0)
    bestMatch <- commonProteins(gsProteinsOfClassJ, clustering)
    
    common <- bestMatch$common
    clusterSize <- table(clustering)[bestMatch$cluster]
    
    truePositive <- common
    falsePositive <- clusterSize - common
    falseNegative <- length(gsProteinsOfClassJ) - common
    
    precision <- truePositive / (truePositive + falsePositive)
    recall <- truePositive / (truePositive + falseNegative)
    fmeasure <- 2 * (recall * precision) / (precision + recall)
    
    fmeasureCluster$fmeasure <- fmeasure
    fmeasureCluster$cluster <- bestMatch$cluster
    fmeasureCluster$common <- common
    
    return(fmeasureCluster)
}

# proteins: Vector with all proteins of class j
# Clustering: tclust result
# return(bestMatch): Vector of (cluster, common). 
#  cluster is the cluster with the best match to proteins_vector of class j
#  common is # of matching proteins between proteins_vector and cluster
commonProteins <- function(proteins, clustering) {
    bestMatch <- list(cluster = 0, common = 0)
    commonClusters <- data.frame(cluster = numeric(0), common = numeric(0))
    class <- 0
    
    for (i in 1:length(proteins)) {
        class <- clustering[proteins[i]]
        if (length(which(commonClusters$cluster == class)) > 0) { # class is allready in 'commonClusters', increment 'common' value
            row <- which(commonClusters$cluster == class, arr.ind = TRUE)
            commonClusters[row,2] = commonClusters[row,2] + 1
            
        } else { # class have not been found in 'commonClusters'. Append with 'class' with 'common'=1 to commonClusters
            commonClusters <- rbind(commonClusters, data.frame(cluster = class, common = 1))
        }
    }
    
    numberOfCandidates <- length(which(commonClusters$common == max(commonClusters$common), arr.ind = TRUE))
    if (numberOfCandidates == 1) { # One cadidate for best match
        row <- which.max(commonClusters$common)
        bestMatch$cluster <- commonClusters[row, 1]
        bestMatch$common <- commonClusters[row, 2]
        
    } 
    else if (numberOfCandidates > 1) { # Multiple candidates for best match -> find the one with smallest size
        maxCommon <- max(commonClusters$common)
        rowsOfCandidates <- which(commonClusters$common == maxCommon)
        
        smallestClass <- 0
        for (i in 1:length(rowsOfCandidates)) { # check size of candidate clusters
            rowOfCandidateInCommonClusters <- rowsOfCandidates[i]
            gsClass <- commonClusters[rowOfCandidateInCommonClusters,1]
            
            # With multiple candidates we choose the gsCluster with the smallest total size as best match
            sizeOfGsCluster <- table(clustering)[gsClass]
            if (smallestClass == 0 || sizeOfGsCluster < smallestClass) {
                smallestClass <- gsClass
            }
        }
        if (smallestClass == 0) {
            stop("smallestClass is 0. A smallestClass was never found.")
        }
        bestMatch$cluster <- smallestClass
        bestMatch$common <- maxCommon
        
    } else {
        stop("No candidates found for common clusters")
    }
    
    if (length(bestMatch$cluster) > 1) {
        stop("bestMatch contains more than one cluster match. It should only hold one.")
    }
    
    return(bestMatch)
}

# clustering 1-dim: Index is specific protein, clustering[i] is class of protein i.
# j: Protein class
# Return(proteins_vector): Vector of all proteins of class j found in clustering
getProteinsOfClass <- function(clustering, j) {
    proteins <- c()
    for (i in 1:length(clustering)) { # i is a protein. Looping over all proteins #232
        if (clustering[i] == j) { # Found protein of class j
            proteins <- c(proteins, i)
        }
        if (length(proteins) >= table(clustering)[j]) { # We have found all proteins of cluster j
            break;
        }
    }
    
    return(proteins)
}
