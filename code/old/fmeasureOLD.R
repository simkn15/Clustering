########## TODO ##########
# Implement method to count how many proteins has wrong "class"
##
# rename 'common' to something that makes more sense: inCommon, inCommonProteins, inCommonAmount, proteinsInCommon, matchingProteins
## 
# Clean Code
##
##########################

# clustering : tclust result
# gsClustring : gold standard in same format as clustering
getQualityOfClustering <- function(clustering, gsClustering) {
  fMeasureClustering <- 0.0
  
  totalProteins <- length(clustering)
  tableClustering <- table(clustering)
  amountOfClusters <- length(tableClustering)
  
  for (j in 1:amountOfClusters) {
    proteinsInReference <- tableClustering[j]
    proteins_vector <- getProteinsOfClass(clustering, j)
    maxMeasure <- fMeasureClass(proteins_vector, gsClustering)
    
    # Take all proteins of cluster j into account
    fMeasureClustering <- fMeasureClustering + (maxMeasure * proteinsInReference)
  }
  
  # Divide by totalProteins to get the mean F-measure
  fMeasureClustering <- fMeasureClustering / totalProteins
  
  return(fMeasureClustering)
}

fMeasureClass <- function(proteins_vector, gsClustering) {
  fmeasure <- 0.0
  #bestMatch <- data.frame(cluster = numeric(0), common = numeric(0))
  #bestMatch <- rbind(bestMatch, commonProteins(proteins_vector, gsClustering)[1,])
  
  bestMatch <- commonProteins(proteins_vector, gsClustering)[1,]
  common <- bestMatch$common
  gsClusterSize <- table(gsClustering)[bestMatch$cluster]
  
  truePositive <- common
  falsePositive <- length(proteins_vector) - common
  falseNegative <- gsClusterSize - common
  
  precision <- (truePositive / (truePositive + falsePositive))
  recall <- (truePositive / (truePositive + falseNegative))
  fmeasure <- 2 * (recall * precision) / (precision + recall)
  
  return(fmeasure)
}

# proteins_vector: Vector with all proteins of class j
# gsClustering: gs clustering, sorted according to tclust result
# return(bestMatch): Vector of (cluster, common). 
#  cluster is the gs cluster with the best match to proteins_vector of class j
#  common is # of matching proteins between proteins_vector and gs cluster
commonProteins <- function(proteins_vector, gsClustering) {
  bestMatch <- data.frame(cluster = numeric(0), common = numeric(0))
  commonClusters <- data.frame(cluster = numeric(0), common = numeric(0))
  class <- 0
  
  for (i in 1:length(proteins_vector)) {
    class <- gsClustering[proteins_vector[i]]
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
    bestMatch <- rbind(bestMatch, commonClusters[row,])
    
  } else if (numberOfCandidates > 1) { # Multiple candidates for best match -> find the one with smallest size
    maxCommon <- max(commonClusters$common)
    rowsOfCandidates <- which(commonClusters$common == maxCommon)
    
    smallestClass <- 0
    for (i in 1:length(rowsOfCandidates)) { # check size of candidate clusters
      rowOfCandidateInCommonClusters <- rowsOfCandidates[i]
      gsClass <- commonClusters[rowOfCandidateInCommonClusters,1]
      
      # With multiple candidates we choose the gsCluster with the smallest total size as best match
      sizeOfGsCluster <- table(gsClustering)[gsClass]
      if (smallestClass == 0 || sizeOfGsCluster < smallestClass) {
        smallestClass <- gsClass
      }
    }
    if (smallestClass == 0) {
      stop("smallestClass is 0. A smallestClass was never found.")
    }
    bestMatch <- rbind(bestMatch, data.frame(cluster = smallestClass, common = maxCommon))
    
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
  proteins_vector <- c()
  tableClustering <- table(clustering)
  for (i in 1:length(clustering)) { # i is a protein. Looping over all proteins #232
    if (length(proteins_vector) >= tableClustering[j]) { # We have found all proteins of cluster j
      break;
    }
    if (clustering[i] == j) { # Found protein of class j
      proteins_vector <- append(proteins_vector, i)
    }
  }
  
  return(proteins_vector)
}
