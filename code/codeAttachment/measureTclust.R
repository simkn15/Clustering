setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(transclustr)
library(ggplot2)
library(plyr)
library(cluster)
library(gtools)
library(ClusterR)
library(stringi)

source("fmeasure.R")
source("utilities.R")
source("randomization.R")

####################################################################################################
# Runs tclust() on the simMatrix within range minTheshold to maxThreshold
#
# simMatrix: Similarity matrix
# proteins: Names of all proteins
# randomSim: TRUE if the similarity matrix has been randomized
# minThreshold: Starting threshold
# maxThreshold: Last threshold to run
# step: Incrementation value for threshold
# writeToFile: Boolean stating if measures should be output to a file
# fileName: Name of the output file
#
# Returns bestClustering which holds information about the best found clustering within the range according to F-measure.
# Randomized simimilarity matrix: Returns bestClustering with 0's. F-measure is not calculated.
####################################################################################################
doTclustWithinRange <- function(simMatrix, proteins, randomSim = FALSE, minThreshold, maxThreshold, step = 1, writeToFile = FALSE, fileName = "output.txt") {
    if (writeToFile) {
        file <- file(fileName, "w")
    }
    
    bestClustering <- list(threshold = 0, fmeasure = 0.0, clusters = 0, cost = 0.0, common = 0)

    if (!randomSim) {
        write(paste0("Threshold F-measure Clusters Common Cost"), file, append = TRUE)
    }
    if (randomSim) { # Creating headers
        write(paste0("Threshold Cost Clusters"), file, append = TRUE)
    }
    
    for (threshold in seq(minThreshold, maxThreshold, step)) {
        tclustResult <- tclust(simmatrix = simMatrix, convert_dissimilarity_to_similarity = FALSE, threshold = threshold)
        tclustResultDataFrame <- data.frame(protein = proteins, cluster = tclustResult$clusters[[1]])
        # ## Add 1 to all classes to match first class = 1 in gold standard
        tclustResultDataFrame$cluster = tclustResultDataFrame$cluster + 1
        
        clusters <- length(table(tclustResultDataFrame$cluster))
        cost <- tclustResult$costs[[1]]
        if (!randomSim) {
            fmeasure <- getQualityOfClustering(tclustResultDataFrame$cluster, gsClusteringDataFrame)
            if (fmeasure$fmeasure > bestClustering$fmeasure) {
                bestClustering$fmeasure <- fmeasure$fmeasure
                bestClustering$threshold <- threshold
                bestClustering$clusters <- clusters
                bestClustering$cost <- cost
                bestClustering$common <- fmeasure$common
            }
        }
        
        if (writeToFile) {
            if (!randomSim) { string <- paste(threshold, round(fmeasure$fmeasure, 7), clusters, fmeasure$common, cost) }
            if (randomSim) { string <- paste(threshold, round(cost), clusters) }
            print(string)
            write(string, file, append = TRUE)
        }
    }

    if (writeToFile) { close(file) }
    
    return(bestClustering)
}