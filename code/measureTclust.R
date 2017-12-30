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
source("randomV2.R")

####################################################################################################
# Cluster
####################################################################################################
doTclustWithinRange <- function(simMatrix, proteins, randomSim = FALSE, minThreshold, maxThreshold, step = 1, writeToFile = FALSE, fileName = "output.txt") {
    if (writeToFile) {
        file <- file(fileName, "w")
    }
    
    start.time <- Sys.time()
    if (randomSim == FALSE) { bestClustering <- list(threshold = 0, fmeasure = 0.0, clusters = 0, cost = 0.0, common = 0) }
    if (randomSim) {write(paste0("Threshold Cost Clusters"), file, append = TRUE)} # Creating headers
    
    for (threshold in seq(minThreshold, maxThreshold, step)) {
        tclustResult <- tclust(simmatrix = simMatrix, convert_dissimilarity_to_similarity = FALSE, threshold = threshold)
        tclustResultDataFrame <- data.frame(protein = proteins, cluster = tclustResult$clusters[[1]])
        # ## Add 1 to all classes to match first class = 1 in gold standard
        tclustResultDataFrame$cluster = tclustResultDataFrame$cluster + 1
        
        clusters <- length(table(tclustResultDataFrame$cluster))
        cost <- tclustResult$costs[[1]]
        if (randomSim == FALSE) {
            fmeasure <- getQualityOfClustering(tclustResultDataFrame$cluster, gsClusteringDataFrame)
            if (fmeasure$fmeasure > bestClustering$fmeasure) {
                bestClustering$fmeasure <- fmeasure$fmeasure
                bestClustering$threshold <- threshold
                bestClustering$clusters <- clusters
                bestClustering$cost <- cost
                bestClustering$common <- fmeasure$common
            }
            string <- paste0("Threshold = ", threshold, " | F-measure = ", round(fmeasure$fmeasure, 7), " | clusters = ", clusters, " | common = ", fmeasure$common, " | cost = ", cost)
        }
        
        if (randomSim == TRUE) { string <- paste(threshold, round(cost), clusters) }
        
        print(string)
        if (writeToFile) { write(string, file, append = TRUE) }
    }
    if (randomSim == FALSE) {
        stringOne <- paste("Best Clustering:")
        print(stringOne)
        stringTwo <- paste0("Threshold = ", bestClustering$threshold, " | F-measure = ", round(bestClustering$fmeasure, 7), " | clusters = ", bestClustering$clusters, " | common = ", bestClustering$common, " | cost = ", bestClustering$cost)
        print(stringTwo)
    }
    
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    time.taken
    timeTaken <- paste0("Total time for clustering: ", time.taken, " hours")
    if (writeToFile) {
        if (randomSim == FALSE) {
            write(stringOne, file, append = TRUE)
            write(stringTwo, file, append = TRUE)   
        }
        # write(timeTaken, file, append = TRUE)
        close(file)
    }    
}

####################################################################################################
# Read in data
# Build similarity matrix
# Run doTclustWithinRange()
# Update parameters according to desired run
####################################################################################################
randomSim <- TRUE
readSmallData <- FALSE
readBigData <- !readSmallData

if (readSmallData && !readBigData) {
    table <- read.table("./data/brown/sfld_brown_et_al_amidohydrolases_protein_similarities_for_beh.txt", sep = "", skip = 5)
    df_table <- as.data.frame(table)
    proteins = levels(df_table[,1])
    simMatrix <- buildSimilarityMatrix(proteins, df_table)
    
    gold_df_before_split <- as.data.frame(readLines("./data/brown/sfld_brown_et_al_amidohydrolases_families_gold_standard.txt"))
    gold <- as.data.frame(stri_split_fixed(gold_df_before_split[,1], "\t", simplify = TRUE))
    colnames(gold) <- c("proteins", "class")
    gold_sorted <- gold[match(proteins,gold$proteins),]
    gsClusteringDataFrame <- as.numeric(as.factor(gold_sorted$class))
}
if (readBigData && !readSmallData) {
    table <- read.table("./data/big/simBig.txt", sep = "")
    df_table <- as.data.frame(table)
    proteins = levels(df_table[,1])
    simMatrix <- buildSimilarityMatrixFromBlast(proteins, df_table)
    
    gold_df_before_split <- as.data.frame(readLines(".data/big/gold_standard_3rd_column.txt"))
    gold <- as.data.frame(stri_split_fixed(gold_df_before_split[,1], "\t", simplify = TRUE))
    gold[,2] <- NULL # Deleting second column, as it is not currently needed
    colnames(gold) <- c("proteins", "class")
    gold_sorted <- gold[match(proteins,gold$proteins),]
    gsClusteringDataFrame <- as.numeric(as.factor(gold_sorted$class))
}

if (!randomSim) {
    step <- 1
    writeToFile <- TRUE
    fileName <- "measureSmallDataOriginal.txt"
    
    if (readSmallData) {
        fileName <- paste0("measureSmallDataOriginal.txt")
    }
    if (readBigData) {
        fileName <- paste0("measureBigDataOriginal.txt")
    }
    
    minThreshold <- round(min(simMatrix)) - 1
    maxThreshold <- round(max(simMatrix)) + 1
    doTclustWithinRange(simMatrix, proteins, randomSim, minThreshold, maxThreshold, step, writeToFile, fileName)
}

if (randomSim) {
    step <- 1
    writeToFile <- TRUE
    
    # minDimensions <- 5
    # maxDimensions <- 100
    # stepDimensions <- 5
    # dimensions <- c(80, 85, 90, 95, 100)
    # seeds <- c(50)
    # for (seed in seeds) {
    #     for (dim in dimensions) {
    #         # for (dimensions in seq(minDimensions, maxDimensions, stepDimensions)) {
    #         if (readSmallData) {
    #             fileName <- paste0("AP2outputRandomSmallData-S", seed, "-K", dim, ".txt")
    #         }
    #         if (readBigData) {
    #             fileName <- paste0("AP2outputRandomBigData-S", seed, "-K", dim, ".txt")
    #         }
    # 
    #         simMatrixRandom <- buildRandomSimMatrixAp2(proteins, simMatrix, dim, seed)
    #         minThreshold <- round(min(simMatrixRandom)) - 1 # -1 to make sure to get all singletons on first run
    #         maxThreshold <- round(max(simMatrixRandom)) + 1 # +1 to make sure to get 1 cluster with all proteins
    #         doTclustWithinRange(simMatrixRandom, proteins, randomSim, minThreshold, maxThreshold, step, writeToFile, fileName)
    #     }
    # }
    
    dimensions <- c(seq(15, 100, 10))
    seeds <- c(7, 21, 42, 50)
    for (seed in seeds) {
        for (dim in dimensions) {
            if (readSmallData) {
                fileName <- paste0("AP3outputRandomSmallData-S", seed, "-K", dim, ".txt")
            }
            if (readBigData) {
                fileName <- paste0("AP3outputRandomBigData-S", seed, "-K", dim, ".txt")
            }

            simMatrixRandom <- buildRandomSimMatrixAp3(proteins, simMatrix, dim, seed)
            minThreshold <- round(min(simMatrixRandom)) - 1 # -1 to make sure to get all singletons on first run
            maxThreshold <- round(max(simMatrixRandom)) + 1 # +1 to make sure to get 1 cluster with all proteins
            doTclustWithinRange(simMatrixRandom, proteins, randomSim, minThreshold, maxThreshold, step, writeToFile, fileName)
        }
    }

    for (seed in seeds) {
        for (dim in dimensions) {
            # for (dimensions in seq(minDimensions, maxDimensions, stepDimensions)) {
            if (readSmallData) {
                fileName <- paste0("AP4outputRandomSmallData-S", seed, "-K", dim, ".txt")
            }
            if (readBigData) {
                fileName <- paste0("AP4outputRandomBigData-S", seed, "-K", dim, ".txt")
            }

            simMatrixRandom <- buildRandomSimMatrixAp4(proteins, simMatrix, dim, seed)
            minThreshold <- round(min(simMatrixRandom)) - 1 # -1 to make sure to get all singletons on first run
            maxThreshold <- round(max(simMatrixRandom)) + 1 # +1 to make sure to get 1 cluster with all proteins
            doTclustWithinRange(simMatrixRandom, proteins, randomSim, minThreshold, maxThreshold, step, writeToFile, fileName)
        }
    }
}



