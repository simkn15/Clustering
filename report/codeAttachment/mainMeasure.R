setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(transclustr)
library(ggplot2)
library(plyr)
library(cluster)
library(gtools)
library(ClusterR)
library(stringi)

source("utilities.R")
source("randomization.R")
source("measureTclust.R")

####################################################################################################
# Read in data
# Build similarity matrix
# Update parameters according to desired run
# Run doTclustWithinRange()
####################################################################################################
randomSim <- FALSE
readSmallData <- TRUE
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
    gold[,2] <- NULL # Deleting second column
    colnames(gold) <- c("proteins", "class")
    gold_sorted <- gold[match(proteins,gold$proteins),]
    gsClusteringDataFrame <- as.numeric(as.factor(gold_sorted$class))
}

if (!randomSim) {
    step <- 1
    writeToFile <- TRUE
    
    if (readSmallData) {
        fileName <- paste0("measureSmallDataOriginalFixedColumns.txt")
    }
    if (readBigData) {
        fileName <- paste0("measureBigDataOriginalFixedColumns.txt")
    }
    
    minThreshold <- round(min(simMatrix)) - 1
    maxThreshold <- round(max(simMatrix)) + 1
    doTclustWithinRange(simMatrix, proteins, randomSim, minThreshold, maxThreshold, step, writeToFile, fileName)
}

if (randomSim) {
    step <- 1
    writeToFile <- TRUE
    
    minDimensions <- 5
    maxDimensions <- 100
    stepDimensions <- 5
    dimensions <- c(80, 85, 90, 95, 100)
    seeds <- c(50)
    for (seed in seeds) {
        for (dim in dimensions) {
            if (readSmallData) {
                fileName <- paste0("AP2outputRandomSmallData-S", seed, "-K", dim, ".txt")
            }
            if (readBigData) {
                fileName <- paste0("AP2outputRandomBigData-S", seed, "-K", dim, ".txt")
            }
            
            simMatrixRandom <- buildRandomSimMatrixAp2(proteins, simMatrix, dim, seed)
            minThreshold <- round(min(simMatrixRandom)) - 1 # -1 to make sure to get all singletons on first run
            maxThreshold <- round(max(simMatrixRandom)) + 1 # +1 to make sure to get 1 cluster with all proteins
            doTclustWithinRange(simMatrixRandom, proteins, randomSim, minThreshold, maxThreshold, step, writeToFile, fileName)
        }
    }
}



