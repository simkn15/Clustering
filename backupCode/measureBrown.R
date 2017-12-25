# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# library(transclustr)
# library(ggplot2)
# library(plyr)
# library(cluster)
# library(gtools)
# library(ClusterR)
# library(stringi)

#round(x, digits = 0)


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

##########
# Read in the gold standard data
##########
gold_df_before_split <- as.data.frame(readLines("sfld_brown_et_al_amidohydrolases_families_gold_standard.txt"))
gold <- as.data.frame(stri_split_fixed(gold_df_before_split[,1], "\t", simplify = TRUE))
# Set names on columns
colnames(gold) <- c("proteins", "class")
# Sort the gold standard to match the sorting of the dataset
gold_sorted <- gold[match(proteins,gold$proteins),]
# Convert classes into unique numeric values
gsClusteringDataFrame <- as.numeric(as.factor(gold_sorted$class))

####################################################################################################
# Cluster
####################################################################################################
writeToFile <- TRUE
simRandom <- FALSE
if (writeToFile && !simRandom) { file <- file("outputWithCommonAttr.txt", "w") }
if (writeToFile && simRandom) { file <- file("outputRANDOM.txt", "w") }
if (simRandom) {
    test_dist = as.vector(simMatrix)
    simMatrix <- matrix(sample(x=test_dist, replace = TRUE, size = length(test_dist)), nrow = length(proteins), ncol = length(proteins))
}

start.time <- Sys.time()
bestClustering <- list(threshold = 0, fmeasure = 0.0, clusters = 0, cost = 0.0, common = 0)
for (threshold in seq(2, 2, 1)) {
    tclustResult <- tclust(simmatrix = simMatrix, convert_dissimilarity_to_similarity = FALSE, threshold = threshold)
    tclustResultDataFrame <- data.frame(protein = proteins, cluster = tclustResult$clusters[[1]])
    # ## Add 1 to all classes to match first class = 1 in gold standard
    tclustResultDataFrame$cluster = tclustResultDataFrame$cluster + 1
    
    # fmeasure <- getQualityOfClustering(clusteringResultDataFrame, gsClusteringDataFrame)
    fmeasure <- getQualityOfClustering(tclustResultDataFrame$cluster, gsClusteringDataFrame)
    
    clusters <- length(table(tclustResultDataFrame$cluster))
    cost <- tclustResult$costs[[1]]
    stringOne <- paste0( "Threshold = ", threshold, " | F-measure = ", round(fmeasure$fmeasure, 7), " | clusters = ", clusters, " | common = ", fmeasure$common, " | cost = ", cost)
    print(stringOne)
    if (writeToFile) { write(stringOne, file, append = TRUE) }

    if (fmeasure$fmeasure > bestClustering$fmeasure) {
        bestClustering$fmeasure <- fmeasure$fmeasure
        bestClustering$threshold <- threshold
        bestClustering$clusters <- clusters
        bestClustering$cost <- cost
        bestClustering$common <- fmeasure$common
    }
}
stringTwo <- paste("Best Clustering:")
print(stringTwo)
stringThree <- paste0("Threshold = ", bestClustering$threshold, " | F-measure = ", round(bestClustering$fmeasure, 7), " | clusters = ", bestClustering$clusters, " | common = ", bestClustering$common, " | cost = ", bestClustering$cost)
print(stringThree)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
timeTaken <- paste0("Total time for clustering: ", time.taken, " seconds")
if (writeToFile) {
    write(stringTwo, file, append = TRUE)
    write(stringThree, file, append = TRUE)
    write(timeTaken, file, append = TRUE)
    close(file)
}

