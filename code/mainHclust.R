setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("fmeasure.R")
source("utilities.R")
source("hclustDivisive.R")

####################################################################################################
# Read in data and build similarity matrix
# hclust
####################################################################################################
# Adjust variables to desired bounds for binary search
# The binary search uses a upperBound and a lowerBound to find the proper threshold in the current iteration. 
step <- 30 # Default step for threshold. Binary search is skipped if 'step <- 1'. 30 is used for the binary split
minSplit <- 2 # Minimum number of splits. Minimum bound will increase if splits are < 2
maxSplit <- 10 # Maximum number of splits. Upper bound will decrease if splits are >= maxSplit.
readSmallDataSet <- TRUE
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

hc <- hclustDivisive(simMatrix, proteins, step, minSplit, maxSplit, minThreshold, maxThreshold, TRUE, FALSE)
plot(hc, xlab = "protein", ylab = "height")
plot(hc, hang = -1, xlab = "Protein", ylab = "Height")


# # Can't plot without making a hc object and reassign
# hc <- hclust(dist(USArrests), "ave")
# hc$merge <- mergeMatrix
# hc$height <- mergeHeights
# hc$order <- order
# hc$labels <- labels
# plot(hc, xlab = "protein", ylab = "height")
# plot(hc, hang = -1, xlab = "Protein", ylab = "Height")

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
