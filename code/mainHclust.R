# Main file for running hierarchical clustering methods.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utilities.R")
source("hclustDivisive.R")
source("hcGap.R")
source("plots.R")

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

# hc <- hclustDivisive(simMatrix, proteins, step, minSplit, maxSplit, minThreshold, maxThreshold, binarySearch = FALSE)
# hc <- hclustDivisiveGap(simMatrix, proteins, step, minSplit, maxSplit, minThreshold, maxThreshold, binarySearch = FALSE, GAP = TRUE, dimensions = 5, seed = 42, randomAp = 4)
# plot(hc$hc, xlab = "protein", ylab = "height")
# plot(hc$hc, hang = -1, xlab = "Protein", ylab = "Height")


# g <- plotCostDifference(hc4$gap)
# ggdraw(g)
# pdf("xxxxxplotCostGapSmallDataStep5K5Seed42NoBinary.pdf")
# print(ggdraw(g))
# dev.off()

# # Can't plot without making a hc object and reassign
# hc <- hclust(dist(USArrests), "ave")
# hc$merge <- mergeMatrix
# hc$height <- mergeHeights
# hc$order <- order
# hc$labels <- labels
# plot(hc, xlab = "protein", ylab = "height")
# plot(hc, hang = -1, xlab = "Protein", ylab = "Height")

# hclust original dataset
# hc <- hclust((max(as.dist(simMatrix)) + 1) - as.dist(simMatrix))
# plot(hc, xlab = "protein", ylab = "threshold")
# plot(hc, hang = -1)



