# Main file for running hierarchical clustering methods.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utilities.R")
source("hclustDivisive.R")
source("hcGap.R")
library(reshape)

plotCostDifference <- function(gap) {
    thresholds <- unique(gap$threshold)
    
    df_plot <- data.frame(threshold = integer(0), costDiff = integer(0), costDiffScaledAverage = integer(0), costDiffScaledPairs = integer(0), foldChange = integer(0))
    for (i in 1:length(thresholds)) {
        t <- thresholds[i]
        rows <- gap[gap$threshold == t,]
        sumCostOriginal <- 0
        sumCostRandom <- 0
        sumNumberOfProteins <- 0
        for (row in 1:nrow(rows)) {
            sumCostOriginal <- sumCostOriginal + rows[row,]$costOriginal
            sumCostRandom <- sumCostRandom + rows[row,]$costRandom
            sumNumberOfProteins <- sumNumberOfProteins + rows[row,]$numberOfProteins
        }
        costDiff <- abs(sumCostRandom - sumCostOriginal)
        costDiffScaledAverage <- costDiff / sumNumberOfProteins
        costDiffScaledPairs <- costDiff / ((sumNumberOfProteins * (sumNumberOfProteins - 1)) / 2)
        if (sumCostOriginal == 0 || sumCostRandom == 0) {
            foldChange <- 0
        }
        else {
            foldChange <- foldchange(sumCostRandom, sumCostOriginal)
        }
        
        # Change variable names
        df <- data.frame(threshold = t, actual = costDiff, avg = costDiffScaledAverage, avgPairs = costDiffScaledPairs, foldChange = foldChange)
        df_plot <- rbind(df_plot, df)
    }
    
    df_plot <- melt(df_plot, id="threshold")
    negativeFoldChange <- which(df_plot[,3] < 0)
    for (i in 1:nrow(df_plot)) {
        if (df_plot[i,]$value == 0) {
            df_plot[i,]$value <- 0.001
        }
        if (df_plot[i,]$value < 0 && df_plot[i,]$variable == "foldChange") {
            df_plot[i,]$value <- abs(df_plot[i,]$value)
        }
    }
    
    title <- paste0("Cost difference : actual vs. randomized")
    g <- ggplot() + 
        labs(col = "Variable:") +
        scale_x_continuous(name = "Threshold", breaks = seq(0, max(df_plot[,1]), 20)) +
        scale_y_continuous(name = "Cost", trans='log2', labels=function(n){format(n, scientific = FALSE, digits = 1)}, breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) + 
        geom_line(data = df_plot, aes(x = threshold, y = value, colour = variable)) +
        geom_point(data = df_plot[negativeFoldChange,], aes(x = threshold, y = value)) +
        theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
        ggtitle(title)
    g <- add_sub(g, "Black point on foldChange is conversion from negative to positive.")
    return(g)
}

####################################################################################################
# Read in data and build similarity matrix
# hclust
####################################################################################################
# Adjust variables to desired bounds for binary search
# The binary search uses a upperBound and a lowerBound to find the proper threshold in the current iteration.
# randomApproaches <- c(3, 4)
# seeds <- c(7, 21, 42)
# minSplit <- 2 # Minimum number of splits. Minimum bound will increase if splits are < 2
# maxSplit <- 10 # Maximum number of splits. Upper bound will decrease if splits are >= maxSplit.
# 
# table <- read.table("./data/brown/sfld_brown_et_al_amidohydrolases_protein_similarities_for_beh.txt", sep = "", skip = 5)
# df_table <- as.data.frame(table)
# proteins = levels(df_table[,1])
# simMatrix <- buildSimilarityMatrix(proteins, df_table)
# minThreshold <- min(simMatrix) - 1 # The starting threshold. -1 such that the first run returns 1 cluster.
# maxThreshold <- round(max(simMatrix)) + 1 # Threshold where all clusters will be singletons
# 
# steps <- c(1, 5)
# for (step in steps) {
#     for (seed in seeds) {
#         for (ap in randomApproaches) {
#             print(paste0("Doing ap ", ap, " no binary on small data"))
#             hc <- hclustDivisiveGap(simMatrix, proteins, step, minSplit, maxSplit, minThreshold, maxThreshold, binarySearch = FALSE, GAP = TRUE, dimensions = 5, seed = seed, randomAp = ap)
#             g <- plotCostDifference(hc$gap)
#             ggdraw(g)
#             fileName <- paste0("plotCostGapSmallDataAp", ap, "Step", step, "K5Seed", seed, "NoBinary.pdf")
#             workspaceName <- paste0("hcGapSmallDataAp", ap, "Step", step, "K5Seed", seed, "NoBinary.RData")
#             pdf(fileName)
#             print(ggdraw(g))
#             dev.off()
#             save.image(workspaceName)
#         }
#     }
# }
# 
# steps <- c(30)
# for (step in steps) {
#     for (seed in seeds) {
#         for (ap in randomApproaches) {
#             print(paste0("Doing ap ", ap, " binary on small data"))
#             hc <- hclustDivisiveGap(simMatrix, proteins, step, minSplit, maxSplit, minThreshold, maxThreshold, binarySearch = TRUE, GAP = TRUE, dimensions = 5, seed = seed, randomAp = ap)
#             g <- plotCostDifference(hc$gap)
#             ggdraw(g)
#             fileName <- paste0("plotCostGapSmallDataAp", ap, "Step", step, "KS5eed", seed, "Binary.pdf")
#             workspaceName <- paste0("hcGapSmallDataAp", ap, "Step", step, "K5Seed", seed, "Binary.RData")
#             pdf(fileName)
#             print(ggdraw(g))
#             dev.off()
#             save.image(workspaceName)
#         }
#     }
# }
# 
# 
# table <- read.table("./data/big/simBig.txt", sep = "")
# df_table <- as.data.frame(table)
# proteins = levels(df_table[,1])
# simMatrix <- buildSimilarityMatrixFromBlast(proteins, df_table)
# minThreshold <- min(simMatrix) - 1
# maxThreshold <- round(max(simMatrix)) + 1
# 
# steps <- c(1, 5)
# for (step in steps) {
#     for (seed in seeds) {
#         for (ap in randomApproaches) {
#             print(paste0("Doing ap ", ap, " no binary on big data"))
#             hc <- hclustDivisiveGap(simMatrix, proteins, step, minSplit, maxSplit, minThreshold, maxThreshold, binarySearch = FALSE, GAP = TRUE, dimensions = 5, seed = seed, randomAp = ap)
#             g <- plotCostDifference(hc$gap)
#             ggdraw(g)
#             fileName <- paste0("plotCostGapBigDataAp", ap, "Step", step, "K5Seed", seed, "NoBinary.pdf")
#             workspaceName <- paste0("hcGapBigDataAp", ap, "Step", step, "K5Seed", seed, "NoBinary.RData")
#             pdf(fileName)
#             print(ggdraw(g))
#             dev.off()
#             save.image(workspaceName)
#         }
#     }
# }
# 
# 
# steps <- c(30)
# for (step in steps) {
#     for (seed in seeds) {
#         for (ap in randomApproaches) {
#             print(paste0("Doing ap ", ap, " binary on big data"))
#             hc <- hclustDivisiveGap(simMatrix, proteins, step, minSplit, maxSplit, minThreshold, maxThreshold, binarySearch = TRUE, GAP = TRUE, dimensions = 5, seed = seed, randomAp = ap)
#             g <- plotCostDifference(hc$gap)
#             ggdraw(g)
#             fileName <- paste0("plotCostGapBigDataAp", ap, "Step", step, "K5Seed", seed, "Binary.pdf")
#             workspaceName <- paste0("hcGapBigDataAp", ap, "Step", step, "K5Seed", seed, "Binary.RData")
#             pdf(fileName)
#             print(ggdraw(g))
#             dev.off()
#             save.image(workspaceName)
#         }
#     }
# }






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
# hc3 <- hclustDivisiveGap(simMatrix, proteins, step, minSplit, maxSplit, minThreshold, maxThreshold, binarySearch = FALSE, GAP = TRUE, dimensions = 5, seed = 42, randomAp = 3)
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

# Fails
# hct <- list(merge = mergeMatrix, height = height, order = order, labels = labels)
# plot(hct, xlab = "Protein", ylab = "Threshold")
# plot(hct, hang = -1)

# hclust original dataset
# hc <- hclust((max(as.dist(simMatrix)) + 1) - as.dist(simMatrix))
# plot(hc, xlab = "protein", ylab = "threshold")
# plot(hc, hang = -1)

# Approx optimal threshold
# > mean(as.vector(as.dist(simMatrix)))
#[1] 55.91576

# hc <- hclust(dist(USArrests), "ave")
# plot(hc)
# plot(hc, hang = -1)

########### Code used to summarizing plots for the THC result section
# Load workspace, run line code.
gap <- hc$gap
# Run the next line for each loaded workspace
gap[,2] <- gap[,2] + hc$gap[,2]

# Make column 2 the average for the 9 results
gap[,2] <- gap[,2] / 9
g <- plotCostDifference(gap)
pdf("THCAp3SmallData.pdf")
print(ggdraw(g))
dev.off()





