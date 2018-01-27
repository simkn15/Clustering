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
    
    title <- paste0("Cost difference : orignal vs. randomized")
    g <- ggplot() + 
        xlab("Threshold") + 
        scale_y_continuous(name = "Cost", trans='log2', labels=function(n){format(n, scientific = FALSE, digits = 2)}) + 
        scale_x_continuous(breaks = seq(0, max(df_plot[,1]), 20)) +
        geom_line(data = df_plot, aes(x = threshold, y = value, colour = variable)) +
        geom_point(data = df_plot[negativeFoldChange,], aes(x = threshold, y = value)) +
        ggtitle(title) +
        theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
        labs(col = "Variable:")
    g <- add_sub(g, "Black point on foldChange is conversion from negative to positive.")
    return(g)
}

####################################################################################################
# Read in data and build similarity matrix
# hclust
####################################################################################################
# Adjust variables to desired bounds for binary search
# The binary search uses a upperBound and a lowerBound to find the proper threshold in the current iteration.
randomApproaches <- c(4, 3)
seeds <- c(7, 21, 42)
minSplit <- 2 # Minimum number of splits. Minimum bound will increase if splits are < 2
maxSplit <- 10 # Maximum number of splits. Upper bound will decrease if splits are >= maxSplit.
DEBUG <- FALSE

table <- read.table("./data/brown/sfld_brown_et_al_amidohydrolases_protein_similarities_for_beh.txt", sep = "", skip = 5)
df_table <- as.data.frame(table)
proteins = levels(df_table[,1])
simMatrix <- buildSimilarityMatrix(proteins, df_table)
minThreshold <- min(simMatrix) - 1 # The starting threshold. -1 such that the first run returns 1 cluster.
maxThreshold <- round(max(simMatrix)) + 1 # Threshold where all clusters will be singletons

steps <- c(1, 5)
for (ap in randomApproaches) {
    for (step in steps) {
        for (seed in seeds) {
            hc <- hclustDivisiveGap(simMatrix, proteins, step, minSplit, maxSplit, minThreshold, maxThreshold, binarySearch = FALSE, GAP = TRUE, dimensions = 5, seed = seed, randomAp = ap)
            g <- plotCostDifference(hc$gap)
            ggdraw(g)
            fileName <- paste0("plotCostGapSmallDataStep", step, "K5Seed", seed, "NoBinaryAp", ap, ".pdf")
            workspaceName <- paste0("hcGapSmallDataStep", step, "K5Seed", seed, "NoBinaryAp", ap, ".RData")
            pdf(fileName)
            print(ggdraw(g))
            dev.off()
            save.image(workspaceName)
        }
    }
}


step <- c(20, 30)
for (ap in randomApproaches) {
    for (step in steps) {
        for (seed in seeds) {
            hc <- hclustDivisiveGap(simMatrix, proteins, step, minSplit, maxSplit, minThreshold, maxThreshold, binarySearch = TRUE, GAP = TRUE, dimensions = 5, seed = seed, randomAp = ap)
            g <- plotCostDifference(hc$gap)
            ggdraw(g)
            fileName <- paste0("plotCostGapSmallDataStep", step, "K5Seed", seed, "BinaryAp", ap, ".pdf")
            workspaceName <- paste0("hcGapSmallDataStep", step, "K5Seed", seed, "BinaryAp", ap, ".RData")
            pdf(fileName)
            print(ggdraw(g))
            dev.off()
            save.image(workspaceName)
        }
    }  
}


table <- read.table("./data/big/simBig.txt", sep = "")
df_table <- as.data.frame(table)
proteins = levels(df_table[,1])
simMatrix <- buildSimilarityMatrixFromBlast(proteins, df_table)
minThreshold <- min(simMatrix) - 1
maxThreshold <- round(max(simMatrix)) + 1

steps <- c(1, 5)
for (ap in randomApproaches) {
    for (step in steps) {
        for (seed in seeds) {
            hc <- hclustDivisiveGap(simMatrix, proteins, step, minSplit, maxSplit, minThreshold, maxThreshold, binarySearch = FALSE, GAP = TRUE, dimensions = 5, seed = seed, randomAp = ap)
            g <- plotCostDifference(hc$gap)
            ggdraw(g)
            fileName <- paste0("plotCostGapBigDataStep", step, "K5Seed", seed, "NoBinaryAp", ap, ".pdf")
            workspaceName <- paste0("hcGapBigDataStep", step, "K5Seed", seed, "NoBinaryAp", ap, ".RData")
            pdf(fileName)
            print(ggdraw(g))
            dev.off()
            save.image(workspaceName)
        }
    }
}


steps <- c(30)
for (ap in randomApproaches) {
    for (step in steps) {
        for (seed in seeds) {
            hc <- hclustDivisiveGap(simMatrix, proteins, step, minSplit, maxSplit, minThreshold, maxThreshold, binarySearch = TRUE, GAP = TRUE, dimensions = 5, seed = seed, randomAp = ap)
            g <- plotCostDifference(hc$gap)
            ggdraw(g)
            fileName <- paste0("plotCostGapBigDataStep", step, "K5Seed", seed, "BinaryAp", ap, ".pdf")
            workspaceName <- paste0("hcGapBigDataStep", step, "K5Seed", seed, "BinaryAp", ap, ".RData")
            pdf(fileName)
            print(ggdraw(g))
            dev.off()
            save.image(workspaceName)
        }
    }
}






# step <- 5 # Default step for threshold. Binary search is skipped if 'step <- 1'. 30 is used for the binary split
# minSplit <- 5 # Minimum number of splits. Minimum bound will increase if splits are < 2
# maxSplit <- 10 # Maximum number of splits. Upper bound will decrease if splits are >= maxSplit.
# readSmallDataSet <- TRUE
# readBigDataSet <- !readSmallDataSet
# DEBUG <- FALSE
# 
# if (readSmallDataSet) {
#     table <- read.table("./data/brown/sfld_brown_et_al_amidohydrolases_protein_similarities_for_beh.txt", sep = "", skip = 5)
#     df_table <- as.data.frame(table)
#     proteins = levels(df_table[,1])
#     simMatrix <- buildSimilarityMatrix(proteins, df_table)
# 
#     minThreshold <- min(simMatrix) - 1 # The starting threshold. -1 such that the first run returns 1 cluster.
#     maxThreshold <- round(max(simMatrix)) + 1 # Threshold where all clusters will be singletons
# 
#     if (DEBUG) { fileName <- "hcSmallData" }
# }
# if (readBigDataSet) {
#     table <- read.table("./data/big/simBig.txt", sep = "")
#     df_table <- as.data.frame(table)
#     proteins = levels(df_table[,1])
#     simMatrix <- buildSimilarityMatrixFromBlast(proteins, df_table)
# 
#     minThreshold <- min(simMatrix) - 1
#     maxThreshold <- round(max(simMatrix)) + 1
# 
#     if (DEBUG) { fileName <- "hcBigData" }
# }
# 
# # hc <- hclustDivisive(simMatrix, proteins, step, minSplit, maxSplit, minThreshold, maxThreshold, binarySearch = FALSE)
# hc <- hclustDivisiveGap(simMatrix, proteins, step, minSplit, maxSplit, minThreshold, maxThreshold, binarySearch = TRUE, GAP = TRUE, dimensions = 5)
# # plot(hc$hc, xlab = "protein", ylab = "height")
# plot(hc$hc, hang = -1, xlab = "Protein", ylab = "Height")
# g <- plotCostDifference(hc$gap)
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
# plot(hc, hang = -1, xlab = "Protein", ylab = "Threshold")

# Approx optimal threshold
# > mean(as.vector(as.dist(simMatrix)))
#[1] 55.91576

# hc <- hclust(dist(USArrests), "ave")
# plot(hc)
# plot(hc, hang = -1)
