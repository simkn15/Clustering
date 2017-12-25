setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(transclustr)
library(ggplot2)
library(plyr)
library(cluster)
library(gtools)
library(ClusterR)
library(stringi)

clusteringWithTclust <- function(simMatrix, proteinLabels, threshold) {
    tclustResult <- tclust(simmatrix = simMatrix, convert_dissimilarity_to_similarity = FALSE, threshold = threshold)
    tclustResultDataFrame <- data.frame(protein = proteinLabels, cluster = tclustResult$clusters[[1]])
    tclustResultDataFrame$cluster = tclustResultDataFrame$cluster + 1
    
    return(tclustResultDataFrame)
}

getProteinLabelsFromClustering <- function(tclustRes) {
    amountOfClasses <- length(table(tclustRes$cluster))
    
    # Initialize each cluster object
    proteinLabels <- list()
    for (i in 1:amountOfClasses) {
        proteinLabels[[i]] <- list(cid = "", proteins = c(), startIndex = 0, endIndex = 0, height = 0, parent = "")
    }
    
    # Collect labels for all clusters
    for (i in 1:length(tclustRes$cluster)) {
        labelCurrentProtein <- tclustRes$protein[i] # change 'tclustRes$protein[i]' to 'i', to get protein number instead of label
        classCurrentProtein <- tclustRes$cluster[i]
        proteinLabels[[classCurrentProtein]]$proteins <- c(proteinLabels[[classCurrentProtein]]$proteins, c(as.character(labelCurrentProtein)))
    }
    
    return(proteinLabels)
}

# Used for the small data set
buildSimilarityMatrix <- function(proteins, df_table) {
    if (length(proteins) < 2) { stop("Cannot build similarity matrix with less than 2 proteins")}
    # Create an empty matrix
    simMatrix = matrix(0, nrow = length(proteins), ncol = length(proteins))

    # Set row and column names for readability
    row.names(simMatrix) = proteins
    colnames(simMatrix) = proteins
    
    # Extract the pairwise similarities
    for(x in 1:(length(proteins)-1)){
        # Select all rows containing p1
        p1 = proteins[x]
        presec = df_table[df_table[,1] == p1 | df_table[,2] == p1, ]
        
        for(y in (x+1):length(proteins)){
            p2 = proteins[y]
            # We have to check for both directions p1 -> p2 and p2 <- p1. Rule is, to be more conservative, we take the minimum of both values
            if(length(presec[presec[,1] == p2, 3])==0 | length(presec[presec[,2] == p2, 3])==0){
                #Allright, one value was missing, we take 0 as fallback (since there is no meaningful minimum)
                simMatrix[x,y] = simMatrix[y,x] = 0
            }else{
                #Take the minimum of both directions
                simMatrix[x,y] = simMatrix[y,x] = min(presec[presec[,1] == p2, 3], presec[presec[,2] == p2, 3])
            }
        }
    }
    # Set the diagonal to the highest observed value (not actually necessary, since self-similarity is never queried)
    # diag(simMatrix) = max(simMatrix)
    
    return(simMatrix)
}

# Used for the big data set
buildSimilarityMatrixFromBlast <- function(proteins, df_table) {
    if (length(proteins) < 2) { stop("Cannot build similarity matrix with less than 2 proteins")}
    # Create an empty matrix
    simMatrix = matrix(0, nrow = length(proteins), ncol = length(proteins))
    
    # Set row and column names for readability
    row.names(simMatrix) = proteins
    colnames(simMatrix) = proteins
    
    # Extract the pairwise similarities
    for(x in 1:(length(proteins)-1)){
        # Select all rows containing p1
        p1 = proteins[x]
        presec = df_table[df_table[,1] == p1 | df_table[,2] == p1, ]
        
        for(y in (x+1):length(proteins)){
            p2 = proteins[y]
            # We have to check for both directions p1 -> p2 and p2 <- p1. Rule is, to be more conservative, we take the minimum of both values
            d1.set <- presec[presec[,1] == p2, 3]
            d2.set <- presec[presec[,2] == p2, 3]
            d1.length <- length(presec[presec[,1] == p2, 3])
            d2.length <- length(presec[presec[,2] == p2, 3])
            if (d1.length == 0 | d2.length == 0){
                #Allright, one value was missing, we take 0 as fallback (since there is no meaningful minimum)
                simMatrix[x,y] = simMatrix[y,x] = 0
            }else{
                #Take the minimum of both directions
                d1.maxValue <- 0
                d2.maxValue <- 0
                for (i in 1:d1.length) { # find d1.maxValue
                    if (d1.set[i] > d1.maxValue) {
                        d1.maxValue <- d1.set[i]
                    }
                }
                for (i in 1:d2.length) { # find d2.maxValue
                    if (d2.set[i] > d2.maxValue) {
                        d2.maxValue <- d2.set[i]
                    }
                }
                
                simMatrix[x,y] = simMatrix[y,x] = min(d1.maxValue, d2.maxValue)
            }
        }
    }
    # Set the diagonal to the highest observed value (not actually necessary, since self-similarity is never queried)
    # diag(simMatrix) = max(simMatrix)
    
    return(simMatrix)
}

plotCostsSmall <- function() {
    dimensions <- c(seq(10,100,10))
    dimensionsAsChar <- as.character(dimensions)
    df <- data.frame(threshold = integer(0), cost = integer(0), dimension = character())
    names <- c("threshold", "cost")
    # Read original data
    table <- read.table(paste0("./measure/measureSmallDataOriginal.txt"), sep = "", fill = TRUE)
    df_table <- as.data.frame(table)
    rows <- nrow(df_table) - 3 # Skip last 3 lines
    columns <- c(3,19)
    df_table <- df_table[1:rows,]
    colThreshold <- sort(df_table[,3])
    df_table <- as.data.frame(cbind(colThreshold, df_table[,19]))
    colnames(df_table) <- names
    # df_table <- df_table[order(df_table$threshold),]
    df_table$dimension <- rep("2(Original data)", rows)
    
    df <- rbind(df, df_table)
    for (i in 1:length(dimensions)) {
        table <- read.table(paste0("./randomSmall/outputRandomSmallData-S42-K", dimensions[i], ".txt"), sep = "")
        df_table <- as.data.frame(table)
        df_table <- as.data.frame(cbind(df_table[,3], df_table[,19]))
        colnames(df_table) <- names
        
        df_table$dimension <- rep(dimensionsAsChar[i], nrow(df_table))
        
        df <- rbind(df, df_table)
    }
    
    title <- paste0("Small Data Set")
    g <- ggplot(df, aes(x = threshold, y = cost, color = dimension)) + ggtitle(title) + geom_line()
    return(g)
}

plotCostsBig <- function() {
    dimensions <- c(seq(10,100,10))
    # dimensions <- c(5) # use vector of size one for one vs. one plot
    dimensionsAsChar <- as.character(dimensions)
    df <- data.frame(threshold = integer(0), cost = integer(0), dimension = character())
    names <- c("threshold", "cost")
    for (i in 1:length(dimensions)) {
        table <- read.table(paste0("./randomBig/outputRandomBigData-S42-K", dimensions[i], ".txt"), sep = "")
        df_table <- as.data.frame(table)
        df_table <- as.data.frame(cbind(df_table[,3], df_table[,7]))
        colnames(df_table) <- names
        
        df_table$dimension <- rep(dimensionsAsChar[i], nrow(df_table))
        
        df <- rbind(df, df_table)
    }
    
    # Read original data
    table <- read.table(paste0("./measure/measureBigDataOriginal.txt"), sep = "", fill = TRUE)
    df_table <- as.data.frame(table)
    rows <- nrow(df_table) - 3 # Skip last 3 lines
    columns <- c(3,19)
    df_table <- df_table[1:rows,]
    colThreshold <- sort(df_table[,3])
    df_table <- as.data.frame(cbind(colThreshold, df_table[,19]))
    colnames(df_table) <- names
    # df_table <- df_table[order(df_table$threshold),]
    df_table$dimension <- rep("2(Original data)", rows)
    
    df <- rbind(df, df_table)
    
    title <- paste0("Big Data Set")
    g <- ggplot(df, aes(x = threshold, y = cost, color = dimension)) + ggtitle(title) + geom_line()
    return(g)
}

plotCostsSmallOneVsOneWithGap <- function(dimensions = c(5)) {
    dimensionsAsChar <- as.character(dimensions)
    df_plot <- data.frame(threshold = integer(0), cost = integer(0), dimension = character())
    names <- c("threshold", "cost")
    # Read original data
    table <- read.table(paste0("./measure/measureSmallDataOriginal.txt"), sep = "", fill = TRUE)
    df_original <- as.data.frame(table)
    rows <- nrow(df_original) - 3 # Skip last 3 lines
    columns <- c(3,19)
    df_original <- df_original[1:rows,]
    colThreshold <- sort(df_original[,3])
    df_original <- as.data.frame(cbind(colThreshold, df_original[,19]))
    colnames(df_original) <- names
    # df_table <- df_table[order(df_table$threshold),]
    df_original$dimension <- rep("2(Original data)", rows)
    
    df_plot <- rbind(df_plot, df_original)
    
    table <- read.table(paste0("./randomSmall/outputRandomSmallData-S42-K", dimensions[1], ".txt"), sep = "")
    df_random <- as.data.frame(table)
    df_random <- as.data.frame(cbind(df_random[,3], df_random[,19]))
    colnames(df_random) <- names
    
    df_random$dimension <- rep(dimensionsAsChar[1], nrow(df_random))
    
    df_plot <- rbind(df_plot, df_random)
    
    # Find biggest GAP original vs. random
    rows <- min(nrow(df_original), nrow(df_random))
    maxGapThreshold <- 0
    maxGapCost <- 0
    for (row in 1:rows) {
        gap <- abs(df_original[row,2] - df_random[row,2])
        if (gap > maxGapCost) {
            maxGapThreshold <- df_original[row,1]
            maxGapCost <- gap
        }
    }
    
    title <- paste0("Small Data Set : Biggest gap at threshold ", maxGapThreshold)
    g <- ggplot(df_plot, aes(x = threshold, y = cost, color = dimension)) + ggtitle(title) + geom_line() + geom_vline(xintercept = maxGapThreshold, linetype="dotted")
    g
    return(g)
}

############# One vs. One plotting of costs with biggest gap and foldChange
plotCostsBigOneVsOneWithGap <- function(dimensions = c(5)) {
    dimensionsAsChar <- as.character(dimensions)
    df_plot <- data.frame(threshold = integer(0), cost = integer(0), dimension = character())
    names <- c("threshold", "cost")
    # Read original data
    table <- read.table(paste0("./measure/measureSmallDataOriginal.txt"), sep = "", fill = TRUE)
    df_original <- as.data.frame(table)
    rows <- nrow(df_original) - 3 # Skip last 3 lines
    columns <- c(3,19)
    df_original <- df_original[1:rows,]
    colThreshold <- sort(df_original[,3])
    df_original <- as.data.frame(cbind(colThreshold, df_original[,19]))
    colnames(df_original) <- names
    # df_table <- df_table[order(df_table$threshold),]
    df_original$dimension <- rep("2(Original data)", rows)
    
    df_plot <- rbind(df_plot, df_original)
    
    table <- read.table(paste0("./randomSmall/outputRandomSmallData-S42-K", dimensions[1], ".txt"), sep = "")
    df_random <- as.data.frame(table)
    df_random <- as.data.frame(cbind(df_random[,3], df_random[,19]))
    colnames(df_random) <- names
    
    df_random$dimension <- rep(dimensionsAsChar[1], nrow(df_random))
    
    df_plot <- rbind(df_plot, df_random)
    
    # Find biggest GAP original vs. random
    rows <- min(nrow(df_original), nrow(df_random))
    maxGapThreshold <- 0
    maxGapCost <- 0
    fold <- 0
    for (row in 1:rows) {
        costOriginal <- df_original[row,2]
        costRandom <- df_random[row,2]
        gap <- abs(costRandom - costOriginal)
        if (gap > maxGapCost) {
            fold <- foldchange(costRandom, costOriginal)
            maxGapThreshold <- df_original[row,1]
            maxGapCost <- gap
        }
    }
    
    title <- paste0("Big Data Set : Biggest gap threshold = ", maxGapThreshold, " | foldChange = ", round(fold, digits = 4))
    g <- ggplot(df_plot, aes(x = threshold, y = cost, color = dimension)) + ggtitle(title) + geom_line() + geom_vline(xintercept = maxGapThreshold, linetype="dotted")
    g
    return(g)
}
# 
# dimensions = c(seq(5, 100, 5))
# for (d in dimensions) {
#     dd <- c(d)
#     dataPlot <- plotCostsBigOneVsOneWithGap(dd)
#     dataPlot
#     fileName <- paste0("costsPlotSmallDataK", d, ".pdf")
#     pdf(fileName)
#     print(dataPlot)
#     dev.off()
# }
############# One vs. One plotting of costs with biggest gap and foldChange


# smallDataPlot <- plotCostsSmall()
# smallDataPlot
# bigDataPlot <- plotCostsBig()
# bigDataPlot
# # 
# pdf("costsPlotSmallDataStep10.pdf")
# print(smallDataPlot)
# dev.off()
# # 
# pdf("costsPlotBigDataStep10.pdf")
# print(bigDataPlot)
# dev.off()








