setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(transclustr)
library(ggplot2)
library(plyr)
library(cluster)
library(gtools)
library(ClusterR)
library(stringi)
library(RColorBrewer)

source("randomization.R")

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

# Returns the similarity subset of simMatrix where only 'proteins' are occuring
getSimMatrixForSpecificProteins <- function(simMatrix, proteins) {
    proteins <- as.character(proteins) # Precaution that no error will occure
    return(simMatrix[proteins, proteins])
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
# Directions can have multiple similarities
buildSimilarityMatrixFromBlast <- function(proteins, df_table) {
    if (length(proteins) < 2) { stop("Cannot build similarity matrix with less than 2 proteins")}
    # Create an empty matrix
    simMatrix = matrix(0, nrow = length(proteins), ncol = length(proteins))
    
    # Set row and column names for readability
    row.names(simMatrix) = proteins
    colnames(simMatrix) = proteins
    
    # Extract the pairwise similarities
    for(x in 1:(length(proteins)-1)) {
        # Select all rows containing p1
        p1 = proteins[x]
        presec = df_table[df_table[,1] == p1 | df_table[,2] == p1, ]
        
        for(y in (x+1):length(proteins)) {
            p2 = proteins[y]
            d1.set <- presec[presec[,1] == p2, 3]
            d2.set <- presec[presec[,2] == p2, 3]
            d1.length <- length(d1.set)
            d2.length <- length(d2.set)
            if (d1.length == 0 | d2.length == 0) {
                #Allright, one value was missing, we take 0 as fallback (since there is no meaningful minimum)
                simMatrix[x,y] = simMatrix[y,x] = 0
            } else{
                #Take the maximum of both directions
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
                # Take the minimum of both directions
                simMatrix[x,y] = simMatrix[y,x] = min(d1.maxValue, d2.maxValue)
            }
        }
    }
    # Set the diagonal to the highest observed value (not actually necessary, since self-similarity is never queried)
    # diag(simMatrix) = max(simMatrix)
    
    return(simMatrix)
}

##############################################################
# Below functions are all for plotting costs
##############################################################
plotCostsSmallOneVsOneWithGap <- function(dimensions = c(5)) { # Outdated
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
    
    table <- read.table(paste0("./randomSmallAp2/outputRandomSmallData-S42-K", dimensions[1], ".txt"), sep = "")
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
    g <- ggplot(df_plot, aes(x = threshold, y = cost, color = dimension)) +
        ggtitle(title) + 
        geom_line() + 
        geom_vline(xintercept = maxGapThreshold, linetype="dotted") +
        scale_y_continuous(labels=function(n){format(n, scientific = FALSE)})
    return(g)
}

############# One vs. One plotting of costs with biggest gap and foldChange
plotCostsBigOneVsOneWithGap <- function(dimensions = 10, seed = 7, ranApproach = 2) {
    # ranApproach <- 2
    # seed <- 7
    # dimensions <- c(10)
    dimensionsAsChar <- as.character(dimensions)
    df_plot <- data.frame(threshold = integer(0), cost = integer(0), dimension = character())
    names <- c("threshold", "cost")
    # Read original data
    table <- read.table(paste0("./measure/measureBigDataOriginal.txt"), sep = "", fill = TRUE)
    df_original <- as.data.frame(table)
    rows <- nrow(df_original) - 2 # Skip last 3 lines
    df_original <- df_original[1:rows,]
    df_original <- as.data.frame(cbind(df_original[,3], df_original[,19]))
    colnames(df_original) <- names
    df_original$dimension <- rep("2(Original data)", rows)
    df_plot <- rbind(df_plot, df_original)
    
    # Read random data
    path <- paste0("./randomBig/ap", ranApproach, "/seed", seed, "/")
    fileName <- paste0("AP", ranApproach, "outputRandomBigData-S", seed, "-K", dimensions, ".txt")
    table <- read.table(paste0(path, fileName), header = TRUE, sep = "")
    df_random <- as.data.frame(table)
    df_random[,3] <- NULL # Removing column 'clusters'
    colnames(df_random) <- names
    
    df_random$dimension <- rep(dimensionsAsChar, nrow(df_random))
    
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
    g <- ggplot(df_plot, aes(x = threshold, y = cost, color = dimension)) + 
        ggtitle(title) + geom_line() + 
        geom_vline(xintercept = maxGapThreshold, linetype="dotted") +
        scale_y_continuous(name = "Cost", labels=function(n){format(n, scientific = FALSE)}) +
        labs(color = "Dimensions")
    return(g)
}

plotHistogramRandomVsOriginalSimilarities <- function(dimension = 10, seed = 42, randomApproaches = 2) {
    colNames <- c("Similarity", "Dataset")
    df_plot <- data.frame(similarity = integer(0), dataset = character())
    # Read in Original data
    table <- read.table("./data/big/simBig.txt", sep = "")
    df_table <- as.data.frame(table)
    proteins = levels(df_table[,1])
    simMatrix <- buildSimilarityMatrixFromBlast(proteins, df_table) # Big data
    
    rows <- length(proteins)^2
    df_plot <- as.data.frame(cbind(as.vector(simMatrix), rep("Original", rows)))
    names(df_plot) <-  colNames
    
    # Read in Random
    if (randomApproaches == 2) {
        simRandom <- buildRandomSimMatrixAp2(proteins, simMatrix, dimension, seed)
        df_randomSim <- as.data.frame(cbind(as.vector(simRandom), rep("Random", rows)))
        names(df_randomSim) <- colNames
        df_plot <- rbind(df_plot, df_randomSim)
    }
    if (randomApproaches == 3) {
        simRandom <- buildRandomSimMatrixAp3(proteins, simMatrix, dimension, seed)
        df_randomSim <- as.data.frame(cbind(as.vector(simRandom), rep("Random", rows)))
        names(df_randomSim) <- colNames
        df_plot <- rbind(df_plot, df_randomSim)
    }
    if (randomApproaches == 4) {
        simRandom <- buildRandomSimMatrixAp4(proteins, simMatrix, dimension, seed)
        df_randomSim <- as.data.frame(cbind(as.vector(simRandom), rep("Random", rows)))
        names(df_randomSim) <- colNames
        df_plot <- rbind(df_plot, df_randomSim)
    }
    df_plot[,1] <- as.numeric(as.character(df_plot[,1]))

    g <- ggplot(data = df_plot, aes(x = Similarity, fill = Dataset)) + 
        geom_histogram(binwidth = 10,
                       position = "identity", colour = "black", boundary = 0) +
        scale_x_continuous(name = "Similarity", breaks = seq(0, max(df_plot$Similarity)+20, 20)) +
        scale_y_log10(name = "Frequency", labels=function(n){format(n, scientific = FALSE)}) +
        ggtitle("Frequency of similarities") +
        theme_bw() +
        theme(axis.line = element_line(size=1, colour = "black"),
              panel.grid.major = element_line(colour = "#d3d3d3"),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), panel.background = element_blank(),
              axis.text.x=element_text(colour="black", size = 9, angle = 90),
              axis.text.y=element_text(colour="black", size = 9),
              legend.position = "none") +
        facet_grid(. ~Dataset , scales = "free")
    
    return(g)
}

plotRandomVsOriginalAndPrintToFile <- function() {
    ranApproaches <- c(2, 3, 4)
    seeds <- c(7, 21, 42, 50)
    dimensions = c(2, 3, 4, 5, 10)
    for (ra in ranApproaches) {
        for (s in seeds) {
            for (d in dimensions) {
                costPlot <- plotCostsBigOneVsOneWithGap(d, s, ra)
                simPlot <- plotHistogramRandomVsOriginalSimilarities(d, s, ra)
                fileName <- paste0("plotBigDataCostAndSim-AP", ra, "-S", s, "-K", d, ".pdf")
                pdf(fileName)
                print(costPlot)
                print(simPlot)
                dev.off()
            }    
        }
    }
}

# Plotting of max cost for randomized big datasets. All approaches, seeds and dimensions. Returns a plot for each approach. 
plotCostsFromOutputFileBig <- function() {
    approaches <- c(2, 3, 4)
    seeds <- c(7, 21, 42, 50)
    dimensions <- seq(5, 100, 5)
    
    for (ap in approaches) {
        df_plot <- data.frame(dimension = integer(0), cost = integer(0), seed = integer(0))
        for (s in seeds) {
            for (dim in dimensions) {
                path <- paste0("./randomBig/ap", ap, "/seed", s, "/")
                fileName <- paste0("AP", ap, "outputRandomBigData-S", s, "-K", dim, ".txt")
                file <- paste0(path, fileName)
                print(file)
                df <- read.table(file, sep = "", skip = 1)
                df_plot <- rbind(df_plot, data.frame(dimension = dim, cost = max(df[,2]), seed = s))
            }
        }
        title <- paste0("Maximum cost for the randomized big dataset | Approach ", ap)
        df_plot[,3] <- as.factor(df_plot[,3])
        g <- ggplot(data = df_plot, aes(x = dimension, y = cost, colour = seed)) + 
            geom_line(linetype = "solid") +
            scale_x_continuous(breaks = dimensions) + 
            scale_y_continuous(breaks = seq(1000000, 4000000, 100000)) +
            ggtitle("Big dataset")
        
        plotFileName <- paste0("plotCostsRandomBigData-Ap", ap, ".pdf")
        pdf(plotFileName)
        print(g)
        dev.off()
    }
}

# Plotting of max cost for randomized big datasets. All approaches, seeds and dimensions. Returns a plot for each approach. 
plotCostsFromOutputFileSmall <- function() {
    approaches <- c(2)
    seeds <- c(42)
    dimensions <- seq(5, 100, 5)
    
    for (ap in approaches) {
        df_plot <- data.frame(dimension = integer(0), cost = integer(0), seed = integer(0))
        for (s in seeds) {
            for (dim in dimensions) {
                path <- paste0("./randomSmallAp2/")
                fileName <- paste0("outputRandomSmallData-S", s, "-K", dim, ".txt")
                file <- paste0(path, fileName)
                print(file)
                df <- read.table(file, sep = "", skip = 1)
                df_plot <- rbind(df_plot, data.frame(dimension = dim, cost = max(df[,19]), seed = s))
            }
        }
        title <- paste0("Maximum cost for the randomized small dataset | Approach ", ap)
        df_plot[,3] <- as.factor(df_plot[,3])
        g <- ggplot(data = df_plot, aes(x = dimension, y = cost)) + 
            geom_line(linetype = "solid") +
            scale_x_continuous(breaks = dimensions) + 
            scale_y_continuous(breaks = seq(200000, 500000, 20000)) +
            ggtitle("Small dataset")
        
        plotFileName <- paste0("plotCostsRandomSmallData-Ap", ap, ".pdf")
        pdf(plotFileName)
        print(g)
        dev.off()
    }
}

# Used to make the cost plots for the actual datasets. Makes a pdf for each dataset.
plotCostsFromMeasureFile <- function() {
    names <- c("Threshold", "F-measure", "Clusters", "Common", "Cost")
    table <- read.table(paste0("./measure/measureSmallDataOriginalFixedColumns.txt"), sep = "", fill = TRUE, skip = 1)
    names(table) <- names
    g <- ggplot(data = table, aes(x = Threshold, y = Cost)) + 
        geom_line(linetype = "solid") +
        scale_x_continuous(breaks = seq(0, 330, 20)) +
        scale_y_continuous(breaks = seq(5000, 100000, 5000)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        ggtitle("Small dataset")
    fileName <- paste0("plotCostSmallDataOriginal.pdf")
    pdf(fileName)
    print(g)
    dev.off()
    
    table <- read.table(paste0("./measure/measureBigDataOriginalFixedColumns.txt"), sep = "", fill = TRUE, skip = 1)
    names(table) <- names
    g <- ggplot(data = table, aes(x = Threshold, y = Cost)) + 
        geom_line(linetype = "solid") +
        scale_x_continuous(breaks = seq(0, 300, 10)) +
        scale_y_continuous(breaks = seq(50000, 1000000, 50000)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        ggtitle("Big dataset")
    fileName <- paste0("plotCostBigDataOriginal.pdf")
    pdf(fileName)
    print(g)
    dev.off()
}


# Method for plotting gap cost difference returned by hclustDivisiveGap
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