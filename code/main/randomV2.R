setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(cluster)
library(fpc)
library(cowplot)
library(MASS)
library(stringi)

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

plot_clustering <- function(X, clusters, title, bb){
    plot_data = data.frame(
        x= X[,1],
        y= X[,2],
        result = clusters
    )
    
    g <- ggplot() + ggtitle(title) + 
        theme_classic() + 
        #coord_fixed(xlim = c(-2000, 1500), ylim = c(-500,1500)) +
        #coord_fixed(xlim = c(-190, 150), ylim = c(-190,150)) +
        coord_fixed(xlim = c(-320, 700), ylim = c(-320,700)) +
        xlab("x") + ylab("y") + 
        geom_point(data = plot_data, aes(x = x, y = y, color = factor(result))) + guides(color=FALSE) + 
        geom_path(data=as.data.frame(bb), aes(x=V1, y=V2), color="blue")
    return(g)
}

calculateIndexInMatrix <- function(IndexInOrder, amountOfProteins) {
    index <- IndexInOrder
    row <- index %% amountOfProteins
    col <- ceiling(index / amountOfProteins)
    if (row == 0) { row <- amountOfProteins }
    if (col == 0) { col <- amountOfProteins }
    
    return(c(row, col))
}


buildRandomSimMatrixAp2 <- function(proteins, simMatrix, k = 10, seed = 42) {
    set.seed(seed)
    simMatrixAsDist <- as.dist((max(simMatrix) + 1) - simMatrix) # Why was +1 needed ?
    
    mds <- isoMDS(simMatrixAsDist, k = k)
    mdsPoints <- mds$points
    Q = cov(mdsPoints)
    
    eval = eigen(Q)$values
    evec = eigen(Q)$vectors
    
    mdsPoints.projected = t(t(evec) %*% t(mdsPoints))
    
    R <- matrix( nrow = nrow(mdsPoints.projected), ncol = 0)
    for (i in 1:ncol(mdsPoints.projected)) {
        column <- runif(nrow(mdsPoints), min = min(mdsPoints.projected[,i], max = max(mdsPoints.projected[,i])))
        R <- cbind(R, column)
    }
    
    distRandom <- dist(R)
    simMatrix <- as.matrix((max(distRandom) + 1) - distRandom)
    row.names(simMatrix) = proteins
    colnames(simMatrix) = proteins
    
    return(simMatrix)
}

buildRandomSimMatrixAp3 <- function(proteins, simMatrix, k = 10, seed = 42) {
    set.seed(seed)
    amountOfZeroes <- table(simMatrix)[1]
    amountOfProteins <- length(proteins)
    
    simMatrixAsDist <- as.dist((max(simMatrix) + 1) - simMatrix) # Why was +1 needed ?
    
    mds <- isoMDS(simMatrixAsDist, k = k)
    mdsPoints <- mds$points
    Q = cov(mdsPoints)
    
    eval = eigen(Q)$values
    evec = eigen(Q)$vectors
    
    mdsPoints.projected = t(t(evec) %*% t(mdsPoints))
    
    R <- matrix( nrow = nrow(mdsPoints.projected), ncol = 0)
    for (i in 1:ncol(mdsPoints.projected)) {
        column <- runif(nrow(mdsPoints), min = min(mdsPoints.projected[,i], max = max(mdsPoints.projected[,i])))
        R <- cbind(R, column)
    }
    
    distRandom <- dist(R)
    simMatrix <- as.matrix((max(distRandom) + 1) - distRandom)
    
    # Find limit for setting value to 0
    tableSimMatrix <- table(simMatrix)
    tableHeaderValues <- sort(unique(as.vector(simMatrix)))
    sum <- 0
    limit <- 0
    for (i in 1:length(tableSimMatrix)) {
        sum <- sum + tableSimMatrix[i]
        if (sum >= amountOfZeroes) { # Found limit
            limit <- tableHeaderValues[i]
            break
        }
    }
    
    # Set values which are <= to limit to 0, until same amount of zeroes as original simMatrix
    countZeroes <- 0
    for (i in 1:amountOfProteins) {
        for (j in 1:amountOfProteins) {
            if (simMatrix[i,j] <= limit) {
                simMatrix[i,j] <- 0
                countZeroes <- countZeroes + 1
            }
            if (countZeroes >= amountOfZeroes) {
                break
            }
        }
        if (countZeroes >= amountOfZeroes) {
            break
        }
    }
    
    row.names(simMatrix) = proteins
    colnames(simMatrix) = proteins
    
    return(simMatrix)
}

buildRandomSimMatrixAp4 <- function(proteins, simMatrix, k = 10, seed = 42) {
    set.seed(seed)
    amountOfProteins <- length(proteins)
    
    simMatrixAsDist <- as.dist((max(simMatrix) + 1) - simMatrix) # Why was +1 needed ?
    
    mds <- isoMDS(simMatrixAsDist, k = k)
    mdsPoints <- mds$points
    Q = cov(mdsPoints)
    
    eval = eigen(Q)$values
    evec = eigen(Q)$vectors
    
    mdsPoints.projected = t(t(evec) %*% t(mdsPoints))
    
    R <- matrix( nrow = nrow(mdsPoints.projected), ncol = 0)
    for (i in 1:ncol(mdsPoints.projected)) {
        column <- runif(nrow(mdsPoints), min = min(mdsPoints.projected[,i], max = max(mdsPoints.projected[,i])))
        R <- cbind(R, column)
    }
    
    distRandom <- dist(R)
    simMatrixRandom <- as.matrix((max(distRandom) + 1) - distRandom)

    # Find ordering of both matrices
    orderOriginal <- order(simMatrix)
    orderRandom <- order(simMatrixRandom)

    # Set ordered values into same ordering of random
    for (i in 1:length(orderOriginal)) {
        posInOriginal <- calculateIndexInMatrix(orderOriginal[i], amountOfProteins)
        posInRandom <- calculateIndexInMatrix(orderRandom[i], amountOfProteins)

        simMatrixRandom[posInRandom[1], posInRandom[2]] <- simMatrix[posInOriginal[1], posInOriginal[2]]
    }

    row.names(simMatrixRandom) = proteins
    colnames(simMatrixRandom) = proteins
    
    return(simMatrixRandom)
}

#simMatrix <- buildRandomSimMatrix(5)
