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

# Why do we need to project ?
buildRandomSimMatrix <- function(proteins, simMatrix, k = 10, seed = 42) {
    set.seed(42)
    # table <- read.table("../data/brown/sfld_brown_et_al_amidohydrolases_protein_similarities_for_beh.txt", sep = "", skip = 5)
    # #table <- read.table("simBig.txt", sep = "", skip = 5)
    # df_table <- as.data.frame(table)
    # proteins = levels(df_table[,1])
    
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

#simMatrix <- buildRandomSimMatrix(5)
