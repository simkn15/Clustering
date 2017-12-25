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

library(ggplot2)
library(cluster)
library(fpc)
library(cowplot)
library(MASS) #library(stringi)
set.seed(2)

source("utilities.R")
####################################################################################################
# Read in data and build similarity matrix
####################################################################################################
# gold_df_before_split <- as.data.frame(readLines("../data/brown/sfld_brown_et_al_amidohydrolases_families_gold_standard.txt"))
# gold <- as.data.frame(stri_split_fixed(gold_df_before_split[,1], "\t", simplify = TRUE))

table <- read.table("../data/brown/sfld_brown_et_al_amidohydrolases_protein_similarities_for_beh.txt", sep = "", skip = 5)
# table <- read.table("simBig.txt", sep = "", skip = 5)
df_table <- as.data.frame(table)

# Get the involved proteins
proteins = levels(df_table[,1])

# Build similarity matrix
simMatrix <- buildSimilarityMatrix(proteins, df_table)
# simMatrix <- buildSimilarityMatrixFromBlast(proteins, df_table)

# multidimensional scaling R -> Google
# mydata <- (max(as.dist(simMatrix)) + 1) - as.dist(simMatrix)
d <- as.dist((max(simMatrix)+.001) - simMatrix)

head(table(d))

#cmd <- cmdscale(d, eig=TRUE, k=2)
cmd <- isoMDS(d, k = 2)

head(cmd$points)

### This should make a similarity matrix
d_new = dist(cmd$points)

mm <- as.matrix((max(d_new) + 1) - d_new)
###
# d_diff = d - d_new
# head(d_diff)
# 
# max(abs(as.vector(d_new) - as.vector(d)))

## Plotting with points with cluster colors according to gold standard
plt_data = data.frame(x=cmd$points[,1], y = cmd$points[,2], class = gold_sorted[,2]) # Load gold data from measureBrown.R

ggplot(data = plt_data) + geom_point(aes(x=x, y=y, color = class)) +coord_fixed(xlim = c(-330, 330), ylim = c(0,200)) 
##

#Create a cluster with center x,y with num point and varianze v.x, v.y
# create.cluster <- function(x, y, num, var.x, var.y){
#     ex <- rnorm(num, 0, var.x)
#     ey <- rnorm(num, 0, var.y)
#     return (cbind(ex+x, ey+y))
# }
# 
# v = 0.45
# c1 <- create.cluster(4, 5, 50, v,v)
# c2 <- create.cluster(6, 7, 80, v,v)
# c3 <- create.cluster(10, 9, 50, v,v)
# c4 <- create.cluster(9, 7, 80, v,v)

#X <- scale(rbind(c1, c2, c3, c4), scale = FALSE, center =TRUE)
X <- cmd$points

plot_clustering <- function(X, clusters, title, bb){
    plot_data = data.frame(
        x= X[,1],
        y= X[,2],
        result = clusters
    )
    
    g <- ggplot() + ggtitle(title) + 
        theme_classic() + 
        #coord_fixed(xlim = c(-2000, 1500), ylim = c(-500,1500)) +
        coord_fixed(xlim = c(-330, 330), ylim = c(0,50)) +
        xlab("x") + ylab("y") + 
        geom_point(data = plot_data, aes(x = x, y = y, color = factor(result))) + guides(color=FALSE) + 
        geom_path(data=as.data.frame(bb), aes(x=V1, y=V2), color="blue")
    return(g)
}

kclust = kmeans(X, centers = 4, nstart=10)

Q = cov(X)

eval = eigen(Q)$values
evec = eigen(Q)$vectors

X.projected = t(t(evec) %*% t(X))

bb <- rbind(c(min(X.projected[,1]),min(X.projected[,2])),
            c(max(X.projected[,1]),min(X.projected[,2])),
            c(max(X.projected[,1]),max(X.projected[,2])),
            c(min(X.projected[,1]),max(X.projected[,2])),
            c(min(X.projected[,1]),min(X.projected[,2])))

bb.projected <- t(solve(t(evec)) %*% t(bb))

plot_data = plot_clustering(X, kclust$cluster, "Original Data k=21", bb.projected)
plot_data

r.x <- runif(nrow(X), min = min(X.projected[,1]), max = max(X.projected[,1]))
r.y <- runif(nrow(X), min = min(X.projected[,2]), max = max(X.projected[,2]))
R <- cbind(r.x, r.y)

R.projected <- t(solve(t(evec)) %*% t(R))

kclust = kmeans(R.projected, centers = 4, nstart=10)
plt_random = plot_clustering(R.projected, kclust$cluster, "Random Data k=21", bb.projected)
plt_random

null = c()
for(i in c(1:25)){
    r.x <- runif(nrow(X), min = min(X.projected[,1]), max = max(X.projected[,1]))
    r.y <- runif(nrow(X), min = min(X.projected[,2]), max = max(X.projected[,2]))
    R <- cbind(r.x, r.y)
    R.projected <- t(solve(t(evec)) %*% t(R))
    null = rbind(null,sapply(c(2:21), function(k) kmeans(R.projected, k, nstart=100, iter.max = 25)$tot.withinss))
}


plt_data <- data.frame(
    x = c(2:21),
    orig = sapply(c(2:21), function(k) kmeans(X, k, nstart=100)$tot.withinss),
    null = colMeans(null),
    sd = apply(null,2,sd)
)

g <- ggplot(data = plt_data, aes(x = x)) + ggtitle("Within Sum of Squares") + 
    theme_classic() + 
    theme(legend.position="bottom") +
    xlab("Number of Clusters k") + ylab("Sum of Squares") + 
    geom_errorbar(aes(y= null, ymin=null-sd, ymax=null+sd, colour="null"), width=.2) +
    geom_line(aes(y = orig, colour="orig")) + 
    geom_line(aes(y = null, colour="null")) +
    scale_color_discrete(name="", labels = c("Randomized Data", "Actual Data"))
plt_diff <- g


plt_gap_data <- data.frame(
    x = c(2:21),
    gap = log(plt_data$null) - log(plt_data$orig),
    sd = log(plt_data$null+plt_data$sd) - log(plt_data$orig) - log(plt_data$null) + log(plt_data$orig)
)

g <- ggplot(data = plt_gap_data, aes(x = x, y=gap)) + ggtitle("GAP Statistic") + 
    theme_classic() + 
    xlab("Number of Clusters k") + ylab("GAP Statistic") + 
    geom_errorbar(aes(ymin=gap-sd, ymax=gap+sd, colour="Gap"), width=.2) +
    geom_line(aes(colour="Gap")) + guides(color=FALSE)
#scale_color_discrete(name="", labels = c("Randomized Data", "Actual Data"))
plt_gap <- g
plt_gap


combined <- plot_grid(plot_data, plt_random, plt_diff, plt_gap, labels = "AUTO", ncol=2)
combined

# factor = 3.7
# pdf("gap.pdf", width=2*factor, height=2*factor)
# print(combined)
# dev.off()
