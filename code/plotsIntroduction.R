#########################
# Code used for making plots for the Introduction section of the report
#########################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(transclustr)
library(ggplot2)

plot_clustering <- function(X, clusters, title){
    plot_data = data.frame(
        x = X[,1],
        y = X[,2],
        class = factor(clusters)
    )
    
    g <- ggplot() + ggtitle(title) + 
        theme(legend.position = "none") + 
        xlab("x") + ylab("y") + 
        geom_point(data = plot_data, aes(x = x, y = y, color = class))
    return(g)
}

# Variables
SD <- 1
CLUSTER_SIZE <- 50
THRESHOLD = 15

################################################################################
# Read test data
################################################################################
table <- read.table("./data/g2/g2-2-10.txt", sep = "", skip = 5)
df_g2_2_10 <- as.data.frame(table)
names(df_g2_2_10) <- c("x","y")

table <- read.table("./data/g2/g2-2-30.txt", sep = "", skip = 5)
df_g2_2_30 <- as.data.frame(table)
names(df_g2_2_30) <- c("x","y")

table <- read.table("./data/g2/g2-2-50.txt", sep = "", skip = 5)
df_g2_2_50 <- as.data.frame(table)
names(df_g2_2_50) <- c("x","y")

table <- read.table("./data/g2/g2-2-70.txt", sep = "", skip = 5)
df_g2_2_70 <- as.data.frame(table)
names(df_g2_2_70) <- c("x","y")

################################################################################
# Draw the test data
################################################################################
plot_g2_2_10 <- ggplot(data=df_g2_2_10,aes(x=x,y=y)) +
    geom_point() +
    coord_fixed()

plot_g2_2_30 <- ggplot(data=df_g2_2_30,aes(x=x,y=y)) +
    geom_point() +
    coord_fixed()

plot_g2_2_50 <- ggplot(data=df_g2_2_50,aes(x=x,y=y)) +
    geom_point() +
    coord_fixed()

plot_g2_2_70 <- ggplot(data=df_g2_2_70,aes(x=x,y=y)) +
    geom_point() +
    coord_fixed()

pdf("g2-2-all.pdf")
print(plot_g2_2_10)
print(plot_g2_2_30)
print(plot_g2_2_50)
print(plot_g2_2_70)
dev.off()

pdf("g2-2-10.pdf")
print(plot_g2_2_10)
dev.off()

pdf("g2-2-30.pdf")
print(plot_g2_2_30)
dev.off()

pdf("g2-2-50.pdf")
print(plot_g2_2_50)
dev.off()

pdf("g2-2-70.pdf")
print(plot_g2_2_70)
dev.off()

combined <- plot_grid(plot_g2_2_10, plot_g2_2_30, plot_g2_2_50, plot_g2_2_70, labels = "AUTO", ncol=2)
combined
pdf("g2-2-combined.pdf")
print(combined)
dev.off()

################################################################################
# Create a dist object from the test data
################################################################################
dist_set <- dist(df_g2_2_70)

################################################################################
# Cluster using tclust
################################################################################
tclust_res <- tclust(dist_set, threshold = 525)

################################################################################
# Draw the result
################################################################################
set_res <- as.data.frame(cbind(df_g2_2_70,tclust_res$clusters[[1]]))
names(set_res) <- c("x","y","cluster")
find_hull <- function(df) df[chull(df$x, df$y), ]
hulls <- ddply(set_res, "cluster", find_hull)
print(ggplot(data=as.data.frame(set_res),aes(x=x,y=y,colour=factor(cluster),fill=factor(cluster))) +
          geom_point() +
          theme(legend.position = "none") +
          coord_fixed() +
          geom_polygon(data = hulls,alpha = 0.25) +
          ggtitle(paste("Threshold:",THRESHOLD," Cost:",tclust_res$costs[[1]]))
)

k <- 2
kmeans.res = kmeans(df_g2_2_70, algorithm = "Lloyd", centers = 10, nstart=1)
# title <- paste0("k-means: Dataset = ", dataset.name, ", k = ", clusters)
plt <- plot_clustering(df_g2_2_70, kmeans.res$cluster, title = "")
plt

pdf("g2-2-70-k10.pdf")
print(plt)
dev.off()
