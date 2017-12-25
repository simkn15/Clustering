library(transclustr)
library(ggplot2)
library(plyr)

# Variables
SD <- 1
CLUSTER_SIZE <- 50
THRESHOLD = 15

################################################################################
# Genereate some test data
################################################################################
a <- cbind(rnorm(CLUSTER_SIZE, mean = 10, sd = SD), rnorm(n  =CLUSTER_SIZE,  mean = 10,sd = SD))
b <- cbind(rnorm(CLUSTER_SIZE, mean = 0,  sd = SD),  rnorm(n =CLUSTER_SIZE, mean = 10,sd = SD))
c <- cbind(rnorm(CLUSTER_SIZE, mean = 10, sd = SD), rnorm(n  =CLUSTER_SIZE,  mean = 0, sd = SD))
d <- cbind(rnorm(CLUSTER_SIZE, mean = 0,  sd = SD),  rnorm(n =CLUSTER_SIZE, mean = 0, sd = SD))
set <- rbind(a,b,c,d)
df_set <- as.data.frame(set)
names(df_set) <- c("x","y")

################################################################################
# Draw the test data
################################################################################
ggplot(data=df_set,aes(x=x,y=y)) +
   geom_point() +
   coord_fixed() +
   ggtitle("Test Data")

################################################################################
# Create a dist object from the test data
################################################################################
dist_set <- dist(set)

################################################################################
# Cluster using tclust
################################################################################
start.time <- Sys.time()
tclust_res<- tclust(dist_set, threshold = THRESHOLD)
end.time <- Sys.time()
time.taken <- end.time - start.time

################################################################################
# Draw the result
################################################################################
set_res <- as.data.frame(cbind(set,tclust_res$clusters[[1]]))
names(set_res) <- c("x","y","cluster")
find_hull <- function(df) df[chull(df$x, df$y), ]
hulls <- ddply(set_res, "cluster", find_hull)
print(ggplot(data=as.data.frame(set_res),aes(x=x,y=y,colour=factor(cluster),fill=factor(cluster))) +
         geom_point() +
         theme(legend.position = "none") +
         coord_fixed() +
         geom_polygon(data = hulls,alpha = 0.25) +
         ggtitle(paste("Threshold:",THRESHOLD," Cost:",tclust_res$costs[[1]],"Time:",time.taken,"s"))
)
