library(ggplot2)
library(cluster)
library(fpc)
library(cowplot)
set.seed(2)

#Create a cluster with center x,y with num point and varianze v.x, v.y
create.cluster <- function(x, y, num, var.x, var.y){
    ex <- rnorm(num, 0, var.x)
    ey <- rnorm(num, 0, var.y)
    return (cbind(ex+x, ey+y))
}

v = 0.45
c1 <- create.cluster(4, 5, 50, v,v)
c2 <- create.cluster(6, 7, 80, v,v)
c3 <- create.cluster(10, 9, 50, v,v)
c4 <- create.cluster(9, 7, 80, v,v)

X <- scale(rbind(c1, c2, c3, c4), scale = FALSE, center =TRUE)

plot_clustering <- function(X, clusters, title, bb){
    plot_data = data.frame(
        x=X[,1],
        y= X[,2],
        result = clusters
    )
    
    g <- ggplot() + ggtitle(title) + 
        theme_classic() + 
        coord_fixed(xlim = c(-5, 5), ylim = c(-4.5,4.5)) +
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

plot_data = plot_clustering(X, kclust$cluster, "Original Data k=4", bb.projected)
plot_data

r.x <- runif(nrow(X), min = min(X.projected[,1]), max = max(X.projected[,1]))
r.y <- runif(nrow(X), min = min(X.projected[,2]), max = max(X.projected[,2]))
R <- cbind(r.x, r.y)
R.projected <- t(solve(t(evec)) %*% t(R))

kclust = kmeans(R.projected, centers = 4, nstart=10)
plt_random = plot_clustering(R.projected, kclust$cluster, "Random Data k=4", bb.projected)
plt_random

null = c()
for(i in c(1:25)){
    r.x <- runif(nrow(X), min = min(X.projected[,1]), max = max(X.projected[,1]))
    r.y <- runif(nrow(X), min = min(X.projected[,2]), max = max(X.projected[,2]))
    R <- cbind(r.x, r.y)
    R.projected <- t(solve(t(evec)) %*% t(R))
    null = rbind(null,sapply(c(2:16), function(k) kmeans(R.projected, k, nstart=100, iter.max = 25)$tot.withinss))
}


plt_data <- data.frame(
    x = c(2:16),
    orig = sapply(c(2:16), function(k) kmeans(X, k, nstart=100)$tot.withinss),
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
    x = c(2:16),
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

factor = 3.7
pdf("./gap.pdf", width=2*factor, height=2*factor)
print(combined)
dev.off()
