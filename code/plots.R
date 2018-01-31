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