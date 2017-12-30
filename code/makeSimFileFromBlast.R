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

table <- read.table("../data/big/blast_original.txt", sep = "")
df_table <- as.data.frame(table)
df_table <- df_table[,c(1,2,11)]
table <- df_table
for (i in 1:nrow(df_table)) {
    if (table[i,3] != 0) {
        table[i,3] <- -log(table[i,3])
    }
}

#proteins <- levels(table[,1])

#simMatrix <- buildSimilarityMatrix(proteins, table)
require(gdata)
file <- file("./data/big/simBig.txt", "w")
# for (i in 1:nrow(table)) {
#     string <- paste0(table[i,])
#     
# }
write.fwf(table, file, rownames = FALSE, colnames = FALSE, width = 15)
close(file)

# con <- file("simBig.txt") 
# length(readLines(con))
# 
# con <- file("../data/big/blast_original.txt") 
# length(readLines(con))



