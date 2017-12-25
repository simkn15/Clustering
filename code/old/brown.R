library(transclustr)
library(ggplot2)
library(plyr)

# Variables
SD <- 1
CLUSTER_SIZE <- 50
THRESHOLD = 50

################################################################################
# Genereate some test data
################################################################################
table <- read.table("sfld_brown_et_al_amidohydrolases_protein_similarities_for_beh.txt", sep = "", skip = 5)
df_table <- as.data.frame(table)

# Get the involved proteins
proteins = levels(df_table[,1])

# Create an empty matrix
sim_matrix = matrix(0,nrow = length(proteins), ncol = length(proteins))

# Set row and column names for readability
row.names(sim_matrix) = proteins
colnames(sim_matrix) = proteins

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
      sim_matrix[x,y] = sim_matrix[y,x] = 0
    }else{
      #Take the minimum of both directions
      sim_matrix[x,y] = sim_matrix[y,x] = min(presec[presec[,1] == p2, 3], presec[presec[,2] == p2, 3])
    }
  }
}
# Set the diagonal to the highest observed value (not actually necessary, since self-similarity is never queried)
diag(sim_matrix) = max(sim_matrix)

################################################################################
# Cluster using tclust
################################################################################
start.time <- Sys.time()
tclust_res<- tclust(simmatrix = sim_matrix, convert_dissimilarity_to_similarity = FALSE, threshold = THRESHOLD)
end.time <- Sys.time()
time.taken <- end.time - start.time

################################################################################
# Draw the result by
################################################################################
#Perform a multidimensional scaling to transfer objects in a two-dimensional Euclidean space
# Not important at the moment
