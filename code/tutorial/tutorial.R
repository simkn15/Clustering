#Read data from file
file <- read.table("data.txt", header = FALSE, sep=",")
#file <- scan("/Users/slk/Dropbox/sdu/clustering ISA 5th semester/code/data.txt", sep <- ",")
print(file)
M = matrix(file, nrow <- 3, ncol <- 3, byrow <- TRUE)
print
#Write a function performing matrix multiplication
size = dim(M)
mult <- function(a, b) {
  res <- matrix(, nrow = 0, ncol = 0)
  if (is.matrix(a) && is.matrix(b) && dim(a)[2] == dim(b)[1]) {
    res <- matrix(1, nrow = dim(a)[1], ncol = dim(b)[2])
    for (i in 1:dim(a)[1]) {
      for (j in 1:dim(b)[2]) {
        res[i,j] <- sum(a[i,] * b[,j])
      }
    }
  }
  res
}
v <- c(1,2,3,1,2,3,1,2,3)
m1 <- matrix(v, nrow = 3, ncol = 3, byrow = TRUE)
m2 <- matrix(rev(v), nrow = 3, ncol = 3, byrow = TRUE)
print(mult(m1, m2))
#Write a function calculating the determinant of a square matrix. Again, check against the built-in capabilities.

#Create large random matrices and
#(a) compare your program's performance (in terms of runtime) versus the built-in capabilities
#     for growing matrix sizes
#(b) create a plot depicting the runtime vs. matrix-size
