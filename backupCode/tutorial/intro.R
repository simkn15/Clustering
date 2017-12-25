myString <- "Hello, World!"
print(myString)

apple <- c('red','green',"yellow")
print(apple)

list <- list(c(2,5,3),21.3,sin)
print(list)

M = matrix( c('a','a','b','c','b','a'), nrow = 2, ncol = 3, byrow = TRUE)
print(M)

a <- array(c('green','yellow'), dim = c(3,3,2))
print(a)

apple_colors <- c('green','green','yellow','red','red','red','green')
factor_apple <- factor(apple_colors)
print(factor_apple)
print(nlevels(factor_apple))

BMI <- data.frame(
  gender = c("Male", "Male", "Female"),
  height = c(152, 171.5, 165),
  weight = c(81,93,78),
  Age = c(42, 38, 26)
)
print(BMI)

var1 = c(0,1,2,3)
var2 <- c("learn","R")
c(TRUE,1) -> var3
print(var1)
print(var2)
print(var3)
cat ("var1 is ", var1, "\n")
cat ("var2 is ", var2, "\n")
cat ("var3 is ", var3, "\n")

print(ls())

v <- c(2,4,6)
t <- c(8,3,4)
print(v^t)

v <- c(2,5.5,6,9)
t <- c(8,2.5,14,9)
print(v == t)

v <- c(3,0,TRUE,2+2i)
t <- c(4,0,FALSE,2+3i)
print(v|t) #v&t

v <- 2:8
print(v)
v1 <- 8
v2 <- 12
t <- 1:10
print(v1 %in% t)
print(v2 %in% t)
M <- matrix(c(2,6,5,1,10,4), nrow = 2, ncol = 3, byrow = TRUE)
print(M)
t <- M %*% t(M)
print(t)

new.function <- function(a) {
  for(i in 1:a) {
    b <- i^2
    print(b)
  }
}
printer <- function(value) {
  i <- 0
  while (i < value) {
    print(i)
    i <- i + 1
  }
}

printYear <- function() {
  for (year in c(2010,2011,2012,2013,2014,2015)){
    print(paste("The year is", year))
  }
}

print1 <- function(value) {
  for (i in 1:10 ) {
    print(i)
  }
}

printer(4)
print1()
printYear()

a <- "Hello"
b <- 'How'
c <- "are you? "
print(paste(a,b,c))
print(paste(a,b,c, sep = "-"))
print(paste(a, b, c, sep = "", collapse = ""))
print(nchar(a))
print(toupper(c))
res <- substring(a, 2, 3)
print(res)

print(seq(5, 9, by <- 0.4))

v <- c(3,8,4,5,0,11,-9,304)
print(sort(v, decreasing <- TRUE))

# Elements are arranged sequentially by row.
M <- matrix(c(3:14), nrow = 4, byrow = TRUE)
print(M)

# Elements are arranged sequentially by column.
N <- matrix(c(3:14), nrow = 4, byrow = FALSE)
print(N)

# Define the column and row names.
rownames = c("row1", "row2", "row3", "row4")
colnames = c("col1", "col2", "col3")

P <- matrix(c(3:14), nrow = 4, byrow = TRUE, dimnames = list(rownames, colnames))
print(P)

v <- c(1,2,3,4)
P[,1] <- P[,1]*v
print(P)
