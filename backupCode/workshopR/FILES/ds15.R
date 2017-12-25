# Antag nu, at vi ?nsker at udf?re matrixoperationer som transponering, matrixinversion eller matrixmultiplikation
# M?ske vi ogs? gerne vil konstruere en identitetsmatrix?

# Husk at en m ? n-matrix og en n ? p-matrix kan multipliceres til en m ? p-matrix

# Vi initialiserer en 3x2-matrix A og en 3x2-matrix B
A <- matrix(c(9, -8, -2, 7, 4, 4), 3, 2)
A
B <- matrix(c(0, -11, 5, -7, 3, 6), 3, 2)
B
# Vi transponerer matricen B ved at skrive t(B) - bem?rk at dimensionerne af den transponerede B-matrix er 2x3
t(B)
# Vi udf?rer matrixmultiplikation ved at skrive %*%, s? her multipliceres en 3x2-matrix og en 2x3-matrix til en 3x3-matrix
A %*% t(B)
# Vi initialiserer en 3x3-matrix X
X <- matrix(c(25, -2, 4, -2, 4, 1, 4, 1, 9), 3, 3)
X
# Funktionen diag, der returnerer diagonalelementerne, omsluttes af funktionen sum for at l?gge diagonalelementerne sammen
sum(diag(X))
# Summen af diagonalelementerne kaldes ogs? matricens spor (engelsk: trace)
# Vi kan ogs? beregne determinanten
det(X)
# Vi kan ogs? beregne den inverse matrix
XI <- solve(X)
XI
# Vi kan kontrollere resultatet
XI %*% X
# Dette giver (numerisk) en identitetsmatrix
round(XI %*% X)

# En n ? n-identitetsmatrix kan konstrueres p? f?lgende m?de ved at skrive diag(n)
diag(3)
