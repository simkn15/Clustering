# Antag nu, at vi ønsker at udvælge en enkelt række eller søjle fra en matrix

# Vi initialiserer nu en matrix med 3 rækker og 4 søjler
mat <- matrix(1:12, 3, 4)
mat
# Normalt, når man udvælger en række eller søjle fra en matrix, så resulterer det i en dimensionsløs vektor
# Her benyttes normal indeksering til at udvælge matricens første række
mat[1, ]
# Her benyttes normal indeksering til at udvælge matricens tredje søjle
mat[, 3]
# Når man tilføjer argumentet drop = FALSE bevarer R derimod dimensionerne - så hvis en række udvælges fås en 1 × n-matrix
mat[1, , drop = FALSE]
# Tilsvarende, når man vælger en søjle med drop = FALSE, fås en n × 1-matrix
mat[, 3, drop = FALSE]
