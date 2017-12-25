# Antag nu, at vi ønsker at udvælge søjler fra en dataramme efter deres navn

# For at udvælge en enkelt søjle, brug et af disse udtryk:
# dfrm[["navn"]] - returnerer søjlen kaldet navn
# dfrm$navn - samme som forrige, bare en alternativ syntaks

# For at udvælge en eller flere søjler og pakke dem ind i en dataramme, brug disse udtryk:
# dfrm["navn"] - udvælger en søjle og pakker den ind i en dataramme
# dfrm[c("navn_1", "navn_2", ..., "navn_k")] - udvælger flere søjler og pakker dem ind i en dataramme
# Man kan også bruge matrixlignende indeksering til at udvælge en eller flere søjler:
# dfrm[, "navn"]
# dfrm[, c("navn_1", "navn_2", ..., "navn_k")]

# Samme udfordring som i forrige fil med returtyperne ved brug af matrixlignende indeksering!

# Prøv selv at udvælge søjler fra en dataramme efter deres navne!
