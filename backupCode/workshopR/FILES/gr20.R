# Beskrivelse af fremgangsmåde:
# Vi bruger ppoints-funktionen til at generere en række af punkter mellem 0 og 1
# Vi transformerer disse punkter til kvartiler ved brug af quantile-funktionen for den antagede fordeling
# Vi sorterer vores stikprøvedata.
# Vi plotter de sorterede data mod de beregnede kvartiler
# Vi bruger abline til at plotte den diagonale linje

# Rateparameter
lambda <- 1/10
# Vi genererer en tilfældig stikprøve fra eksponentialfordelingen med middelværdi 10 (eller ækvivalent rate 1/10)
y <- rexp(1000, rate = lambda)
# Vi plotter de sorterede data mod de beregnede kvartiler som beskrevet i fremgangsmåden ovenfor
plot(qexp(ppoints(y), rate = lambda), sort(y), main = "QQ-plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
# Kvartilfunktionen for eksponentialfordelingen er qexp, som tager et rateargument

# Vi tilføjer den diagonale linje
abline(a = 0, b = 1)