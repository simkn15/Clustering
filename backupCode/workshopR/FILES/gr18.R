# x er en vektor af diskrete værdier, i dette eksempel genereret fra Poissonfordelingen
x <- rpois(100, 5)
# Vi kunne godt bruge hist-funktionen til at fremstille et histogram af diskrete data
# Funktionen er imidlertid mere rettet mod kontinuerte data, så vi bruger i stedet plot-funktionen med type = "h"
plot(table(x), type = "h", lwd = 5, ylab = "Frekvens")
# Bemærk at lwd = 5 gør søjlerne lidt bredere

# Hvis vi foretrækker et histogram med relative hyppigheder, så kan hvert table-element skaleres med længden af x
plot(table(x) / length(x), type = "h", lwd = 5, ylab = "Frekvens")