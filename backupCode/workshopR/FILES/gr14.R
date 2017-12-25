# Funktionen abline kan bruges til at tegne en vertikal eller horisontal linje:
# abline(v = x) tegner en vertikal linje ved x
# abline(h = y) tegner en horisontal linje ved y

# Antag nu, at vi har en vektor med en stikprøve bestående af 100 uafhængige tal fra den uniforme fordeling på (0, 1)
samp <- runif(100, 0, 1)
# Vi plotter stikprøvens punkter
plot(samp, ylim = c(-0.2, 1.2))
# Vi beregner stikprøvens middelværdi
m <- mean(samp)
# Vi indtegner en linje gennem middelværdien
abline(h = m)
# Vi vil gerne illustrere stikprøvens standardafvigelse grafisk - så vi beregner den og indtegner prikkede
# linjer ved plus-minus 1 standardafvigelse og plus-minus 2 standardafvigelser fra stikprøvens middelværdi
stdevs <- m + c(-2, -1, +1, +2) * sd(samp)
abline(h = stdevs, lty = "dotted")