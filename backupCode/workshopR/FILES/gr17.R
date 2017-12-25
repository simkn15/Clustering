# Antag nu, at vi ønsker et glattere tæthedsestimat end vi får ud af histogramkonstruktionen

# Vi kan bedre visualisere en underliggende fordeling, hvis vi tilføjer en approksimeret tæthedsfunktion

# Vi genererer en stikprøve fra gammafordelingen
samp <- rgamma(500, 2, 2)
# Vi plotter histogrammet med sandsynlighedstætheder
hist(samp, 20, freq = FALSE, ylab = "Tæthed", main = "Histogram med tæthedsestimat")
# Vi plotter den estimerede tæthedsfunktion
lines(density(samp))

# Funktionen density approksimerer formen af tæthedsfunktionen ikke-parametrisk
# Hvis vi kender den faktiske underliggende fordeling, så kan vi med fordel plotte tæthedsfunktionen i stedet