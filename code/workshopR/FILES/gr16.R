library(MASS)

# Datasættet Cars93 indeholder 27 variable, som beskriver 93 bilmodeller fra og med 1993:
# En numerisk variabel er MPG.city (miles per gallon i byen)

# Et histogram fremstilles
hist(Cars93$MPG.city, main = "Miles per gallon i byen (1993)", xlab = "Miles per gallon", ylab = "Frekvens")
# Funktionen hist afgør antallet af søjler og standardalgoritmen vælger 7 søjler for dette eksempel
# Desværre er det svært at genkende fordelingens form med kun 7 søjler, så jeg foreslår 20 søjler:
hist(Cars93$MPG.city, 20, main = "Miles per gallon i byen (1993)", xlab = "Miles per gallon", ylab = "Frekvens")