# Antag nu, at vi har parrede observationer: (x_1, y_1), (x_2, y_2), ..., (x_n, y_n)
# Hvis data lagres i to parallelle vektorer, x og y, så kan disse bruges som argumenter i plot-funktionen:
# plot(x, y)
# Hvis data lagres i en dataramme med to søjler, dfrm, så kan datarammen plottes:
# plot(dfrm)

# Det indbyggede datasæt cars indeholder to søjler, speed (afbildes på x-aksen) og dist (afbildes på y-aksen)

# Fremstil et scatterplot
plot(cars)
# Funktionens formål er, at fremstille et plot af (x, y)-parrene i grafikvinduet, så den returnerer ikke noget

# For datarammer indeholdende mere end to søjler fremstilles flere scatterplot

# For at få et scatterplot skal data være numeriske (bemærk: plot er en polymorf funktion!)