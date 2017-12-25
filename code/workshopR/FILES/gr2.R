# Ved kald af plot:
# - Brug main-argumentet til en titel
# - Brug xlab-argumentet til en etiket p? x-aksen
# - Brug ylab-argumentet til en etiket p? y-aksen

# Tilf?j en titel og etiketter, s? plottet bliver mere interessant og lettere at fortolke
plot(cars,
     main = "Biler: Hastighed versus bremselængde (1920)",
     xlab = "Hastighed (miles i timen)",
     ylab = "Bremselængde (fod)")

# Alternativt kan en titel og etiketter f?rst udelades og efterf?lgende tilf?jes med et kald til title-funktionen
plot(cars, ann = FALSE)
title(main = "Biler: Hastighed versus bremsel?ngde (1920)", xlab = "Hastighed (miles i timen)", ylab = "Bremsel?ngde (fod)")
# Titlen fort?ller os noget om plottet og etiketterne indeholder nu m?leenheder

# Ovenst?ende er klare forbedringer i forhold til plottet, som vi fremstillede i forrige fil