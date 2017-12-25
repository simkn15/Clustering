library(faraway)

# Antag nu, at vi plotter par af datapunkter og gerne vil tilføje en linje, der illustrerer deres lineære regression

# Vi modellerer strongx-datasættet fra faraway-pakken
data(strongx)
# Lav et modelobjekt
model <- lm(crossx ~ energy, data = strongx)
# Plot (x, y)-parrene
plot(crossx ~ energy, data = strongx, xlab = "Energi", ylab = "Tværsnit")
# Alternativ syntaks
plot(strongx$energy, strongx$crossx, xlab = "Energi", ylab = "Tværsnit")
# Brug abline-funktionen til at indtegne den tilpassede regressionslinje
abline(model)