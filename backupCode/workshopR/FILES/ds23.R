library(MASS)

# Antag nu, at vi ønsker en lettere måde til at udvælge rækker og søjler fra en dataramme eller en matrix

# Indeksering kan være besværlig at bruge i visse tilfælde
# Funktionen subset er en mere bekvem og læsbar måde at udvælge rækker og søjler

# Find modelnavnet for biler, der kan overstige 30 miles per gallon (MPG) i byen
subset(Cars93, select = Model, subset = (MPG.city > 30))
# Find modelnavnet og prisintervallet for firecylindrede biler produceret i USA
subset(Cars93, select = c(Model, Min.Price, Max.Price), subset = (Cylinders == 4 & Origin == "USA"))
# Find producentens navn og modellens navn for alle biler, hvor miles per gallon (MPG) på motorvejen er større end medianen
subset(Cars93, select = c(Manufacturer, Model), subset = c(MPG.highway > median(MPG.highway)))
