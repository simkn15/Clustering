library(MASS)

# Såkaldte "conditioning plots" er en anden måde at undersøge og illustrere effekten af en faktor på data

# Datasættet Cars93 indeholder 27 variable, som beskriver 93 bilmodeller fra og med 1993:
# To numeriske variable er MPG.city (miles per gallon i byen) og Horsepower (motorens hestekræfter)
# En kategorisk variabel er Origin (oprindelse), som kan være USA eller non-USA, alt efter hvor modellen blev bygget

# Er der forskellige sammenhænge for USA-modeller og ikke-USA-modeller?
coplot(Horsepower ~ MPG.city | Origin, data = Cars93, xlab = "Miles per gallon i byen", ylab = "Hestekræfter")

# Hvis virkelig vi behøver et monster med 300 hestekræfter, så er vi nødt til at købe en bil, der er bygget i USA
# Hvis vi ønsker højt antal miles per gallon i byen har vi flere valgmuligheder blandt ikke-USA-modeller