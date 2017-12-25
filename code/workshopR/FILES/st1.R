# Antag nu, at vi har en stikprøve fra en population
# Givet denne stikprøve ønsker vi at bestemme konfidensintervallet for middelværdien

# Lav en stikprøve
x <- rnorm(50, mean = 100, sd = 15)
# Brug t.test-funktionen på din vektor (giver en masse output, heriblandt et konfidensinterval)
t.test(x)
# I dette eksempel er konfidensintervallet approksimativt 95,42 < mu < 103,92, hvilket nogle gange skrives som (95,42; 103,92)
# Vi kan øge konfidensniveauet til 99% ved at skrive conf.level = 0.99
t.test(x, conf.level = 0.99)
# Dette ændrer bredden af konfidensintervallet, som nu approksimativt er 93,99 < mu < 105,34
