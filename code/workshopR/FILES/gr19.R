library(MASS)

# Cars93-datasættet indeholder en prissøjle og vi vil undersøge om den er normalfordelt ved at fremstille et QQ-plot
qqnorm(Cars93$Price, main="QQ-plot: Price")
# En diagonal linje tilføjes til QQ-plottet
qqline(Cars93$Price)
# Hvis data havde en perfekt normalfordeling, så ville punkterne ligge præcist på den diagonale linje
# Mange punkter er tæt på linjen, specielt i den midterste del, men punkterne i halerne er temmelig langt væk!

# For mange punkter over linjen indikerer en generel venstreskævhed, hvorfor vi laver en logaritmisk transformation
qqnorm(log(Cars93$Price), main="QQ-plot: log(Price)")
qqline(log(Cars93$Price))
# Det ser ud som om, at log(Price) er approksimativt normalfordelt