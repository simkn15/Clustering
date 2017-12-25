# Antag nu, at vi ønsker at plotte en tæthedsfunktion

# Vi definerer en vektor x over definitionsmængden
x <- seq(from = -3, to = +3, length.out = 100)
# Vi anvender fordelingens tæthedsfunktion på x og plotter resultatet (her standardnormalfordelingen)
plot(x, dnorm(x))

# Vi vil i det følgende fremstille et 2x2-plot af fire tæthedsfunktioner

# Tæthedsfunktionens definitionsmængde
x <- seq(from = 0, to = 6, length.out = 100)
# Lav et 2x2-plotningsområde
par(mfrow = c(2, 2))
# Plot en uniform tæthedsfunktion
plot(x, dunif(x, min = 2, max = 4), main = "Uniform", type = "l", ylim = c(0, 0.6))
# Plot en normal tæthedsfunktion
plot(x, dnorm(x, mean = 3, sd = 1), main = "Normal", type = "l", ylim = c(0, 0.6))
# Plot en eksponentialtæthedsfunktion
plot(x, dexp(x, rate = 1/2), main = "Eksponential", type = "l", ylim = c(0, 0.6))
# Plot en gammatæthedsfunktion
plot(x, dgamma(x, shape = 2, rate = 1), main = "Gamma", type = "l", ylim = c(0, 0.6))

# Det følgende eksempel er en standardnormalfordeling, hvor området 1 <= z <= 2 er skraveret

# Vi fremstiller plottet ved først at plotte tæthedsfunktionen og dernæst skravere området med polygon-funktionen
# polygon-funktionen tegner en række linjer omkring det markerede område og udfylder det

# Først tegner vi tæthedskurven
x <- seq(from = -3, to = +3, length.out = 100)
y <- dnorm(x)
plot(x, y, main="Standardnormalfordeling", type = "l", xlab = "Kvartil", ylab = "Tæthed")
abline(h = 0)
# Vi definerer nu det interessante område med en række linjestykker, som er definerede ud fra (x, y)-koordinater
# Bemærk at polygon-funktionen vil forbinde det første og det sidste (x, y)-punkt for at lukke polygonen
# Her følger polygonen tæthedskurven, hvor 1 <= z <= 2
region.x <- x[1 <= x & x <= 2]
region.y <- y[1 <= x & x <= 2]
# Vi tilføjer begyndelses- og slutstykker
region.x <- c(region.x[1], region.x, tail(region.x, 1))
region.y <- c(0, region.y, 0)
# Til sidst kalder vi polygon-funktionen for at plotte kanten af området og udfylde det
polygon(region.x, region.y, density = 10)
# Som standard udfylder polygon ikke området, men density = 10 udfylder det med tynde linjer i en 45 graders vinkel
# Vi kan også vælge i stedet at udfylde området med en farve
polygon(region.x, region.y, density = -1, col = "red")
