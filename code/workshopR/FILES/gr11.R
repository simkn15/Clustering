# Her plottes et diagram med tre søjler, som farves henholdsvis rød, hvid og blå (her bruges eksplicitte farver)
barplot(c(3, 5, 4), col = c("red", "white", "blue"))

# Vi arbejder videre med det indbyggede airquality-datasæt
attach(airquality)
# I den forrige fil beregnede vi gruppernes middelværdier ved hjælp af mean-funktionen
heights <- tapply(Temp, Month, mean)

# En typisk effekt er, at skravere søjlerne efter deres rang (lavere søjler er lysere, højere søjler er mørkere)

# For at skravere søjlediagrammet konverterer vi søjlernes rang til relativ højde, udtrykt som en værdi mellem 0 og 1
rel.hts <- rank(heights) / length(heights)
# Vi kan bruge gray-funktionen til at generere en vektor af gråtoner (eller rainbow-funktionen for en vektor af farver)
# Som argument tager denne funktion en vektor af numeriske værdier mellem 0 og 1 og returnerer en gråtone for hvert element:
# 0.0 = kulsort, 1.0 = kridhvid
# Bemærk at vi her inverterer de relative højder, så højere søjler er mørke, ikke lyse!
grays <- gray(1 - rel.hts)
# Vi fremstiller et skraveret søjlediagram
barplot(heights, col = grays)
# Igen kan vi pynte lidt på diagrammet
rel.hts <- (heights - min(heights)) / (max(heights) - min(heights))
grays <- gray(1 - rel.hts)
barplot(heights,
        col = grays,
        ylim = c(50, 90), xpd = FALSE,
        main = "Middeltemperatur per måned",
        names.arg = c("Maj", "Jun", "Jul", "Aug", "Sep"),
        ylab = "Temperatur (grader Fahrenheit)")