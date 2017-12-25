# Funktionen barplot producerer et simpelt s?jlediagram
# Bem?rk at denne funktion antager, at h?jderne af s?jlerne er lagret i en vektor, hvilket ikke altid er tilf?ldet

# Ofte har vi en vektor af numeriske data og en parallel faktor, der grupperer data:
# M?ske vi gerne vil fremstille et s?jlediagram over gruppernes middelv?rdier eller gruppernes totaler?

# For eksempel indeholder det indbyggede airquality-datas?t en numerisk Temp-s?jle og en Month-s?jle

# Vi kan fremstille et s?jlediagram over middeltemperaturer per m?ned i to skridt

# F?rst beregner vi middelv?rdierne (h?jderne af s?jlerne)
heights <- tapply(airquality$Temp, airquality$Month, mean)
# Vi kan nu fremstille et s?jlediagram ud fra disse
barplot(heights)
# Vi kan med fordel tilf?je en titel og etiketter til s?jlediagrammet
barplot(heights,
        main = "Middeltemperatur per m?ned",
        names.arg = c("Maj", "Jun", "Jul", "Aug", "Sep"),
        ylab = "Temperatur (grader Fahrenheit)")
