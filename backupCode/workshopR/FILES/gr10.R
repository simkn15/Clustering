library(gplots)

# Vi arbejder videre med det indbyggede airquality-datas?t
attach(airquality)

# S?jlediagrammer viser ofte punktestimater ved hj?lp af s?jlernes h?jder, men indeholder sj?ldent konfidensintervaller

# I den forrige fil beregnede vi gruppernes middelv?rdier ved hj?lp af mean-funktionen
heights <- tapply(Temp, Month, mean)
# Vi kan beregne konfidensintervaller for disse middelv?rdier p? en tilsvarende m?de
# Funktionen t.test returnerer en liste af v?rdier, heriblandt elementet conf.int, som er en vektor af to elementer:
# conf.int[1] indeholder den nedre gr?nse og conf.int[2] indeholder den ?vre gr?nse
lower <- tapply(Temp, Month, function(v) t.test(v)$conf.int[1])
upper <- tapply(Temp, Month, function(v) t.test(v)$conf.int[2])
# Vi kan nu fremstille et s?jlediagram med konfidensintervaller ved hj?lp af funktionen barplot2
barplot2(heights, plot.ci = TRUE, ci.l = lower, ci.u = upper)
# Vi kan eksempelvis trimme s?jlerne (xpd), tilf?je en titel (main) og give s?jlerne etiketter (names.arg)
barplot2(heights, plot.ci = TRUE, ci.l = lower, ci.u = upper,
         ylim = c(50, 90), xpd = FALSE,
         main = "Middeltemperatur per m?ned",
         names.arg = c("Maj", "Jun", "Jul", "Aug", "Sep"),
         ylab = "Temperatur (grader Fahrenheit)")
