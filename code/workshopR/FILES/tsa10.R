library(zoo)
library(xts)
library(tseries)

# Hent historiske finansielle data fra Yahoo Finance
ibm <- get.hist.quote(instrument = "ibm", quote = "Close", start = "1970-01-01", end = "1974-12-31")
# Vi kan beregne den gennemsnitlige pris per måned, hvis vi bruger funktionerne apply.monthly og mean sammen
apply.monthly(as.xts(ibm), mean)
# En mere interessant anvendelse er beregning af volatilitet per kalendermåned
# Volatilitet er målt som standardafvigelsen af daglige log-afkast, der beregnes som diff(log(ibm))
diff(log(ibm))
# Vi beregner deres standardafvigelse, måned for måned
apply.monthly(as.xts(diff(log(ibm))), sd)
# Vi skalerer det daglige antal for at estimere den annualiserede volatilitet (der er i gennemsnit 251 handelsdage per år)
sqrt(251) * apply.monthly(as.xts(diff(log(ibm))), sd)
# Vi plotter resultatet (figur 5)
plot(sqrt(251) * apply.monthly(as.xts(diff(log(ibm))), sd),
     main="IBM: Månedlig volatilitet",
     ylab="Annualiseret standardafvigelse")

# Se også funktionerne apply.daily, apply.weekly, apply.quarterly og apply.yearly