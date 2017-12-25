library(zoo)
library(xts)
library(tseries)

# Hent historiske finansielle data fra Yahoo Finance
ibm <- get.hist.quote(instrument = "ibm", quote = "Close", start = "1970-01-01", end = "1979-12-31")
# Se de 6 ældste observationer
head(ibm)
# Se de 6 nyeste observationer
tail(ibm)
# Se de 20 nyeste observationer
tail(ibm, 20)
# Funktionerne first og last fra xts-pakken bruger kalenderperioder frem for antal observationer
# Vi kan bruge first og last til at udvælge data efter antal dage, uger, måneder eller år
first(as.xts(ibm), "3 weeks")
last(as.xts(ibm), "month")
# Bemærk at en uge er defineret af kalenderen, ikke blot vilkårlige syv på hinanden følgende dage
# Brug eventuelt den indbyggede hjælpefunktion, hvis du vil vide mere om funktionerne first og last
help(first.xts)
help(last.xts)