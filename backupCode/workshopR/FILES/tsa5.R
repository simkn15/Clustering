library(zoo)
library(tseries)

# Hent historiske finansielle data fra Yahoo Finance
ibm <- get.hist.quote(instrument = "ibm", quote = "Close", start = "1970-01-01", end = "1979-12-31")
# Forbrugerprisindeksdata (CPI-data) indlæses
cpi <- read.zoo("CPIAUCNS.csv", header = TRUE, sep = ",", format = "%Y-%m-%d")
# Bemærk at de to tidsrækker har forskellige tidsstempler, da den ene er daglige data og den anden er månedlige data
# og CPI-dataene er tidsstemplede for den første dag i hver måned, selv hvis denne dag er en helligdag eller weekend
merge(ibm, cpi)
# Funktionen merge finder som standard foreningsmængden af alle datoer
# Outputtet indeholder alle datoer fra begge input og manglende observationer udfyldes med NA-værdier
# Disse NA-værdier kan erstattes med den nyeste observation ved at bruge na.locf-funktionen
na.locf(merge(ibm, cpi))
# NA'erne er erstattet og funktionen eliminerede den første observation, da der ingen IBM-aktiepris var for denne dag
# Vi kan få fællesmængden af alle datoerne ved at sætte all = FALSE
merge(ibm, cpi, all = FALSE)
# Outputtet er begrænset til observationerne, som er fælles for begge filer

# Bemærk at fællesmængden begynder 1. april, ikke 1. januar (da 1/1, 1/2 og 1/3 alle er helligdage eller weekender)
# Der er ingen IBM-aktiepris for disse datoer og dermed ingen fællesmængde med CPI-data! Dette håndterer vi om lidt!

# locf står for "last observation carried forward"