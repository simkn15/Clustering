library(zoo)
library(tseries)

# Hent historiske finansielle data fra Yahoo Finance
ibm <- get.hist.quote(instrument = "ibm", quote = "Close", start = "1970-01-01", end = "1979-12-31")
# Forbrugerprisindeksdata (CPI-data) indlæses
cpi <- read.zoo("CPIAUCNS.csv", header = TRUE, sep = ",", format = "%Y-%m-%d")
# Så vidt R ved, så har vi kun en observation for den første dag i hver måned og ingen observationer for de andre dage
# Dog ved vi, at hver CPI-værdi gælder for de efterfølgende dage i måneden
# Vi laver et objekt af længde nul med hver dag i årtiet, men ingen data
dates <- seq(from = as.Date("1970-01-01"), to = as.Date("1979-12-31"), by = 1)
empty <- zoo(, dates)
# Så tager vi foreningsmængden af CPI-data og objektet af længde nul, hvilket giver et datasæt fyldt med NA-værdier
filled.cpi <- merge(cpi, empty, all = TRUE)
filled.cpi
# Den resulterende tidsrække indeholder hver kalenderdag med NA'er, hvor der ikke var en observation
# Et mere almindeligt behov er, at erstatte hver NA med den nyeste observation fra og med denne dato
# Funktionen na.locf fra zoo-pakken gør præcist dette!
filled.cpi <- na.locf(merge(cpi, empty, all = TRUE))
filled.cpi
# Januars værdi (37,8) videreføres indtil 1. februar, hvor den erstattes af februars værdi (38,0)
# Hver dag har den nyeste CPI-værdi fra og med denne dato, hvilket kan bruges til at løse problemet fra forrige fil
# Husk: Den daglige pris på IBM-aktien og de månedlige CPI-data havde ingen fællesmængde på visse dage!
# En løsning er, at udvide IBM-data til at inkludere CPI-datoerne og tage fællesmængden (giver månedlige observationer)
filled.ibm <- na.locf(merge(ibm, zoo(, index(cpi))))
merge(filled.ibm, cpi, all = FALSE)
# En anden løsning er, at udfylde CPI-data og så tage fællesmængden med IBM-data (giver daglige observationer)
merge(ibm, filled.cpi, all = FALSE)
