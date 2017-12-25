library(zoo)
library(tseries)
library(quantmod)

# Hent historiske finansielle data fra Yahoo Finance
ibm.infl <- get.hist.quote(instrument = "ibm", quote = c("Close", "AdjClose"), start = "1970-01-01", end = "1979-12-31")
# Normaliser til en startværdi på 100
ibm.infl$Close <- ibm.infl$Close / as.numeric(ibm.infl$Close[1]) * 100
ibm.infl$AdjClose <- ibm.infl$AdjClose / as.numeric(ibm.infl$AdjClose[1]) * 100
# Fremstil et plot (figur 1)
plot(ibm.infl, screens = 1, lty = c("solid", "dotted"),
     main = "IBM: Historisk versus inflationsjusteret",
     xlab = "Dato", ylab = "Relativ pris", ylim = range(coredata(ibm.infl)))
# Tilføj en signaturforklaring
legend("topleft", horiz = TRUE, legend = c("Hist", "Just"), lty = c("dotted", "solid"), bty = "n")

# Hent historiske finansielle data fra Yahoo Finance
aapl <- get.hist.quote(instrument = "aapl", quote = c("Close", "Volume"), start = "2000-01-01", end = "2009-12-31")
goog <- get.hist.quote(instrument = "goog", quote = c("Close", "Volume"), start = "2000-01-01", end = "2009-12-31")
msft <- get.hist.quote(instrument = "msft", quote = c("Close", "Volume"), start = "2000-01-01", end = "2009-12-31")
# Fremstil et plot (figur 2)
plot(aapl$Close, main = "Sammenligning af aktiepriser", ylim = c(0, 800), col = "gray", xlab = "Dato", ylab = "Aktiepris (USD)")
lines(goog$Close, col = "blue")
lines(msft$Close, col = "red")
# Tilføj en signaturforklaring
legend("topleft",  horiz = TRUE, legend = c("Apple", "Google", "Microsoft"), col = c("gray", "blue", "red"), lty = 1, bty = "n")

# Fremstil et diagram ved hjælp af quantmod-pakken, som indeholder indbyggede funktioner til visualisering af aktiedata (figur 3)
getSymbols("AAPL", from = "2000-01-01", to = "2009-12-31")
barChart(AAPL)
# Et lignende diagram i et andet farveskema (figur 4)
candleChart(AAPL, theme = "white")