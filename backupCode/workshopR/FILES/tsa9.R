library(zoo)

# Vi vender nu igen-igen-igen tilbage til et tidligere eksempel

# Priser
prices <- c(132.45, 130.85, 130.00, 129.55, 130.85)
# Datoer
dates <- as.Date(c("2010-01-04", "2010-01-05", "2010-01-06", "2010-01-07", "2010-01-08"))
# Lav et zoo-objekt, som indeholder prisen for IBM-aktien for de første fem dage af 2010
ibm.daily <- zoo(prices, dates)
ibm.daily
# Beregn differencerne for på hinanden følgende priser for IBM-aktien
diff(ibm.daily)
# For at beregne den procentvise ændring i IBM-aktien skal vi dividere den daglige ændring med prisen
# Disse to tidsrækker har forskellige starttider og forskellige længder!
# Heldigvis håndterer R dette for os, når vi dividerer en tidsrække med en anden!
diff(ibm.daily) / ibm.daily
# Vi kan skalere resultatet med 100 for at beregne den procentvise ændring
100 * (diff(ibm.daily) / ibm.daily)
# Hvis vi tager eksempelvis logaritmen eller kvadratroden bevares tidsstemplerne
log(ibm.daily)
# Det er også let at beregne differencen mellem logaritmerne til priserne
diff(log(ibm.daily))