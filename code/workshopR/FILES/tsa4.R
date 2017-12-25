library(zoo)

# Vi vender nu tilbage til et tidligere eksempel

# Priser
prices <- c(132.45, 130.85, 130.00, 129.55, 130.85)
# Datoer
dates <- as.Date(c("2010-01-04", "2010-01-05", "2010-01-06", "2010-01-07", "2010-01-08"))
# Lav et zoo-objekt, som indeholder prisen for IBM-aktien for de første fem dage af 2010
ibm.daily <- zoo(prices, dates)
ibm.daily
# Den i'te observation fra tidsrækken ibm.daily udvælges som ibm.daily[i]
ibm.daily[1]
ibm.daily[2]
# Vi kan også udvælge en observation efter dato (vores indeks er opbygget af Date-objekter)
ibm.daily[as.Date("2010-01-05")]
# Eller vi kan udvælge observationer efter en vektor af Date-objekter
dates <- seq(as.Date("2010-01-04"), as.Date("2010-01-08"), by = 2)
ibm.daily[dates]
# Funktionen window er lettere at bruge til udvælgelsen, hvis datoerne er fortløbende
window(ibm.daily, start = as.Date("2010-01-05"), end = as.Date("2010-01-07"))