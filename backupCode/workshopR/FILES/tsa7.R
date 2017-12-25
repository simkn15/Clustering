library(zoo)

# Vi vender nu igen tilbage til et tidligere eksempel

# Priser
prices <- c(132.45, 130.85, 130.00, 129.55, 130.85)
# Datoer
dates <- as.Date(c("2010-01-04", "2010-01-05", "2010-01-06", "2010-01-07", "2010-01-08"))
# Lav et zoo-objekt, som indeholder prisen for IBM-aktien for de første fem dage af 2010
ibm.daily <- zoo(prices, dates)
ibm.daily
# Data forskydes en dag frem med k = +1 (morgendagens data bliver dagens data)
lag(ibm.daily, k = +1, na.pad = TRUE)
# Data forskydes en dag tilbage med k = -1 (gårsdagens data bliver dagens data)
lag(ibm.daily, k = -1, na.pad = TRUE)
# Vi kan skriver na.pad = TRUE for at udfylde med NA (ellers ville tidsrækken blot blive forkortet)