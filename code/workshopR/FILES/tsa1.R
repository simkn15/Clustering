library(zoo)

# Priser
prices <- c(132.45, 130.85, 130.00, 129.55, 130.85)
# Datoer
dates <- as.Date(c("2010-01-04", "2010-01-05", "2010-01-06", "2010-01-07", "2010-01-08"))
# Lav et zoo-objekt, som indeholder prisen for IBM-aktien for de første fem dage af 2010
ibm.daily <- zoo(prices, dates)

# Priser
prices <- c(131.18, 131.20, 131.17, 131.15, 131.17)
# Antal timer fra midnat, startende fra klokken 9:30 (1 sekund = 0,00027778 time)
seconds <- c(9.5, 9.500278, 9.500556, 9.500833, 9.501111)
# Lav et zoo-objekt, som indeholder prisen for IBM-aktien i intervaller af 1 sekund
ibm.sec <- zoo(prices, seconds)

# Det rene data kan udtrækkes ved hjælp af funktionen coredata, som returnerer en simpel vektor
coredata(ibm.daily)
# Dato- eller tidsdelen kan udtrækkes ved hjælp af funktionen index
index(ibm.daily)