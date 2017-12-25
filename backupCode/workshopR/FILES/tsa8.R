library(zoo)

# Vi vender nu igen-igen tilbage til et tidligere eksempel

# Priser
prices <- c(132.45, 130.85, 130.00, 129.55, 130.85)
# Datoer
dates <- as.Date(c("2010-01-04", "2010-01-05", "2010-01-06", "2010-01-07", "2010-01-08"))
# Lav et zoo-objekt, som indeholder prisen for IBM-aktien for de første fem dage af 2010
ibm.daily <- zoo(prices, dates)
ibm.daily
# Beregn differencerne for på hinanden følgende priser for IBM-aktien: (x_2 - x_1), (x_3 - x_2), (x_4 - x_3), ...
diff(ibm.daily)
# Forskellen med etiketten "2010-01-05" er ændringen i forhold til den foregående dato (2010-01-04)
# Bemærk at denne række er ét element kortere, da R ikke kan beregne ændringen for den første dato

# Forbrugerprisindeksdata (CPI-data) indlæses
cpi <- read.zoo("CPIAUCNS.csv", header = TRUE, sep = ",", format = "%Y-%m-%d")
# Som standard beregner diff på hinanden følgende forskelle, men man kan også beregne forskelle med større afstande
head(cpi, 24)
# Eksempelvis kan vi beregne ændringen fra de foregående 12 måneder, hvilket giver år-til-år-ændringen
diff(cpi, lag = 12)
