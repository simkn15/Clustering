library(Lahman)

# Antag nu, at vi ?nsker at udtage tilf?ldige elementer fra et datas?t

# Modern World Series-data (delm?ngde af BattingPost-data fra BattingLahman-pakken) indl?ses
ws <- subset(BattingPost, round == "WS" & yearID >= 1903)
# Vi kan udv?lge 10 tilf?ldige ?r ved hj?lp af sample-funktionen
sample(ws$yearID, 10)
# Funktionen sample udv?lger normalt elementer uden tilbagel?gning (dvs. de samme elementer udv?lges ikke to gange)

# Nogle statistiske procedurer (specielt bootstrap) kr?ver udv?lgelse med tilbagel?gning
sample(ws$yearID, 10, replace = TRUE)
# Dette betyder, at de samme elementer kan optr?de flere gange i stikpr?ven