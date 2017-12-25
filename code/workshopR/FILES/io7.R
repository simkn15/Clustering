# Antag nu, at vi ønsker at vælge en fil interaktivt
myData <- file.choose()
# Eller at vi ønsker at indlæse tabellariske data fra udklipsholderen
myData <- read.table("clipboard", header = TRUE)
