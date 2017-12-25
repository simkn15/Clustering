# Antag nu, at vi ønsker at indlæse data fra en kommasepareret fil (CSV-fil)

# CSV-filformatet er populært, fordi mange programmer kan importere og eksportere data i dette format:
# R, Excel, andre regneark-programmer, mange databaseadministrationsværktøjer og statistiske pakker

# Det er en fil med tabellariske data, hvor hver linje i filen indeholder en række data med elementer separeret af komma
# Her er en meget simpel CSV-fil med tre rækker og tre søjler (første linje er overskrifter, som også er kommaseparerede)
# label,lbound,ubound
# low,0,0.674
# mid,0.674,1.64
# high,1.64,2.33

# Funktionen read.csv indlæser data og laver en dataramme, som er den sædvanlige R-repræsentation for tabellariske data
# Denne funktion antager at din fil indeholder søjleoverskrifter, medmindre vi fortæller den andet
tbl <- read.csv("tabledata.csv")
tbl
# Bemærk at read.csv tog søjlenavnene fra den første linje og brugte dem til datarammen
# Hvis ikke filen havde indeholdt overskrifter, så ville vi skrive header = FALSE
tbl <- read.csv("tabledata2.csv", header = FALSE)
tbl
# Bemærk at read.csv automatisk fortolker ikke-numeriske data som en faktor (kategorisk variabel)
# Lad os se på strukturen af tbl
str(tbl)
# Nogle gange ønsker vi virkelig at data skal fortolkes som strenge, ikke faktorer
# I dette tilfælde, sæt as.is-parameteret til TRUE
tbl <- read.csv("tabledata.csv", as.is = TRUE)
str(tbl)
