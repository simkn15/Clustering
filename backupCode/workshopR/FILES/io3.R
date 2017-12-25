# Antag nu, at vi gerne vil indl?se en tekstfil, der indeholder en tabel af data

# Tabellariske datafiler er meget almindelige - det er tekstfiler med et simpelt format:
# - Hver linje indeholder en observation
# - Inden for hver linje er felter separeret af et skilletegn, f.eks. mellemrum, tabulator, kolon eller komma
# - Hver linje indeholder det samme antal felter

# Her ses et eksempel p?, hvordan en tabellarisk datafil kan v?re formateret:
# Fisher R.A. 1890 1962
# Pearson Karl 1857 1936
# Cox Gertrude 1900 1978
# Yates Frank 1902 1994
# Smith Kirstine 1878 1939
# J?rgensen Bent 1954 2015

# Den indbyggede funktion read.table kan l?se denne fil (som standard antager den, at datafelter er separeret af mellemrum)
dfrm <- read.table("statisticians.txt")
dfrm
# Husk at s?tte din arbejdsmappe!

# Hvis vores fil havde brugt et kolon som skilletegn kunne den have v?ret indl?st p? f?lgende m?de:
# dfrm <- read.table("statisticians.txt", sep = ":")

# Det kan ikke ses p? det printede output, men read.table fortolkede fornavn og efternavn som faktorer, ikke strenge
class(dfrm$V1)
# For at undg? at read.table fortolker karakterstrenge som faktorer kan vi s?tte stringsAsFactors-parameteret til FALSE
dfrm <- read.table("statisticians.txt", stringsAsFactor = FALSE)
# Nu er problemet l?st
class(dfrm$V1)

# Hvis et felt indeholder strengen "NA" antager read.table at v?rdien mangler og konverterer den til NA
# Vores datafil kunne m?ske have benyttet en anden streng eller et andet signal til at indikere manglende v?rdier?
# I dette tilf?lde kan vi benytte na.strings-parameteret (SAS-konventionen er f.eks. at benytte punktum):
# dfrm <- read.table("filename.txt", na.strings = ".")

# Vi kan ogs? indl?se en fil med s?jleoverskrifter
dfrm <- read.table("statisticians2.txt", header = TRUE, stringsAsFactor = FALSE)
dfrm
