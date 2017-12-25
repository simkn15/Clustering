# Antag nu, at vi ?nsker at skrive data til en kommasepareret fil (CSV-fil)

# Indl?s noget data
tbl <- read.csv("tabledata.csv")
tbl
# Funktionen write.table skriver tabellariske data til en ASCII-fil i CSV-format
# Hver datar?kke laver en ny linje i filen med data separeret af komma
write.csv(tbl, file = "tabledata_modified.csv", row.names = TRUE)
# Filen oprettes altid i vores nuv?rende arbejdsmappe!

# Bem?rk at funktionen som standard skriver s?jleoverskrifter til filen (brug col.names = FALSE for at ?ndre dette)
# Hvis ikke vi skriver row.names = FALSE, s? f?r vi r?kkenumre i hver linje

# Funktionen write.table kan bruges til at gemme tabellariske data i andre formater