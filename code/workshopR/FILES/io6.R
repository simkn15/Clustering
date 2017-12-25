# Antag nu, at vi ?nsker at indl?se tabellariske eller CSV-data direkte fra internettet til vores arbejdsrum

# B?de read.csv og read.table kan indl?se data direkte fra internettet, hvis du angiver en URL i stedet for et filnavn
armeben <- read.table("http://imada.sdu.dk/~chjoe11/ST520/F17/armeben.txt", header = TRUE)
armeben
