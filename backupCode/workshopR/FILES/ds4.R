# I R kaldes hver mulig værdi af en kategorisk variabel for et niveau og en vektor af niveauer kaldes en faktor

# Antag nu, at vi har en vektor med karakterstrenge eller heltal, som vi ønsker at R skal behandle som en faktor

# R-funktionen factor koder din vektor af diskrete værdier som en faktor:
# f <- factor(v)
# Hvis din vektor kun indeholder en delmængde af de mulige værdier, så tilføj et andet argument med de mulige niveauer
# f <- factor(v, levels)

# Konvertering af dine kategoriske data til en faktor kræver for det meste blot et kald til factor-funktionen
f <- factor(c("Sejr", "Sejr", "Nederlag", "Uafgjort", "Sejr", "Nederlag"))
f
# Bemærk at factor-funktionen identificerer de forskellige niveauer i de kategoriske data og pakker dem i en faktor

# Hvis din vektor kun indeholder en delmængde af de mulige niveauer, så vil R have et ufuldstændigt billede
# Antag at wday indeholder strengværdier, der giver ugedagen, hvor data blev observeret
wday <- c("Ons", "Tor", "Man", "Ons", "Tor", "Tor", "Tor", "Tir", "Tor", "Tir")
f <- factor(wday)
f
# R tror at mandag, onsdag, tirsdag og torsdag er de eneste mulige niveauer (fredag er ikke listet)
# Laboratoriepersonalet har tilsyneladende ingen observationer fra fredage, så R ved ikke at det er en mulig værdi
# Derfor lister vi nu eksplicit de mulige niveauer for wday
f <- factor(wday, c("Man", "Tir", "Ons", "Tor", "Fre"))
f
# Nu forstår R, at f er en faktor med fem forskellige niveauer - og den kender den korrekte rækkefølge

# Oprindeligt sorterede den dem alfabetisk efter ugedagenes navne på dansk

# I mange situationer er det ikke nødvendigt at kalde factor-funktionen eksplicit
# Når en R-funktion skal bruge en faktor konverterer den ofte automatisk data til en faktor