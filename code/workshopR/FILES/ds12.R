# Antag nu, at en liste indeholder NULL-værdier, som vi ønsker at fjerne

# At finde og fjerne NULL-elementer kan være overraskende drilsk - her beskrives en fremgangsmåde:
# 1. R kalder sapply for at anvende funktionen is.null på hvert element i listen
# 2. sapply returnerer en vektor af logiske værdier, som er sande, hvis det tilhørende listeelement er NULL
# 3. R udvælger værdier fra listen udfra denne vektor
# 4. R tildeler NULL til de valgte elementer, hvilket fjerner dem fra listen

# Nysgerrige deltagere vil måske undre sig over, hvordan en liste kan indeholde NULL-elementer, når elementer fjernes ved
# at sætte dem til NULL, men det er dog muligt at oprette lister med NULL-elementer som det illustreres i dette eksempel:
mylist <- list("Hej", NULL, "dig")
mylist
# Her fjernes NULL-elementer fra mylist
mylist[sapply(mylist, is.null)] <- NULL
mylist
