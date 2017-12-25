# Bemærk: Denne fil indeholder udelukkende beskrivelser af fremgangsmåder og producerer intet output!

# Antag nu, at vi ønsker at fjerne elementer fra en liste ud fra en betingelse
# Måske vi ønsker at fjerne elementer, som er større end eller mindre end en bestemt værdi?

# Husk:
# 1. En liste kan indekseres af en logisk vektor - når vektorelementet er TRUE, så vælges det tilhørende listeelement
# 2. Vi kan fjerne listeelementer ved at tildele værdien NULL til dem

# Her fjernes elementer fra en liste, mylist, hvis elementernes værdier er 0:
# mylist[mylist == 0] <- NULL
# Bemærk at vi konstruerede en logisk vektor, som identificerede de uønskede værdier (mylist == 0)
# Så udvalgte vi hereffter disse elementer fra listen og tildelte NULL til dem

# Dette udtryk vil fjerne NA-værdier fra listen:
# mylist[is.na(mylist)] <- NULL

# Så langt, så godt - problemerne opstår, når ikke vi let kan opbygge den logiske vektor!
# Dette sker ofte, når man bruger en funktion, der ikke kan håndtere en liste (eksempelvis abs-funktionen)
# Antag at vi ønsker at fjerne listeelementer, hvis absolutte værdier er mindre end 1
# mylist[abs(mylist) < 1] <- NULL vil give følgende fejl: Error in abs(mylist) : non-numeric argument to function
# Den simpleste løsning er, at udflade listen til en vektor ved at kalde unlist og så teste på vektoren:
# mylist[abs(unlist(mylist)) < 1] <- NULL
# En mere elegant løsning ville være at bruge lapply-funktionen til at anvende abs på hvert element i listen:
# mylist[lapply(mylist, abs) < 1] <- NULL

# OPGAVE: Afprøv dette på en liste kaldet mylist indeholdende tallet 0, en NA-værdi og forskellige positive og negative tal

# Lister kan også indeholde komplekse objekter, ikke blot atomare værdier
# Antag at mymodels er en liste af lineære modeller, som er lavet med lm-funktionen
# Dette udtryk vil fjerne enhver model, hvis R^2-værdi er mindre end 0,30:
# mymodels[sapply(mymodels, function(m) summary(m)$r.squared < 0.3)] <- NULL
