# Antag nu, at vi ønsker at udvælge søjler fra en dataramme efter deres position

# For at udvælge en enkelt søjle, brug denne listeoperator:
# dfrm[[n]]
# Denne returnerer én søjle - mere specifikt den n'te søjle af datarammen dfrm

# For at udvælge en eller flere søjler og pakke dem ind i en dataramme, brug disse udtryk:
# dfrm[n] - returnerer en dataramme bestående af kun den n'te søjle af dfrm
# dfrm[c(n_1, n_2, ..., n_k)] - returnerer en dataramme opbygget af søjlerne på positionerne n_1, n_2, ..., n_k i dfrm
# Man kan også bruge matrixlignende indeksering til at udvælge en eller flere søjler:
# dfrm[, n]
# dfrm[, c(n_1, n_2, ..., n_k)]

# Datasættet suburbs indlæses
suburbs <- read.csv("suburbs.csv")
suburbs
# Brug simpel listenotation til at udvælge præcis én søjle, eksempelvis den første
suburbs[[1]]
# Første søjle pakket ind i en dataramme
suburbs[1]
# Første og tredje søjle pakket ind i en dataramme
suburbs[c(1, 3)]

# Som nævnt kan matrixnotationen også bruges til at udvælge søjler, men pas på:
# Nogle gange får du en søjle, andre gange en dataramme - alt efter hvor mange indekser du bruger!

# Her returneres en søjle
suburbs[, 1]
# Her returneres en dataramme ved brug af samme matrixlignende syntaks med flere indekser
suburbs[, c(1, 4)]
# Hvis vi også i førstnævnte tilfælde ønsker en dataramme frem for en søjle kan vi tilføje drop = FALSE
suburbs[, 1, drop = FALSE]

# Brug eventuelt i stedet listeoperatorerne beskrevet ovenfor
