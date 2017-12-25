# Antag nu, at vi ønsker at beregne fordelingsfunktionen for en kontinuert stokastisk variabel

# Brug fordelingsfunktionen, som beregner P(X <= x):
# For alle indbyggede sandsynlighedsfordelinger findes en fordelingsfunktion, hvis navn er "p" efterfulgt af fordelingens navn

# Den signifikante forskel fra forrige fil er, at kontinuerte variable ikke har en "sandsynlighed" i et enkelt punkt, P(X = x)
# I stedet har de en tæthed i et punkt!

# Det antages at mænds højder er normalfordelte med en middelværdi på 70 tommer og en standardafvigelse på 3 tommer
# Vi kan bruge pnorm til at beregne sandsynligheden for, at en mand er lavere end 66 tommer under denne antagelse
# Udtrykt mere matematisk vil vi gerne beregne P(X <= 66) givet at X ~ N(70, 3)
pnorm(66, mean = 70, sd = 3)
# Tilsvarende kan vi beregne sandsynligheden for, at en eksponentialfordelt variabel med middelværdi 40 kan være mindre end 20
pexp(20, rate = 1/40)
# Ligesom for diskrete sandsynligheder kan vi få overlevelsesfunktionen, P(X > x), ved at skrive lower.tail = FALSE
# Her beregnes sandsynligheden for, at den føromtalte eksponentialfordelte variabel kan være større end 50
pexp(50, rate = 1/40, lower.tail = FALSE)
# Intervalsandsynligheden for en kontinuert stokastisk variabel beregnes som differencen mellem to kumulative sandsynligheder:
# P(x_1 < X < x_2) = P(X < x_2) - P(X < x_1)
# For den samme eksponentialfordelte variabel beregnes sandsynligheden for, at den kan falde mellem 20 og 50, P(20 < X < 50)
pexp(50, rate = 1/40) - pexp(20, rate = 1/40)
