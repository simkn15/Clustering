# Antag nu, at vi ønsker at beregne den simple eller den kumulative sandsynlighed forbundet med en diskret stokastisk variabel

# For den simple sandsynlighed, P(X = x), brug tæthedsfunktionen:
# For alle indbyggede sandsynlighedsfordelinger findes en tæthedsfunktion, hvis navn er "d" efterfulgt af fordelingens navn
# For en kumulativ sandsynlighed, P(X <= x), brug fordelingsfunktionen:
# For alle indbyggede sandsynlighedsfordelinger findes en fordelingsfunktion, hvis navn er "p" efterfulgt af fordelingens navn

# Lad os antage, at vi har en binomial stokastisk variabel X over 10 forsøg, hvor hvert forsøg har successandsynlighed 1/2
# Så kan vi beregne sandsynligheden for at observere x = 7 ved at kalde dbinom
dbinom(7, size = 10, prob = 0.5)
# Fordelingsfunktionen for binomialfordelingen hedder pbinom og her beregnes den kumulative sandsynlighed for x = 7
pbinom(7, size = 10, prob = 0.5)
# Sandsynligheden for at observere X <= 7 er omkring 0,945

# Komplementet til den kumulative sandsynlighed er overlevelsesfunktionen, P(X > x)
# Alle fordelingsfunktioner lader os finde denne højre halesandsynlighed ved at skrive lower.tail = FALSE
pbinom(7, size = 10, prob = 0.5, lower.tail = FALSE)
# Sandsynligheden for at observere X > 7 er omkring 0,055

# Intervalsandsynligheden, P(x_1 < X <= x_2), er sandsynligheden for at observere X mellem grænserne x_1 og x_2
# Den beregnes blot som differencen mellem de to kumulative sandsynligheder: P(X <= x_2) - P(X <= x_1)
# Her beregnes P(3 < X <= 7) for vores binomiale variabel
pbinom(7, size = 10, prob = 0.5) - pbinom(3, size = 10, prob = 0.5)

# R lader os specificere flere x-værdier og returnerer en vektor af de tilhørende sandsynligheder
# Her beregnes to kumulative sandsynligheder, P(X <= 3) og P(X <= 7), i ét kald til pbinom
pbinom(c(3, 7), size = 10, prob = 0.5)
# Vi kan således beregne intervalsandsynligheder ved hjælp af diff-funktionen
diff(pbinom(c(3, 7), size = 10, prob = 0.5))

# diff-funktionen beregner differencen mellem på hinanden følgende elementer i en vektor