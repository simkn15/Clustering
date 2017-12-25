# Antag nu, at vi givet en sandsnlighed p og en fordeling ønsker at bestemme den tilhørende kvartil for p:
# Dvs. værdien x, sådan at P(X <= x) = p

# Hver indbygget fordeling indeholder en kvartilfunktion, som konverterer sandsynligheder til kvartiler
# Funktionens navn er "q" efterfulgt af fordelingens navn (f.eks. er qnorm kvartilfunktionen for normalfordelingen)
# Kvartilfunktionens første argument er sandsynligheden og de resterende argumenter er fordelingens parametre
qnorm(0.05, mean = 100, sd = 15)
# Et almindeligt eksempel på beregning af kvartiler er, når vi beregner grænserne for et konfidensinterval
# Hvis vi gerne vil have 95% konfidensintervallet (alpha = 0,05) for en standardnormal variabel skal vi bruge
# kvartilerne med sandsynlighederne for alpha / 2 = 0,025 og (1 - alpha) / 2 = 0,975, som beregnes med qnorm:
qnorm(0.025)
qnorm(0.975)
# Det første argument kan imidlertid også være en vektor af sandsynligheder (så returneres en vektor af kvartiler)
qnorm(c(0.025, 0.975))

# Udforsk også funktioner som qbinom, qgeom, qpois, qt, qexp, qgamma og qchisq