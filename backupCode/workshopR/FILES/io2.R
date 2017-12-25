# Antag nu, at vores output indeholder for mange cifre eller for få cifre - og vi gerne vil printe flere eller færre

# R formaterer normalt output med decimaltal til at have 7 cifre
pi
100 * pi
# Dette fungerer for det meste fint, men det kan være irriterende, når man har mange tal og begrænset plads
# Det kan være misvisende, hvis der kun er få betydende cifre i dine tal og R fortsat printer 7
# Funktionen print lader dig variere antallet af printede cifre med digits-parameteret
print(pi, digits = 4)
print(100 * pi, digits = 4)
# Bemærk at print formaterer hele vektorer ad gangen
pnorm(-3:3)
print(pnorm(-3:3), digits = 3)
# Bemærk at print formaterer vektorelementer konsistent:
# Den finder antallet af cifre, som er nødvendige for at formatere det mindste tal
# Dernæst formaterer den alle tal til at have samme bredde (ikke nødvendigvis samme antal cifre)
# Dette er meget brugbart, når man skal formatere en hel tabel
q <- seq(from = 0, to = 3, by = 0.5)
tbl <- data.frame(Quant = q, Lower = pnorm(-q), Upper = pnorm(q))
# Uformateret print
tbl
# Formateret print: færre cifre
print(tbl, digits = 2)
# Du kan også ændre formatet af ALT output ved brug af options-funktionen
pi
options(digits = 15)
pi