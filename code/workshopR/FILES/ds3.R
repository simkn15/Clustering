# Forståelse af "genbrugsreglen" - hvordan håndterer R vektorer af forskellig længde?

# Når vi udfører vektoraritmetik, så udfører R operationerne element-efter-element
# Det fungerer fint, når begge vektorer har samme længde: R parrer elementerne i vektorerne og anvender operationen på disse par

# Når vektorerne har forskellig længde, så bruger R den såkaldte "genbrugsregel":
# Den behandler vektorelementerne i par, startende med det første element i begge vektorer
# På et bestemt tidspunkt bliver elementerne i den korteste vektor opbrugt, men den længste har stadig ubehandlede elementer
# R går tilbage til begyndelsen af den korteste vektor og "genbruger" dens elementer, men fortsætter med den længste
# Den genbruger elementerne i den korteste vektor så ofte som det er nødvendigt indtil operationen er færdig

# Vektoren 1:6 er længere end vektoren 1:3, men R genbruger elementerne i 1:3 (se nedenstående "illustration")
(1:6) + (1:3)

#   1:6  1:3  (1:6)+(1:3)
#  -----------------------
#   1    1    2
#   2    2    4
#   3    3    6
#   4    1    5
#   5    2    7
#   6    3    9

# Også funktioner bruger denne regel
# Funktionen cbind kan lave søjlevektorer ud fra eksempelvis 1:6 og 1:3 (disse vil selvsagt ikke blive lige høje)
cbind(1:6)
cbind(1:3)
cbind(1:6, 1:3)

# Hvis den længste vektors længde ikke er et multiplum af den korteste vektors længde giver R en advarsel
(1:6) + (1:5)

# Når man forstår denne regel vil man indse, at operationer mellem en vektor og en skalar er anvendelser af den
# Her genbruges 10 gentagne gange ind til vektoradditionen er udført
(1:6) + 10
