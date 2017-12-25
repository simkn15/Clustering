# Antag nu, at vi har flere grupper af data med en vektor for hver gruppe
# Vi vil gerne kombinere vektorerne og lave en parallel faktor, der identificerer hver v?rdis oprindelige gruppe

# Eksempel:
# I en sp?rgeskemaunders?gelse er b?de 1. ?rsstuderende, 2. ?rsstuderende og 3. ?rsstuderende blevet spurgt om f?lgende:
# "Hvor stor en procentdel af tiden f?ler du dig fortr?stningsfuld p? universitetet?"

# F?rste ?r
fy <- c(0.60, 0.35, 0.44, 0.62, 0.60)
# Andet ?r
sy <- c(0.70, 0.61, 0.63, 0.87, 0.85, 0.70, 0.64)
# Tredje ?r
ty <- c(0.76, 0.71, 0.92, 0.87)
# Lav en liste, som indeholder vektorerne - og brug stack-funktionen til at kombinere listen til en dataramme med to s?jler
comb <- stack(list(first = fy, second = sy, third = ty))
# Datarammens s?jler kaldes values og ind, hvor f?rste s?jle indeholder data og anden s?jle indeholder den parallelle faktor
comb
# Vi kan nu udf?re eksempelvis en variansanalyse p? de to s?jler (denne funktion kr?ver data i ovenst?ende format)
# ANOVA-funktionen, aov, kr?ver en vektor med unders?gelsens resultater og en parallel faktor, som identificerer gruppen
aov(values ~ ind, data = comb)
