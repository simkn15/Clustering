library(isdals)

# Antag nu, at vi har en stikprøve fra hver af to populationer og ønsker at vide om de kunne have samme middelværdi

# Vi indlæser datasættet lucerne
data(lucerne)

# Udfør en t-test ved at kalde funktionen t.test
# Som standard antager t.test, at dine data ikke er parrede - hvis observationerne er parrede skriver vi paired = TRUE
# I begge tilfælde vil t.test beregne en p-værdi - og ifølge konventionen:
# - Hvis p < 0,05 er middelværdierne sandsynligvis forskellige
# - Hvis p > 0,05 er der ingen sådan evidens
# Hvis en af stikprøvernes størrelser er lille (dvs. mindre end 20 datapunkter) kræves normalfordeling af populationerne
# Hvis de to populationer har samme varians, så skriv var.equal = TRUE for at opnå en mindre konservativ test
t.test(lucerne$seeds.exp, lucerne$seeds.bent, paired = TRUE)
