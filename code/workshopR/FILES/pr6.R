# Antag nu, at vi ønsker at generere en tilfældig række af simulerede møntkast eller Bernoulli-forsøg

# Vi kan bruge funktionen sample med tilbagelægning ved at skrive replace = TRUE

# Her genereres en række af 10 simulerede møntkast (P = plat, K = krone)
sample(c("P", "K"), 10, replace = TRUE)
# Her genereres en række af 20 Bernoulli-forsøg - tilfældige succeser eller fiaskoer (TRUE = succes, FALSE = fiasko)
sample(c(TRUE, FALSE), 20, replace = TRUE)
# Som standard vil sample vælge et uniformt tilfældigt element, så TRUE og FALSE har hver sandsynlighed 0,5
# Med et Bernoulli-forsøg er successandsynligheden p ikke nødvendigvis 0,5
# Her har TRUE sandsynlighed 0,8 og FALSE sandsynlighed 0,2
sample(c(TRUE, FALSE), 20, replace = TRUE, prob = c(0.8, 0.2))
# For det særlige tilfælde med binære værdier kan vi bruge rbinom, tilfældighedsgeneratoren for binomiale variate
rbinom(10, 1, 0.8)
