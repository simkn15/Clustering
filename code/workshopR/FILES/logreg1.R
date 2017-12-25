library(faraway)

# Et almindeligt problem inden for statistisk modellering er, at prædiktere et udfald med binære værdier:
# - Vil en behandling være effektiv eller ej?
# - Vil priserne stige eller falde?
# - Hvem vil vinde spillet, hold A eller hold B?
# Logistisk regression er brugbart til at modellere disse situationer!
# Det giver ikke bare et "thumbs up"- eller "thumbs down"-svar, men beregner en sandsynlighed for hvert af de to mulige udfald

# For at udføre logistisk regression kan man kalde glm-funktionen med family = binomial

# Datasættet pima fra faraway-pakken indlæses
data(pima)
# Vis de første observationer
head(pima)
# Faraway giver et eksempel på prædiktion af en variabel med binære værdier: test fra datasættet pima er sand, hvis patienten
# er testet positiv for diabetes. Prædiktorerne er diastolisk blodtryk og BMI (Faraway laver dog i stedet lineær regression!)
b <- factor(pima$test)
m <- glm(b ~ diastolic + bmi, family = binomial, data = pima)
# Opsummeringen af den resulterende model viser, at p-værdierne for diastolic- og bmi-variablene er 0,805 og (numerisk) 0
summary(m)
# Vi kan med andre ord konkludere, at kun bmi-variablen er signifikant, hvorfor modellen kan reduceres
m.reduced <- glm(b ~ bmi, family = binomial, data = pima)
# Lad os bruge modellen til at beregne sandsynligheden for, at en med gennemsnitligt BMI (32,0) testes positiv for diabetes
newdata <- data.frame(bmi = 32.0)
predict(m.reduced, type = "response", newdata = newdata)
# Ifølge denne model er sandsynligheden omkring 33,3%
newdata <- data.frame(bmi = quantile(pima$bmi, 0.90))
predict(m.reduced, type = "response", newdata = newdata)
# Den samme beregning for nogen i den 90. percentil giver en sandsynlighed på 54,9%

# I kaldet til predict-funktionen skriver vi type = "response" (så den returnerer en sandsynlighed frem for log-odds)