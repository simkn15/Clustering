# Antag nu, at vores data er organiseret efter søjler og vi ønsker at samle det i en dataramme

# Vi har i dette eksempel to numeriske prædiktorvariable, en kategorisk prædiktorvariabel og en responsvariabel

pred1 <- c(-2.7528917,  -0.3626909,  -1.0416039,   1.2666820,   0.7806372,
           -1.0832624,  -2.0883305,  -0.7063653,  -0.8394022,  -0.4966884)
pred2 <- c(-1.40784130,  0.31286963, -0.69685664, -1.27511434, -0.27292745,
            0.73383339,  0.96816822, -0.84476203,  0.31530793, -0.08030948)

pred3 <- c("AM", "AM", "PM", "PM", "AM", "AM", "PM", "PM", "PM", "AM")

resp  <- c(12.57715, 21.02418, 18.94694, 18.98153, 19.59455,
           20.71605, 22.70062, 18.40691, 21.00930, 19.31253)

# En dataramme er en samling søjler, som hver svarer til en observeret variabel
# Hvis data allerede er organiseret i søjler er det let at opbygge en dataramme

# Funktionen data.frame kan konstruere en dataramme udfra vektorer, hvor hver vektor er en observeret variabel
dfrm <- data.frame(pred1, pred2, pred3, resp)
dfrm
# Funktionen data.frame bruger søjlenavnene fra programmets variable, men vi kan eksplicit angive søjlenavne
dfrm <- data.frame(p1 = pred1, p2 = pred2, p3 = pred3, r = resp)
dfrm

# Måske er vores data organiseret i vektorer, som er samlet i en liste frem for at være individuelle programvariable
mylist <- list(p1 = pred1, p2 = pred2, p3 = pred3, r = resp)
# I dette tilfælde kan funktionen as.data.frame bruges til at lave en dataramme ud fra listen af vektorer
as.data.frame(mylist)
