# Antag nu, at vores data er organiseret efter rækker og vi ønsker at samle det i en dataramme

# Data kommer ofte som en samling observationer - hvor vi får rækker én ad gangen, frem for søjler én ad gangen

# Hver observation er en tupel, der indeholder flere værdier - én for hver observeret variabel
# I det følgende lagres hver tupel som en dataramme med én række (hvis man kun har numeriske data kan vektorer bruges)
obs <- list()
obs[[1]] <- data.frame(pred1 = -2.7528917, pred2 = -1.40784130, pred3 = "AM", resp = 12.57715)
obs[[2]] <- data.frame(pred1 = -0.3626909, pred2 = 0.31286963, pred3 = "AM", resp = 21.02418)
obs[[3]] <- data.frame(pred1 = -1.0416039, pred2 = -0.69685664, pred3 = "PM", resp = 18.94694)
obs[[4]] <- data.frame(pred1 = 1.2666820, pred2 = -1.27511434, pred3 = "PM", resp = 18.98153)
obs[[5]] <- data.frame(pred1 = 0.7806372, pred2 = -0.27292745, pred3 = "AM", resp = 19.59455)
obs[[6]] <- data.frame(pred1 = -1.0832624, pred2 = 0.73383339, pred3 = "AM", resp = 20.71605)
obs[[7]] <- data.frame(pred1 = -2.0883305, pred2 = 0.96816822, pred3 = "PM", resp = 22.70062)
obs[[8]] <- data.frame(pred1 = -0.7063653, pred2 = -0.84476203, pred3 = "PM", resp = 18.40691)
obs[[9]] <- data.frame(pred1 = -0.8394022, pred2 = 0.31530793, pred3 = "PM", resp = 21.00930)
obs[[10]] <- data.frame(pred1 = -0.4966884, pred2 = -0.08030948, pred3 = "AM", resp = 19.31253)

# Vi binder disse rækker sammen til en dataramme, hvilket er hvad rbind-funktionen gør

# Her bindes de to første observationer sammen
rbind(obs[[1]], obs[[2]])
# Vi vil gerne binde alle observationer sammen, ikke blot de første to - så vi bruger do.call-funktionen
do.call(rbind, obs)

# Hvis obs er en liste af lister frem for en liste af datarammer med én række:
# Vi kan transformere rækkerne til datarammer ved brug af Map-funktionen og så bruge ovenstående fremgangsmåde:
# dfrm <- do.call(rbind, Map(as.data.frame, obs))
