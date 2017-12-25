# Antag nu, at vi har en lille datamængde, som er for lille til at retfærdiggøre at oprette en inputfil
# Vi vil gerne blot indtaste data direkte til vores arbejdsrum

# For meget små datasæt kan vi med fordel indtaste data ved hjælp af c() konstruktøren til vektorer
scores <- c(61, 66, 90, 88, 100)
# Alternativ kan man lave en tom dataramme og bruge det indbyggede regnearklignende redigeringsværktøj til at udfylde den
scores <- data.frame()
scores <- edit(scores)
