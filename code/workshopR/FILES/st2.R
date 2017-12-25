# Antag nu, at vi har en stikprøve af værdier fra en population bestående af succeser og fiaskoer
# Baseret på stikprøvedata vil vi gerne konstruere et konfidensinterval for populationens succesdel

# Antag at et nyhedsbrev indeholder et afsnit, der foregiver at identificere aktier, som med stor sandsynlighed vil stige
# Der meldes eksempelvis, at en bestemt aktie fulgte et mønster - og at den steg 6 ud af 9 gange dette mønster forekom
# Forfatterne konkluderer, at sandsynligheden for at aktien stiger igen derfor er 6/9 eller 66,7%

# Ved brug af prop.test kan vi opnå et konfidensinteral for den sande proportion af gangene aktien stiger efter mønsteret
# Her er antallet af observationer n = 9 og antallet af succeser er x = 6
# Outputtet indeholder et 95% konfidensinterval: (0,309; 0,910)
prop.test(6, 9)

# Det er måske ikke så klogt at sige, at sandsynligheden for en stigning er 66,7%
