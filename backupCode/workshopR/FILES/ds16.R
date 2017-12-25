# Antag nu, at vi ønsker at tildele beskrivende navne til rækkerne eller søjlerne i en matrix

# Betragt denne matrix af korrelationer mellem priserne for IBM-, Microsoft- og Google-aktien
tech.corr <- matrix(c(1.000, 0.556, 0.390, 0.556, 1.000, 0.444, 0.390, 0.444, 1.000), 3, 3)
tech.corr
# Vi tilføjer række- og søjlenavne, som øger læsbarheden af outputtet
rownames(tech.corr) <- c("IBM", "MSFT", "GOOG")
colnames(tech.corr) <- c("IBM", "MSFT", "GOOG")
tech.corr
# Nu ved læseren hvilke rækker og søjler, der hører til hvilke aktier

# En anden fordel er, at vi nu kan referere til matrixelementer med disse navne
# Hvad er korrelationen mellem IBM og GOOG?
tech.corr["IBM", "GOOG"]
