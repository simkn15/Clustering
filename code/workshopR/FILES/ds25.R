# Antag nu, at vores dataramme indeholder NA-værdier, hvilket skaber et problem for os

# Her har vi en dataramme med NA-værdier
x <- c(-0.9714511, NA, 0.3367627, 1.7520504, 0.4918786)
y <- c(-0.4578746, 3.1663282, NA, 0.7406335, 1.4543427)
dfrm <- data.frame(x, y)
dfrm
# Her vil cumsum fejle, fordi inputtet indeholder NA-værdier
cumsum(dfrm)
# En løsning er, at fjerne alle rækker, der indeholder NA-værdier - det er hvad funktionen na.omit gør
na.omit(dfrm)
# Så kan cumsum fortsætte summeringen
cumsum(na.omit(dfrm))
# Dette virker også for vektorer og matricer, men ikke for lister

# Pas på: Hvis man smider for mange observationer væk, så kan det resultere i meningsløse resultater
# na.omit vli fjerne hele rækker, ikke bare NA-værdierne, hvilket kan eliminere meget brugbar information
