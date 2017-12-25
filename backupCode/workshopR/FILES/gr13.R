# Bemærk: Denne fil indeholder udelukkende beskrivelser af fremgangsmåder og producerer intet output!

# Hvis du ønsker at ændre type, bredde eller farve på en linje kan du tilføje yderligere parametre:

# Brug lty-parameteret til at vælge linjetype:
# lty = "solid" eller lty = 1 (standard)
# lty = "dashed" eller lty = 2
# lty = "dotted" eller lty = 3
# lty = "dotdash" eller lty = 4
# lty = "longdash" eller lty = 5
# lty = "twodash" eller lty = 6
# lty = "blank" eller lty = 0 (undertrykker tegning)

# Dette kald til plot-funktionen tegner eksempelvis en stiplet linje:
# plot(x, y, type = "l", lty = "dashed")

# Brug lwd-parameteret til at vælge linjens bredde eller tykkelse (som standard har de tykkelse 1):
# plot(x, y, type = "l", lwd = 2)

# Brug col-parameteret til at vælge linjens farve (som standard er de sorte):
# plot(x, y, type = "l", col = "red")

# Kaldet til plot-funktionen initialiserer grafikvinduet og tegner de første linjer med blåt
# plot(x, y1, type = "l", col = "blue")
# De efterfølgende kald til lines-funktionen tegner yderligere linjer, først med rødt og så med gult
# lines(x, y2, col = "red")
# lines(x, y3, col = "yellow")