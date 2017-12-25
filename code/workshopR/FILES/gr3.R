# Kald plot med type = "n" for at initialisere grafikrammen uden at vise data
plot(cars,
     main = "Biler: Hastighed versus bremsel?ngde (1920)",
     xlab = "Hastighed (miles i timen)",
     ylab = "Bremsel?ngde (fod)",
     type = "n")
# Kald grid-funktionen for at tegne gitteret
grid()
# Kald lavniveau-grafikfunktioner som points og lines for at tilf?je grafikken ovenp? gitteret
points(cars)
