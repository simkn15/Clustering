# Antag nu, at vi har parrede observationer i to vektorer, x og y, og en parallel faktor f, som indikerer deres grupper
# Vi vil gerne fremstille et scatterplot af x og y, som skelner mellem deres grupper

# Det indbyggede iris-datas?t indeholder parrede m?l for kronbladl?ngde og -bredde og
# hver m?ling har endvidere en art-egenskab, som indikerer arten af den m?lte blomst!

# At plotte flere grupper i ?t scatterplot skaber ikke-informativt rod, hvis ikke vi kan skelne grupperne fra hinanden
with(iris, plot(Petal.Length, Petal.Width, xlab = "Kronbladl?ngde", ylab = "Kronbladbredde"))
# Vi kan bruge pch-argumentet i plot-funktionen og dermed plotte hvert punkt med et tegn svarende til dens gruppe
with(iris, plot(Petal.Length, Petal.Width, xlab = "Kronbladl?ngde", ylab = "Kronbladbredde", pch = as.integer(Species)))
