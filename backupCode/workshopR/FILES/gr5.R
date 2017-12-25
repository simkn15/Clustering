# Vi bygger nu videre på eksemplet fra forrige fil
with(iris, plot(Petal.Length, Petal.Width, xlab = "Kronbladlængde", ylab = "Kronbladbredde", pch = as.integer(Species)))
# Vi tilføjer en signaturforklaring
legend("topleft", horiz = TRUE, legend = c("setosa", "versicolor", "virginica"), pch = 1:3, bty = "n")

# Vi kan også undgå at skrive arternes navne i R-koden, men i stedet bruge art-faktoren f
f <- factor(iris$Species)
with(iris, plot(Petal.Length, Petal.Width, xlab = "Kronbladlængde", ylab = "Kronbladbredde", pch = as.integer(f)))
legend("topleft", horiz = TRUE, legend = as.character(levels(f)), pch = 1:length(levels(f)), bty = "n")
# Dette har følgende fordele:
# 1. Man undgår at glemme et artsnavn eller skrive et artsnavn forkert
# 2. Hvis nogen tilføjer endnu en art, så tilpasser koden sig automatisk