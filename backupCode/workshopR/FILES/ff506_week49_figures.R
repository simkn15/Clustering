# Eksempel 1
data <- read.table("amphetamine.txt", header = TRUE)
plot(data$Dose,
     data$FoodConsumption,
     xlab = "Amfetamindosis (mg/kg)",
     ylab = "Fortæring af mad (g/kg)",
     col = "dark red",
     pch = 16)

# Eksempel 2
data <- read.table("rice.txt", header = TRUE)
plot(data$StrawSi,
     data$RiceAs,
     xlab = "Silicium i strå (g/kg tørstof)",
     ylab = expression(paste("Arsen i poleret ris (", mu, "g/kg tørstof)")),
     col = "dark red",
     pch = 16)

# Eksempel 3
data <- read.table("snakes.txt", header = TRUE)
plot(data$Length,
     data$Weight,
     xlab = "Længde (cm)",
     ylab = "Vægt (g)",
     col = "dark red",
     pch = 16)

result <- lm(data$Weight ~ data$Length)
abline(result)

summary(result)

cor(data$Weight, data$Length)

pre <- predict(result)
segments(data$Length, data$Weight, data$Length, pre, col="dark red")

library(calibrate)

textxy(data$Length, data$Weight, signif(residuals(result), 5))