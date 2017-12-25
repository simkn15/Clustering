# Det indbyggede pressure-datasæt plottes som et scatterplot
plot(pressure, xlab = "Temperatur (grader Celsius)", ylab = "Tryk (mmHg)", main = "Trykdata: Kviksølvs damptryk")
# Vi kan i stedet tegne en linje
plot(pressure, xlab = "Temperatur (grader Celsius)", ylab = "Tryk (mmHg)", main = "Trykdata: Kviksølvs damptryk", type = "l")
