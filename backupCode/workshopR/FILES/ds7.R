# Antag nu, at vi ?nsker at udv?lge listeelementer efter position

# Vi har en liste med 4 elementer
years <- list(1960, 1964, 1976, 1994)
years
# Vi kan tilg? enkeltelementer ved at bruge dobbelte firkantede parenteser
years[[1]]
# Vi kan udtr?kke dellister med enkelte firkantede parenteser
years[c(1, 2)]
# Bem?rk at der er en vigtig forskel p? years[[1]] og years[1], da f?rstn?vnte er et element og sidstn?vnte er en liste
class(years[[1]])
class(years[1])
