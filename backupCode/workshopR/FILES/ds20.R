# Antag nu, at vi ønsker at tilføje en eller flere nye rækker til en dataramme

# Datasættet suburbs indlæses
suburbs <- read.csv("suburbs.csv")
suburbs
# Vi laver en ny dataramme med én række med nye data
newRow <- data.frame(city = "West Dundee", county = "Kane", state = "IL", pop = 5428)
# Herefter kan rbind-funktionen bruges til at tilføje denne dataramme med én række til vores eksisterende dataramme
suburbs <- rbind(suburbs, newRow)
# Havde vi i stedet ønsket at tilføje en ny søjle, så kunne vi have brugt cbind
# Bemærk: Det er vigtigt at den nye række bruger samme søjlenavne som den eksisterende dataramme!

# Det havde også været muligt at tilføje flere rækker, da rbind tillader flere argumenter
suburbs <- rbind(suburbs,
                 data.frame(city = "West Dundee", county = "Kane", state = "IL", pop = 5428),
                 data.frame(city = "East Dundee", county = "Kane", state = "IL", pop = 2955))