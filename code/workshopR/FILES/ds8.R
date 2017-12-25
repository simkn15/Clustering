# Antag nu, at vi ?nsker at udv?lge listeelementer efter navn

# Hvert element i en liste kan have et navn - og hvis det er navngivet kan elementet udv?lges efter dets navn

# Denne tildeling laver en liste af fire navngivne heltal
years <- list(Kennedy = 1960, Johnson = 1964, Carter = 1976, Clinton = 1994)
# De f?lgende to udtryk returnerer samme v?rdi (elementet som er navngivet "Kennedy")
years[["Kennedy"]]
years$Kennedy
# De f?lgende to udtryk returnerer dellister af years-listen
years[c("Kennedy", "Johnson")]
years["Carter"]
