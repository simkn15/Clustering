# Antag nu, at vi ?nsker at inds?tte et eller flere elementer i en vektor

# Funktionen append inds?tter data i en vektor og after-parameteret bruges til at angive inds?ttelsesstedet for nye elementer
# Hvis after-parameteret udelades, s? inds?ttes nye elementer til sidst (samme resultat som i forrige fil)

# Her inds?tter 99 i midten (efter 5)
append(1:10, 99, after = 5)
# Hvis man skriver after = 0 inds?ttes nye elementer forrest
append(1:10, 99, after = 0)
