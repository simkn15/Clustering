# Antag nu, at vi ?nsker at lave en matrix og initialisere den ud fra givne v?rdier

# Her laver vi en matrix ved f?rst at samle data i en vektor og dern?st forme vektoren som en matrix med 2 r?kker og 3 s?jler
dat <- c(1.1, 1.2, 2.1, 2.2, 3.1, 3.2)
mat <- matrix(dat, 2, 3)
mat
# Dette kaldes en 2x3-matrix
# Det f?rste argument er data, det andet argument er antal r?kker og det tredje argument er antal s?jler
# Bem?rk: Matricen udfyldes s?jle efter s?jle, ikke r?kke efter r?kke

# Vi kan ogs? lave en s?kaldt 2x3-nulmatrix (en 2x3-matrix udfyldt kun med 0'er)
matrix(0, 2, 3)
# Da f?rste argument er en enkelt v?rdi bruger R blot genbrugsreglen og gentager v?rdien indtil matricen er udfyldt
# Vi kunne ogs? have udfyldt den med NA
matrix(NA, 2, 3)
# Man kan ogs? skrive det hele i ?n linje, men det bliver m?ske lidt sv?rere at l?se
mat <- matrix(c(1.1, 1.2, 1.3, 2.1, 2.2, 2.3), 2, 3)
# Nogle gange ville man m?ske v?lge at indtaste selve data i en rektangul?r form, der afsl?rer matricens struktur
dat <- c(1.1, 1.2, 1.3,
         2.1, 2.2, 2.3)
# Udtrykt p? denne m?de kan l?seren hurtigt identificere de to r?kker og de tre s?jler
# Ved at skrive byrow = TRUE udfyldes matricen r?kke efter r?kke (frem for s?jle efter s?jle)
mat <- matrix(dat, 2, 3, byrow = TRUE)
# Man kan eventuelt inds?tte data direkte i udtrykket
mat <- matrix(c(1.1, 1.2, 1.3,
                2.1, 2.2, 2.3),
              2, 3, byrow = TRUE)
# Man kan ogs? lave en vektor om til en matrix ved blot at tildele dimensioner til vektoren
dat <- c(1.1, 1.2, 1.3, 2.1, 2.2, 2.3)
dim(dat) <- c(2, 3)
dat
