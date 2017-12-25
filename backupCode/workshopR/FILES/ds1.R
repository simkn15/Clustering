# Antag nu, at vi ?nsker at tilf?je yderligere elementer til en vektor

v <- c(1, 2, 3)
# Tilf?j en enkelt v?rdi til v
v <- c(v, 4)
v

w <- c(5, 6, 7, 8)
# Tilf?j en hel vektor til v
v <- c(v, w)
v

# Lav en vektor med tre elementer
v <- c(1, 2, 3)
# Tildel en v?rdi til det tiende element
v[10] <- 10
# R udviser automatisk vektoren
v
