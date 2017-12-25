# Antag nu, at vi ønsker at udflade alle elementerne i en liste til en vektor

# Dette er en liste af tal
iq.scores <- list(89.73383, 116.5565, 113.0454)
# De fleste statistiske funktioner virker kun på vektorer, ikke på lister, så vi kan ikke direkte beregne middelværdien
mean(iq.scores)
# I stedet må vi udflade listen til en vektor ved hjælp af unlist-funktionen og beregne middelværdien af den resulterende vektor
mean(unlist(iq.scores))
