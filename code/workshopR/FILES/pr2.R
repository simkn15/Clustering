# Antag nu, at vi ønsker at generere alle kombinationer af n elementer, udtaget k ad gangen

# Vi kan bruge combn(1:5, 3) til at generere alle kombinationer af tallene 1 til 5, udtaget 3 ad gangen
combn(1:5, 3)
# Funktionen er ikke begrænset til tal - vi kan også generere kombinationer af strenge
# Her genereres alle kombinationer af 5 behandlinger, udtaget 3 ad gangen
combn(c("B1", "B2", "B3", "B4", "B5"), 3)
# Bemærk: Når antallet af elementer vokser kan antallet af kombinationer eksplodere, specielt hvis k ikke er tæt på 1 eller n