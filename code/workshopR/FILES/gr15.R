library(MASS)

# Antag nu, at vores datasæt indeholder en numerisk variabel og en faktor (kategorisk variabel)
# Vi vil gerne fremstille flere boksplot af den numeriske variabel, opdelt efter faktorniveauer

# UScereal-datasættet indeholder mange variable om morgensmadsprodukter

# En variabel er sukkermængde per portion og en anden er hyldeposition (hvor man tæller fra gulvet)
# Leverandører af morgenmadsprodukter kan forhandle hyldeposition med henblik på at øge salgspotentialet
# Hvor mon de placerer morgenmadsprodukter med højt sukkerindhold?
# Dette kan vi undersøge ved at fremstille et boksplot per hylde!
boxplot(sugars ~ shelf, data = UScereal)
# Vi tilføjer etiketter på akserne, så det resulterende boksplot bliver lettere at fortolke
boxplot(sugars ~ shelf, data = UScereal,
        main = "Sukkerindhold efter hylde",
        xlab = "Hylde", ylab = "Sukker (gram per portion)")

# Boksplottene viser, at hylde 2 har flest sukkerholdige morgenmadsprodukter
# Måske denne hylde er i øjenhøjde for små børn, som kan påvirke forældrenes valg?