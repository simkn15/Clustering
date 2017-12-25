# Antag nu, at vi gerne vil beregne antallet af kombinationer af n elementer, udtaget k ad gangen

# Et almindeligt problem, når man skal beregne sandsynligheder for diskrete variable, er at tælle antallet af kombinationer:
# Antallet af forskellige delmængder af størrelse k, som kan dannes ud fra en mængde med n elementer
# Tallet er givet ved n! / r!(n - r)!, men det er meget mere bekvemt at bruge choose-funktionen (specielt for store n og k)

# Find antallet af måder, hvorpå 3 elementer kan udtages af en mængde med 5 elementer
choose(5, 3)
# Find antallet af måder, hvorpå 3 elementer kan udtages af en mængde med 50 elementer
choose(50, 3)
# Find antallet af måder, hvorpå 30 elementer kan udtages af en mængde med 50 elementer
choose(50, 30)

# Disse tal er også kendt som binomialkoefficienter