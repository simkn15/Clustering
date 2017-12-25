# Antag nu, at vi ?nsker at generere tilf?ldige tal og tillige ?nsker at kunne reproducere samme tilf?ldige tal hver gang

# Vi kalder funktionen set.seed for at initialisere tilf?ldighedsgeneratoren til en kendt tilstand
set.seed(42)
# Generer 15 uniformt tilf?ldige tal
runif(15)
# Geninitialiser tilf?ldighedsgeneratoren til samme kendte tilstand
set.seed(42)
# Generer de samme 15 uniformt "tilf?ldige" tal
runif(15)

# P? denne m?de f?s de samme resultater fra k?rsel til k?rsel!