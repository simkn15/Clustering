# Antag nu, at vi ?nsker at generere tilf?ldige tal

# Det simple tilf?lde, hvor et uniformt tilf?ldigt tal mellem 0 og 1 ?nskes genereret, h?ndteres af runif-funktionen
runif(1)
# Argumentet til runif-funktionen er antallet af tilf?ldige v?rdier, som ?nskes genereret
# Vi kan s?ledes let generere en vektor med 10 s?danne v?rdier
runif(10)

# R kan dog ogs? generere tilf?ldige variate fra andre fordelinger
# For en given fordeling er navnet p? tilf?ldighedsgeneratoren "r" efterfulgt af fordelingens forkortede navn (f.eks. rnorm)
# Dette eksempel genererer en tilf?ldig v?rdi fra standardnormalfordelingen
rnorm(1)

# Der er tilf?ldighedsgeneratorer for alle indbyggede fordelinger

# En uniform variat mellem -3 og 3
runif(1, min = -3, max = 3)
# En standardnormal variat
rnorm(1)
# En normal variat med middelv?rdi 100 og standardafvigelse 15
rnorm(1, mean = 100, sd = 15)
# En binomial variat
rbinom(1, size = 10, prob = 0.5)
# En Poisson-variat
rpois(1, lambda = 10)
# En eksponentialvariat
rexp(1, rate = 0.1)
# En gammavariat
rgamma(1, shape = 2, rate = 0.1)

# Som med runif er det f?rste argument antallet af tilf?ldige variate, som ?nskes genereret
# Efterf?lgende argumenter er fordelingens parametre, f.eks.
# - middelv?rdi (mean) og standardafvigelse (sd) for normalfordelingen
# - antal fors?g (size) og successandsynlighed i hvert fors?g (prob) for binomialfordelingen

# Ovenst?ende eksempler bruger simple skalarer som parametre, men parametre kan ogs? v?re vektorer
# I dette tilf?lde vil R bladre igennem vektoren, n?r den genererer tilf?ldige v?rdier
rnorm(3, mean = c(-10, 0, +10), sd = 1)
# Her genereres tre normale tilf?ldige v?rdier fra fordelingerne med middelv?rdier p? henholdsvis -10, 0 og +10
# Dette er smart, hvis man eksempelvis arbejder med hierarkiske modeller, hvor selve parametrene er tilf?ldige
# I det f?lgende eksempel beregnes 100 udtr?k af en normal variat, hvis middelv?rdi selv er tilf?ldigt fordelt
# Hyperparametre: mu = 0 og sigma = 0,2
means <- rnorm(100, mean = 0, sd = 0.2)
rnorm(100, mean = means, sd = 1)
