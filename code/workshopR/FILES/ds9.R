# Antag nu, at vi ?nsker at lave en liste, der sammenk?der navne og v?rdier

# Funktionen list muligg?r at navngive listeelementer og skabe en sammenh?ng mellem hvert navn og den tilh?rende v?rdi
mylist <- list(far.left = 0.023, left = 0.159, mid = 0.500, right = 0.841, far.right = 0.977)
mylist
# Alternativ syntaks
mylist <- list()
mylist$far.left <- 0.023
mylist$left <- 0.159
mylist$mid <- 0.500
mylist$right <- 0.841
mylist$far.right <- 0.977

# Nogle gange har man to parallelle vektorer: en vektor med navne og en vektor med tilh?rende v?rdier
names <- c("far.left", "left", "mid", "right", "far.right")
values <- pnorm(-2:2)
# Man kan sammenk?de navnene og v?rdierne ved at lave en tom liste og udfylde den med vektoriseret tildeling
mylist <- list()
mylist[names] <- values
mylist
?pnorm
