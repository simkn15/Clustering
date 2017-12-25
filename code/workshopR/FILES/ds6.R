# Antag nu, at vi gerne vil lave en liste

# Lister kan være meget simple, som denne liste af tre tal
mylist <- list(0.5, 0.841, 0.977)
mylist
# Når R printer listen, så identificerer den hvert element med positionen ([[1]], [[2]], [[3]]) og printer værdien herunder

# Lister kan, i modsætning til vektorer, indeholde elementer af forskellige typer
mylist <- list(3.14, "Buh", c(1, 1, 2, 3), mean)
mylist

# Samme liste kan opbygges ved at lave en tom liste og udfylde den
mylist <- list()
mylist[[1]] <- 3.14
mylist[[2]] <- "Buh"
mylist[[3]] <- c(1, 1, 2, 3)
mylist[[4]] <- mean

# Listeelementer kan navngives
mylist <- list(mid = 0.5, right = 0.841, far.right = 0.977)
mylist

