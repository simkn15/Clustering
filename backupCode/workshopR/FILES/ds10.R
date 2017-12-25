# Antag nu, at vi gerne vil fjerne et element fra en liste

# Vi har fra en tidligere fil en liste af fire navngivne heltal
years <- list(Kennedy = 1960, Johnson = 1964, Carter = 1976, Clinton = 1994)
years
# Vi kan fjerne et element ved at tildele v?rdien NULL til det - her fjernes alts? Johnson
years[["Johnson"]] <- NULL
years
# Vi kan ogs? fjerne flere elementer p? denne m?de
years[c("Carter", "Clinton")] <- NULL
years
