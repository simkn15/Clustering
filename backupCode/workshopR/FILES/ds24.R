# Antag nu, at vi har konverteret en matrix eller en liste til en dataramme
# R giver navne til søjlerne, men navnene er i bedste fald ikke særlig informative og i værste fald bizarre

# Søjlerne i en dataramme skal have navne
# Hvis vi konverterer nedenstående matrix til en dataramme vil R kalde søjlerne V1, V2 og V3
mat <- matrix(c(-0.818, -0.819, 0.385, -2.155, -0.667, -0.946, 1.531, -0.535, -0.494, -0.205, -0.611, -0.316), 4, 3)
as.data.frame(mat)
# Hvis matricen havde søjlenavnet defineret ville R have brugt disse navne i stedet

# Derimod giver det nogle meget mærkværdige navne, hvis man konverterer en liste til en dataramme
mylist <- list(c(-0.284, -1.114, -1.097, -0.873), c(-1.673, 0.929, 0.306, 0.778), c(0.323, 0.368, 0.067, -0.080))
as.data.frame(mylist)
# Igen, hvis listeelementerne havde navne ville R have brugt dem
# Heldigvis kan vi overskrive disse navne med vores egne navne
dfrm <- as.data.frame(mylist)
colnames(dfrm) <- c("before", "treatment", "after")
dfrm
