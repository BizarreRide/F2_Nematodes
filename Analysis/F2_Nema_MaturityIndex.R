

# take the values of a data frame and multiply them by a coefficient in column B of the mastertable, if their colnames matches the name of column A in the mastertable

master <- read.csv("Data/F2_Nema_MastertableFam.csv", sep=";")

MI.matrix <- matrix(NA,30,25)

for (i in 1:25){
  for (j in 1:30){
    if (colnames(fam)[i] == master[j,2]) {
      MI.matrix[,i] <- fam[,i]*master[j,4]/
    }
  }
}


if (colnames(fam)[i] == master[j,2] & master[j,5] == "MI") {
  print(fam)
}


