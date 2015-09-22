#################
# F2_Eearthworms
# Contrast matrix for Post-Hoc Multicomparisons
# Given two way interaction terms
# Create a contrast matrix to perform pairwisecomparisons only within groups
# Quentin Schorpp
# 18.05.2015
#################


##############################################################

# Define Parameters
k <- 2 # pairwise comparison
n <- 5 # Nr of factor levels for pairwise comparisons, e.g. levels of Interaction factor 1
groups <- 2 # Nr of groups (e.g. levels of interaction-factor 2) to draw comparisons within

# create interaction factor 
# interaction.factor <- interaction(factor1,factor2)

# create empty contrast matrix
x <- groups*choose(n,k) # row number for matrix, i.e. all comparisons that should be drawn
y <-  length(levels(fety2$agsam))+2 # column number, e.g. all levels of the interaction factor
cm1 <- matrix(0,x,y) # empty contrast matrix with two mor columns for pairwise factor combinations

# Fill in pairwise factor combinations, [first create interaction factor!]
# Use as ID of all combinations of interest
for(i in 1:2) {
  ID <- rbind(t(combn(levels(fety2$agsam)[(n-(n-1)):n],k)),
             t(combn(levels(fety2$agsam)[(2*n-(n-1)):(2*n)],k)))
  cm1[,i] <- ID[,i]
}
cm1 <- data.frame(cm1, stringsAsFactors=FALSE) # turn into fety2 frame
colnames(cm1)[(1+2):(2*n+2)] <- levels(fety2$agsam) # fill in column names

# write 1 or -1 if colnames match the names of ID-columns 1 or 2
for (i in 3:12) {
  for (j in 1:20) {
    if(colnames(cm1)[i]==cm1$X1[j]) {cm1[j,i] = -1}
    if(colnames(cm1)[i]==cm1$X2[j]) {cm1[j,i] = 1}
  }
}


rownames(cm1) <- paste(cm1$X1, cm1$X2, sep=" - ") # create rownames
cm1= cm1[,-c(1,2)] # delete ID-columns 1 and 2
cm1 <- as.matrix(cm1)
class(cm1) <- "numeric"

rm(k,n,groups,i,j,x,y,ID)

