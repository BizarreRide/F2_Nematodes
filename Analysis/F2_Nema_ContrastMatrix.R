#################
# F2_Eearthworms
# Contrast matrix for Post-Hoc Multicomparisons
# Given two way interaction terms
# Create a contrast matrix to perform pairwisecomparisons only within groups
# Quentin Schorpp
# 18.05.2015
#################


##############################################################

data <- indices
data$factor1 <- data$age_class
data$factor2 <- data$samcam
data$interaction1 <- with(data, interaction(factor1,factor2))
data$interaction2 <- with(data, interaction(factor2,factor1))

# Define Parameters
k <- 2 # pairwise comparison
n <- length(unique(data$factor1)) # Nr of factor levels of factor 1, e.g. groups, treatment,....
m <- length(unique(data$factor2)) # Nr Nr of factor levels of factor 2, e.g. time, repeated measurements,...


# create interaction factor 
# interaction.factor <- interaction(factor1,factor2)

# create empty contrast matrix
x <- m*choose(n,k)+n*choose(m,k)  # row number for matrix, i.e. all comparisons that should be drawn
y <-  length(levels(data$interaction1))+2 # column number = all levels of the interaction factor 
                                                         # +2 for index
cm1 <- matrix(0,x,y) # empty contrast matrix with two mor columns for pairwise factor combinations

# Fill in pairwise factor combinations, [first create interaction factor!]
# Use as ID of all combinations of interest

# Pairwise comparisons written out:

IDz1 <- cm1[1:(m*choose(n,k)),1:2]
z1 <- choose(n,k)

for (i in 1:m) {
  ID <- t(combn(levels(data$interaction1)[(i*n-(n-1)):(i*n)],k))
  nam <- paste("ID",i, sep=".")
  assign(nam,ID)
  IDz1[(1+i*z1-z1):(i*z1),] <- ID
}


IDz2 <- cm1[(m*z1+1):x,1:2]
z2 <- choose(m,k)

for (i in 1:n) {
ID <- t(combn(sort(levels(data$interaction1))[(i*m-(m-1)):(i*m)],k))
nam <- paste("ID",i, sep=".")
assign(nam,ID)
IDz2[(1+i*z2-z2):(i*z2),] <- ID
}


for(i in 1:2) {
  cm1[1:(m*choose(n,k)),i] <- IDz1[,i]
  cm1[(m*z1+1):x,i] <- IDz2[,i]
}

cm1 <- data.frame(cm1, stringsAsFactors=FALSE) # turn into fety2 frame
colnames(cm1)[3:y] <- levels(data$interaction1) # fill in column names

# write 1 or -1 if colnames match the names of ID-columns 1 or 2
for (i in 3:y) {
  for (j in 1:x) {
    if(colnames(cm1)[i]==cm1$X1[j]) {cm1[j,i] = -1}
    if(colnames(cm1)[i]==cm1$X2[j]) {cm1[j,i] = 1}
  }
}


rownames(cm1) <- paste(cm1$X1, cm1$X2, sep=" - ") # create rownames
cm1= cm1[,-c(1,2)] # delete ID-columns 1 and 2
cm1 <- as.matrix(cm1)
class(cm1) <- "numeric"

rm(k,n,m,x,y,z1,z2,ID, IDz1, IDz2)
