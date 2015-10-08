
#$$$$$$$$$$$$$$$$$$$$$$$$$
# F2 Nematodes
# Calculate Indices in R
# Quentin Schorpp
# 07.10.2015
#$$$$$$$$$$$$$$$$$$$$$$$$$

# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("data/RMAkeLikeFile.R")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# take the values of a data frame and multiply them by a coefficient in column B of the mastertable, if their colnames matches the name of column A in the mastertable

master <- read.csv("Data/F2_Nema_MastertableFam.csv", sep=";")
master$c.p <- as.factor(master$c.p)

# s,e and b weighings ####
# The following code aims to assign weightings (values 0.8, 3.2,5,..) to certain combination of feeding Type and c-p values
# i.e. Fu2  <-  0.8.
# However, not all these combinations get a weighting!
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# data.frames with weightings

bwt <- data.frame(
  FeedingType=as.factor(c(2,3)),
  c.p=as.factor(c(2,2)),
  bwt=c(0.8,0.8))

swt <- data.frame( 
  FeedingType=as.factor(c(2,2,5,5,8,8)),
  c.p=as.factor(c(3,4,4,5,4,5)),
  swt=c(1.8,3.2,3.2,5,3.2,5))

ewt <- data.frame(
  FeedingType=as.factor(c(2,3)),
  c.p=as.factor(c(2,1)),
  ewt=c(0.8,3.2))

cwt <- data.frame(
  FeedingType=as.factor(2),
  c.p=as.factor(2),
  cwt=0.8)


# assign weightings to master table
master <- master %>% 
  full_join(bwt, by=c("FeedingType","c.p"), nomatch="0") %>%
  full_join(swt, by=c("FeedingType","c.p"), nomatch="0") %>%
  full_join(ewt, by=c("FeedingType","c.p"), nomatch="0") %>%
  full_join(cwt, by=c("FeedingType","c.p"), nomatch="0") 

master[is.na(master)] <- 0

str(master)

# multiply abundances by their weights



# calculate  b, e and s , the sums of the weighted (relative?) abundances.


FaPro <- fam %>% 
  add_rownames("ID") %>% 
  gather(Family, Value, -ID) %>% #convert to long format
  inner_join(master, by = "Family") %>% 
  mutate(b=Value*bwt, s=Value*swt, e=Value*ewt, c=Value*cwt) %>%
  group_by(ID.x) %>% 
  summarise(b=sum(b), s = sum(s),  e=sum(e), c=sum(c)) %>%
  mutate(BI=100*b/(b+e+s), SI=100*(s/(s+b)), EI=(100*(e/(e+b))), CI=100*(c/e))
FaPro



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Maturity Index ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MI.sum <- fam %>% 
  add_rownames("ID") %>% 
  gather(Family, Value, -ID) %>% #convert to long format
  inner_join(master, by = "Family") %>% 
  filter(MI.PPI == "MI") %>% 
  group_by(ID.x) %>% 
  summarise(Total = sum(Value))


MI.matrix <- matrix(NA,30,25)

# c-p weighting of relative abundances
for (i in 1:25){
  for (j in 1:30){
    if (colnames(fam)[i] == master[j,2]) {
      MI.matrix[,i] <- fam[,i]*master[j,4]/MI.sum
    }
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Example from Thierry ####
# calculate summed abundances only for species that count for the MI
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(dplyr)
library(tidyr)
data(varespec, package = "vegan")
attributes <- data.frame(
  Species = c(colnames(varespec), "spec1", "spec2"),
  Attribute = c(rep(c("MI", "PI"), c(14, 30)), "MI", "PI")
)
varespec %>% 
  add_rownames("ID") %>% 
  gather(Species, Value, -ID) %>% #convert to long format
  inner_join(attributes, by = "Species") %>% 
  filter(Attribute == "MI") %>% 
  group_by(ID) %>% 
  summarise(Total = sum(Value))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








fam %>% 
  add_rownames("ID") %>% 
  gather(Family, Value, -ID) %>% #convert to long format
  full_join(master, by = "Family") %>% 
  group_by(ID.x) %>% 
  summarise(EI = sum(Value))

  summarise(EI = 100*(sum(ewt)/(sum(ewt)+sum(bwt))))



data(varespec, package = "vegan")
attributes <- data.frame(
  Species = c(colnames(varespec), "spec1", "spec2"),
  Attribute = c(rep(c("MI", "PI"), c(14, 30)), "MI", "PI"),
  index1=c(rep(c(1,2,3,4),each=11),3,3)
)

attributes$new[attributes$Attribute=="MI" & attributes$index1==c(2,3)] <- 0.8



MI.matrix <- matrix(NA,30,25)

for (i in 1:25){
  for (j in 1:30){
    if (colnames(fam)[i] == master[j,2]) {
      MI.matrix[,i] <- fam[,i]*master[j,4]/
    }
  }
}

for (i in 1:25){
  for (j in 1:30){
    if (colnames(fam)[i] == master[j,2]) {
      if(master[j,5] == "MI") {
        vector1 <- apply(fam,1,sum)
      }
    }
  }}

vector1 - apply(fam,1,sum)
sum(fam)


tfam <- data.frame(t(fam))

for (i in 1:25){
  for (j in 1:30){
 if (rownames(tfam)[i]==rownames(master)[j]){
      tfam$mat<- master[,5]
}}}




# Species data from vegan package:
data(varespec)

# create attributes table
attributes <- matrix(NA, length(varespec), 2)
attributes[,1] <- colnames(varespec)
attributes[,2] <- c(rep("MI",14),rep("PI",30))

# add species to the attribute table
x <- c("spec1","MI")
y <- c("spec2","PI")

attributes <- rbind(attributes, x, y)
row.names(attributes) <- c(1:46)

# calculate rowsums only for species contained in the attributes table 
# and having the entry "MI" in the attributes table
for (i in 1:44){
  for (j in 1:46){
    if ((colnames(fam)[i] == master[j,1]) & (master[j,2] == "MI")) {
      apply(varspec,1,sum)
    }
  }}


library(dplyr)
library(tidyr)

MI.sum <- fam %>% 
  add_rownames("ID") %>% 
  gather(Family, Value, -ID) %>% #convert to long format
  inner_join(master, by = "Family") %>% 
  filter(MI.PPI == "MI") %>% 
  group_by(ID.x) %>% 
  summarise(Total = sum(Value))

MI.sum$ID.x <- as.numeric(MI.sum$ID.x)
MI.sum <- MI.sum[order(MI.sum$ID.x, decreasing = FALSE),]

library(dplyr)
library(tidyr)
data(varespec, package = "vegan")
attributes <- data.frame(
  Species = c(colnames(varespec), "spec1", "spec2"),
  Attribute = c(rep(c("MI", "PI"), c(14, 30)), "MI", "PI")
)
varespec %>% 
  add_rownames("ID") %>% 
  gather(Species, Value, -ID) %>% #convert to long format
  inner_join(attributes, by = "Species") %>% 
  filter(Attribute == "MI") %>% 
  group_by(ID) %>% 
  summarise(Total = sum(Value))


gather(fam, Species,orks,-ID)

fam <- add_rownames(fam,"ID")

fam <- gather(fam, Species,orks,-ID)

filter(fam, Species=="Pratylenchidae")
