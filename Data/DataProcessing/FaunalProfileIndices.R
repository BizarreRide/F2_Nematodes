#§§§§§§§§§§§§§§§§§§§§§§§§
# F2 Nematodes
# Calculate Faunal Profile Indices
# Quentin Schorpp
# 07.10.2015
#$$$$$$$$$$$$$$$$$$$$$$$$$

# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#source("data/RMAkeLikeFile.R")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# take the values of a data frame and multiply them by a coefficient in column B of the mastertable, if their colnames matches the name of column A in the mastertable

master <- read.csv("Data/F2_Nema_MastertableFam.csv", sep=";")
master$c.p <- as.factor(master$c.p)

# Faunal Profile ####
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
# calculate BI, SI, EI and CI 


# specify data first!!!!!!!

FaPro <- data %>% 
  add_rownames("ID") %>% 
  gather(Family, Value, -ID) %>% #convert to long format
  inner_join(master, by = "Family") %>% 
  mutate(b=Value*bwt, s=Value*swt, e=Value*ewt, c=Value*cwt) %>%
  group_by(ID.x) %>% 
  summarise(b=sum(b), s = sum(s),  e=sum(e), c=sum(c)) %>%
  mutate(BI=100*b/(b+e+s), SI=100*(s/(s+b)), EI=(100*(e/(e+b))), CI=100*(c/e))
FaPro <- FaPro[order(as.numeric(FaPro$ID.x)),]

rm(ewt,bwt,swt,cwt)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%