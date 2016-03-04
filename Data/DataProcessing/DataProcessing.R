#§§§§§§§§§§§§§§§§§§§§§§§§§§
# F2 Nematodes
# Data Processing
# Quentin Schorpp
# 06.08.2015
#§§§§§§§§§§§§§§§§§§§§§§§§§


# Load and slice Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Community Data (Family level)

# nema is a data set with.... 
# ... family abundance counted in subset of 100 Individuals
# ... Feedingtypes and Indices calculated with Excel (?) 
nema <- read.delim("Data/spe.txt") # Not really species data
nema$samcam <- as.factor(nema$samcam)
nema$field.ID <- as.factor(nema$field.ID)

# Extract family data
fam.org <- nema[, which(names(nema)=="Aphelenchidae"):which(names(nema)=="Pratylenchidae")]

# Load environmental data from separate file
env.org <- read.delim("Data/env.txt")
env.org$samcam <- as.factor(env.org$samcam)
env.org$field.ID <- as.factor(env.org$field.ID)
# cultivation period length
env.org$cpl <- strptime(env.org$date, "%d.%m.%Y") - strptime(paste("01.01.",env.org$y_est, sep=""), "%d.%m.%Y")
env1 <- droplevels(env.org[16:45,])

# Total count Data from separate file
counts.org <- read.delim("Data/counts.txt")
counts.org$counts <- counts.org$counts.av*10 # upscale to 100 ml Extract from 1oog dry soil => [Ind/100g soil]
counts.env <- env.org[rep(seq_len(nrow(env.org)), each=3),-1] # adjust environmental data

# Averaged over replicates for field.ID (since there are at least 3 counts per sample)
counts3SC <- aggregate(. ~ field.ID + samcam, counts.org, mean)
counts <- counts3SC[16:45,]
  
#counts %>% group_by(samcam, field.ID) %>% summarise_each(funs(mean))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Nematode counts Notes ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# counts is the average number of three counts,
# Nematode aliquots of 10 ml were taken from 100 ml Extract and counted
# Three times per subsample and three subsamples per composite sample

# A - multiply by 10 to get the number of Individuals per Extract.
# B - Calculatehow many gramm dry soil, does an m² with 10 cm depth have?
# C - multiply by the coefficint of B

# B - soil per m²
# bulk density ~ 1.4 g/cm³
# Volume <- 100*100*10
# soil.weight <- 1.4*Volume
# 
# counts$countm2 <- counts$counts.av*10*soil.weight
# str(counts)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




