###########################
# F2 Nematodes
# Data Processing
# Quentin Schorpp
# 06.08.2015
###########################

source("Analysis/Sources/RequiredPackages.R")
source("Data/DataProcessing/DataProcessing.R")    # Load Data; Slice nema, create fam; Upscale counts to 100g soil;
                                                  # Average count data, over three subsamples

env1 <- droplevels(env.org[16:45,])
source("Data/DataProcessing/EnvDataProcessing.R") # Slice Data (defined as env1) into groups (management, soil, etc)
                                                  # select environemtnal variables via PCA (orthognality) => env.fin
                                                  # create variable "intensity" from PCA


source("Data/DataProcessing/AverageData.R")       # Average Data over sampling campaigns (samcam), env.av, fam.av, counts.av
                                                  # Calculate differences for families between samcam => fam.slope

fam1 <- fam.av
source("Data/DataProcessing/FamDatProcessing.R")  # calculate relative and upscaled family abundances => fam.rel, fam.usc
                                                  # select abundant families => fam.fin

data <- fam.av                                       # calculate Nematode Indices; WITHOUT NCR
source("Data/DataProcessing/FaunalProfileIndices.R") 
source("Data/DataProcessing/MaturityIndices.R")
source("Data/DataProcessing/FeedingTypes.R")         # Feeding types coded by numbers 1(bcdef) - 8; 
                                                     # it is possible to change this coding!!!
source("Data/DataProcessing/cpValueAbundances.R")

