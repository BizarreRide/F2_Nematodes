#§§§§§§§§§§§§§§§§§§§§§§§§
# F2 Nematodes
# Family Data Processing
# Quentin Schorpp
# 07.10.2015
#$$$$$$$$$$$$$$$$$$$$$$$$$



# Nematode Families ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# fam is the number of Individuals according to one family, found in a subsample of ca. 100 Individuals
# The subsample is taken from all Individals in 100 ml Extract, extracted from 100g dry soil equivalent


# Family proportions
fam.rel <- round(data/rowSums(data),2)

# Family proportions upscaled to total bundances (counts)
fam.usc <- round(fam.rel*counts[15:45,"counts"],2)

## Selected families ####
fam.pa <- decostand(data, "pa")
# calculate sum per species
fam.sum <- apply(fam.pa,2,sum)
#sort(fam.sum)
# remove species that occur at less than 5 sites
fam.fin <- data[, !fam.sum<5]
rm(fam.sum)

biodiv <- biodiv.fun(round(fam.usc,0))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
