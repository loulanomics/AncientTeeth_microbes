###################################################
### organize shotgun metagenomic data, derived from 
### metaphlan3, for ancient teeth microbe study
### with Ashley Brennaman
### Lou LaMartina, start Mar 16 2022
###################################################


# set working directory
setwd("~/Desktop/Lab/Projects/Others/Brennaman/Metagenome")


# packages
library(indicspecies)
library(ggplot2)
library(RColorBrewer)


# load results of metaphlan (taxonomy & relative abundance from metagenomic data)
# (done by Sterling)
metaphlan <- read.csv("./RData/metaphlan_results.csv")


####################
### extract data ###
####################

# assign OTU IDs
tax <- data.frame(OTU = paste0("OTU", sprintf("%0.3d", 1:nrow(metaphlan))),
                  Clade = metaphlan$clade_name)


# convert clade from metaphlan to data frame
tax.ls <- strsplit(tax$Clade, "\\|")
names(tax.ls) <- tax$OTU

for (i in tax$OTU) {
  
  # designate # of data frame columns
  length(tax.ls[[i]]) <- 7
}


# put in data frame
tax <- data.frame(tax, do.call(rbind, tax.ls)) 


# define taxonomy classifications / assignments
classes <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(tax)[3:9] <- classes


# remove identifiers - a letter and __ at beginning of names,
# e.g., k__Bacteria, p__Firmicutes ...
for (i in classes) {
  tax[[i]] <- gsub("^[[:alpha:]]__", "", tax[[i]])
}


# species also have genus names before an underscore, remove them so it's similar to dada2 layout
tax$Species <- sapply(strsplit(tax$Species, "_"), function(x) paste0(tail(x, -1), collapse = "_"))
tax[tax == ""] <- NA


# same thing, remove "_unclassified" from genus names, denotes unknown species i think
tax$Genus <- gsub("_unclassified", "", tax$Genus)


# note, metaphlan has their own way of calculating relative abundance,
# using genome size (so there are no raw "counts")
relabun <- data.frame(metaphlan[-c(1:2)])
rownames(relabun) <- tax$OTU
relabun <- data.frame(t(relabun))


# load sample metadata
info <- read.csv("../Amplicon/Final/ancientTeeth_sample_info.csv")


# change metagenome file names to sample names
rownames(relabun) <- paste0(c(rep("Soil", 2), rep("MCPFC_", 14)), 
                                sapply(strsplit(rownames(relabun), "_"), '[', 2))


# subset info to samples that were shotgun sequenced
info <- info[info$Sample_ID %in% rownames(relabun),]




####################
### indicspecies ###
####################
# see 03_ancientTeeth_removeSoil.R for ASV version & reasoning

# stat
indic.stat <- multipatt(relabun, c(rep("Soil", 2), rep("Teeth", 14)), control = how(nperm = 999))


# extract info
indic_all.df <- cbind(indic.stat$A[,1:2], indic.stat$B[,1:2], indic.stat$sign[,4:5])
colnames(indic_all.df) <- c("Soil_A", "Teeth_A", "Soil_B", "Teeth_B", "stat", "p.value")
indic_all.df$OTU <- rownames(indic_all.df)
indic_all.df <- merge(indic_all.df, tax, by = "OTU")


# how can you have a p.value of NA? --
# ASVs occur in all samples belonging to all groups. 
# The association with the set of all sites cannot be statistically tested, 
# because there is no external group for comparison.
# calling them "ubiquitous"
ubiq.df <- indic_all.df[is.na(indic_all.df$p.value),]


# remove those
indic.df <- indic_all.df[is.na(indic_all.df$p.value) == FALSE,]


# check i got them all
identical(sort(c(ubiq.df$OTU, indic.df$OTU)), sort(indic_all.df$OTU))


# subset samples
teeth_relabun <- relabun[! rownames(relabun) %in% c("Soil1", "Soil2"),]
teeth_relabun <- teeth_relabun[, colSums(teeth_relabun) > 0]

soil_relabun <- relabun[rownames(relabun) %in% c("Soil1", "Soil2"),]
soil_relabun <- soil_relabun[, colSums(soil_relabun) > 0]



####################
### teeth indicators

# A = 1.0 ASVs with a low p value have high teeth B values
# this means ASVs are significantly teeth-associated when:
#    - they are not found in any soil samples
#    - they are found in most teeth samples
# i'm still going to ignore this because microbe diversity can
# change person to person


# teeth A = 1.0 -> ASVs found only in teeth samples / not found in soil
teeth.df <- subset(indic.df, Teeth_A == 1.0)


# proportion entire dataset
teeth_totalProp <- sum(relabun[colnames(relabun) %in% teeth.df$OTU]) / sum(relabun)


# proportion teeth samples
teeth_inTeeth <- sum(teeth_relabun[ colnames(teeth_relabun) %in% teeth.df$OTU]) / sum(teeth_relabun)


# proportion soil samples (should be zero)
teeth_inSoil <- sum(soil_relabun[ colnames(soil_relabun) %in% teeth.df$OTU]) / sum(soil_relabun)



########################
### non-teeth indicators

# teeth A < 1.0 -> ASVs found in teeth and soil, or just soil
nonTeeth.df <- subset(indic.df, Teeth_A < 1.0)


# proportion entire dataset
nonTeeth_totalProp <- sum(relabun[colnames(relabun) %in% nonTeeth.df$OTU]) / sum(relabun)


# proportion teeth samples
nonTeeth_inTeeth <- sum(teeth_relabun[ colnames(teeth_relabun) %in% nonTeeth.df$OTU]) / sum(teeth_relabun)


#  proportion soil samples
nonTeeth_inSoil <- sum(soil_relabun[ colnames(soil_relabun) %in% nonTeeth.df$OTU]) / sum(soil_relabun)



###################
### ubiquitous asvs

# proportion entire dataset
ubiq_totalProp <- sum(relabun[colnames(relabun) %in% ubiq.df$OTU]) / sum(relabun)


# proportion teeth samples
ubiq_inTeeth <- sum(teeth_relabun[colnames(teeth_relabun) %in% ubiq.df$OTU]) / sum(teeth_relabun)


# what proportion of uniquitous asvs make up soil samples?
ubiq_inSoil <- sum(soil_relabun[colnames(soil_relabun) %in% ubiq.df$OTU]) / sum(soil_relabun)


# check i got them all
identical(sort(c(ubiq.df$OTU, nonTeeth.df$OTU, teeth.df$OTU)), colnames(relabun))
teeth_totalProp + nonTeeth_totalProp + ubiq_totalProp
teeth_inTeeth + nonTeeth_inTeeth + ubiq_inTeeth
teeth_inSoil + nonTeeth_inSoil + ubiq_inSoil


# make data frame
props.df <- data.frame(measure = c("teeth_totalProp", "nonTeeth_totalProp", 
                                   "ubiq_totalProp", "nonTeeth_inSoil", 
                                   "teeth_inTeeth", "nonTeeth_inTeeth",
                                   "teeth_inSoil", "ubiq_inSoil", "ubiq_inTeeth"),
                       indic = c("teeth", "nonTeeth", 
                                 "ubiquitous", "nonTeeth", 
                                 "teeth", "nonTeeth", 
                                 "teeth", "ubiquitous", "ubiquitous"),
                       data = c("total", "total", 
                                "total", "soil", 
                                "teeth", "teeth",
                                "soil", "soil", "teeth"),
                       prop = c(teeth_totalProp, nonTeeth_totalProp, 
                                ubiq_totalProp, nonTeeth_inSoil, 
                                teeth_inTeeth, nonTeeth_inTeeth,
                                teeth_inSoil, ubiq_inSoil, ubiq_inTeeth),
                       y = c(0.75, 0.97, 0.25, 0.75,
                             0.75, 0.97, NA, 0.25, 0.25))

# plot
ggplot(props.df, aes(x = data, y = prop, fill = indic)) +
  geom_bar(stat = "identity") +
  geom_text(aes(x = data, y = y, label = round(prop, 2))) +
  scale_fill_brewer(palette = 1) +
  theme_minimal() +
  labs(y = "proportion of reads", x = "sample set", fill = "indicator groups")




############
### save ###
############

teeth_relabun <- data.frame(Sample_ID = rownames(teeth_relabun), teeth_relabun)
write.csv(teeth_relabun, "./RData/ancientTeeth_OTU_relAbund_noSoil.csv", na = "", row.names = F)

write.csv(indic_all.df, "./RData/ancientTeeth_OTU_indicspecies_results.csv", na = "", row.names = F)

teeth_tax <- tax[tax$OTU %in% colnames(teeth_relabun),]
write.csv(teeth_tax, "./RData/ancientTeeth_OTU_taxonomy_noSoil.csv", na = "", row.names = F)

relabun <- data.frame(Sample_ID = rownames(relabun), relabun)
write.csv(relabun, "./RData/ancientTeeth_OTU_relAbund_wSoil.csv", na = "", row.names = F)
write.csv(tax, "./RData/ancientTeeth_OTU_taxonomy_wSoil.csv", na = "", row.names = F)
