###############################################
### identify & remove ASVs associated with soil
### for ancient teeth microbe study
### with Ashley Brennaman
### Lou LaMartina, started Feb 3, 2022
###############################################


# set working directory
setwd("~/Desktop/Lab/Projects/Others/Brennaman/Final/Amplicon")


# packages
library(ggplot2)
library(vegan)
library(reshape2)
library(ape)
library(ggrepel)
library(indicspecies)
library(dada2)


# load ASV counts
counts <- read.csv("./RData/Brennaman_all_ASV_counts.csv")
rownames(counts) <- counts$Sample_ID
counts <- counts[-1]


# load taxonomy
tax <- read.csv("./RData/Brennaman_all_ASV_taxonomy.csv")


# convert to relative abundance
relabun <- counts / rowSums(counts)


# load sample info
info <- read.csv("./RData/Brennaman_sample_info.csv")



####################
### indicspecies ###
####################

# A = 1.0 -> ASVs found only samples in that group
# B = 1.0 -> found in all samples belonging to that group

# suspect that soil ASVs will be found in teeth samples,
# but teeth ASVs will not be found in soil samples

# soil A = probability that the ASV belongs to "soil"
#          ASV not found in any teeth samples
# soil B = probability of finding ASV in soil samples
#          ASV found in all soil samples

# teeth A = probability that the ASV belongs to "teeth"
#           ASV not found in any soil samples
# teeth B = probability of finding ASV in teeth samples
#           ASV found in all teeth samples

# therefore, teeth A is most important stat
# A = 1.0 -> ASVs found only in teeth samples


# stat
indic.stat <- multipatt(relabun, as.vector(info$Sample_type), control = how(nperm = 999))


# extract data
indic_all.df <- cbind(indic.stat$A[,1:2], indic.stat$B[,1:2], indic.stat$sign[,4:5])
colnames(indic_all.df) <- c("Soil_A", "Teeth_A", "Soil_B", "Teeth_B", "stat", "p.value")
indic_all.df$ASV <- rownames(indic_all.df)
indic_all.df <- merge(indic_all.df, tax, by = "ASV")


# how can you have a p.value of NA? --
# ASVs occur in all samples belonging to all groups. 
# The association with the set of all sites cannot be statistically tested, 
# because there is no external group for comparison.
# calling them "ubiquitous"
ubiq.df <- indic_all.df[is.na(indic_all.df$p.value),]


# remove those
indic.df <- indic_all.df[is.na(indic_all.df$p.value) == FALSE,]


# check i got them all
identical(sort(c(ubiq.df$ASV, indic.df$ASV)), sort(indic_all.df$ASV))


# subset samples
teeth_smp_counts <- counts[! rownames(counts) %in% c("MCPFC_S3","MCPFC_S4","MCPFC_S5","MCPFC_S6"),]
ncol(teeth_smp_counts[rowSums(teeth_smp_counts) > 0, colSums(teeth_smp_counts) > 0]) # [1] 21113

soil_smp_counts <- counts[rownames(counts) %in% c("MCPFC_S3","MCPFC_S4","MCPFC_S5","MCPFC_S6"),]
ncol(soil_smp_counts[rowSums(soil_smp_counts) > 0, colSums(soil_smp_counts) > 0]) # [1] 2909



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
teeth_totalProp <- sum(counts[colnames(counts) %in% teeth.df$ASV]) / sum(counts)


# proportion teeth samples
teeth_inTeeth <- sum(teeth_smp_counts[ colnames(teeth_smp_counts) %in% teeth.df$ASV]) / sum(teeth_smp_counts)


# proportion soil samples (should be zero)
teeth_inSoil <- sum(soil_smp_counts[ colnames(soil_smp_counts) %in% teeth.df$ASV]) / sum(soil_smp_counts)



########################
### non-teeth indicators

# teeth A < 1.0 -> ASVs found in teeth and soil, or just soil
nonTeeth.df <- subset(indic.df, Teeth_A < 1.0)


# proportion entire dataset
nonTeeth_totalProp <- sum(counts[colnames(counts) %in% nonTeeth.df$ASV]) / sum(counts)


# proportion teeth samples
nonTeeth_inTeeth <- sum(teeth_smp_counts[ colnames(teeth_smp_counts) %in% nonTeeth.df$ASV]) / sum(teeth_smp_counts)


#  proportion soil samples
nonTeeth_inSoil <- sum(soil_smp_counts[ colnames(soil_smp_counts) %in% nonTeeth.df$ASV]) / sum(soil_smp_counts)



###################
### ubiquitous asvs

# proportion entire dataset
ubiq_totalProp <- sum(counts[colnames(counts) %in% ubiq.df$ASV]) / sum(counts)


# proportion teeth samples
ubiq_inTeeth <- sum(teeth_smp_counts[colnames(teeth_smp_counts) %in% ubiq.df$ASV]) / sum(teeth_smp_counts)


# what proportion of uniquitous asvs make up soil samples?
ubiq_inSoil <- sum(soil_smp_counts[colnames(soil_smp_counts) %in% ubiq.df$ASV]) / sum(soil_smp_counts)


# check i got them all
identical(sort(c(ubiq.df$ASV, nonTeeth.df$ASV, teeth.df$ASV)), colnames(counts))
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

# keep teeth-associated ASVs
teeth_counts <- teeth_smp_counts[colnames(teeth_smp_counts) %in% teeth.df$ASV]
teeth_counts <- teeth_counts[rowSums(teeth_counts) > 0, colSums(teeth_counts) > 0]
teeth_relabun <- teeth_counts / rowSums(teeth_counts)

teeth_tax <- subset(tax, ASV %in% teeth.df$ASV)

teeth_relabun.t <- data.frame(t(teeth_relabun))
identical(teeth_tax$ASV, rownames(teeth_relabun.t))
teeth_data <- cbind(teeth_tax[c(9,1,3:8)], teeth_relabun.t)

identical(teeth_tax$ASV, colnames(teeth_counts))
teeth_counts.m <- as.matrix(teeth_counts)
colnames(teeth_counts.m) <- teeth_tax$FASTA


# FASTA
uniquesToFasta(teeth_counts.m,
               ids = paste0(teeth_tax$ASV, "__",
                            teeth_tax$Kingdom, "__",
                            teeth_tax$Phylum, "__",
                            teeth_tax$Class, "__",
                            teeth_tax$Order, "__",
                            teeth_tax$Family, "__",
                            teeth_tax$Genus),
               fout = "./RData/Brennaman_teeth_ASVs.fasta")


# write
write.csv(data.frame(Sample = rownames(teeth_counts), teeth_counts), 
          "./Final/Brennaman_teeth_ASV_counts.csv", row.names = F, na = "")

write.csv(teeth_tax, "./Final/Brennaman_teeth_ASV_taxonomy.csv", row.names = F, na = "")
