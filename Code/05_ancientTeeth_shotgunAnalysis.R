##################################################
### teeth-associated microbe assessment in shotgun
### metagenomic data, derived from metaphlan3,
### for ancient teeth microbe study
### with Ashley Brennaman
### Lou LaMartina, start Mar 16 2022
############################################

# notes:
#   "OTU" refers to clades assigned with metaphlan from metagenomic data.
#   like OTUs, each clade is unique and represents a group of many unique sequences.

# "ASV" refers to ampicon sequence variants assigned by dada2.
#   unlike OTUs, these are single, unique sequences.

# therefore:
# ASV : dada2 : amplicon sequencing :: OTU : metaphlan :: shotgun sequencing


# set working directory
setwd("~/Desktop/Lab/Projects/Others/Brennaman/Final/Metagenome")


# packages
library(indicspecies)
library(ggplot2)


# load results of metaphlan (taxonomy & relative abundance from metagenomic data)
# done by Sterling
metaphlan <- read.csv("./RData/metaphlan_results.csv")




################
### OTU data ###
################

############
### taxonomy

# assign OTU IDs
tax_otu <- data.frame(OTU = paste0("OTU", sprintf("%0.3d", 1:nrow(metaphlan))),
                  Clade = metaphlan$clade_name)


# convert clade from metaphlan to data frame
tax_otu.ls <- strsplit(tax_otu$Clade, "\\|")
names(tax_otu.ls) <- tax_otu$OTU

for (i in tax_otu$OTU) {
  
  # designate # of data frame columns
  length(tax_otu.ls[[i]]) <- 7
}


# put in data frame
tax_otu <- data.frame(tax_otu, do.call(rbind, tax_otu.ls)) 


# define taxonomy classifications / assignments
classes <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(tax_otu)[3:9] <- classes


# remove identifiers - a letter and __ at beginning of names,
# e.g., k__Bacteria, p__Firmicutes ...
for (i in classes) {
  tax_otu[[i]] <- gsub("^[[:alpha:]]__", "", tax_otu[[i]])
}


# species also have genus names before an underscore, remove them so it's similar to dada2 layout
tax_otu$Species <- sapply(strsplit(tax_otu$Species, "_"), function(x) paste0(tail(x, -1), collapse = "_"))
tax_otu[tax_otu == ""] <- NA


# same thing, remove "_unclassified" from genus names, denotes unknown species i think
tax_otu$Genus <- gsub("_unclassified", "", tax_otu$Genus)



##########
### counts

# note, metaphlan has their own way of calculating relative abundance,
# using genome size
relabun_otu <- data.frame(metaphlan[-c(1:2)])
rownames(relabun_otu) <- tax_otu$OTU
relabun_otu <- data.frame(t(relabun_otu))



################
### ASV data ###
################

# load ASV taxonomy
tax_asv <- read.csv("../Amplicon/Final/ancientTeeth_teeth_ASVs_taxonomy.csv")


# load ASV counts
counts_asv <- read.csv("../Amplicon/Final/ancientTeeth_teeth_ASVs_counts.csv")
rownames(counts_asv) <- counts_asv$Sample
counts_asv <- counts_asv[-1]


# convert to relative abundance
relabun_asv <- counts_asv / rowSums(counts_asv)




################
### metadata ###
################

# load sample metadata
info <- read.csv("../Amplicon/Final/ancientTeeth_sample_info.csv")


# change metagenome file names to sample names
rownames(relabun_otu) <- paste0(c(rep("Soil", 2), rep("MCPFC_", 14)), 
                                sapply(strsplit(rownames(relabun_otu), "_"), '[', 2))


# subset asv counts & info to samples that were shotgun sequenced
counts_asv <- counts_asv[rownames(counts_asv) %in% rownames(relabun_otu),]
relabun_asv <- relabun_asv[rownames(relabun_asv) %in% rownames(relabun_otu),]
info <- info[info$Sample_ID %in% rownames(relabun_otu),]



####################
### indicspecies ###
####################
# see 03_ancientTeeth_removeSoil.R for ASV version

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
indic.stat <- multipatt(relabun_otu, c(rep("Soil", 2), rep("Teeth", 14)), control = how(nperm = 999))


# extract info
indic_all.df <- cbind(indic.stat$A[,1:2], indic.stat$B[,1:2], indic.stat$sign[,4:5])
colnames(indic_all.df) <- c("Soil_A", "Teeth_A", "Soil_B", "Teeth_B", "stat", "p.value")
indic_all.df$OTU <- rownames(indic_all.df)
indic_all.df <- merge(indic_all.df, tax_otu, by = "OTU")


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
teeth_relabun_otu <- relabun_otu[! rownames(relabun_otu) %in% c("Soil1", "Soil2"),]
teeth_relabun_otu <- teeth_relabun_otu[, colSums(teeth_relabun_otu) > 0]

soil_relabun_otu <- relabun_otu[rownames(relabun_otu) %in% c("Soil1", "Soil2"),]
soil_relabun_otu <- soil_relabun_otu[, colSums(soil_relabun_otu) > 0]



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
teeth_totalProp <- sum(relabun_otu[colnames(relabun_otu) %in% teeth.df$OTU]) / sum(relabun_otu)


# proportion teeth samples
teeth_inTeeth <- sum(teeth_relabun_otu[ colnames(teeth_relabun_otu) %in% teeth.df$OTU]) / sum(teeth_relabun_otu)


# proportion soil samples (should be zero)
teeth_inSoil <- sum(soil_relabun_otu[ colnames(soil_relabun_otu) %in% teeth.df$OTU]) / sum(soil_relabun_otu)



########################
### non-teeth indicators

# teeth A < 1.0 -> ASVs found in teeth and soil, or just soil
nonTeeth.df <- subset(indic.df, Teeth_A < 1.0)


# proportion entire dataset
nonTeeth_totalProp <- sum(relabun_otu[colnames(relabun_otu) %in% nonTeeth.df$OTU]) / sum(relabun_otu)


# proportion teeth samples
nonTeeth_inTeeth <- sum(teeth_relabun_otu[ colnames(teeth_relabun_otu) %in% nonTeeth.df$OTU]) / sum(teeth_relabun_otu)


#  proportion soil samples
nonTeeth_inSoil <- sum(soil_relabun_otu[ colnames(soil_relabun_otu) %in% nonTeeth.df$OTU]) / sum(soil_relabun_otu)



###################
### ubiquitous asvs

# proportion entire dataset
ubiq_totalProp <- sum(relabun_otu[colnames(relabun_otu) %in% ubiq.df$OTU]) / sum(relabun_otu)


# proportion teeth samples
ubiq_inTeeth <- sum(teeth_relabun_otu[colnames(teeth_relabun_otu) %in% ubiq.df$OTU]) / sum(teeth_relabun_otu)


# what proportion of uniquitous asvs make up soil samples?
ubiq_inSoil <- sum(soil_relabun_otu[colnames(soil_relabun_otu) %in% ubiq.df$OTU]) / sum(soil_relabun_otu)


# check i got them all
identical(sort(c(ubiq.df$OTU, nonTeeth.df$OTU, teeth.df$OTU)), colnames(relabun_otu))
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




#########################
### OTU ~ ASV compare ###
#########################

# glom asvs to phylum
phylum_asv <- data.frame(ASV = colnames(relabun_asv), t(relabun_asv))
phylum_asv <- merge(tax_asv[1:3], phylum_asv, by = "ASV")
phylum_asv$Phylum <- paste0(phylum_asv$Kingdom, "__", phylum_asv$Phylum)
phylum_asv <- aggregate(. ~ Phylum, sum, data = phylum_asv[-c(1:2)])
rownames(phylum_asv) <- phylum_asv$Phylum
phylum_asv[phylum_asv == ""] <- "unclassified"
phylum_asv <- phylum_asv[-1]


# glom otus to phylum
phylum_otu <- data.frame(OTU = colnames(teeth_relabun_otu[-1]), t(teeth_relabun_otu[-1]))
phylum_otu <- merge(tax_otu[c(1,3,4)], phylum_otu, by = "OTU")
phylum_otu$Phylum <- paste0(phylum_otu$Kingdom, "__", phylum_otu$Phylum)
phylum_otu <- aggregate(. ~ Phylum, sum, data = phylum_otu[-c(1:2)])
rownames(phylum_otu) <- phylum_otu$Phylum
phylum_otu[phylum_otu == ""] <- "unclassified"
phylum_otu <- phylum_otu[-1]


# rank abundance of phyla in shotgun OTUs
toptax <- names(sort(rowSums(phylum_otu), decreasing = T))



############################
### OTU phylum sample counts

# sum & sort phyla
phylum_otu_sums <- data.frame(Phylum = toptax, Sum = sort(rowSums(phylum_otu), decreasing = T))


# add numbers for ordering
phylum_otu_sums$Top <- paste(sprintf("%0.2d", 1:length(toptax)), phylum_otu_sums$Phylum)


# add dataset variable
phylum_otu_sums$Data <- "Shotgun"



############################
### ASV phylum sample counts

# add them to tax table, call everyone else "other"
phylum_tax_asv <- unique(data.frame(Phylum = paste0(tax_asv$Kingdom, "__", tax_asv$Phylum), Top = NA))
phylum_tax_asv$Top[phylum_tax_asv$Phylum %in% toptax] <- 
  as.character(phylum_tax_asv$Phylum[phylum_tax_asv$Phylum %in% toptax])
phylum_tax_asv$Top[is.na(phylum_tax_asv$Top)] <- "Other"


# sort by abundance
toporder_asv <- phylum_tax_asv[phylum_tax_asv$Top %in% toptax,]
toporder_asv <- toporder_asv[match(toptax, toporder_asv$Top),]
phylum_tax_asv <- phylum_tax_asv[! phylum_tax_asv$Top %in% toptax,]
phylum_tax_asv <- rbind(toporder_asv, phylum_tax_asv)


# add numbers for ordering
phylum_tax_asv$Top <- paste(c(sprintf("%0.2d", 1:length(toptax)),
                              rep(length(toptax) + 1, nrow(phylum_tax_asv) - length(toptax))), phylum_tax_asv$Top)


# sum phyla
phylum_asv_sums <- data.frame(Phylum = names(sort(rowSums(phylum_asv), decreasing = T)),
                              Sum = sort(rowSums(phylum_asv), decreasing = T))


# add taxonomy names
phylum_asv_sums <- merge(phylum_tax_asv, phylum_asv_sums, by = "Phylum")


# add dataset variable
phylum_asv_sums$Data <- "Amplicon"
phyla <- rbind(phylum_otu_sums, phylum_asv_sums[c(1,3,2,4)])
phyla <- aggregate(Sum ~ Top + Data, sum, data = phyla)


# plot
data.plot <-
  ggplot(phyla, aes(x = Data, y = Sum, fill = Top)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c(brewer.pal(11, "Paired"), "grey90")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 8, color = "black", face = "bold"),
        axis.title.x = element_text(size = 8, color = "black", face = "bold"),
        legend.title = element_text(size = 8, color = "black", face = "bold"),
        legend.text = element_text(size = 8, color = "black", face = "italic"),
        panel.border = element_blank(),
        axis.line = element_blank()) +
  guides(fill = guide_legend(keyheight = 1.5, keywidth = 0.5, units = "in", ncol = 1)) +
  labs(y = "Taxa proportion", x = "Dataset", fill = "Kingdom__Phylum")
data.plot

ggsave("./Plots/phyla compare.pdf", plot = data.plot, device = "pdf", width = 6, height = 6, units = "in")



###############
### stats & plot

# spearman rho
cor.test(phylum_otu_sums$Sum, phylum_asv_sums$Sum[phylum_asv_sums$Top != "12 Other"], method = "spearman") 
# Spearman's rank correlation rho
# 
# data:  phylum_otu_sums$Sum and phylum_asv_sums$Sum[phylum_asv_sums$Top != "12 Other"]
# S = 166, p-value = 0.4682
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2454545


# new data frame for plot
ranks <- phyla

# make multiplier variable for second y axis
mult = max(phylum_otu_sums$Sum) / max(phylum_asv_sums$Sum)
ranks$Sum[ranks$Data == "Amplicon"] <- ranks$Sum[ranks$Data == "Amplicon"] * mult


# plot
rank.plot <-
  ggplot(ranks, aes(x = Top, y = Sum, fill = Data)) +
  geom_bar(stat="identity", position = position_dodge()) +
  scale_y_continuous(sec.axis = sec_axis(~ . / mult,
                                         name = "Amplicon sum relative abundance")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 8, color = "black", face = "bold"),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 8, color = "black", face = "bold"),
        legend.text = element_text(size = 8, color = "black"),
        panel.border = element_rect(size = .75, color = "grey80", fill = NA),
        axis.ticks = element_line(size = 0.25),
        axis.line = element_line(size = 0.25)) +
  labs(x = "Shotgun data phyla", y = "Shotgun sum relative abundance")
rank.plot

ggsave("./Plots/ranks.pdf", plot = rank.plot, device = "pdf", width = 6, height = 3, units = "in")





############
### save ###
############

# keep teeth-associated OTUs
teeth_relabun_otu <- teeth_relabun_otu[,colnames(teeth_relabun_otu) %in% teeth.df$OTU]
teeth_relabun_otu <- teeth_relabun_otu[rowSums(teeth_relabun_otu) > 0, colSums(teeth_relabun_otu) > 0]
teeth_relabun_otu <- data.frame(Sample_ID = rownames(teeth_relabun_otu), teeth_relabun_otu)
teeth_tax_otu <- subset(tax_otu, OTU %in% teeth.df$OTU)


# write
write.csv(teeth_relabun_otu, "./RData/ancientTeeth_teeth_OTU_relabun.csv", row.names = F, na = "")
write.csv(teeth_tax_otu, "./RData/ancientTeeth_teeth_OTU_taxonomy.csv", row.names = F, na = "")

