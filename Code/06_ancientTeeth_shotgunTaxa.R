####################################################
### teeth-associated microbe assessment in shotgun
### metagenomic data for ancient teeth microbe study
### with Ashley Brennaman
### Lou LaMartina, start Mar 16 2022
####################################################


# set working directory
setwd("~/Desktop/Lab/Projects/Others/Brennaman/Metagenome")


# packages
library(reshape2)
library(ggplot2)
library(RColorBrewer)


#################
### load data ###
#################
# see 05_ancientTeeth_shotgunDataPrep.R

# OTU relative abundances
relabun <- read.csv("./RData/ancientTeeth_OTU_relAbund_noSoil.csv")
rownames(relabun) <- relabun$Sample_ID
relabun <- relabun[-1]


# taxonomy
tax <- read.csv("./RData/ancientTeeth_OTU_taxonomy_noSoil.csv")


# sample metadata
info <- read.csv("../Amplicon/Final/ancientTeeth_sample_info.csv")


# subset info to samples that were shotgun sequenced
info <- info[info$Sample_ID %in% rownames(relabun),]




################
### bar plot ###
################

# get top 7 taxa of each cemetary area
toptax <- c()
for (i in unique(info$Cem_area)) {
  toptax[[i]] <- names(sort(colSums(relabun[rownames(relabun) %in%
                                                    info$Sample_ID[info$Cem_area == i],]),
                            decreasing = T)[1:7])
}
toptax <- unique(unlist(toptax))


# add them to tax table, call everyone else "other"
tax$Top <- NA
tax$Top[tax$OTU %in% toptax] <- tax$OTU[tax$OTU %in% toptax]
tax$Top[is.na(tax$Top)] <- "Other"


# sort by abundance
tax <- tax[match(c(toptax, tax$OTU[! tax$OTU %in% toptax]), tax$OTU),]


# add numbers for ordering
tax$Top <- paste(c(sprintf("%0.2d", 1:length(toptax)),
                         rep(length(toptax) + 1, nrow(tax) - length(toptax))), tax$Top)


# add info & taxa
top_relabun <- data.frame(Sample_ID = rownames(relabun), relabun)
top_relabun <- merge(top_relabun, info[c("Sample_ID", "Cem_area", "Lot_number")], by = "Sample_ID")
top.m <- melt(top_relabun, id.vars = c("Sample_ID", "Cem_area", "Lot_number"),
                  value.name = "Relabun", variable.name = "OTU")
top.m <- merge(top.m, tax[c("OTU", "Top")], by = "OTU")


# glom to top taxa
top.m <- aggregate(Relabun ~ Sample_ID + Cem_area + Top + Lot_number, sum, data = top.m)


# labels
taxlabels <- c(tax$Top[tax$OTU %in% toptax], "15 Other")
names(taxlabels) <- c(tax$Clade[tax$OTU %in% toptax], "Other")
newtax <- c(
  "k__Bacteria" =
    "Kingdom: Bacteria",
  "k__Bacteria|p__Firmicutes" =
    "Phylum: Firmicutes",
  "k__Bacteria|p__Firmicutes|c__Clostridia" =
    "Class: Clostridia",
  "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales" =
    "Phylum: Firmicutes\nOrder: Clostridiales",
  "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales_Family_XIII_Incertae_Sedis" =
    "Phylum: Firmicutes\nOrder: Clostridiales\nFamily: VIII incertae sedis",
  "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales_Family_XIII_Incertae_Sedis|g__Clostridiales_Family_XIII_Incertae_Sedis_unclassified" =
    "Phylum: Firmicutes\nOrder: Clostridiales\nFamily: VIII incertae sedis/unclassified",
  "k__Archaea" =
    "Kingdom: Archaea",
  "k__Archaea|p__Euryarchaeota" = 
    "Phylum: Euryarchaeota",
  "k__Archaea|p__Euryarchaeota|c__Methanobacteria" =
    "Phylum: Euryarchaeota\nClass: Methanobacteria",
  "k__Bacteria|p__Actinobacteria" =
    "Phylum: Actinobacteria",
  "k__Bacteria|p__Actinobacteria|c__Actinobacteria" =
    "Phylum: Actinobacteria\nClass: Actinobacteria",
  "k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales" =
    "Phylum: Actinobacteria\nOrder: Actinomycetales",
  "k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae" =
    "Phylum: Actinobacteria\nFamily: Actinomycetaceae",
  "k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces" =
    "Phylum: Actinobacteria\nGenus: Actinomyces",
  "Other" = "Other")
names(taxlabels) <- newtax

top.m$idlab <- sapply(strsplit(as.character(top.m$Sample_ID), "_"), '[', 2)
top.m$idlab <- paste(sprintf("%0.2d", as.numeric(top.m$idlab)), " -")


# plot
bar.plot <-
  ggplot(top.m, aes(x = idlab, y = Relabun, fill = Top)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(labels = names(taxlabels),
    values = c(brewer.pal(11, "Paired"),
               "gold3", "darkslategray1", "darkslategray4", "grey90")) +
  scale_x_discrete(limits = rev(sort(unique(top.m$idlab)))) +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 5, color = "black", margin = unit(c(0,-0.59,0,0), "cm")),
        axis.title.y = element_text(size = 8, color = "black", face = "bold"),
        axis.title.x = element_text(size = 8, color = "black", face = "bold"),
        legend.title = element_text(size = 8, color = "black", face = "bold"),
        legend.text = element_text(size = 6, color = "black", face = "italic"),
        legend.box.margin = margin(c(-20, 0, -20, -20)),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_blank()) +
  guides(fill = guide_legend(keyheight = 1.75, keywidth = 0.5, units = "in", ncol = 1)) +
  labs(y = "Taxa proportion", x = "Sample ID", fill = "Lowest\ntaxonomic\nassignment")
bar.plot

ggsave("./Plots/barplot.pdf", plot = bar.plot, device = "pdf", width = 6, height = 7, units = "in")
