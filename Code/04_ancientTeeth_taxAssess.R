##################################################
### teeth-associated taxonomy composition analysis
### for ancient teeth microbe study
### with Ashley Brennaman
### Lou LaMartina, started Feb 3, 2022
##################################################


# set working directory
setwd("~/Desktop/Lab/Projects/Others/Brennaman/Amplicon")


# packages
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(vegan)


# load ASV counts
counts <- read.csv("./Final/ancientTeeth_teeth_ASV_counts.csv")
rownames(counts) <- counts$Sample
counts <- counts[-1]


# load taxonomy
tax <- read.csv("./Final/ancientTeeth_teeth_ASV_taxonomy.csv")

info <- read.csv("./Final/ancientTeeth_sample_info.csv")


# convert to relative abundance
relabun <- counts / rowSums(counts)




#############################
### lowest classification ###
#############################

# make sure columns are in this order:
tax <- tax[c("ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "FASTA")]


# create vector of classifications
classes <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")


# change blank classes to NA
tax[tax == ""] <- NA



##################################
### create list with ASVs that are 
###   unclassified at each level

# empty list
unknownASVs.ls <- list()

# for each taxonomic class/column in tax data frame,
for (i in classes){
  
  # in the "kingdom" column,
  if (i == "Kingdom") {
    
    # if it's blank, obtain its ASV names -> have list of all ASVs 
    # that are not classified to kingdom level.
    unknownASVs.ls[[i]] <- tax[is.na(tax[[i]]), "ASV"]
    
    # in the other columns,
  } else if (i != "Kingdom") {
    
    # if they are blank, obtain their ASV names, but don't include ones found previously.
    unknownASVs.ls[[i]] <- setdiff(tax[is.na(tax[[i]]), "ASV"], 
                                   tax[is.na(tax[[match(i, classes)]]), "ASV"])
  }
}



##############################
### make their new names their 
###   lowest classification

# empty list
classes.ls <- list()

# for each column in tax data frame,
for (i in classes) {
  
  # in kingdom column,
  if (i == "Kingdom") {
    
    # if they are blank, call them "unclassified."
    classes.ls[[i]] <- data.frame(ASV = subset(tax, ASV %in% unknownASVs.ls[[i]])$ASV,
                                  Name = "Unclassified")  
    
    # for all the others,
  } else if (i != "Kingdom") {
    
    # call them by their lowest level of classification = the column before it.
    classes.ls[[i]] <- data.frame(ASV = subset(tax, ASV %in% unknownASVs.ls[[i]])$ASV,
                                  Name = paste0(classes[match(i, classes) - 1], "__", 
                                                subset(tax, ASV %in% unknownASVs.ls[[i]])[[match(i, classes)]]))
  }
}



####################
### add to tax table
###   "Name" column

# combine into single data frame 
classes.df <- dplyr::bind_rows(classes.ls)


# for those classified down to genus, give them genus names
classes.df <- rbind(classes.df,
                    data.frame(ASV = tax[is.na(tax$Genus) == F, "ASV"],
                               Name = paste0("Genus__", tax[is.na(tax$Genus) == F, "Genus"])))


# remove characters R doesn't like
classes.df$Name <- gsub("-", "_", classes.df$Name)
classes.df$Name <- gsub("\\(", "_", classes.df$Name)
classes.df$Name <- gsub("\\)", "", classes.df$Name)


# add to taxa table
tax <- merge(classes.df, tax, by = "ASV")



############################
### counts of lowest classes

# get ASVs for each class
classASVs.ls <- list()

for (i in unique(tax$Name)) {
  classASVs.ls[[i]] <- as.character(tax$ASV[tax$Name == i])
}

length(unique(unique(tax$Name)))
# [1] 1470 unique taxa


# sum those ASVs in each sample
counts_class.ls <- list()

for( i in unique(tax$Name)) {
  counts_class.ls[[i]] <- rowSums(counts[colnames(counts) %in% unlist(classASVs.ls[i])])
}
  

# make new counts table - "Name" as column instead of ASVs
counts_class <- data.frame(do.call(cbind, counts_class.ls))


# convert to relative abundance
relabun_class <-  counts_class / rowSums(counts_class)




################
### bar plot ###
################

# get top taxa of each cemetary area
toptax <- c()
for (i in unique(info$Cem_area)) {
  toptax[[i]] <- names(sort(colSums(relabun_class[rownames(relabun_class) %in% 
                                   info$Sample_ID[info$Cem_area == i],]), 
             decreasing = T)[1:7])
}
toptax <- unique(unlist(toptax))


# add them to tax table, call everyone else "other"
tax_class <- unique(data.frame(Name = tax$Name, Top = NA))
tax_class$Top[tax_class$Name %in% toptax] <- as.character(tax_class$Name[tax_class$Name %in% toptax])
tax_class$Top[is.na(tax_class$Top)] <- "Other"


# sort by abundance
toporder <- tax_class[tax_class$Top %in% toptax,]
toporder <- toporder[match(toptax, toporder$Top),]
tax_class <- tax_class[! tax_class$Top %in% toptax,]
tax_class <- rbind(toporder, tax_class)


# add numbers for ordering
tax_class$Top <- paste(c(sprintf("%0.2d", 1:length(toptax)),
                                 rep(length(toptax) + 1, nrow(tax_class) - length(toptax))), tax_class$Top)


# add info & taxa
relabun_class_cem <- data.frame(Sample_ID = rownames(relabun_class), relabun_class)
relabun_class_cem <- merge(relabun_class_cem, info[c("Sample_ID", "Cem_area", "Lot_number")], by = "Sample_ID")
rel_cem.m <- melt(relabun_class_cem, id.vars = c("Sample_ID", "Cem_area", "Lot_number"),
                  value.name = "Relabun", variable.name = "Name")
rel_cem.m <- merge(rel_cem.m, tax_class, by = "Name")


# glom to top taxa
rel_cem.m <- aggregate(Relabun ~ Sample_ID + Cem_area + Top + Lot_number, sum, data = rel_cem.m)


# labels
taxlabels <- unique(tax[tax$Name %in% toptax, 2:8])
taxlabels <- merge(taxlabels, tax_class[1:15,], by = "Name")
taxlabels <- c(
"01 Genus__Actinomyces" = "Phylum: Actinobacteria
Genus: Actinomyces",
"02 Genus__Methanobrevibacter" = "Phylum: Euryarchaeota
Genus: Methanobrevibacter",
"03 Family__Family_XIII" = "Phylum: Firmicutes
Order: Clostridiales
Family XIII",
"04 Genus__Parvimonas" = "Phylum: Firmicutes
Genus: Parvimonas",
"05 Genus__Family_XIII_AD3011_group" = "Phylum: Firmicutes
Order: Clostridiales
Family XIII (AD3011)",
"06 Class__Subgroup_6" = "Phylum: Acidobacteria
Class: Subgroup 6",
"07 Genus__Pseudoramibacter" = "Phylum: Firmicutes
Genus: Pseudoramibacter",
"08 Genus__Anaerosalibacter" = "Phylum: Firmicutes
Genus: Anaerosalibacter",
"09 Genus__Paenibacillus" = "Phylum: Firmicutes
Genus: Paenibacillus",
"10 Genus__Massilia" = "Phylum: Proteobacteria
Genus: Massilia",
"11 Genus__Flexilinea" = "Phylum: Chloroflexi
Genus: Flexilinea",
"12 Genus__Actinophytocola" = "Phylum: Actinobacteria
Genus: Actinophytocola",
"13 Genus__Fretibacterium" = "Phylum: Synergistetes
Genus: Fretibacterium",
"14 Genus__Desulfobulbus" = "Phylum: Proteobacteria
Genus: Desulfobulbus",
"15 Other" = "Other")

rel_cem.m$idlab <- sapply(strsplit(as.character(rel_cem.m$Sample_ID), "_"), '[', 2)
rel_cem.m$idlab <- paste(sprintf("%0.2d", as.numeric(rel_cem.m$idlab)), " -")


# plot
bar.plot <-
  ggplot(rel_cem.m, aes(x = idlab, y = Relabun, fill = Top)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(labels = taxlabels, 
                    values = c(brewer.pal(11, "Paired"), 
                               "gold3", "darkslategray1", "darkslategray4", "grey90")) +
  scale_x_discrete(limits = rev(sort(unique(rel_cem.m$idlab)))) +
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
  guides(fill = guide_legend(keyheight = 1.5, keywidth = 0.5, units = "in", ncol = 1)) +
  labs(y = "Taxa proportion", x = "Sample ID", fill = "Lowest\ntaxonomic\nassignment")
bar.plot

#ggsave("./Plots/barplot.pdf", plot = bar.plot, device = "pdf", width = 6, height = 7, units = "in")




#######################################
### similarity of duplicate samples ###
#######################################

# duplicates = same person / human remains subject

# biologial duplicates = split the calculus samples, extracted DNA separately
bio_dups <- list("Lot_9243" = c("MCPFC_39", "MCPFC_81"),
                 "Lot_10299" = c("MCPFC_59", "MCPFC_82"),
                 "Lot_10537" = c("MCPFC_78", "MCPFC_83"))


# technical duplicates = single extract, twice PCR
tech_dups <- list("Lot_10091" = c("MCPFC_20", "MCPFC_84"),
                  "Lot_10371" = c("MCPFC_40", "MCPFC_85"),
                  "Lot_3073" = c("MCPFC_60", "MCPFC_86"),
                  "Lot_10561" = c("MCPFC_80", "MCPFC_87"))


# combine
dups <- rbind(data.frame(Lot_number = names(bio_dups), 
                         do.call(rbind, bio_dups), Dup = "Bio"),
              data.frame(Lot_number = names(tech_dups), 
                         do.call(rbind, tech_dups), Dup = "Tech"))


# ignore warning message
dups <- melt(dups, id.vars = c("Lot_number", "Dup"), value.name = "Sample_ID")[-3]
dups <- dups[order(dups$Sample_ID),]


# get relative abundances of those
relabun_dups <- relabun[rownames(relabun) %in% dups$Sample_ID,]


# analysis of similarity
dup.bray <- scores(vegdist(relabun_dups, method = "bray"))
identical(rownames(dup.bray), dups$Sample_ID)
dup.sim <- anosim(dup.bray, dups$Lot_number)
summary(dup.sim)
# Call:
#   anosim(x = dup.bray, grouping = dups$Lot_number) 
# Dissimilarity: user supplied square matrix 
# 
# ANOSIM statistic R:     1 
# Significance: 0.001 
# 
# Permutation: free
# Number of permutations: 999
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
#   0.255 0.371 0.449 0.538 
# 
# Dissimilarity ranks between and within classes:
#   0%   25%  50%   75% 100%  N
# Between    8 28.75 49.5 70.25   91 84
# Lot_10091  2  2.00  2.0  2.00    2  1
# Lot_10299  7  7.00  7.0  7.00    7  1
# Lot_10371  3  3.00  3.0  3.00    3  1
# Lot_10537  4  4.00  4.0  4.00    4  1
# Lot_10561  6  6.00  6.0  6.00    6  1
# Lot_3073   1  1.00  1.0  1.00    1  1
# Lot_9243   5  5.00  5.0  5.00    5  1

# R stat: ratio between dissimilarities between sites within a group
# and the dissimilarities between sites that are in different groups. 
# The closer this value is to 1, the more the sites within a group are 
# similar to each other and dissimilar to sites in other groups.
# The significance of the R-statistic is determined by permuting 
# the membership of sites in groups.

