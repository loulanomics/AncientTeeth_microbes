########################################
### processing 16S rRNA V4 amplicon data
### for ancient teeth microbe study
### with Ashley Brennaman
### Lou LaMartina, started Feb 2, 2022
########################################


# set working directory
setwd("~/Desktop/Lab/Projects/Others/Brennaman/Final/Amplicon")


# packages
library(dada2)
library(phyloseq)
library(decontam)
library(ggplot2)




###################################
### prepare data for processing ###
###################################

# set file paths
path <- "./fastq"
pathF <- "./fastq/fastqF"
pathR <- "./fastq/fastqR"


# set paths for filtered reads
filtered_pathF <- "./fastq/fastqF/Filtered"
filtered_pathR <- "./fastq/fastqR/Filtered"


# sort forward and reverse reads into file paths
fastqFs <- sort(list.files(pathF, pattern = "_R1_", full.names = TRUE))
fastqRs <- sort(list.files(pathR, pattern = "_R2_", full.names = TRUE))


# extract file names
Sample_names <- gsub("_L001_R1_001.fastq.gz", "", gsub("trim_sort_", "", basename(fastqFs)))




################################
### Inspect sequence quality ###
################################

# visualize quality of reads
qualityF <- plotQualityProfile(fastqFs[1:4])
qualityR <- plotQualityProfile(fastqRs[1:4])


# save quality profiles
ggsave("./Plots/qualityF.pdf", plot = qualityF, device = "pdf", width = 12, height = 8, units = "in")
ggsave("./Plots/qualityR.pdf", plot = qualityR, device = "pdf", width = 12, height = 8, units = "in")


# check: if there are not the same number of F and R files, stop.
if(length(fastqFs) != length(fastqRs)) 
  stop("Forward and reverse files do not match.")


# inspect sequence distribution - this takes a few mins
fastq_lengths1.ls <- list()
for(i in 1:length(fastqFs)){
  fastq_lengths1.ls[[i]] <- data.frame(table(nchar(getSequences(fastqFs[i]))))
}

fastq_lengths1 <- do.call(rbind, fastq_lengths1.ls)
colnames(fastq_lengths1) <- c("SeqLength", "Frequency")
fastq_lengths1$SeqLength <- as.numeric(as.character(fastq_lengths1$SeqLength))


# visualize
dist1 <-
  ggplot(fastq_lengths1, aes(x = SeqLength, y = Frequency)) +
  geom_bar(stat = "identity", width = 0.5, fill = "black")
dist1

ggsave("./Plots/seq_distribution_beforeFilt.pdf", plot = dist1, device = "pdf",
       width = 10, height = 4, units = "in")




#######################
### Quality control ###
#######################

# give filtered files new names and paths
filteredFs <- file.path(filtered_pathF, paste0(Sample_names, "_F_filt.fastq.gz"))
filteredRs <- file.path(filtered_pathR, paste0(Sample_names, "_R_filt.fastq.gz"))


# filter based on quality and read length
start = Sys.time()
filtered_out <- filterAndTrim(fastqFs, filteredFs, fastqRs, filteredRs, matchIDs = TRUE,
                              maxEE = 2, maxN = 0, truncQ = 10,
                              rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = TRUE)
Sys.time() - start
# 20.03538 mins


# inspect how many reads were filtered out of each sample
filtered_out
(1 - (filtered_out[,2] / filtered_out[,1])) * 100
mean((1 - (filtered_out[,2] / filtered_out[,1])) * 100)
# [1] 16.85279 % removed


# plot quality profiles of filtered reads
filtF.plot <- plotQualityProfile(filteredFs[1:4]); filtF.plot
filtR.plot <- plotQualityProfile(filteredRs[1:4]); filtR.plot


# save quality profiles
ggsave("./Plots/filt_qualityF.pdf", plot = filtF.plot, device = "pdf", width = 12, height = 8, units = "in")
ggsave("./Plots/filt_qualityR.pdf", plot = filtR.plot, device = "pdf", width = 12, height = 8, units = "in")


# set sample names to the ID only
names(filteredFs) <- Sample_names
names(filteredRs) <- Sample_names


# inspect sequence distribution
fastq_lengths2.ls <- list()
for(i in 1:length(filteredFs)){
  fastq_lengths2.ls[[i]] <- data.frame(table(nchar(getSequences(filteredFs[i]))))
}

fastq_lengths2 <- do.call(rbind, fastq_lengths2.ls)
colnames(fastq_lengths2) <- c("SeqLength", "Frequency")
fastq_lengths2$SeqLength <- as.numeric(as.character(fastq_lengths2$SeqLength))


# visualize
dist2 <-
  ggplot(fastq_lengths2, aes(x = SeqLength, y = Frequency)) +
  geom_bar(stat = "identity", width = 0.5, fill = "black")
dist2

ggsave("./Plots/seq_distribution_afterFilt.pdf", plot = dist2, device = "pdf", 
       width = 10, height = 4, units = "in")




############################
### Learning error rates ###
############################

# learn and visualize error rates of F reads
start = Sys.time()
errorF <- learnErrors(filteredFs, multithread = TRUE)
errF.plot <- plotErrors(errorF, nominalQ = TRUE, err_in = TRUE, err_out = TRUE); errF.plot


# learn and visualize error rates of R reads
errorR <- learnErrors(filteredRs, multithread = TRUE)
errR.plot <- plotErrors(errorR, nominalQ = TRUE, err_in = TRUE, err_out = TRUE); errR.plot
Sys.time() - start
# 1.004478 hours


# save error plots
ggsave("./Plots/errorF.pdf", plot = errF.plot, device = "pdf", width = 12, height = 8, units = "in")
ggsave("./Plots/errorR.pdf", plot = errR.plot, device = "pdf", width = 12, height = 8, units = "in")




################################
### Merging paired-end reads ###
################################

# create list of merged reads
Mergers <- vector("list", length(Sample_names))
names(Mergers) <- Sample_names


# sample inference and merging paired-end reads
start = Sys.time()
for(i in Sample_names) {
  cat("\nProcessing:", i, "(", match(i, Sample_names), ") :", format(Sys.time(), "%H:%M %p"), "\n")
  derepF <- derepFastq(filteredFs[[i]])
  dadaF <- dada(derepF, err = errorF, multithread = TRUE)
  derepR <- derepFastq(filteredRs[[i]])
  dadaR <- dada(derepR, err = errorR, multithread = TRUE)
  Merger <- mergePairs(dadaF, derepF, dadaR, derepR)
  Mergers[[i]] <- Merger
}
Sys.time() - start
# 1.006945 hours


# removing dereps to save memory
rm(derepF, derepR)


# contruct a sequence table
Sequence_table <- makeSequenceTable(Mergers)


# dimensions: num rows x num columns
dim(Sequence_table)
# 94 29593




##################################
### Quality control: processed ###
##################################


########
### trim

# inspect sequence distribution
seq_distribution <- data.frame(table(nchar(getSequences(Sequence_table)))) 
colnames(seq_distribution) <- c("SeqLength", "Frequency")
seq_distribution$SeqLength <- as.numeric(as.character(seq_distribution$SeqLength))


# visualize and save
dist3 <- 
  ggplot(seq_distribution, aes(x = SeqLength, y = Frequency)) +
  geom_bar(stat = "identity", width = 0.5, fill = "black")
dist3

ggsave("./Plots/seq_distribution_afterMerge.pdf", plot = dist3, device = "pdf", width = 10, height = 4, units = "in")


# remove reads of non target length, 5% above and below the median 
median(nchar(getSequences(Sequence_table)))
# [1] 253
min_len <- floor(median(nchar(getSequences(Sequence_table))) * 0.95); min_len
# [1] 240
max_len <- ceiling(median(nchar(getSequences(Sequence_table))) * 1.05); max_len
# [1] 266


# modify sequence table with new guidelines
Seq_table_trimmed <- Sequence_table[ ,nchar(colnames(Sequence_table)) 
                                     %in% seq(min_len, max_len)]


dim(Seq_table_trimmed)
# 94 22897



###################
### Remove chimeras

# removing chimeras with denovo screening
start = Sys.time()
Sequence_table_nochim <- removeBimeraDenovo(Seq_table_trimmed, method = "consensus",
                                            multithread = TRUE, verbose = TRUE)
Sys.time() - start
#  1.916502 mins


# how many unique sequences were moved?
ncol(Sequence_table)               # [1] 29593
ncol(Sequence_table_nochim)        # [1] 22017


# what percentage of reads were identified as chimeras?
(1 - sum(Sequence_table_nochim) / sum(Sequence_table)) * 100
# 27.80531




#######################
### Assign taxonomy ###
#######################

Sequence_taxa_table <- assignTaxonomy(Sequence_table_nochim,
                                      "~/Desktop/Lab/Projects/Misc/silva_nr_v132_train_set.fa.gz",
                                      multithread = TRUE)


# identify non-bacterial/archael ASVs
taxa.df <- data.frame(Sequence_taxa_table)
taxa.df <- data.frame(FASTA = rownames(taxa.df), taxa.df)
nonspecifics <- rbind(subset(taxa.df, Kingdom == "Eukaryota"),
                      subset(taxa.df, Order == "Chloroplast"),
                      subset(taxa.df, Family == "Mitochondria"))
nonspecifics$Source <- "not bacteria"




######################################
### identify negative control ASVs ###
######################################

# create info
counts_nochim.df <- data.frame(Sequence_table_nochim)
rownames(counts_nochim.df) <- gsub("-", "_", rownames(counts_nochim.df))
info <- data.frame(Sample_name = rownames(counts_nochim.df))
rownames(info) <- info$Sample_name


# create phyloseq object
NoMock_object <- phyloseq(otu_table(as.matrix(counts_nochim.df), taxa_are_rows = FALSE),
                          tax_table(as.matrix(taxa.df)),
                          sample_data(info))


# add logical variable to sample info stating if it's negative control
sample_data(NoMock_object)$NTC <- FALSE
sample_data(NoMock_object)$NTC[sample_data(NoMock_object)$Sample_name %in% 
                                 c("MCPFC_Blank_S92", "MCPFC_ExtC_S86")] <- TRUE
sample_data(NoMock_object)


# run check for contaminants using prevalence method and negative control
Contam_prev_NTC <- isContaminant(NoMock_object, method = "prevalence", neg = "NTC")


# how many are contaminants?
table(Contam_prev_NTC$contaminant)
# FALSE  TRUE 
# 21994    23 


# what taxa are they?
Contam_prev_NTC$FASTA <- rownames(Contam_prev_NTC)
NTC_contaminants <- subset(merge(Contam_prev_NTC[6:7], 
                                 taxa.df, by = "FASTA"), contaminant == TRUE)[-2]
NTC_contaminants$Source <- "NTC"




####################################
### identify mock community ASVs ###
####################################

# save as fasta
uniquesToFasta(Sequence_table_nochim,
               ids = paste0("NOCHIM", sprintf("%06d", 1:ncol(Sequence_table_nochim))), 
               fout = "./RData/Brennaman_wMock.fasta")


# in terminal:
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# $ makeblastdb -dbtype nucl -in zymo_mock_all.fasta -input_type fasta
#
# $ blastn -db zymo_mock_all.fasta -query Brennaman_wMock.fasta -task blastn -perc_identity 100 -outfmt 6 -out mock_align_nochim.txt
#
# $ echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n$(cat mock_align_nochim.txt)" > mock_align_nochim.txt
#
# $ sed 's/\t/,/g' mock_align_nochim.txt > mock_align_nochim.csv
#
# $ rm zymo_mock_all.fasta.*
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# load alignment results
blastn_nochim <- read.csv("./RData/mock_align_nochim.csv")


# keep only those full-length
blastn_nochim <- subset(blastn_nochim, length >= min_len)
length(unique(blastn_nochim$qseqid))
# [1] 9


# create vector of fastas associated with SEQ IDs
fastas_nochim.df <- data.frame(FASTA = colnames(Sequence_table_nochim))
fastas_nochim.df$qseqid <- paste0("NOCHIM", sprintf("%06d", 1:ncol(Sequence_table_nochim)))


# add FASTA seqs
blastn_nochim <- merge(blastn_nochim, fastas_nochim.df, by = "qseqid")


# keep the mock ASVs
Mock_ASVs <- subset(taxa.df, FASTA %in% blastn_nochim$FASTA)
Mock_ASVs$Source <- "Mock"




##############################
### contamination analysis ###
##############################

# combine
contaminants <- rbind(nonspecifics, NTC_contaminants, Mock_ASVs)
length(unique(contaminants$FASTA))
# [1] 349


# calculate their prevalence in the dataset
Proportion_dataset <- (colSums(Sequence_table_nochim[, colnames(Sequence_table_nochim)
                                                     %in% contaminants$FASTA]) / sum(Sequence_table_nochim))
Proportion_dataset <- data.frame(Proportion_dataset, FASTA = names(Proportion_dataset))


# add that
contaminants <- merge(Proportion_dataset, contaminants, by = "FASTA")



#########################
### organize and save ###
#########################

# remove contaminants
counts.df <- counts_nochim.df[! colnames(counts_nochim.df) %in% contaminants$FASTA]
taxa.df <- subset(taxa.df, ! FASTA %in% contaminants$FASTA)
identical(rownames(taxa.df), colnames(counts.df))


# simplify ASV names
taxa.df <- data.frame(ASV = paste0("ASV", sprintf("%05d", 1:ncol(counts.df))), taxa.df)


# create FASTA
uniquesToFasta(as.matrix(counts.df),
               ids = paste0(taxa.df$ASV, "__",
                            taxa.df$Kingdom, "__",
                            taxa.df$Phylum, "__",
                            taxa.df$Class, "__",
                            taxa.df$Order, "__",
                            taxa.df$Family, "__",
                            taxa.df$Genus),
               fout = "./RData/Brennaman_v4_wSoil.fasta")


# change ASVs in counts
identical(rownames(taxa.df), colnames(counts.df))
colnames(counts.df) <- taxa.df$ASV


# remove mock & blanks samples, empty ASVs
counts.df <- counts.df[! rownames(counts.df) %in% c("MCPFC_Blank_S92", "MCPFC_ExtC_S86","MCPFC_PosC_S91"),]
counts.df <- counts.df[rowSums(counts.df) > 0, colSums(counts.df) > 0]


# these should have been removed in contamination analysis?
empties <- subset(taxa.df, ! ASV %in% colnames(counts.df))[-1]
empties <- data.frame(FASTA = empties$FASTA,
                      Proportion_dataset = (colSums(Sequence_table_nochim[, colnames(Sequence_table_nochim)
                                                                         %in% empties$FASTA]) / sum(Sequence_table_nochim)),
                      empties[-1], Source = "Unidentified / empty")
contaminants <- rbind(contaminants, empties)


# also change sample names
info <- read.csv("./RData/Brennaman_sample_info.csv")
info$File_name <- gsub("-", "_", info$File_name)
identical(info$File_name, rownames(counts.df))
rownames(counts.df) <- info$Sample_ID
counts.df <- data.frame(Sample_ID = rownames(counts.df), counts.df)


# save
write.csv(contaminants, "./RData/Brennaman_contaminats.csv", na = "", row.names = F)
write.csv(counts.df, "./RData/Brennaman_all_ASV_counts.csv", row.names = F, na = "")
write.csv(taxa.df[c(1,3:8,2)], "./RData/Brennaman_all_ASV_taxaonomy.csv", row.names = F, na = "")



