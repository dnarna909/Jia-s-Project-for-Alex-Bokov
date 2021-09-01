rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 80000000)

## Import RNAseq information file 
library(xlsx)
cols <- read.xlsx2("E:/2018-10-30 Maggie data/RNAseq Data/updated.sampleinformation_6.xlsx", 1, header=TRUE)
cols <- cols[, c("ProbeIDs", "Order", "tissue", "treatment1", "treatment2","pre.post", "subject", "Sample_id", "Clamp_ID") ]
cols$number <- 1:length(rownames(cols))
cols$subject <- factor(cols$subject, levels = c("DM", "obese", "lean"))
cols$treatment1 <- factor(cols$treatment1, levels = c("eritoran", "placebo"))
cols$treatment2 <- factor(cols$treatment2, levels = c("lipid", "saline", "none"))
cols$pre.post <- factor(cols$pre.post, levels = c("post", "pre"))

## Import clamp data
clamp <- read.xlsx2("E:/2021-03-26 Bulk Fat RNAseq/Subjects with Clamp Infor and All Infor.xlsx", 1, header=TRUE)
clamp <- clamp[, -which(names(clamp) %in% c("fat", "muscle", "monocyte") )]
# change date
library("zoo")
clamp$DATE <- as.Date(as.numeric(clamp$DATE), origin="1899-12-30")

## merge clamp data and RNAseq information
cols_clamp <- merge(cols, clamp, 
                    by = sort(intersect(colnames(cols), colnames(clamp))) )
cols_clamp$subject <- factor(cols_clamp$subject, levels = c("DM", "obese", "lean"))
cols_clamp$treatment1 <- factor(cols_clamp$treatment1, levels = c("eritoran", "placebo"))
rownames(cols_clamp) <- cols_clamp$ProbeIDs

## Import counts from RNAseq
# Prepare original counts data for all the package analysis
#counts <- read.delim("E:/2018-10-30 Maggie data/RNAseq Data/allSamples_normData.txt", row.names = 1)
#hist(unlist(counts[1,]),breaks = 50)
#counts <- read.delim("E:/2018-10-30 Maggie data/RNAseq Data/Norm_logRPKM_allSamples_Genes_ReadCount.txt", row.names = 1)
#hist(unlist(counts[1,]),breaks = 50)
# normalized readcount
#counts <- read.delim("E:/2018-10-30 Maggie data/RNAseq Data/Norm_allSamples_Genes_ReadCount.txt", row.names = 1)
#hist(unlist(counts[1,]),breaks = 50)

# row readcount
txt <- read.delim("E:/2018-10-30 Maggie data/RNAseq Data/allSamples_Genes_ReadCount.txt", row.names = 1)
hist(unlist(txt[1,]),breaks = 50)
txt <- txt[-which(rownames(txt) %in% c("psiTPTE22", "tAKR", "no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique") ), ]
hist(unlist(txt[1,]),breaks = 50)
counts <- txt
total.reads <- colSums(counts)

# RPKM 
# RPKM is not comparable across different samples.
txt <- read.delim("E:/2018-10-30 Maggie data/RNAseq Data/RPKM_allSamples_Genes_ReadCount.txt", row.names = 1)
hist(unlist(txt[1,]),breaks = 50)
txt <- txt[-which(rownames(txt) %in% c("psiTPTE22", "tAKR", "no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique") ), ]
RPKM <- txt

# calculate gene length
per_million_scaling_factor <- total.reads/1000000
RPM <- data.frame(matrix(nrow = nrow(counts), ncol = ncol(counts)))
rownames(RPM) <- rownames(counts)
colnames(RPM) <- colnames(counts)
for (i in 1:ncol(RPM)) { RPM[, i] <- counts[, i]/per_million_scaling_factor[i] }

# calculate gene length from reads and RPKM, however, the gene length will be 0 if the gene reads is 0
genelength <- data.frame(matrix(nrow = nrow(counts), ncol = ncol(counts)))
rownames(genelength) <- rownames(counts)
colnames(genelength) <- colnames(counts)
for (i in 1:ncol(genelength)) { genelength[, i] <- RPM[, i]/RPKM[, i] }
genelength[is.na(genelength)] = 0
genelength2 <- as.data.frame((rowSums(genelength)/ncol(genelength))*1000)
colnames(genelength2) <- "calcualted.length"
head(genelength2)

#import gene length
geneLength <- readRDS("E:/Gene length/USCS_human/hg38.ncbiRefSeq.geneLength.rds")#19768
geneLength <- readRDS("E:/Gene length/USCS_human/hg38.refGene.geneLength.rds")#19866

geneLength$Symbol <- rownames(geneLength)
geneLength <- as.data.frame(geneLength[which(rownames(geneLength) %in% rownames(counts)), ])

df <- merge(genelength2, geneLength, by = 0, all.x = TRUE)



# calculate TPM
RPK <- data.frame(matrix(nrow = nrow(counts), ncol = ncol(counts)))
rownames(RPK) <- rownames(counts)
colnames(RPK) <- colnames(counts)
for (i in 1:ncol(RPK)) { RPK[, i] <- counts[, i]/genelength[, i] }
RPK[is.na(RPK)] = 0
total.RPK <- colSums(RPK)
per_million_scaling_factor <- total.RPK/1000000
TPM <- data.frame(matrix(nrow = nrow(counts), ncol = ncol(counts)))
rownames(TPM) <- rownames(counts)
colnames(TPM) <- colnames(counts)
for (i in 1:ncol(TPM)) { TPM[, i] <- RPM[, i]/per_million_scaling_factor[i] }
hist(unlist(TPM[1,]),breaks = 50)

# Sample means should be equal in TPM
colSums(RPKM)[1:2]
colSums(TPM)[1:2]
colMeans(RPKM)[1:2]
colMeans(TPM)[1:2]

# test
txt <- read.delim("E:/2018-08-25 bulk islet aging data analysis/bulk RNA in aged islets/TPM_allSamples_Genes_ReadCount.txt", row.names =1)
colMeans(txt)[1:2]

# match to RNAseq
counts <- counts[, which(names(counts) %in% cols_clamp$ProbeIDs )]

# Save files
save(counts, cols_clamp, file="E:/2021-03-26 Bulk Fat RNAseq/RNAseq Data/counts_data.Rdata")





rm(list=ls())# Clean workspace
if(!is.null(dev.list())) dev.off()# Clear plots

# if using normalized/RPKM data:
# counts2 <- read.delim("RPKM_allSamples_Genes_ReadCount.txt", row.names = 1)

# for look at the number for each group
# length(rownames(subset(cols, tissue %in% "monocyte" & treatment2 == "saline" & treatment1 == "placebo" & pre.post == "pre" & subject == "obese")))

