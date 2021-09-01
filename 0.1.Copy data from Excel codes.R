rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
memory.limit(size = 80000000)

library(xlsx)
library(ggplot2)
library(ggpmisc)
library(readxl)


## Extracting data from clamp sheet
subject_infor <- read_excel("E:/2021-03-26 Bulk Fat RNAseq/Subjects with Clamp Infor.xlsx",
                            sheet = "Subjects.infor", col_types = c("numeric", "text"))
col.location <- read_excel("E:/2021-03-26 Bulk Fat RNAseq/Subjects with Clamp Infor.xlsx",
                            sheet = "col.location")

for ( i in 1:nrow(subject_infor)) {
  i=i
  x <- read_excel("E:/2021-03-26 Bulk Fat RNAseq/Erit-2 Clamp Calculations Dr M with lean mass (003)nm updated with 206 215_2.xlsx", 
                  sheet = subject_infor$Sheet_Name[i], col_names = FALSE)
  for (b in 1:nrow(col.location)) {
    b =b
    subject_infor[i, col.location$col[b]] <- x[col.location$nrow[b], col.location$ncol_number[b]]
  }
}

## organize the data extracted from clamp sheet
# change date
library("zoo")
subject_infor$DATE <- as.Date(as.numeric(subject_infor$DATE), origin="1899-12-30")

# Consistent gender
table(subject_infor$Gender)
Female = c("F", "Female", "FEMALE"  )    
Male = c("M",   "male",   "Male",   "MALE" )
Gender.new <- ifelse(subject_infor$Gender %in% Female, "F", ifelse(subject_infor$Gender %in% Male, "M", "NA"))
subject_infor$Gender <- Gender.new

# calculate ibm
subject_infor$ibm <- as.numeric(subject_infor$LBM)/as.numeric(subject_infor$WT)
subject_infor <- as.data.frame(subject_infor)

# save data
saveRDS(subject_infor, file = "E:/2021-03-26 Bulk Fat RNAseq/Subjects basic infor from Clamp.Rds")
write.csv(subject_infor, file = "E:/2021-03-26 Bulk Fat RNAseq/Subjects basic infor from Clamp.csv", row.names = FALSE)
write.xlsx(subject_infor, file = "E:/2021-03-26 Bulk Fat RNAseq/Subjects basic infor from Clamp.xlsx", row.names = FALSE)


## Import clamp data 
clamp <- read.xlsx2("E:/2021-03-26 Bulk Fat RNAseq/Subjects with Clamp Infor.xlsx", sheetName = "Subjects", header=TRUE)

# Organize FFA parameters
groups = c( "FFA_basal", "FFA_20", "FFA_80")
Corr.Result_p <- as.data.frame(clamp[, groups])
rownames(Corr.Result_p) <- clamp$Subject_ID
for ( col in colnames(Corr.Result_p)) { Corr.Result_p <- Corr.Result_p[!grepl("ERROR", Corr.Result_p[, col]),] }

# calculate AUC
library(bayestestR)
for (g in 1:nrow(Corr.Result_p)) {
  Corr.Result_p[g, "AUC"] <- area_under_curve(c(0,20,80), as.numeric(Corr.Result_p[g, groups]), method = "trapezoid")
  Corr.Result_p[g, "Total"] <- (as.numeric(Corr.Result_p[g, groups[1]]) * (80-0)) 
}
Corr.Result_p[["FFA.revAUC"]] <- Corr.Result_p[["Total"]] - Corr.Result_p[["AUC"]]
hist(Corr.Result_p[["FFA.revAUC"]])
Corr.Result_p[["Subject_ID"]] <- rownames(Corr.Result_p)

# add FFA.revAUC to clamp data
clamp <- merge(clamp, Corr.Result_p[, c("Subject_ID", "FFA.revAUC")], by = "Subject_ID", all.x = TRUE)

## Combine all information
all.infor <- merge(subject_infor, clamp, by= intersect(colnames(subject_infor), colnames(clamp)) , all = TRUE)

# save data
saveRDS(all.infor, file = "E:/2021-03-26 Bulk Fat RNAseq/Subjects with Clamp Infor and All Infor.RDS")
write.csv(all.infor, file = "E:/2021-03-26 Bulk Fat RNAseq/Subjects with Clamp Infor and All Infor.csv", row.names = FALSE)
write.xlsx(all.infor, "E:/2021-03-26 Bulk Fat RNAseq/Subjects with Clamp Infor and All Infor.xlsx", row.names = FALSE)


