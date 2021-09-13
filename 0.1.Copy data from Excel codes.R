#rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
#gc() #free up memrory and report the memory usage.
memory.limit(size = 80000000)

library(xlsx)
library(ggplot2)
library(ggpmisc)
library(readxl)
library(dplyr)
library(forcats)
library(bayestestR)

# set data locations in one place
inputdata <- c(subjects="Subjects with Clamp Infor.xlsx"
               ,dat0="Erit-2 Clamp Calculations Dr M with lean mass (003)nm updated with 206 215_2.xlsx");

## Extracting data from clamp sheet
subject_infor <- read_excel(inputdata['subjects'],
                            sheet = "Subjects.infor", col_types = c("numeric", "text"))
col.location <- read_excel(inputdata['subjects'],
                            sheet = "col.location")

dat0 <- NULL;
for ( i in unique(subject_infor$Sheet_Name)) {
  dat0 <- rbind(dat0
               ,bind_cols(apply(col.location,1,function(thiscol){
                   read_excel(inputdata['dat0'],sheet=i,col_names=F
                                     ,range=thiscol['location']) %>%
                   ifelse(nrow(.)==0,as_tibble(NA),.) %>%
                   setNames(thiscol['col'])})))
};
dat0 <- cbind(subject_infor,dat0);


## organize the data extracted from clamp sheet
# change date
dat0$DATE <- as.Date(dat0$DATE);

# Consistent gender
#+ fix_gender
dat0$Gender <- fct_collapse(dat0$Gender
                            ,`F` = c('F','Female','FEMALE','female')
                            ,`M` =c('M','Male','MALE','male'))

# calculate ibm
dat0$ibm <- with(dat0,LBM/WT)

# save data
#+ save_data
saveRDS(dat0,file='Subjects basic infor from Clamp.Rds')
write.csv(dat0, file = "Subjects basic infor from Clamp.csv", row.names = FALSE)
write.xlsx2(dat0, file = "Subjects basic infor from Clamp.xlsx", row.names = FALSE)


## Import clamp data
#+ clamp_data
clamp <- read_xlsx(inputdata['subjects'],sheet='Subjects');

# Organize FFA parameters
groups = c( "FFA_basal", "FFA_20", "FFA_80")
#+ AUC
clamp$AUC <- apply(clamp,1,function(xx) area_under_curve(c(0,20,80),as.numeric(xx[groups])));
clamp$Total <- clamp[[groups[1]]]*80;
clamp$FFA.revAUC <- with(clamp,Total - AUC);


## Combine all information
all.infor <- merge(dat0, clamp, by= intersect(colnames(dat0), colnames(clamp)) , all = TRUE)

# save data
saveRDS(all.infor, file = "Subjects with Clamp Infor and All Infor.RDS")
write.csv(all.infor, file = "Subjects with Clamp Infor and All Infor.csv", row.names = FALSE)
write.xlsx(all.infor, "Subjects with Clamp Infor and All Infor.xlsx", row.names = FALSE)


