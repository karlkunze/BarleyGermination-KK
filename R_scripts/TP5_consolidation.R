#TP5 consolidation
path<-"C:/Users/Karl/git/BarleyGermination-KK"
setwd(path)
library(tidyverse)
library(readxl)
DHTP5_tr<-read_xlsx("data/Winter_DH/TP5/DHtp5_all_ter_edit.xlsx")
DHTP5_kk<-read_xlsx("data/Winter_DH/TP5/DHtp5_all_kk.xlsx")


#removing notes, which in
str(DHTP5_kk)
str(DHTP5_tr)
TP5_all<-rbind(DHTP5_kk,DHTP5_tr)
TP5_all[!is.na(TP5_all$notes),]

#plot 12359 is indicated to be naked
colnames(TP5_all)[18]<-"note"
table(TP5_all$PLOT)
View(TP5_all)

TP5_all<-TP5_all[!is.na(TP5_all$PLOT),]

TP5_all
TP5_all$PLOT<-gsub("a","",(TP5_all$PLOT))


TP5_all$Check<-NA;TP5_all$Check<-"B"
TP5_all$Mold_kern<-NA
TP5_all$Cooler<-NA
TP5_all$Year<-NA;TP5_all$Year<-"2020"
TP5_all$Timepoint<-NA;TP5_all$Timepoint<-"5"
TP5_all$PM_at_plating<-NA;TP5_all$PM_at_plating<-"152"
TP5_all$Experiment<-NA;TP5_all$Experiment<-"Winter DH Germination"
#from schedule
colnames(TP5_all)[c(1:3,19,5:18,4,20:25)]
TP5_all<-TP5_all[,c(1:3,19,5:18,4,20:25)]

TP5_all[TP5_all$Entry=="Check 1-Flavia",]$Check<-"C"
TP5_all[TP5_all$Entry=="Check 2-Scala",]$Check<-"C"
TP5_all[TP5_all$Entry=="Check 3-DH130910",]$Check<-"C"
TP5_all[TP5_all$Entry=="Check 4-SY Tepee",]$Check<-"C"
TP5_all[TP5_all$Entry=="Check 5-Wintmalt",]$Check<-"C"
TP5_all[TP5_all$Entry=="Check 6-Charles",]$Check<-"C"

View(TP5_all)

TP5_all[TP5_all$Entry=="Check 1-Flavia",]$Entry<-"Flavia"
TP5_all[TP5_all$Entry=="Check 2-Scala",]$Entry<-"Scala"
TP5_all[TP5_all$Entry=="Check 3-DH130910",]$Entry<-"DH130910"
TP5_all[TP5_all$Entry=="Check 4-SY Tepee",]$Entry<-"SY Tepee"
TP5_all[TP5_all$Entry=="Check 5-Wintmalt",]$Entry<-"Wintmalt"
TP5_all[TP5_all$Entry=="Check 6-Charles",]$Entry<-"Charles"


"



table(TP5_all$Entry)
table(Morex2019$Chrom,useNA = "ifany")
TP5_all$Day1Germ<-as.numeric(TP5_all$Day1Germ)
TP5_all$Day2Germ<-as.numeric(TP5_all$Day2Germ)
TP5_all$Day3Germ<-as.numeric(TP5_all$Day3Germ)
TP5_all$Day4Germ<-as.numeric(TP5_all$Day4Germ)
TP5_all$Day5Germ<-as.numeric(TP5_all$Day5Germ)
TP5_all$kernels_remaining<-as.numeric(TP5_all$kernels_remaining)
TP5_all<-TP5_all[with(TP5_all,order(PLOT,Location)),]
colnames(TP5_all)

TP5_all<-add_column(TP5_all,1:nrow(TP5_all),.before="Entry");colnames(TP5_all)[1]<-"Order"
add_column(TP5_all,1:nrow(TP5_all),.before="Entry")
colnames(TP4)[1:10]
write.csv(TP5_all,"data/Winter_DH/TP5/TP5_all.csv")
