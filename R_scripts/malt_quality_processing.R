#file to read in data for malting quality
#need to be aware of positioning on the excel file, as that is how we have info on TMC checks
library(openxlsx)
library(dplyr)
library(stringr)


#MQ_list in case only plot numbers are presented
MQ_list<-as.data.frame(read.xlsx("data/Phenotype_Data/2021/MQ_samples_schedule.xlsx",
                                 sheet = "WinterSamples"))%>%filter(rep==1)



SWE_0128<-as.data.frame(read.xlsx("data/Malting_quality/SWE/SWE 01-28-22 21CYGGS1216-endspringtp2_WMB6011-6051-as.xlsx",
                           sheet = "Data",startRow = 9,detectDates = TRUE))%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=4,PLOT=5,Date=6,DP=7,AA=8)%>%select(1:8)

t<-SWE_0128 
total<-as.data.frame(table(t$Date))[1,2]
t[t$PLOT=="Tradition Malt Check"&!is.na(t$PLOT),c("Entry")]<-"TMC"

t[1:24,]$Set<-"Set 1"; t[25:48,]$Set<-"Set 2";t[49:72,]$Set<-"Set 3";t$Date<-t$Date[3];t<-t%>%slice(1:72)
t$Check<-"E";t[t$Entry%in%c("R","TMC"),]$Check<-"C";t$Check<-as.factor(t$Check)
#specific to split spring and winter
SWE_0128<-t
SWE_0128[1:56,]$Treatment<-"SpringTP2-4";SWE_0128[57:72,]$Treatment<-"WinterTP1-1"

SWE_0128<-SWE_0128%>%filter(Treatment=="WinterTP1-1"|Check=="C")%>%arrange(Extract_number)
SWE_0128$ID<-as.numeric(paste0(1,format(SWE_0128$Date, "%m"),format(SWE_0128$Date, "%d"),sprintf('%02d', SWE_0128$Extract_number)))
SWE_0128[SWE_0128$Check=="C",]$PLOT<-SWE_0128[SWE_0128$Check=="C",]$ID
#SWE_0128$ID<-paste0(SWE_0128$ID,"-T",sapply(str_split(SWE_0128$Treatment,"T"), "[[" , 2))
SWE_0128$ID<-paste0(SWE_0128$ID,"-TP1-1")
SWE_0128[SWE_0128$Entry=="R",]$Treatment<-"Rahr"

#SWE 1-31-2022

SWE_0131<-as.data.frame(read.xlsx("data/Malting_quality/SWE/SWE 01-31-22 21CYGGS6060-6290.xlsx",
                                  sheet = "Data",startRow = 9,detectDates = TRUE))%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=5,PLOT=4,Date=6,DP=7,AA=8)%>%select(1:8)%>%slice(1:72)
total<-as.data.frame(table(SWE_0131$Date))[1,2]
SWE_0131
t<-SWE_0131
t[1:24,]$Set<-"Set 1"; t[25:48,]$Set<-"Set 2";t[49:72,]$Set<-"Set 3";t$Date<-t$Date[3]
t[t$PLOT=="Tradition Malt Check",]$Entry<-"TMC"
t[t$PLOT=="R",]$Entry<-"R"
t$Check<-"E";t[t$Entry%in%c("R","TMC"),]$Check<-"C";t$Check<-as.factor(t$Check)
t[t$Set=="Set 1",]$Treatment<-"WinterTP1-1"
t[t$Set=="Set 2",]$Treatment<-"WinterTP1-1"
t[t$Set=="Set 3",]$Treatment<-"WinterTP1-2"


t$ID<-as.numeric(paste0(1,format(t$Date, "%m"),format(t$Date, "%d"),sprintf('%02d', t$Extract_number)))
t[t$Check=="C",]$PLOT<-t[t$Check=="C",]$ID
t$ID<-paste0(t$ID,"-T",sapply(str_split(t$Treatment,"T"), "[[" , 2))
t[t$Entry=="R",]$Treatment<-"Rahr"
SWE_0131<-t
SWE_0131
#

SWE_0201<-as.data.frame(read.xlsx("data/Malting_quality/SWE/SWE 02-01-22  21CYGGS6300-7018.xlsx",
                                  sheet = "Data",startRow = 9,detectDates = TRUE))%>%rename(Set=1,Extract_number=2,Entry=4,Treatment=5,PLOT=3,Date=6,DP=7,AA=8)%>%select(1:8)%>%slice(1:72)
total<-as.data.frame(table(SWE_0201$Date))[1,2]
total

t<-SWE_0201

t$PLOT
t[1:24,]$Set<-"Set 1"; t[25:48,]$Set<-"Set 2";t$Date<-t$Date[4]
t[t$PLOT=="Tradition Malt Check",]$Entry<-"TMC"
t[t$PLOT=="R",]$Entry<-"R"
t$Check<-"E";t[t$Entry%in%c("R","TMC"),]$Check<-"C";t$Check<-as.factor(t$Check)
table(t$Treatment)
t[t$Set=="Set 1",]$Treatment<-"WinterTP1-2"
t[t$Set=="Set 2",]$Treatment<-"WinterTP1-2"
t$ID<-as.numeric(paste0(1,format(t$Date, "%m"),format(t$Date, "%d"),sprintf('%02d', t$Extract_number)))
t[t$Check=="C",]$PLOT<-t[t$Check=="C",]$ID
t$ID<-paste0(t$ID,"-T",sapply(str_split(t$Treatment,"T"), "[[" , 2))
t[t$Entry=="R",]$Treatment<-"Rahr"
SWE_0201<-t

# 02-02-2022
SWE_0202<-as.data.frame(read.xlsx("data/Malting_quality/SWE/SWE 02-02-22 21CYGGW7019-7194.xlsx",
                                  sheet = "Data",startRow = 9,detectDates = TRUE))%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=4,PLOT=5,Date=6,DP=7,AA=8)%>%select(1:8)%>%slice(1:72)
total<-as.data.frame(table(SWE_0202$Date))[1,2]
total
t<-SWE_0202

t$PLOT
t[1:24,]$Set<-"Set 1"; t[25:48,]$Set<-"Set 2";t$Date<-t$Date[4]

t[t$PLOT=="Tradition Malt Check",]$Entry<-"TMC"
t[t$PLOT=="R",]$Entry<-"R"
t$Check<-"E";t[t$Entry%in%c("R","TMC"),]$Check<-"C";t$Check<-as.factor(t$Check)
table(t$Treatment)
t[1:4,]$Treatment<-"WinterTP1-2"
t[5:48,]$Treatment<-"WinterTP1-3"
t$ID<-as.numeric(paste0(1,format(t$Date, "%m"),format(t$Date, "%d"),sprintf('%02d', t$Extract_number)))
t[t$Check=="C",]$PLOT<-t[t$Check=="C",]$ID
t$ID<-paste0(t$ID,"-T",sapply(str_split(t$Treatment,"T"), "[[" , 2))
t[t$Entry=="R",]$Treatment<-"Rahr"
SWE_0202<-t
# 02-03-2022

# 02-03-2022
SWE_0203<-as.data.frame(read.xlsx("data/Malting_quality/SWE/SWE 02-03-22 21CYGGW7176-7335.xlsx",
                                  sheet = "Data",startRow = 9,detectDates = TRUE))%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=4,PLOT=5,Date=6,DP=7,AA=8)%>%select(1:8)%>%slice(1:72)
total<-as.data.frame(table(SWE_0203$Date))[1,2]
t<-SWE_0203
#Entries are missing, will use the MQsamples excel file for reference

MQ_list

t[1:24,]$Set<-"Set 1"; t[25:48,]$Set<-"Set 2";t$Date<-t$Date[4]
t$PLOT
t[is.na(t$PLOT),]$PLOT<-"Rahr"
t[t$PLOT=="Tradition Malt Check",]$Entry<-"TMC"

t[t$PLOT=="Rahr",]$Entry<-"R"
t[t$PLOT=="Rahr",]$PLOT<-"R"

t$Check<-"E";t[t$Entry%in%c("R","TMC"),]$Check<-"C";t$Check<-as.factor(t$Check)
table(t$Treatment)
t[1:24,]$Treatment<-"WinterTP1-3"
t[25:48,]$Treatment<-"WinterTP1-4"
t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))

t$ID<-as.numeric(paste0(1,format(t$Date, "%m"),format(t$Date, "%d"),sprintf('%02d', t$Extract_number)))
t[t$Check=="C",]$PLOT<-t[t$Check=="C",]$ID
t$ID<-paste0(t$ID,"-T",sapply(str_split(t$Treatment,"T"), "[[" , 2))
t[t$Entry=="R",]$Treatment<-"Rahr"
SWE_0203<-t
# 2-04-2022
SWE_0204<-as.data.frame(read.xlsx("data/Malting_quality/SWE/SWE 02-04-22 21CYGGW7337-7419.xlsx",
                                  sheet = "Data",startRow = 9,detectDates = TRUE))%>%rename(Set=1,Extract_number=2,Entry=4,Treatment=3,PLOT=5,Date=6,DP=7,AA=8)%>%select(1:8)%>%slice(1:72)
t<-SWE_0204
total<-as.data.frame(table(t$Date))[1,2]
t[1:24,]$Set<-"Set 1";t$Date<-t$Date[4]
t$PLOT

t[t$PLOT=="Traditional Malt Check",]$Entry<-"TMC"

t[t$PLOT=="Rahr",]$Entry<-"R"
t[t$PLOT=="Rahr",]$PLOT<-"R"

t$Check<-"E";t[t$Entry%in%c("R","TMC"),]$Check<-"C";t$Check<-as.factor(t$Check)
table(t$Treatment)
t[1:24,]$Treatment<-"WinterTP1-4"
#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
t$ID<-as.numeric(paste0(1,format(t$Date, "%m"),format(t$Date, "%d"),sprintf('%02d', t$Extract_number)))
t[t$Check=="C",]$PLOT<-t[t$Check=="C",]$ID
t$ID<-paste0(t$ID,"-T",sapply(str_split(t$Treatment,"T"), "[[" , 2))
t[t$Entry=="R",]$Treatment<-"Rahr"
SWE_0204<-t

# 02-07-2022
SWE_0207<-as.data.frame(read.xlsx("data/Malting_quality/SWE/SWE 02-07-22 21CYGGW7439-6095.xlsx",
                                  sheet = "Data",startRow = 9,detectDates = TRUE))%>%rename(Set=1,Extract_number=2,Entry=4,Treatment=3,PLOT=5,Date=6,DP=7,AA=8)%>%select(1:8)
t<-SWE_0207
t

total<-as.data.frame(table(t$Date))[1,2]
total
t[1:24,]$Set<-"Set 1";t[25:48,]$Set<-"Set 2";t$Date<-t$Date[4]
t[t$Extract_number==36,c("Entry","Treatment","PLOT")]<-c("TMC","TMC","TMC")
t$Treatment
t[t$Treatment=="RAHR",]$Entry<-"R"

t$Check<-"E";t[t$Entry%in%c("R","TMC"),]$Check<-"C";t$Check<-as.factor(t$Check)
table(t$Treatment)
t<-t%>%arrange(Extract_number)
t[1:18,]$Treatment<-"WinterTP1-4"
t[19:48,]$Treatment<-"WinterTP2-1"
#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
t$ID<-as.numeric(paste0(1,format(t$Date, "%m"),format(t$Date, "%d"),sprintf('%02d', t$Extract_number)))
t[t$Check=="C",]$PLOT<-t[t$Check=="C",]$ID
t$ID<-paste0(t$ID,"-T",sapply(str_split(t$Treatment,"T"), "[[" , 2))
t[t$Entry=="R",]$Treatment<-"Rahr"
SWE_0207<-t
#View(SWE_0207)
#02-10-22
SWE_0207$Treatment
SWE_0210<-as.data.frame(read.xlsx("data/Malting_quality/SWE/SWE 02-10-22 21CYGGW7421-6245.xlsx",
                                  sheet = "Data",startRow = 9,detectDates = TRUE))%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=4,PLOT=5,Date=6,DP=7,AA=8)%>%select(1:8)
t<-SWE_0210
t

total<-as.data.frame(table(t$Date))[1,2]
total
t[1:24,]$Set<-"Set 1";t[25:48,]$Set<-"Set 2";t$Date<-t$Date[4]

t$PLOT
t[t$PLOT=="Rahr",]$Entry<-"R"
t[t$PLOT=="TMC",]$Entry<-"TMC"
t$Check<-"E";t[t$Entry%in%c("R","TMC"),]$Check<-"C";t$Check<-as.factor(t$Check)
table(t$Treatment)
#View(t)
t<-t%>%arrange(Extract_number)
t[1:2,]$Treatment<-"WinterTP2-1"
t[3:6,]$Treatment<-"WinterTP1-4"
t[7:44,]$Treatment<-"WinterTP2-1"
t[45:48,]$Treatment<-"WinterTP2-2"
#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
t$ID<-as.numeric(paste0(1,format(t$Date, "%m"),format(t$Date, "%d"),sprintf('%02d', t$Extract_number)))
t[t$Check=="C",]$PLOT<-t[t$Check=="C",]$ID
t$ID<-paste0(t$ID,"-T",sapply(str_split(t$Treatment,"T"), "[[" , 2))
t[t$Entry=="R",]$Treatment<-"Rahr"
SWE_0210<-t
#02
SWE_0211<-as.data.frame(read.xlsx("data/Malting_quality/SWE/SWE 02-11-22 21CYGGW6246-6347.xlsx",
                                  sheet = "Data",startRow = 9,detectDates = TRUE))%>%rename(Set=1,Extract_number=2,Entry=4,Treatment=3,PLOT=5,Date=6,DP=7,AA=8)%>%select(1:8)
t<-SWE_0211
t

total<-as.data.frame(table(t$Date))[1,2]
total
t[1:24,]$Set<-"Set 1";#t[25:48,]$Set<-"Set 2"
t$Date<-t$Date[4]
t[13,]$PLOT
t[13,]$Entry<-"TMC"
t[13,]$Treatment<-"WinterTP2-2"
t[t$Treatment=="Rahr",]$Entry<-"R"
t[t$Treatment=="Rahr",]$PLOT<-"R"
t$Check<-"E";t[t$Entry%in%c("R","TMC"),]$Check<-"C";t$Check<-as.factor(t$Check)
table(t$Treatment)
#View(t)
t<-t%>%arrange(Extract_number)
t$Treatment<-"WinterTP2-2"
t$Treatment
#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
t$ID<-as.numeric(paste0(1,format(t$Date, "%m"),format(t$Date, "%d"),sprintf('%02d', t$Extract_number)))
t[t$Check=="C",]$PLOT<-t[t$Check=="C",]$ID
t$ID<-paste0(t$ID,"-T",sapply(str_split(t$Treatment,"T"), "[[" , 2))
t[t$Entry=="R",]$Treatment<-"Rahr"
SWE_0211<-t

#
SWE_0215<-as.data.frame(read.xlsx("data/Malting_quality/SWE/SWE 02-15-22 21CYGGW6342-7030.xlsx",
                                  sheet = "Data",startRow = 9,detectDates = TRUE))%>%rename(Set=1,Extract_number=2,Entry=4,Treatment=3,PLOT=5,Date=6,DP=7,AA=8)%>%select(1:8)
t<-SWE_0215
t

total<-as.data.frame(table(t$Date))[1,2]
total
t[1:24,]$Set<-"Set 1";t[25:48,]$Set<-"Set 2"
t$Date<-t$Date[4]
t[25:26,]$PLOT<-"Rahr"
t[t$PLOT=="Rahr",]$Entry<-"R"
t[t$PLOT=="Traditional Malt Check",]$Entry<-"TMC"
t$Check<-"E";t[t$Entry%in%c("R","TMC"),]$Check<-"C";t$Check<-as.factor(t$Check)
table(t$Treatment)
#View(t)
t<-t%>%arrange(Extract_number)
t[1:41,]$Treatment<-"WinterTP2-2"
t[42:48,]$Treatment<-"WinterTP2-3"
#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
t$ID<-as.numeric(paste0(1,format(t$Date, "%m"),format(t$Date, "%d"),sprintf('%02d', t$Extract_number)))
t[t$Check=="C",]$PLOT<-t[t$Check=="C",]$ID
t$ID<-paste0(t$ID,"-T",sapply(str_split(t$Treatment,"T"), "[[" , 2))
t[t$Entry=="R",]$Treatment<-"Rahr"
SWE_0215<-t
SWE_0215
##2 17 2022

SWE_0217<-as.data.frame(read.xlsx("data/Malting_quality/SWE/SWE 02-17-22 21CYGGW7035-7377.xlsx",
                                  sheet = "Data",startRow = 9,detectDates = TRUE))%>%rename(Set=1,Extract_number=2,Entry=4,Treatment=3,PLOT=5,Date=6,DP=7,AA=8)%>%select(1:8)
t<-SWE_0217
t

total<-as.data.frame(table(t$Date))[1,2]
total
t[1:24,]$Set<-"Set 1"#;t[25:48,]$Set<-"Set 2"
t$Date<-t$Date[4]
t$PLOT
t[t$PLOT=="Rahr",]$Entry<-"R"
t[t$PLOT=="TMC",]$Entry<-"TMC"
t$Check<-"E";t[t$Entry%in%c("R","TMC"),]$Check<-"C";t$Check<-as.factor(t$Check)
table(t$Treatment)
#View(t)
t<-t%>%arrange(Extract_number)
t[1:18,]$Treatment<-"WinterTP2-3"
t[19:24,]$Treatment<-"WinterTP2-4"

#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
t$ID<-as.numeric(paste0(1,format(t$Date, "%m"),format(t$Date, "%d"),sprintf('%02d', t$Extract_number)))
t[t$Check=="C",]$PLOT<-t[t$Check=="C",]$ID
t$ID<-paste0(t$ID,"-T",sapply(str_split(t$Treatment,"T"), "[[" , 2))
t[t$Entry=="R",]$Treatment<-"Rahr"
SWE_0217<-t

##


SWE_0218<-as.data.frame(read.xlsx("data/Malting_quality/SWE/SWE 02-18-22 21CYGGW7380-7471.xlsx",
                                  sheet = "Data",startRow = 9,detectDates = TRUE))%>%rename(Set=1,Extract_number=2,Entry=4,Treatment=3,PLOT=5,Date=6,DP=7,AA=8)%>%select(1:8)
t<-SWE_0218
t

total<-as.data.frame(table(t$Date))[1,2]
total
t[1:24,]$Set<-"Set 1"#;t[25:48,]$Set<-"Set 2"
t$Date<-t$Date[4]
t$PLOT
t[t$PLOT=="Rahr",]$Entry<-"R"
t[t$PLOT=="Traditional Malt Check",]$Entry<-"TMC"
t$Check<-"E";t[t$Entry%in%c("R","TMC"),]$Check<-"C";t$Check<-as.factor(t$Check)
table(t$Treatment)
#View(t)
t<-t%>%arrange(Extract_number)
t[1:24,]$Treatment<-"WinterTP2-4"


#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
t$ID<-as.numeric(paste0(1,format(t$Date, "%m"),format(t$Date, "%d"),sprintf('%02d', t$Extract_number)))
t[t$Check=="C",]$PLOT<-t[t$Check=="C",]$ID
t$ID<-paste0(t$ID,"-T",sapply(str_split(t$Treatment,"T"), "[[" , 2))
t[t$Entry=="R",]$Treatment<-"Rahr"
SWE_0218<-t
###
#02-28-2022
SWE_0228<-as.data.frame(read.xlsx("data/Malting_quality/SWE/SWE 02-28-22 21CYGGW7099-7274.xlsx",
                                  sheet = "Data",startRow = 9,detectDates = TRUE))%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=4,PLOT=5,Date=6,DP=7,AA=8)%>%select(1:8)
t<-SWE_0228
t
#View(t)
total<-as.data.frame(table(t$Date))[1,2]
total
t[1:24,]$Set<-"Set 1";t[25:48,]$Set<-"Set 2"
t$Date<-t$Date[4]
t$PLOT
t[t$PLOT=="Rahr",]$PLOT<-"R"
t[t$PLOT=="TMC",]$Entry<-"TMC"
t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
t$Check<-"E";t[t$Entry%in%c("R","TMC"),]$Check<-"C";t$Check<-as.factor(t$Check)
table(t$Treatment)

#0228_iii
#View(t)
t<-t%>%arrange(Extract_number)
#View(t)
t[t$Entry=="R",]$Treatment<-"WinterTP2-3"


t$ID<-as.numeric(paste0(1,format(t$Date, "%m"),format(t$Date, "%d"),sprintf('%02d', t$Extract_number)))
t[t$Check=="C",]$PLOT<-t[t$Check=="C",]$ID
t$ID<-paste0(t$ID,"-T",sapply(str_split(t$Treatment,"T"), "[[" , 2))
t[t$Entry=="R",]$Treatment<-"Rahr"
SWE_0228<-t
###

SWE_0302<-as.data.frame(read.xlsx("data/Malting_quality/SWE/SWE 03-02-22 21CYGGW7049-7358.xlsx",
                                  sheet = "Data",startRow = 9,detectDates = TRUE))%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=4,PLOT=5,Date=6,DP=7,AA=8)%>%select(1:8)
t<-SWE_0302%>%slice(1:37)
t
#View(t)
total<-as.data.frame(table(t$Date))[1,2]
total
t[1:24,]$Set<-"Set 1";t[25:37,]$Set<-"Set 2"
t$Date<-t$Date[4]
t$PLOT
t[t$PLOT=="Rahr",]$PLOT<-"R"
t[t$PLOT=="TMC",]$Entry<-"TMC"
t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
t$Check<-"E";t[t$Entry%in%c("R","TMC"),]$Check<-"C";t$Check<-as.factor(t$Check)
table(t$Treatment)

#0228_iii
#View(t)
t<-t%>%arrange(Extract_number)
#View(t)
t[t$Entry=="R",]$Treatment<-"WinterTP2-3"


t$ID<-as.numeric(paste0(1,format(t$Date, "%m"),format(t$Date, "%d"),sprintf('%02d', t$Extract_number)))
t[t$Check=="C",]$PLOT<-t[t$Check=="C",]$ID
t$ID<-paste0(t$ID,"-T",sapply(str_split(t$Treatment,"T"), "[[" , 2))
t[t$Entry=="R",]$Treatment<-"Rahr"
SWE_0302<-t
View(SWE_0302)
library(tidyverse)
SWE=plyr::rbind.fill(SWE_0128,SWE_0131,SWE_0201,SWE_0202,SWE_0203,SWE_0204,SWE_0207,SWE_0210,SWE_0211,SWE_0215,SWE_0217,SWE_0218,
                 SWE_0228,SWE_0302)
SWE[SWE$Entry%in%c("Tradition Malt Check"),]$Entry<-"TMC"
SWE$order_column<-1:nrow(SWE)
y<-SWE%>%filter(!Treatment=="Rahr")%>%arrange(Treatment,Date,Extract_number)
y$unique_MQ_ID<-1:nrow(y)+3000;y<-y%>%select(Treatment,Date,Extract_number,PLOT,Entry,order_column,unique_MQ_ID,ID)%>%select(ID,unique_MQ_ID)

SWE<-SWE%>%left_join(y,by=c("ID"))#%>%arrange(unique_MQ_ID)
View(SWE)
save(SWE,file = "data/Malting_quality/SWE/WMB21_SWE.Rdata")

#Congress_mash
MQ_list<-as.data.frame(read.xlsx("data/Phenotype_Data/2021/MQ_samples_schedule.xlsx",
                                 sheet = "WinterSamples"))%>%filter(rep==1)




#c(SWE$Entry,SWE$unique_MQ_ID)

CM_0128<-as.data.frame(read.xlsx("data/Malting_quality/CM/CM 01-28-22 21CYGGS1208-endTP2Springas.xlsx",
                                  sheet = "Data",startRow = 8,detectDates = TRUE))%>%select(1:9)%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=4,PLOT=5,Date=6,ME=7,BG=8,FAN=9)%>%select(1:9)%>%
drop_na(Date)
Date1<-w$Date[3]
w<-CM_0128[c(56,58:72),]
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$Entry<-"TMC"
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-"TMC-1"

SP_0128<-as.data.frame(read.xlsx("data/Malting_quality/CM/CM 01-28-22 21CYGGS1208-endTP2Springas.xlsx",
                                 sheet = "SP3",detectDates = TRUE))%>%select(1,2,4,5,6,7,8)%>%rename(SP_order=1,SP_ID=2,Date=3,nm1=4,abs_nm1=5,nm2=6,abs_nm2=7)%>%filter(SP_ID%in%c("NaCl","TMC")|SP_ID>6000)%>%mutate(delta=abs_nm1-abs_nm2)%>%mutate(Date=Date1)

SP_0128[SP_0128$SP_ID=="6052",]$SP_ID<-"6051"
s<-SP_0128
s[s$SP_ID=="TMC",]
s[s$SP_ID=="TMC",]$SP_ID<-c("TMC-1","TMC-1")
value<-s%>%group_by(SP_ID) %>%summarize(mean = mean(delta, na.rm = TRUE),sd=sd(delta))
ID<-s%>%arrange%>%select(SP_ID)%>%unique()%>%mutate(order=1:n_distinct(SP_ID))%>%select(2,1)
sp1=ID%>%full_join(value,by="SP_ID")%>%mutate(sd=round(sd,4))%>%rename(PLOT=SP_ID)#%>%filter(!PLOT=="6052")
sp1
cm<-w%>%full_join(sp1,by="PLOT")%>%rename(sp_mean=mean,sp_sd=sd)%>%select(-order)
cm$Date<-Date1
cm$Treatment<-"WinterTP1-1"
cm$Set<-"Set 3"
cm[cm$PLOT%in%c("NaCl"),]$Entry<-"NaCl"

cm[cm$PLOT%in%c("NaCl"),]$Extract_number<-99

#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
cm$Check<-"E";cm[cm$Entry%in%c("NaCl","TMC","TMC-1","TMC-2"),]$Check<-"C";cm$Check<-as.factor(cm$Check)

cm$ID<-as.numeric(paste0(1,format(cm$Date, "%m"),format(cm$Date, "%d"),sprintf('%02d', cm$Extract_number)))
cm$PLOT
cm[cm$Check=="C",]$PLOT<-cm[cm$Check=="C",]$ID
cm$ID<-paste0(cm$ID,"-T",sapply(str_split(cm$Treatment,"T"), "[[" , 2))

CM_0128<-cm
# Congress mash 1-31

CM_0128<-as.data.frame(read.xlsx("data/Malting_quality/CM/CM 01-28-22 21CYGGS1208-endTP2Springas.xlsx",
                                 sheet = "Data",startRow = 8,detectDates = TRUE))%>%select(1:9)%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=4,PLOT=5,Date=6,ME=7,BG=8,FAN=9)%>%select(1:9)%>%
  drop_na(Date)
Date1<-w$Date[3]
w<-CM_0128[c(56,58:72),]
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$Entry<-"TMC"
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-"TMC-1"

SP_0128<-as.data.frame(read.xlsx("data/Malting_quality/CM/CM 01-28-22 21CYGGS1208-endTP2Springas.xlsx",
                                 sheet = "SP3",detectDates = TRUE))%>%select(1,2,4,5,6,7,8)%>%rename(SP_order=1,SP_ID=2,Date=3,nm1=4,abs_nm1=5,nm2=6,abs_nm2=7)%>%filter(SP_ID%in%c("NaCl","TMC")|SP_ID>6000)%>%mutate(delta=abs_nm1-abs_nm2)%>%mutate(Date=Date1)

SP_0128[SP_0128$SP_ID=="6052",]$SP_ID<-"6051"
s<-SP_0128
s[s$SP_ID=="TMC",]
s[s$SP_ID=="TMC",]$SP_ID<-c("TMC-1","TMC-1")
value<-s%>%group_by(SP_ID) %>%summarize(mean = mean(delta, na.rm = TRUE),sd=sd(delta))
ID<-s%>%arrange%>%select(SP_ID)%>%unique()%>%mutate(order=1:n_distinct(SP_ID))%>%select(2,1)
sp1=ID%>%full_join(value,by="SP_ID")%>%mutate(sd=round(sd,4))%>%rename(PLOT=SP_ID)#%>%filter(!PLOT=="6052")
sp1
cm<-w%>%full_join(sp1,by="PLOT")%>%rename(sp_mean=mean,sp_sd=sd)%>%select(-order)
cm$Date<-Date1
cm$Treatment<-"WinterTP1-1"
cm$Set<-"Set 3"
cm[cm$PLOT%in%c("NaCl"),]$Entry<-"NaCl"

cm[cm$PLOT%in%c("NaCl"),]$Extract_number<-99

#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
cm$Check<-"E";cm[cm$Entry%in%c("NaCl","TMC","TMC-1","TMC-2"),]$Check<-"C";cm$Check<-as.factor(cm$Check)

cm$ID<-as.numeric(paste0(1,format(cm$Date, "%m"),format(cm$Date, "%d"),sprintf('%02d', cm$Extract_number)))
cm$PLOT
cm[cm$Check=="C",]$PLOT<-cm[cm$Check=="C",]$ID
cm$ID<-paste0(cm$ID,"-T",sapply(str_split(cm$Treatment,"T"), "[[" , 2))

CM_0128<-cm