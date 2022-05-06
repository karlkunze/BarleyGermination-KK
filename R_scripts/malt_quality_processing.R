#file to read in data for malting quality
#need to be aware of positioning on the excel file, as that is how we have info on TMC checks
library(openxlsx)
library(dplyr)
library(stringr)
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
View(SWE_0201)
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

MQ_list<-as.data.frame(read.xlsx("data/Phenotype_Data/2021/MQ_samples_schedule.xlsx",
                                 sheet = "WinterSamples"))%>%filter(rep==1)

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
                                  sheet = "Data",startRow = 9,detectDates = TRUE))%>%rename(Set=1,Extract_number=2,Entry=4,Treatment=3,PLOT=5,Date=6,DP=7,AA=8)%>%select(1:8)%>%slice(1:48)
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

