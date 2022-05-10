
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)
#Congress_mash
MQ_list<-as.data.frame(read.xlsx("data/Phenotype_Data/2021/MQ_samples_schedule.xlsx",
                                 sheet = "WinterSamples"))%>%filter(rep==1)




#c(SWE$Entry,SWE$unique_MQ_ID)

CM_0128<-as.data.frame(read.xlsx("data/Malting_quality/CM/CM 01-28-22 21CYGGS1208-endTP2Springas.xlsx",
                                 sheet = "Data",startRow = 8,detectDates = TRUE))%>%select(1:9)%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=4,PLOT=5,Date=6,ME=7,BG=8,FAN=9)%>%select(1:9)%>%
  drop_na(Date)
w<-CM_0128[c(56,58:72),]
Date1<-w$Date[3]
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

cm[cm$PLOT%in%c("NaCl"),]$Extract_number<-97

#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
cm$Check<-"E";cm[cm$Entry%in%c("NaCl","TMC","TMC-1","TMC-2"),]$Check<-"C";cm$Check<-as.factor(cm$Check)

cm$ID<-as.numeric(paste0(1,format(cm$Date, "%m"),format(cm$Date, "%d"),sprintf('%02d', cm$Extract_number)))
cm$PLOT
cm[cm$Check=="C",]$PLOT<-cm[cm$Check=="C",]$ID
cm$ID<-paste0(cm$ID,"-T",sapply(str_split(cm$Treatment,"T"), "[[" , 2))

CM_0128<-cm
# Congress mash 1-31

CM_0131<-as.data.frame(read.xlsx("data/Malting_quality/CM/CM 01-31-22 21CYGGS6136-6440as.xlsx",
                                 sheet = "Data",startRow = 8,detectDates = TRUE))%>%select(1:9)%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=5,PLOT=4,Date=6,ME=7,BG=8,FAN=9)%>%select(1:9)%>%
drop_na(Date)
w<-CM_0131
w$Date<-as.Date("2022-01-31")
Date1<-w$Date[3]
CM_0131

w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$Entry<-"TMC"
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-"TMC"
w<-w%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))%>%arrange(Extract_number)

w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-c("TMC-1","TMC-2","TMC-3","TMC-4","TMC-5","TMC-6","TMC-7")

SP_0131_34<-as.data.frame(read.xlsx("data/Malting_quality/CM/CM 01-31-22 21CYGGS6136-6440as.xlsx",
                                 sheet = "SP 3-4",detectDates = TRUE))%>%select(1,2,4,5,6)%>%rename(SP_order=1,SP_ID=2,Date=3,sp_mean=4,sp_sd=5)%>%mutate(Date=Date1,mutate(SP_order=as.numeric(SP_order)),sp_mean=round(as.numeric(sp_mean),4),sp_sd=round(as.numeric(sp_sd),4)) %>%filter(!row_number() %in%c(26))

s<-SP_0131_34
s[s$SP_ID%in%c(6356),]$SP_order<-16
s[s$SP_ID=="TMC",]
s[s$SP_ID=="TMC",]$SP_ID<-c("TMC-4","TMC-5","TMC-6","TMC-7")
s[s$SP_ID%in%c("NaCl", "0.5% NaCl"),]$SP_ID<-"NaCl"
s$Group<-"Group"
s[1:25,]$Group<-"Group_3"
s[26:nrow(s),]$Group<-"Group_4"

s_34<-s%>%mutate(SP_order=as.numeric(SP_order))%>%arrange(Group,SP_order)%>%mutate(SP_order=1:length(SP_order))%>%rename(PLOT=SP_ID,mean=sp_mean,sd=sp_sd)

####
SP_0131_12<-as.data.frame(read.xlsx("data/Malting_quality/CM/CM 01-31-22 21CYGGS6136-6440as.xlsx",
                                 sheet = "SP 1-2",detectDates = TRUE))%>%filter()%>%select(1,2,4,5,6,7,8)%>%
  rename(SP_order=1,SP_ID=2,Date=3,nm1=4,abs_nm1=5,nm2=6,abs_nm2=7)%>%filter(!SP_order%in%c("#"))%>%mutate(abs_nm1=as.numeric(abs_nm1),abs_nm2=as.numeric(abs_nm2),delta=abs_nm2-abs_nm1)%>%mutate(Date=Date1)


s<-SP_0131_12
s
s[s$SP_ID=="TMC",]
s[s$SP_ID=="TMC",]$SP_ID<-c("TMC-1","TMC-1","TMC-1","TMC-1","TMC-2","TMC-2","TMC-2","TMC-2","TMC-3","TMC-3","TMC-3")
s$Group<-"Group"
s[1:99,]$Group<-"Group_1"
s[100:nrow(s),]$Group<-"Group_2"
value<-s%>%group_by(Group,SP_ID) %>%summarize(mean = mean(delta, na.rm = TRUE),sd=sd(delta))%>%ungroup
#View(value)
ID<-s%>%arrange%>%select(SP_ID)%>%unique()%>%mutate(SP_order=1:n_distinct(SP_ID))%>%select(2,1)
sp1=ID%>%full_join(value,by="SP_ID","Group")%>%mutate(sd=round(sd,4))%>%rename(PLOT=SP_ID)%>%arrange(Group,SP_order)
s_131_12<-sp1
#View(s_131_12)
sp<-plyr::rbind.fill(s_131_12,s_34)%>%arrange(Group,SP_order)%>%mutate(Date=Date1,SP_order=1:length(SP_order))

####
sp
cm<-w%>%full_join(sp,by=c("PLOT","Date"))%>%rename(sp_mean=mean,sp_sd=sd)%>%select(-SP_order)%>%arrange(Extract_number,Group)
cm$Date<-Date1
cm$Extract_number
cm
cm[cm$PLOT%in%c("NaCl"),]$Extract_number<-97:100
cm[cm$Group=="Group_1",]$Set<-"Set 1";cm[cm$Group=="Group_2",]$Set<-"Set 2";cm[cm$Group=="Group_3",]$Set<-"Set 3";cm[cm$Group=="Group_4",]$Set<-"Set 4"
#need to think of a way to deal with number ordering for NaCls
cm[1:48,]$Treatment<-"WinterTP1-1"
cm[48:nrow(cm),]$Treatment<-"WinterTP1-2"
View(cm)
cm[cm$PLOT%in%c("NaCl"),]$Entry<-"NaCl"

cm[cm$PLOT%in%c("NaCl"),]$Extract_number<-99

#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
cm$Check<-"E";cm[cm$Entry%in%c("NaCl","TMC","TMC-1","TMC-2"),]$Check<-"C";cm$Check<-as.factor(cm$Check)

cm$ID<-as.numeric(paste0(1,format(cm$Date, "%m"),format(cm$Date, "%d"),sprintf('%02d', cm$Extract_number)))
cm$PLOT
cm[cm$Check=="C",]$PLOT<-cm[cm$Check=="C",]$ID
cm$ID<-paste0(cm$ID,"-T",sapply(str_split(cm$Treatment,"T"), "[[" , 2))

CM_<-cm