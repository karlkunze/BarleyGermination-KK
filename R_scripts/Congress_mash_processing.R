
library(dplyr)
library(tidyr)
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