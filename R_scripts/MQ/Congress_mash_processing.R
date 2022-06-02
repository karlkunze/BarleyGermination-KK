
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)
#Congress_mash
MQ_list<-as.data.frame(read.xlsx("data/Phenotype_Data/2021/MQ_samples_schedule.xlsx",
                                 sheet = "WinterSamples"))%>%filter(rep==1)




#c(SWE$Entry,SWE$unique_MQ_ID)

CM_0128<-as.data.frame(read.xlsx("data/MQ/CM/CM 01-28-22 21CYGGS1208-endTP2Springas.xlsx",
                                 sheet = "Data",startRow = 8,detectDates = TRUE))%>%select(1:9)%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=4,PLOT=5,Date=6,ME=7,BG=8,FAN=9)%>%select(1:9)%>%
  drop_na(Date)
nrow(CM_0128)
w<-CM_0128[c(58:72),]
w$PLOT
Date1<-w$Date[3]
Date1

# w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT
# w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$Entry<-"TMC"
# w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-"TMC-1"

SP_0128<-as.data.frame(read.xlsx("data/MQ/CM/CM 01-28-22 21CYGGS1208-endTP2Springas.xlsx",
                                 sheet = "SP3",detectDates = TRUE))%>%select(1,2,4,5,6,7,8)%>%rename(SP_order=1,SP_ID=2,Date=3,nm1=4,abs_nm1=5,nm2=6,abs_nm2=7)%>%filter(SP_ID%in%c("NaCl")|SP_ID>6000)%>%mutate(delta=abs_nm1-abs_nm2)%>%mutate(Date=Date1)

SP_0128[SP_0128$SP_ID=="6052",]$SP_ID<-"6051"
s<-SP_0128

# s[s$SP_ID=="TMC",]
# s[s$SP_ID=="TMC",]$SP_ID<-c("TMC-1","TMC-1")
value<-s%>%group_by(SP_ID) %>%summarize(mean = mean(delta, na.rm = TRUE),sd=sd(delta))
s
ID<-s%>%arrange%>%select(SP_ID)%>%unique()%>%mutate(order=1:n_distinct(SP_ID))%>%select(2,1)
sp1=ID%>%full_join(value,by="SP_ID")%>%mutate(sd=round(sd,4))%>%rename(PLOT=SP_ID)#%>%filter(!PLOT=="6052")
sp1$Group<-"Group_1"

sp1<-sp1%>%filter(!PLOT%in%c("TMC"))
cm<-w%>%full_join(sp1,by="PLOT")%>%rename(sp_mean=mean,sp_sd=sd)%>%select(-order)
cm$Date<-Date1
cm$Treatment
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
rm(cm,sp1,s,w)

# Congress mash 1-31

CM_0131<-as.data.frame(read.xlsx("data/MQ/CM/CM 01-31-22 21CYGGS6136-6440as.xlsx",
                                 sheet = "Data",startRow = 8,detectDates = TRUE))%>%select(1:9)%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=5,PLOT=4,Date=6,ME=7,BG=8,FAN=9)%>%select(1:9)%>%
drop_na(Date)
w<-CM_0131
w$Date<-as.Date("2022-01-31")
Date1<-w$Date[3]
w$PLOT

w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$Entry<-"TMC"
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-"TMC"
w<-w%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))%>%arrange(Extract_number)

w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-c("TMC-1","TMC-2","TMC-3","TMC-4","TMC-5","TMC-6","TMC-7")

SP_0131_34<-as.data.frame(read.xlsx("data/MQ/CM/CM 01-31-22 21CYGGS6136-6440as.xlsx",
                                 sheet = "SP 3-4",detectDates = TRUE))%>%select(1,2,4,5,6)%>%rename(SP_order=1,SP_ID=2,Date=3,sp_mean=4,sp_sd=5)%>%as.data.frame()

SP_0131_34
s<-SP_0131_34%>%filter(!SP_order%in%c("#"))%>%mutate(Date=Date1,SP_order=as.numeric(SP_order),sp_mean=round(as.numeric(sp_mean),4),sp_sd=round(as.numeric(sp_sd),4))

s
s[s$SP_ID%in%c(6356),]$SP_order<-16
s[s$SP_ID=="TMC",]
s[s$SP_ID=="TMC",]$SP_ID<-c("TMC-4","TMC-5","TMC-6","TMC-7")
s[s$SP_ID%in%c("NaCl", "0.5% NaCl"),]$SP_ID<-"NaCl"
s$Group<-"Group"
s[1:25,]$Group<-"Group_3"
s[26:nrow(s),]$Group<-"Group_4"
s
s_34<-s%>%mutate(SP_order=as.numeric(SP_order))%>%arrange(Group,SP_order)%>%mutate(SP_order=1:length(SP_order))%>%rename(PLOT=SP_ID,mean=sp_mean,sd=sp_sd)

####
SP_0131_12<-as.data.frame(read.xlsx("data/MQ/CM/CM 01-31-22 21CYGGS6136-6440as.xlsx",
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
cm$PLOT
cm[46:48,]
cm[1:48,]$Treatment<-"WinterTP1-1"
cm[49:nrow(cm),]$Treatment<-"WinterTP1-2"

cm[cm$PLOT%in%c("NaCl"),]$Entry<-"NaCl"

#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
cm$Check<-"E";cm[cm$Entry%in%c("NaCl","TMC","TMC-1","TMC-2"),]$Check<-"C";cm$Check<-as.factor(cm$Check)

cm$ID<-as.numeric(paste0(1,format(cm$Date, "%m"),format(cm$Date, "%d"),sprintf('%02d', cm$Extract_number)))

cm[cm$Check=="C",]$PLOT<-cm[cm$Check=="C",]$ID
cm$ID<-paste0(cm$ID,"-T",sapply(str_split(cm$Treatment,"T"), "[[" , 2))
cm$Entry
View(cm)
CM_0131<-cm

### Congress mash 02 01 2022

CM_0201<-as.data.frame(read.xlsx("data/MQ/CM/CM 02-01-22 21CYGGS6441-7126as.xlsx",
                                 sheet = "Data",startRow = 8,detectDates = TRUE))%>%select(1:10)%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=4,PLOT=5,Date=6,ME=7,BG=8,FAN=9,mean_sp=10)%>%select(1:10)%>%
  drop_na(Date)
w<-CM_0201%>%arrange(Extract_number)

Date1<-w$Date[3]

nrow(w)
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$Entry<-"TMC"
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-"TMC"
#w<-w%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))%>%arrange(Extract_number)
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-c("TMC-1","TMC-2","TMC-3")
cm<-w
cm[1:24,]$Set<-"Set 1"
cm[25:48,]$Set<-"Set 2"
cm$Group<-"Group"
cm[cm$Set=="Set 1",]$Group<-"Group_1"
cm[cm$Set=="Set 2",]$Group<-"Group_2"
#View(cm)
cm[1:17,]
cm[1:16,]$Treatment<-"WinterTP1-2"
cm[17:nrow(cm),]$Treatment<-"WinterTP1-3"
cm$Entry

#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
cm$Check<-"E";cm[cm$Entry%in%c("NaCl","TMC"),]$Check<-"C";cm$Check<-as.factor(cm$Check)

cm$ID<-as.numeric(paste0(1,format(cm$Date, "%m"),format(cm$Date, "%d"),sprintf('%02d', cm$Extract_number)))
cm$PLOT
cm[cm$Check=="C",]$PLOT<-cm[cm$Check=="C",]$ID
cm$ID<-paste0(cm$ID,"-T",sapply(str_split(cm$Treatment,"T"), "[[" , 2))
names(cm)

CM_0201<-cm%>%rename(sp_mean=mean_sp)
CM_0201$Entry
rm(cm,sp1,s,w,ID)
#2-02-2022
#rbind(CM_0128,CM_0131,CM_0201)
table(CM_0201$Treatment)
CM_0202<-as.data.frame(read.xlsx("data/MQ/CM/CM 02-02-22 21CYGGS7129-7314.xlsx",
                                 sheet = "Data",startRow = 8,detectDates = TRUE))%>%select(1:9)%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=4,PLOT=5,Date=6,ME=7,BG=8,FAN=9)%>%select(1:9)%>%
  drop_na(Date)
w<-CM_0202
w$Date<-as.Date("2022-02-02")
Date1<-w$Date[3]
CM_0202

w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$Entry<-"TMC"
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-"TMC"
#w<-w%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))%>%arrange(Extract_number)
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-c("TMC-1","TMC-2","TMC-3")


####
SP_0202<-as.data.frame(read.xlsx("data/MQ/CM/CM 02-02-22 21CYGGS7129-7314.xlsx",
                                    sheet = "SP, %",detectDates = TRUE))%>%filter()%>%select(1,2,4,5,6,7,8)%>%
  rename(SP_order=1,SP_ID=2,Date=3,nm1=4,abs_nm1=5,nm2=6,abs_nm2=7)%>%filter(!SP_order%in%c("#"))%>%mutate(abs_nm1=as.numeric(abs_nm1),abs_nm2=as.numeric(abs_nm2),delta=abs_nm2-abs_nm1)%>%mutate(Date=Date1)


s<-SP_0202%>%mutate(SP_order=1:nrow(SP_0202))
#View(s)
s[s$SP_ID=="TMC",]
s[s$SP_ID=="TMC",]$SP_ID<-c("TMC-1","TMC-1","TMC-1","TMC-1","TMC-2","TMC-2","TMC-2","TMC-2","TMC-3","TMC-3","TMC-3","TMC-3")
s$Group<-"Group"
s[s$SP_ID%in%c("NACL", "NACL(Reblank)"),]$SP_ID<-"NaCl"

s[1:99,]$Group<-"Group_1"
s[100:nrow(s),]$Group<-"Group_2"
value<-s%>%group_by(Group,SP_ID) %>%summarize(mean = mean(delta, na.rm = TRUE),sd=sd(delta))%>%ungroup
#View(value)
ID<-s%>%arrange%>%select(SP_ID)%>%unique()%>%mutate(SP_order=1:n_distinct(SP_ID))%>%select(2,1)
sp=ID%>%full_join(value,by="SP_ID","Group")%>%mutate(sd=round(sd,4))%>%rename(PLOT=SP_ID)%>%arrange(Group,SP_order)
sp$Date<-Date1
SP_0202<-sp

#View(s_131_12)
#sp<-plyr::rbind.fill(s_131_12,s_34)%>%arrange(Group,SP_order)%>%mutate(Date=Date1,SP_order=1:length(SP_order))

####

cm<-w%>%full_join(SP_0202,by=c("PLOT","Date"))%>%rename(sp_mean=mean,sp_sd=sd)%>%select(-SP_order)%>%arrange(Extract_number,Group)
cm$Date<-Date1
cm$Extract_number
#View(cm)
cm[cm$PLOT%in%c("NaCl"),]$Extract_number<-97
cm[cm$Group=="Group_1",]$Set<-"Set 1";cm[cm$Group=="Group_2",]$Set<-"Set 2"#;cm[cm$Group=="Group_3",]$Set<-"Set 3";cm[cm$Group=="Group_4",]$Set<-"Set 4"
#need to think of a way to deal with number ordering for NaCls
cm$PLOT
cm[1:32,]$Treatment<-"WinterTP1-3"
cm[33:nrow(cm),]$Treatment<-"WinterTP1-4"
#View(cm)
cm[cm$PLOT%in%c("NaCl"),]$Entry<-"NaCl"

#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
cm$Check<-"E";cm[cm$Entry%in%c("NaCl","TMC","TMC-1","TMC-2"),]$Check<-"C";cm$Check<-as.factor(cm$Check)

cm$ID<-as.numeric(paste0(1,format(cm$Date, "%m"),format(cm$Date, "%d"),sprintf('%02d', cm$Extract_number)))
cm$PLOT
cm[cm$Check=="C",]$PLOT<-cm[cm$Check=="C",]$ID
cm$ID<-paste0(cm$ID,"-T",sapply(str_split(cm$Treatment,"T"), "[[" , 2))

CM_0202<-cm
rm(cm,s,w,ID)

#02 03 2022




CM_0203<-as.data.frame(read.xlsx("data/MQ/CM/CM 02-03-22 21CYGGS7318-7477as.xlsx",
                                 sheet = "Data",startRow = 8,detectDates = TRUE))%>%select(1:10)%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=4,PLOT=5,Date=6,ME=7,BG=8,FAN=9,sp_mean=10)%>%select(1:9)%>%
  drop_na(Date)
w<-CM_0203
w
Date1<-w$Date[3]

w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$Entry<-"TMC"
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-"TMC"
#w<-w%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))%>%arrange(Extract_number)
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-c("TMC-1","TMC-2","TMC-3")



####
SP_0203<-as.data.frame(read.xlsx("data/MQ/CM/CM 02-03-22 21CYGGS7318-7477as.xlsx",
                                 sheet = "SP, %",detectDates = TRUE))%>%filter()%>%select(1,2,4,5,6,7,8)%>%
  rename(SP_order=1,SP_ID=2,Date=3,nm1=4,abs_nm1=5,nm2=6,abs_nm2=7)%>%filter(!SP_order%in%c("#"))%>%mutate(abs_nm1=as.numeric(abs_nm1),abs_nm2=as.numeric(abs_nm2),delta=abs_nm2-abs_nm1)%>%mutate(Date=Date1)


s<-SP_0203%>%mutate(SP_order=1:nrow(SP_0203))
w[w$PLOT%in%c("7323"),]
s[s$SP_ID%in%c("7223"),]$SP_ID<-"7323"
s
#View(s)
s[s$SP_ID=="TMC",]
s[s$SP_ID=="TMC",]$SP_ID<-c("TMC-1","TMC-1","TMC-2","TMC-2","TMC-3","TMC-3")
s$Group<-"Group"
s[s$SP_ID%in%c("NACL", "NACL(Reblank)","0.5% NaCl"),]$SP_ID<-"NaCl"

s[1:50,]$Group<-"Group_1"
s[51:nrow(s),]$Group<-"Group_2"

value<-s%>%group_by(Group,SP_ID) %>%summarize(mean = mean(delta, na.rm = TRUE),sd=sd(delta))%>%ungroup
#View(value)
ID<-s%>%arrange%>%select(SP_ID)%>%unique()%>%mutate(SP_order=1:n_distinct(SP_ID))%>%select(2,1)
length(ID$SP_ID)
value

sp=ID%>%full_join(value,by="SP_ID","Group")%>%mutate(sd=round(sd,4))%>%rename(PLOT=SP_ID)%>%arrange(Group,SP_order)
sp$Date<-Date1
SP_0203<-sp
SP_0203
#View(s_131_12)
#sp<-plyr::rbind.fill(s_131_12,s_34)%>%arrange(Group,SP_order)%>%mutate(Date=Date1,SP_order=1:length(SP_order))

####
w$Date
SP_0203$P
cm<-w%>%full_join(SP_0203,by=c("PLOT","Date"))%>%rename(sp_mean=mean,sp_sd=sd)%>%select(-SP_order)%>%arrange(Extract_number,Group)
cm$Date<-Date1


cm
cm[cm$Group=="Group_1",]$Set<-"Set 1";cm[cm$Group=="Group_2",]$Set<-"Set 2"#;cm[cm$Group=="Group_3",]$Set<-"Set 3";cm[cm$Group=="Group_4",]$Set<-"Set 4"
#need to think of a way to deal with number ordering for NaCls
cm$PLOT

cm[1:nrow(cm),]$Treatment<-"WinterTP1-4"

cm[cm$PLOT%in%c("NaCl"),]$Entry<-"NaCl"

#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
cm$Check<-"E";cm[cm$Entry%in%c("NaCl","TMC","TMC-1","TMC-2"),]$Check<-"C";cm$Check<-as.factor(cm$Check)

cm$ID<-as.numeric(paste0(1,format(cm$Date, "%m"),format(cm$Date, "%d"),sprintf('%02d', cm$Extract_number)))
cm$PLOT
cm[cm$Check=="C",]$PLOT<-cm[cm$Check=="C",]$ID
cm$ID<-paste0(cm$ID,"-T",sapply(str_split(cm$Treatment,"T"), "[[" , 2))

CM_0203<-cm
cm
rm(cm,sp1,s,w,ID,sp,value)
#2 07 2022




CM_0207<-as.data.frame(read.xlsx("data/MQ/CM/CM 02-07-22 21CYGGW6011-6088as.xlsx",
                                 sheet = "Data",startRow = 8,detectDates = TRUE))%>%select(1:10)%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=5,PLOT=4,Date=6,ME=7,BG=8,FAN=9,sp_mean=10)%>%select(1:9)%>%
  drop_na(Date)
w<-CM_0207
w
Date1<-w$Date[3]

w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$Entry<-"TMC"
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-"TMC"
#w<-w%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))%>%arrange(Extract_number)
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-c("TMC-1")



####
SP_0207<-as.data.frame(read.xlsx("data/MQ/CM/CM 02-07-22 21CYGGW6011-6088as.xlsx",
                                 sheet = "SP, %",detectDates = TRUE))%>%filter()%>%select(1,2,4,5,6,7,8)%>%
  rename(SP_order=1,SP_ID=2,Date=3,nm1=4,abs_nm1=5,nm2=6,abs_nm2=7)%>%filter(!SP_order%in%c("#"))%>%mutate(abs_nm1=as.numeric(abs_nm1),abs_nm2=as.numeric(abs_nm2),delta=abs_nm2-abs_nm1)%>%mutate(Date=Date1)


s<-SP_0207%>%mutate(SP_order=1:nrow(SP_0207))
s
table(s$SP_ID)
#View(s)
s[s$SP_ID=="TMC",]
s[s$SP_ID=="TMC",]$SP_ID<-c("TMC-1")
s$Group<-"Group"
s[s$SP_ID%in%c("NACL", "NACL(Reblank)","0.5% NaCl","NaCl"),]$SP_ID<-"NaCl"

s[1:50,]$Group<-"Group_1"
s[51:nrow(s),]$Group<-"Group_2"

value<-s%>%group_by(Group,SP_ID) %>%summarize(mean = mean(delta, na.rm = TRUE),sd=sd(delta))%>%ungroup
#View(value)
ID<-s%>%arrange%>%select(SP_ID)%>%unique()%>%mutate(SP_order=1:n_distinct(SP_ID))%>%select(2,1)
length(ID$SP_ID)
value

sp=ID%>%full_join(value,by="SP_ID","Group")%>%mutate(sd=round(sd,4))%>%rename(PLOT=SP_ID)%>%arrange(Group,SP_order)
sp$Date<-Date1
SP_0207<-sp
SP_0207
#View(s_131_12)
#sp<-plyr::rbind.fill(s_131_12,s_34)%>%arrange(Group,SP_order)%>%mutate(Date=Date1,SP_order=1:length(SP_order))

####
w$Date
SP_0207
cm<-w%>%full_join(SP_0207,by=c("PLOT","Date"))%>%rename(sp_mean=mean,sp_sd=sd)%>%select(-SP_order)%>%arrange(Extract_number,Group)
cm$Date<-Date1



cm[cm$Group=="Group_1",]$Set<-"Set 1";cm[cm$Group=="Group_2",]$Set<-"Set 2"#;cm[cm$Group=="Group_3",]$Set<-"Set 3";cm[cm$Group=="Group_4",]$Set<-"Set 4"
#need to think of a way to deal with number ordering for NaCls
cm$PLOT

cm[1:nrow(cm),]$Treatment<-"WinterTP2-1"

cm[cm$PLOT%in%c("NaCl"),]$Entry<-"NaCl"

#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
cm$Check<-"E";cm[cm$Entry%in%c("NaCl","TMC","TMC-1","TMC-2"),]$Check<-"C";cm$Check<-as.factor(cm$Check)

cm$ID<-as.numeric(paste0(1,format(cm$Date, "%m"),format(cm$Date, "%d"),sprintf('%02d', cm$Extract_number)))
cm$PLOT
cm[cm$Check=="C",]$PLOT<-cm[cm$Check=="C",]$ID
cm$ID<-paste0(cm$ID,"-T",sapply(str_split(cm$Treatment,"T"), "[[" , 2))

CM_0207<-cm
rm(cm,sp1,s,w,ID,sp,value)






# 2 08 2022
MQ_list[MQ_list$PLOT%in%c(6379),]
CM_0208<-as.data.frame(read.xlsx("data/MQ/CM/CM 02-08-22 21CYGGW6089-6356.xlsx",
                                 sheet = "Data",startRow = 8,detectDates = TRUE))%>%select(1:10)%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=4,PLOT=5,Date=6,ME=7,BG=8,FAN=9,sp_mean=10)%>%select(1:9)%>%
  drop_na(Date)
w[w$PLOT%in%c("6286"),]#mixed up sample with SWE. Did not have enough to do analysis
w<-CM_0208%>%filter(!PLOT%in%c("6286"))
w
Date1<-w$Date[3]
w[w$Entry%in%c("6379"),c("Entry","Treatment","PLOT","Extract_number")]
MQ_list[MQ_list$PLOT%in%c(6379),]
w[w$Entry%in%c("6379"),c("Entry","Treatment","PLOT")]<-c("BS615-45","WinterTP2-2","6379")#this was an out of order sample

w[w$PLOT%in%c("6254"),c("Entry","Treatment","PLOT","Extract_number")]#missing the extract number
MQ_list[MQ_list$PLOT%in%c(6254),]
w[w$PLOT%in%c("6254"),c("Entry")]<-"BS912-128"
w[w$PLOT%in%c("6246"),]$Extract_number<-44

w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$Entry<-"TMC"
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-"TMC"
#w<-w%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))%>%arrange(Extract_number)
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-c("TMC-1","TMC-2","TMC-3","TMC-4","TMC-5")



####
SP_0208<-as.data.frame(read.xlsx("data/MQ/CM/CM 02-08-22 21CYGGW6089-6356.xlsx",
                                 sheet = "SP, %",detectDates = TRUE))%>%filter()%>%select(1,2,4,5,6,7,8)%>%
  rename(SP_order=1,SP_ID=2,Date=3,nm1=4,abs_nm1=5,nm2=6,abs_nm2=7)%>%filter(!SP_order%in%c("#"))%>%mutate(abs_nm1=as.numeric(abs_nm1),abs_nm2=as.numeric(abs_nm2),delta=abs_nm2-abs_nm1)%>%mutate(Date=Date1)


s<-SP_0208%>%mutate(SP_order=1:nrow(SP_0208))%>%filter(!SP_ID%in%c(6286))
s
table(s$SP_ID)
#View(s)
s[s$SP_ID=="TMC",]
s[s$SP_ID=="TMC",]$SP_ID<-c("TMC-1","TMC-1","TMC-2","TMC-2","TMC-3","TMC-3","TMC-4","TMC-4","TMC-5","TMC-5")
s$Group<-"Group"
s[s$SP_ID%in%c("NACL", "NACL(Reblank)","0.5% NaCl","NaCl"),]$SP_ID<-"NaCl"
#View(s)
s[1:50,]$Group<-"Group_1"
s[51:100,]$Group<-"Group_2"
s[101:nrow(s),]$Group<-"Group_3"

value<-s%>%group_by(Group,SP_ID) %>%summarize(mean = mean(delta, na.rm = TRUE),sd=sd(delta))%>%ungroup
#View(value)
ID<-s%>%arrange%>%select(SP_ID)%>%unique()%>%mutate(SP_order=1:n_distinct(SP_ID))%>%select(2,1)
length(ID$SP_ID)
value

sp=ID%>%full_join(value,by="SP_ID","Group")%>%mutate(sd=round(sd,4))%>%rename(PLOT=SP_ID)%>%arrange(Group,SP_order)
sp$Date<-Date1
SP_0208<-sp

#View(s_131_12)
#sp<-plyr::rbind.fill(s_131_12,s_34)%>%arrange(Group,SP_order)%>%mutate(Date=Date1,SP_order=1:length(SP_order))

####
w$Date

cm<-w%>%full_join(SP_0208,by=c("PLOT","Date"))%>%rename(sp_mean=mean,sp_sd=sd)%>%select(-SP_order)%>%arrange(Extract_number,Group)
cm$Date<-Date1


#View(cm)
cm[cm$Group=="Group_1",]$Set<-"Set 1";cm[cm$Group=="Group_2",]$Set<-"Set 2";cm[cm$Group=="Group_3",]$Set<-"Set 3"#;cm[cm$Group=="Group_4",]$Set<-"Set 4"
#need to think of a way to deal with number ordering for NaCls
#View(cm)
#View(cm)
cm$PLOT
cm[1:40,]$PLOT
cm[1:40,]$Treatment<-"WinterTP2-1"
cm[1:40,]$Treatment
cm[41:nrow(cm),]$Treatment<-"WinterTP2-2"
cm[cm$PLOT%in%c("6379"),]$Treatment<-"WinterTP2-2"#these were mislabels as the wrong treatment, the TP1-2 treatment for
#these plot numbers were done on 1-31-2022 if you would like to confirm
cm[cm$PLOT%in%c("NaCl"),]$Entry<-"NaCl"

#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
cm$Check<-"E";cm[cm$Entry%in%c("NaCl","TMC","TMC-1","TMC-2"),]$Check<-"C";cm$Check<-as.factor(cm$Check)

cm$ID<-as.numeric(paste0(1,format(cm$Date, "%m"),format(cm$Date, "%d"),sprintf('%02d', cm$Extract_number)))
cm$PLOT
cm[cm$Check=="C",]$PLOT<-cm[cm$Check=="C",]$ID
cm$ID<-paste0(cm$ID,"-T",sapply(str_split(cm$Treatment,"T"), "[[" , 2))

CM_0208<-cm

rm(cm,sp1,s,w,ID,sp,value)

# 02 14 2022

#MQ_list[MQ_list$PLOT%in%c(6379),]
CM_0214<-as.data.frame(read.xlsx("data/MQ/CM/CM 02-14-22 21CYGGW6380-7009.xlsx",
                                 sheet = "Data",startRow = 8,detectDates = TRUE))%>%select(1:10)%>%rename(Set=1,Extract_number=2,Entry=4,Treatment=3,PLOT=5,Date=6,ME=7,BG=8,FAN=9,sp_mean=10)%>%select(1:9)%>%
  drop_na(Date)

w<-CM_0214
w
Date1<-w$Date[3]
Date1

w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$Entry<-"TMC"
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-"TMC"
#w<-w%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))%>%arrange(Extract_number)
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-c("TMC-1")



####
SP_0214<-as.data.frame(read.xlsx("data/MQ/CM/CM 02-14-22 21CYGGW6380-7009.xlsx",
                                 sheet = "SP",detectDates = TRUE))%>%filter()%>%select(1,2,4,5,6,7,8)%>%
  rename(SP_order=1,SP_ID=2,Date=3,nm1=4,abs_nm1=5,nm2=6,abs_nm2=7)%>%filter(!SP_order%in%c("#"))%>%mutate(abs_nm1=as.numeric(abs_nm1),abs_nm2=as.numeric(abs_nm2),delta=abs_nm2-abs_nm1)%>%mutate(Date=Date1)


s<-SP_0214%>%mutate(SP_order=1:nrow(SP_0214))
s
table(s$SP_ID)
#View(s)
s[s$SP_ID=="TMC",]
s[s$SP_ID=="TMC",]$SP_ID<-c("TMC-1")
s$Group<-"Group"
s[s$SP_ID%in%c("NACL", "NACL(Reblank)","0.5% NaCl","NaCl"),]$SP_ID<-"NaCl"
#View(s)
s[1:nrow(s),]$Group<-"Group_1"
#s[1:50,]$Group<-"Group_1"
#s[51:100,]$Group<-"Group_2"
#s[101:nrow(s),]$Group<-"Group_3"

value<-s%>%group_by(Group,SP_ID) %>%summarize(mean = mean(delta, na.rm = TRUE),sd=sd(delta))%>%ungroup
#View(value)
ID<-s%>%arrange%>%select(SP_ID)%>%unique()%>%mutate(SP_order=1:n_distinct(SP_ID))%>%select(2,1)
length(ID$SP_ID)
value

sp=ID%>%full_join(value,by="SP_ID","Group")%>%mutate(sd=round(sd,4))%>%rename(PLOT=SP_ID)%>%arrange(Group,SP_order)
sp$Date<-Date1
SP_0214<-sp

#View(s_131_12)
#sp<-plyr::rbind.fill(s_131_12,s_34)%>%arrange(Group,SP_order)%>%mutate(Date=Date1,SP_order=1:length(SP_order))

####
w$Date

cm<-w%>%full_join(SP_0214,by=c("PLOT","Date"))%>%rename(sp_mean=mean,sp_sd=sd)%>%select(-SP_order)%>%arrange(Extract_number,Group)
cm$Date<-Date1

#View(cm)
#View(cm)
cm[cm$Group=="Group_1",]$Set
cm[cm$Group=="Group_1",]$Set<-"Set 1"#;cm[cm$Group=="Group_2",]$Set<-"Set 2";cm[cm$Group=="Group_3",]$Set<-"Set 3"#;cm[cm$Group=="Group_4",]$Set<-"Set 4"
#need to think of a way to deal with number ordering for NaCls
#View(cm)
cm$PLOT
#cm[1:40,]$Treatment<-"WinterTP2-1"

cm[1:nrow(cm),]$Treatment<-"WinterTP2-2"

cm[cm$PLOT%in%c("NaCl"),]$Entry<-"NaCl"

#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
cm$Check<-"E";cm[cm$Entry%in%c("NaCl","TMC","TMC-1","TMC-2"),]$Check<-"C";cm$Check<-as.factor(cm$Check)

cm$ID<-as.numeric(paste0(1,format(cm$Date, "%m"),format(cm$Date, "%d"),sprintf('%02d', cm$Extract_number)))
cm$PLOT
cm[cm$Check=="C",]$PLOT<-cm[cm$Check=="C",]$ID
cm$ID<-paste0(cm$ID,"-T",sapply(str_split(cm$Treatment,"T"), "[[" , 2))
#View(cm)
CM_0214<-cm

rm(cm,sp1,s,w,ID,sp,value)

# 02 16 2022


#MQ_list[MQ_list$PLOT%in%c(6379),]
CM_0216<-as.data.frame(read.xlsx("data/MQ/CM/CM 02-16-22 21CYGGW7010-7436.xlsx",
                                 sheet = "Data",startRow = 8,detectDates = TRUE))%>%select(1:10)%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=4,PLOT=5,Date=6,ME=7,BG=8,FAN=9,sp_mean=10)%>%select(1:9)%>%
  drop_na(Date)

w<-CM_0216
w
Date1<-w$Date[3]
Date1

w[w$PLOT%in%c("Traditional Malt Check","TMC"),]$PLOT
w[w$PLOT%in%c("Traditional Malt Check","TMC"),]$Entry<-"TMC"
w[w$PLOT%in%c("Traditional Malt Check","TMC"),]$PLOT<-"TMC"
#w<-w%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))%>%arrange(Extract_number)
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-c("TMC-1","TMC-2","TMC-3")



####
SP_0216<-as.data.frame(read.xlsx("data/MQ/CM/CM 02-16-22 21CYGGW7010-7436.xlsx",
                                 sheet = "SP, %",detectDates = TRUE))%>%filter()%>%select(1,2,4,5,6,7,8)%>%
  rename(SP_order=1,SP_ID=2,Date=3,nm1=4,abs_nm1=5,nm2=6,abs_nm2=7)%>%filter(!SP_order%in%c("#"))%>%mutate(abs_nm1=as.numeric(abs_nm1),abs_nm2=as.numeric(abs_nm2),delta=abs_nm2-abs_nm1)%>%mutate(Date=Date1)


s<-SP_0216%>%mutate(SP_order=1:nrow(SP_0216))
s
s[c(1,51),]$SP_ID<-"NaCl"

#View(s)
s[s$SP_ID=="TMC",]
s[s$SP_ID=="TMC",]$SP_ID<-c("TMC-1","TMC-1","TMC-2","TMC-2","TMC-3","TMC-3")
s$Group<-"Group"
s[s$SP_ID%in%c("NACL", "NACL(Reblank)","0.5% NaCl","NaCl"),]$SP_ID<-"NaCl"
#View(s)
s[1:50,]$Group<-"Group_1"
s[51:nrow(s),]$Group<-"Group_2"
#s[1:50,]$Group<-"Group_1"
#s[51:100,]$Group<-"Group_2"
#s[101:nrow(s),]$Group<-"Group_3"

value<-s%>%group_by(Group,SP_ID) %>%summarize(mean = mean(delta, na.rm = TRUE),sd=sd(delta))%>%ungroup
#View(value)
ID<-s%>%arrange%>%select(SP_ID)%>%unique()%>%mutate(SP_order=1:n_distinct(SP_ID))%>%select(2,1)
length(ID$SP_ID)
value

sp=ID%>%full_join(value,by="SP_ID","Group")%>%mutate(sd=round(sd,4))%>%rename(PLOT=SP_ID)%>%arrange(Group,SP_order)
sp$Date<-Date1
SP_0216<-sp

#View(s_131_12)
#sp<-plyr::rbind.fill(s_131_12,s_34)%>%arrange(Group,SP_order)%>%mutate(Date=Date1,SP_order=1:length(SP_order))

####
w$Date

cm<-w%>%full_join(SP_0216,by=c("PLOT","Date"))%>%rename(sp_mean=mean,sp_sd=sd)%>%select(-SP_order)%>%arrange(Extract_number,Group)
cm$Date<-Date1
cm
#View(cm)
#View(cm)
cm[cm$Group=="Group_1",]$Set
cm[cm$Group=="Group_1",]$Set<-"Set 1";cm[cm$Group=="Group_2",]$Set<-"Set 2"#;cm[cm$Group=="Group_3",]$Set<-"Set 3"#;cm[cm$Group=="Group_4",]$Set<-"Set 4"
#need to think of a way to deal with number ordering for NaCls
#View(cm)
cm$PLOT
cm[1:7,]$Treatment<-"WinterTP2-2"
cm[8:30,]$Treatment<-"WinterTP2-3"
cm[31:nrow(cm),]$Treatment<-"WinterTP2-4"

cm[cm$PLOT%in%c("NaCl"),]$Entry<-"NaCl"

#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
cm$Check<-"E";cm[cm$Entry%in%c("NaCl","TMC","TMC-1","TMC-2"),]$Check<-"C";cm$Check<-as.factor(cm$Check)

cm$ID<-as.numeric(paste0(1,format(cm$Date, "%m"),format(cm$Date, "%d"),sprintf('%02d', cm$Extract_number)))
cm$PLOT
cm[cm$Check=="C",]$PLOT<-cm[cm$Check=="C",]$ID
cm$ID<-paste0(cm$ID,"-T",sapply(str_split(cm$Treatment,"T"), "[[" , 2))
cm$Entry
#View(cm)
CM_0216<-cm

rm(cm,sp1,s,w,ID,sp,value)



## 02 23 2022

#MQ_list[MQ_list$PLOT%in%c(6379),]
CM_0223<-as.data.frame(read.xlsx("data/MQ/CM/CM 02-23-22 21CYGGS7099-7477as.xlsx",
                                 sheet = "Data",startRow = 8,detectDates = TRUE))%>%select(1:10)%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=4,PLOT=5,Date=6,ME=7,BG=8,FAN=9,sp_mean=10)%>%select(1:9)%>%
  drop_na(Date)

w<-CM_0223
w
Date1<-w$Date[3]
Date1

w[w$PLOT%in%c("Traditional Malt Check","TMC"),]$PLOT
w[w$PLOT%in%c("Traditional Malt Check","TMC"),]$Entry<-"TMC"
w[w$PLOT%in%c("Traditional Malt Check","TMC"),]$PLOT<-"TMC"
#w<-w%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))%>%arrange(Extract_number)
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-c("TMC-1","TMC-2")


####
SP_0223<-as.data.frame(read.xlsx("data/MQ/CM/CM 02-23-22 21CYGGS7099-7477as.xlsx",
                                 sheet = "SP, %",detectDates = TRUE))%>%filter()%>%select(1,2,4,5,6,7,8)%>%
  rename(SP_order=1,SP_ID=2,Date=3,nm1=4,abs_nm1=5,nm2=6,abs_nm2=7)%>%filter(!SP_order%in%c("#"))%>%mutate(abs_nm1=as.numeric(abs_nm1),abs_nm2=as.numeric(abs_nm2),delta=abs_nm2-abs_nm1)%>%mutate(Date=Date1)


s<-SP_0223%>%mutate(SP_order=1:nrow(SP_0223))
s

#View(s)
s[s$SP_ID=="TMC",]
s[s$SP_ID=="TMC",]$SP_ID<-c("TMC-1","TMC-1","TMC-2","TMC-2")
s$Group<-"Group"
s[s$SP_ID%in%c("NACL", "NACL(Reblank)","0.5% NaCl","NaCl"),]$SP_ID<-"NaCl"
#View(s)

s[1:nrow(s),]$Group<-"Group_1"
#s[1:50,]$Group<-"Group_1"
#s[51:100,]$Group<-"Group_2"
#s[101:nrow(s),]$Group<-"Group_3"

value<-s%>%group_by(Group,SP_ID) %>%summarize(mean = mean(delta, na.rm = TRUE),sd=sd(delta))%>%ungroup
#View(value)
ID<-s%>%arrange%>%select(SP_ID)%>%unique()%>%mutate(SP_order=1:n_distinct(SP_ID))%>%select(2,1)
length(ID$SP_ID)
value

sp=ID%>%full_join(value,by="SP_ID","Group")%>%mutate(sd=round(sd,4))%>%rename(PLOT=SP_ID)%>%arrange(Group,SP_order)
sp$Date<-Date1
SP_0223<-sp

#View(s_131_12)
#sp<-plyr::rbind.fill(s_131_12,s_34)%>%arrange(Group,SP_order)%>%mutate(Date=Date1,SP_order=1:length(SP_order))

####
w$Date

cm<-w%>%full_join(SP_0223,by=c("PLOT","Date"))%>%rename(sp_mean=mean,sp_sd=sd)%>%select(-SP_order)%>%arrange(Extract_number,Group)
cm$Date<-Date1
cm
#View(cm)
#View(cm)
cm[cm$Group=="Group_1",]$Set
cm[cm$Group=="Group_1",]$Set<-"Set 1"#;cm[cm$Group=="Group_2",]$Set<-"Set 2"#;cm[cm$Group=="Group_3",]$Set<-"Set 3"#;cm[cm$Group=="Group_4",]$Set<-"Set 4"
#need to think of a way to deal with number ordering for NaCls
#View(cm)
cm$PLOT
cm[1:9,]$PLOT
cm[1:9,]$Treatment<-"WinterTP2-3"
cm[10:nrow(cm),]$Treatment<-"WinterTP2-4"


cm[cm$PLOT%in%c("NaCl"),]$Entry<-"NaCl"

#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
cm$Check<-"E";cm[cm$Entry%in%c("NaCl","TMC","TMC-1","TMC-2"),]$Check<-"C";cm$Check<-as.factor(cm$Check)

cm$ID<-as.numeric(paste0(1,format(cm$Date, "%m"),format(cm$Date, "%d"),sprintf('%02d', cm$Extract_number)))
cm$PLOT
cm[cm$Check=="C",]$PLOT<-cm[cm$Check=="C",]$ID
cm$ID<-paste0(cm$ID,"-T",sapply(str_split(cm$Treatment,"T"), "[[" , 2))
#View(cm)
CM_0223<-cm
#View(CM_0223)
rm(cm,sp1,s,w,ID,sp,value)




# 03 01 2022

library(tibble)
#MQ_list[MQ_list$PLOT%in%c(6379),]

CM_0301<-as.data.frame(read.xlsx("data/MQ/CM/CM 03-01-22 21CYGGW7149-7219as.xlsx",
                                 sheet = "Data",startRow = 8,detectDates = TRUE))%>%select(1:10)%>%rename(Set=1,Extract_number=2,Treatment=3,PLOT=4,Date=5,ME=6,BG=7,FAN=8,sp_mean=9)%>%select(1:9)%>%
  drop_na(Date)

w<-CM_0301

Date1<-w$Date[3]
#Date1%>%add_column(dataset, .after = 2)
w$Treatment
w[w$PLOT%in%c("Traditional Malt Check","TMC"),]$PLOT
w[w$PLOT%in%c("Traditional Malt Check","TMC"),]
w[w$PLOT%in%c("Traditional Malt Check","TMC"),]$PLOT<-"TMC"


w<-w%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))%>%arrange(Extract_number)
w
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-c("TMC-1")

####
SP_0301<-as.data.frame(read.xlsx("data/MQ/CM/CM 03-01-22 21CYGGW7149-7219as.xlsx",
                                 sheet = "SP, %",detectDates = TRUE))%>%filter()%>%select(1,2,4,5,6,7,8)%>%
  rename(SP_order=1,SP_ID=2,Date=3,nm1=4,abs_nm1=5,nm2=6,abs_nm2=7)%>%filter(!SP_order%in%c("#"))%>%mutate(abs_nm1=as.numeric(abs_nm1),abs_nm2=as.numeric(abs_nm2),delta=abs_nm2-abs_nm1)%>%mutate(Date=Date1)


s<-SP_0301%>%mutate(SP_order=1:nrow(SP_0301))
s

#View(s)
s[s$SP_ID=="TMC",]
s[s$SP_ID=="TMC",]$SP_ID<-c("TMC-1","TMC-1")
s$Group<-"Group"
s[s$SP_ID%in%c("NACL", "NACL(Reblank)","0.5% NaCl","NaCl","Blank"),]$SP_ID<-"NaCl"
#View(s)

s[1:nrow(s),]$Group<-"Group_1"
#s[1:50,]$Group<-"Group_1"
#s[51:100,]$Group<-"Group_2"
#s[101:nrow(s),]$Group<-"Group_3"
s[s$SP_ID%in%c("7049-WTP2-2","7049-WTP2-3"),]$SP_ID<-"7049"
value<-s%>%group_by(Group,SP_ID) %>%summarize(mean = mean(delta, na.rm = TRUE),sd=sd(delta))%>%ungroup
#View(value)
ID<-s%>%arrange%>%select(SP_ID)%>%unique()%>%mutate(SP_order=1:n_distinct(SP_ID))%>%select(2,1)
length(ID$SP_ID)
#value

sp=ID%>%full_join(value,by="SP_ID","Group")%>%mutate(sd=round(sd,4))%>%rename(PLOT=SP_ID)%>%arrange(Group,SP_order)
sp$Date<-Date1
sp

#View(s_131_12)
#sp<-plyr::rbind.fill(s_131_12,s_34)%>%arrange(Group,SP_order)%>%mutate(Date=Date1,SP_order=1:length(SP_order))

####
w$Date
w
cm<-w%>%select(-sp_mean)%>%full_join(sp,by=c("PLOT","Date"))%>%rename(sp_mean=mean,sp_sd=sd)%>%select(-SP_order)%>%arrange(Extract_number,Group)
cm$Date<-Date1
cm
#View(cm)
#View(cm)
cm[cm$Group=="Group_1",]$Set
cm[cm$Group=="Group_1",]$Set<-"Set 1"#;cm[cm$Group=="Group_2",]$Set<-"Set 2"#;cm[cm$Group=="Group_3",]$Set<-"Set 3"#;cm[cm$Group=="Group_4",]$Set<-"Set 4"
#need to think of a way to deal with number ordering for NaCls
#View(cm)
cm$Treatment
cm$PLOT
cm[1,]
cm[1,]$Treatment<-"WinterTP2-3"
cm[2:nrow(cm),]$Treatment<-"WinterTP2-3"


cm[cm$PLOT%in%c("NaCl"),]$Entry<-"NaCl"

#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
cm$Check<-"E";cm[cm$Entry%in%c("NaCl","TMC","TMC-1","TMC-2"),]$Check<-"C";cm$Check<-as.factor(cm$Check)

cm$ID<-as.numeric(paste0(1,format(cm$Date, "%m"),format(cm$Date, "%d"),sprintf('%02d', cm$Extract_number)))
cm$PLOT
cm[cm$Check=="C",]$PLOT<-cm[cm$Check=="C",]$ID
cm$ID<-paste0(cm$ID,"-T",sapply(str_split(cm$Treatment,"T"), "[[" , 2))
#View(cm)
CM_0301<-cm
#View(CM_0223)
rm(cm,sp1,s,w,ID,sp,value)


#### 03 03 2022
# 03 01 2022

library(tibble)
#MQ_list[MQ_list$PLOT%in%c(6379),]

CM_0303<-as.data.frame(read.xlsx("data/MQ/CM/CM 03-03-22 21CYGGS7220-7476as.xlsx",
                                 sheet = "Data",startRow = 8,detectDates = TRUE))%>%select(1:10)%>%rename(Set=1,Extract_number=2,Treatment=3,PLOT=4,Date=5,ME=6,BG=7,FAN=8,sp_mean=9)%>%select(1:9)%>%
  drop_na(Date)

w<-CM_0303

Date1<-w$Date[3]
#Date1%>%add_column(dataset, .after = 2)
w$Treatment
w[w$PLOT%in%c("Traditional Malt Check","TMC"),]$PLOT
w[w$PLOT%in%c("Traditional Malt Check","TMC"),]
w[w$PLOT%in%c("Traditional Malt Check","TMC"),]$PLOT<-"TMC"

w
w<-w%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))%>%arrange(Extract_number)
w
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]
w[w$PLOT%in%c("Tradition Malt Check","TMC"),]$PLOT<-c("TMC-1")

####
SP_0303<-as.data.frame(read.xlsx("data/MQ/CM/CM 03-03-22 21CYGGS7220-7476as.xlsx",
                                 sheet = "SP, %",detectDates = TRUE))%>%filter()%>%select(1,2,4,5,6,7,8)%>%
  rename(SP_order=1,SP_ID=2,Date=3,nm1=4,abs_nm1=5,nm2=6,abs_nm2=7)%>%filter(!SP_order%in%c("#"))%>%mutate(abs_nm1=as.numeric(abs_nm1),abs_nm2=as.numeric(abs_nm2),delta=abs_nm2-abs_nm1)%>%mutate(Date=Date1)


s<-SP_0303%>%mutate(SP_order=1:nrow(SP_0303))
s

#View(s)
s[s$SP_ID=="TMC",]
s[s$SP_ID=="TMC",]$SP_ID<-c("TMC-1","TMC-1")
s$Group<-"Group"
s[s$SP_ID%in%c("NACL", "NACL(Reblank)","0.5% NaCl","NaCl","Blank"),]$SP_ID<-"NaCl"
#View(s)
s[1:nrow(s),]$Group
s[1:nrow(s),]$Group<-"Group_1"
#s[1:50,]$Group<-"Group_1"
#s[51:100,]$Group<-"Group_2"
#s[101:nrow(s),]$Group<-"Group_3"
value<-s%>%group_by(Group,SP_ID) %>%summarize(mean = mean(delta, na.rm = TRUE),sd=sd(delta))%>%ungroup
#View(value)
ID<-s%>%arrange%>%select(SP_ID)%>%unique()%>%mutate(SP_order=1:n_distinct(SP_ID))%>%select(2,1)
length(ID$SP_ID)
#value

sp=ID%>%full_join(value,by="SP_ID","Group")%>%mutate(sd=round(sd,4))%>%rename(PLOT=SP_ID)%>%arrange(Group,SP_order)
sp$Date<-Date1
sp

#View(s_131_12)
#sp<-plyr::rbind.fill(s_131_12,s_34)%>%arrange(Group,SP_order)%>%mutate(Date=Date1,SP_order=1:length(SP_order))

####
w$Date
w
cm<-w%>%select(-sp_mean)%>%full_join(sp,by=c("PLOT","Date"))%>%rename(sp_mean=mean,sp_sd=sd)%>%select(-SP_order)%>%arrange(Extract_number,Group)
cm$Date<-Date1
cm
#View(cm)
#View(cm)
cm[cm$Group=="Group_1",]$Set
cm[cm$Group=="Group_1",]$Set<-"Set 1"#;cm[cm$Group=="Group_2",]$Set<-"Set 2"#;cm[cm$Group=="Group_3",]$Set<-"Set 3"#;cm[cm$Group=="Group_4",]$Set<-"Set 4"
#need to think of a way to deal with number ordering for NaCls
#View(cm)
cm$PLOT
cm[1:9,]
cm[1:8,]$Treatment<-"WinterTP2-3"
cm[9:nrow(cm),]$Treatment<-"WinterTP2-4"


cm[cm$PLOT%in%c("NaCl"),]$Entry<-"NaCl"

#t<-t%>%mutate(Entry=plyr::mapvalues(PLOT,from=MQ_list$PLOT,to=MQ_list$Entry))
cm$Check<-"E";cm[cm$Entry%in%c("NaCl","TMC","TMC-1","TMC-2"),]$Check<-"C";cm$Check<-as.factor(cm$Check)

cm$ID<-as.numeric(paste0(1,format(cm$Date, "%m"),format(cm$Date, "%d"),sprintf('%02d', cm$Extract_number)))
cm$PLOT
cm[cm$Check=="C",]$PLOT<-cm[cm$Check=="C",]$ID
cm$ID<-paste0(cm$ID,"-T",sapply(str_split(cm$Treatment,"T"), "[[" , 2))
#View(cm)
CM_0303<-cm
#View(CM_0223)
rm(cm,sp1,s,w,ID,sp,value)

CM<-plyr::rbind.fill(CM_0128,CM_0131,CM_0201,CM_0202,CM_0203,CM_0207,CM_0208,CM_0214,CM_0216,CM_0223,CM_0301,CM_0303)
#table(CM$Entry)
#View(CM%>%filter(Treatment%in%c("WinterTP1-1","WinterTP1-2"))%>%arrange(Date,Extract_number))
#View(CM)

CM<-CM%>%filter(!is.na(Extract_number))
CM$ME<-round(as.numeric(CM$ME),4)
CM$ME
#View(CM)
CM[CM$Date%in%c(as.Date("2022-02-08"))&CM$Extract_number%in%c(72),]$Treatment<-"WinterTP2-2"
CM[CM$Date%in%c(as.Date("2022-01-31"))&CM$Extract_number%in%c(52),]$Treatment<-"WinterTP1-1"
CM[68,]$Extract_number<-0.5

CM<-CM[-275,]

#View(CM)

#View(CM)
save(CM,file="data/MQ/CM/WMB21_CM.Rdata")
table(CM$Treatment)
