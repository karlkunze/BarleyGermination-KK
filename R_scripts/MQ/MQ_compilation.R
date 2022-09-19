load("data/MQ/CM/WMB21_CM.Rdata")
load("data/MQ/SWE/WMB21_SWE.Rdata")
load("data/MQ/WMB21_Moistures.Rdata")
load("data/MQ/Leco/Leco.Rdata")
library(dplyr)
##linux
if( .Platform$OS.type == "unix" )
  pathT <- "/home/karl/git/"
if( .Platform$OS.type == "windows" )
  pathT <- "windows"

#some last minute formatting, this change was to accomodate an out of order TMC check in TP1-1 of the the Congress mash
CM[CM$Date=="2022-01-31"&CM$Extract_number==52&CM$Entry=="TMC",c("Treatment","Extract_number")]<-c("WinterTP1-1",0.5)




Experiment_SWE<-SWE%>%filter(!Entry=="TMC")
Experiment_CM<-CM%>%filter(!Entry=="TMC")
Exp_Moisture<-WMB21_Moisture%>%rename(Treatment=Treatment_ID)%>%filter(!PLOT=="TMC")%>%arrange(Date,Order)%>%group_by(PLOT)%>%arrange(Order)%>%mutate(Date_Moisture=Date,Order_Moisture=Order)
Exp_Leco<-Leco%>%filter(!PLOT=="TMC")%>%arrange(Treatment,PLOT)



Check_SWE<-SWE%>%filter(!Treatment=="SpringTP2-4")%>%filter(Entry=="TMC")%>%arrange(Entry,Date,Extract_number)%>%group_by(Treatment)%>%arrange(Date)%>%mutate(numbering = row_number())
Check_CM<-CM%>%filter(Entry=="TMC")%>%arrange(Entry,Date,Extract_number)%>%group_by(Treatment)%>%arrange(Date)%>%mutate(numbering = row_number())%>%ungroup()
Check_Moisture<-WMB21_Moisture%>%rename(Treatment=Treatment_ID)%>%filter(PLOT=="TMC")%>%arrange(Date,Order)%>%group_by(Treatment,PLOT)%>%arrange(Order)%>%
  mutate(numbering=row_number())%>%rename(Entry=PLOT,Date_Moisture=Date,Order_Moisture=Order)
Check_Leco<-Leco%>%filter(PLOT=="TMC")%>%arrange(Treatment,TB)%>%group_by(Treatment)%>%arrange(Treatment,TB)%>%mutate(numbering=row_number())%>%ungroup()%>%dplyr::select(!PLOT)
Check_Leco
#View(Check_Leco)
check<-full_join(Check_CM,Check_SWE,by=c("Entry","Treatment","numbering"))%>%
  full_join(Check_Moisture,by=c("Entry","Treatment","numbering"))%>%full_join(Check_Leco,by=c("Entry","Treatment","numbering"))%>%
  mutate(PLOT=paste0(PLOT.x,"-",PLOT.y))%>%dplyr::select(!c(PLOT.x,PLOT.y))%>%filter(!Treatment=="SpringTP2-4")
#View(check)



Exper<-full_join(Experiment_CM,Experiment_SWE,by=c("Treatment","PLOT","Entry"))%>%full_join(Exp_Moisture,by=c("Treatment","PLOT"))%>%full_join(Exp_Leco,by=c("Treatment","PLOT","Entry"))%>%
  mutate(numbering=1)#%>%filter(!Treatment=="Rahr")

colnames(Exper)
colnames(check)[!colnames(check)%in%colnames(Exper)]
colnames(Exper)[!colnames(Exper)%in%colnames(Exper)]

MQ_WMB21<-Exper%>%plyr::rbind.fill(check)%>%arrange(Treatment,TB)


#
#load("data/MQ/WMB21_Master_Cornell.Rdata")#MQ_WMB21
library(readxl)
if( .Platform$OS.type == "unix" )
  pathT <- "/home/karl/git/"
if( .Platform$OS.type == "windows" )
  pathT <- "windows"

AlaT_markers<-read_excel(paste0(pathT,"NY-winter-barley-analysis/data/genotypes/WMB_DH_AlaAT_KASP_rawdata.xlsx"),sheet = 'Results')%>%
  dplyr::select("Sample Name",Allele_bp)%>%rename(GID=1,AlaAT=2)%>%mutate(GID=toupper(GID))%>%filter(!GID%in%c("UNKNOWN-OMIT"))%>%tidyr::drop_na(AlaAT)



MQ_WMB21$Entry

MQ_WMB21<-MQ_WMB21%>%rename(gid=Entry,N_percent=NitrogenPercent)%>%mutate(MP=N_percent*6.25/(1-Moisture_percent/100),sp_mean=abs(sp_mean),ST=sp_mean/MP)

MQ_WMB21$Timepoint<-substr(MQ_WMB21$Treatment, start = 9, stop =9)
MQ_WMB21$Timepoint_group<-substr(MQ_WMB21$Treatment, start = 11, stop =11)

MQ_WMB21[MQ_WMB21$Timepoint%in%c("1"),]$Timepoint<-"TP1"
MQ_WMB21[MQ_WMB21$Timepoint%in%c("2"),]$Timepoint<-"TP2"

MQ_WMB21[MQ_WMB21$gid=="TMC",]$PLOT<-paste0("TMC-",MQ_WMB21[MQ_WMB21$gid=="TMC",]$Timepoint_group,"-",MQ_WMB21[MQ_WMB21$gid=="TMC",]$TB)
AlaT_markers=AlaT_markers%>%rename(gid=GID)
AlaT_markers
load("data/Genotype_data/AlaT_markers.Rdata")

MQ_WMB21<-MQ_WMB21%>%left_join(AlaT_markers,by="gid")%>%mutate_at(vars(Timepoint_group,Timepoint,gid,AlaAT),  list(factor))
library(lme4)
colnames(MQ_WMB21)


#MQ_WMB21<-MQ_WMB21%>%left_join(AlaAT,by="gid")
#MQ_WMB21$AlaT_Allele<-as.factor(MQ_WMB21$AlaT_Allele)


MQ_WMB21$gid<-toupper(MQ_WMB21$gid)
MQ_WMB21[MQ_WMB21$gid=="KWS SCALA",]$gid<-"SCALA"
MQ_WMB21[MQ_WMB21$gid=="ENDEAVOR",]$gid<-"ENDEAVOR"
MQ_WMB21$TotalProtein<-MQ_WMB21$N_percent*6.25
#Split up into timepoints

save(MQ_WMB21,file="data/MQ/WMB21_Master_Cornell.Rdata")#complete form



