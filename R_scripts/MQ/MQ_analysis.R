load("data/MQ/WMB21_Master_Cornell.Rdata")#MQ_WMB21



MQ_WMB21%>%filter(is.na())
MQ_WMB21<-MQ_WMB21%>%rename(N_percent=NitrogenPercent)%>%mutate(MP=N_percent*6.25/(1-Moisture_percent/100),sp_mean=abs(sp_mean),ST=sp_mean/MP)

MQ_WMB21$Timepoint<-substr(MQ_WMB21$Treatment, start = 9, stop =9)
MQ_WMB21$Timepoint_group<-substr(MQ_WMB21$Treatment, start = 11, stop =11)
MQ_WMB21[MQ_WMB21$Timepoint%in%c("1"),]$Timepoint<-"TP1"
MQ_WMB21[MQ_WMB21$Timepoint%in%c("2"),]$Timepoint<-"TP2"
MQ_WMB21$Timepoint<-as.factor(MQ_WMB21$Timepoint)
View(MQ_WMB21)
#%>%mutate(MaltProtein=NitrogenPercent*6.25,ST=MaltProtein/Moisture_percent)
# MQ$sp_mean