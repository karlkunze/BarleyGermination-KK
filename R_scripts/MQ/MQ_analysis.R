load("data/MQ/WMB21_Master_Cornell.Rdata")#MQ_WMB21
library(dplyr)


MQ_WMB21<-MQ_WMB21%>%rename(gid=Entry,N_percent=NitrogenPercent)%>%mutate(MP=N_percent*6.25/(1-Moisture_percent/100),sp_mean=abs(sp_mean),ST=sp_mean/MP)

MQ_WMB21$Timepoint<-substr(MQ_WMB21$Treatment, start = 9, stop =9)
MQ_WMB21$Timepoint_group<-substr(MQ_WMB21$Treatment, start = 11, stop =11)

MQ_WMB21[MQ_WMB21$Timepoint%in%c("1"),]$Timepoint<-"TP1"
MQ_WMB21[MQ_WMB21$Timepoint%in%c("2"),]$Timepoint<-"TP2"

MQ_WMB21[MQ_WMB21$gid=="TMC",]$PLOT<-paste0("TMC-",MQ_WMB21[MQ_WMB21$gid=="TMC",]$Timepoint_group,"-",MQ_WMB21[MQ_WMB21$gid=="TMC",]$TB)
MQ_WMB21<-MQ_WMB21%>%mutate_at(vars(Timepoint_group,Timepoint,
                                    gid),  list(factor))
library(lme4)
colnames(MQ_WMB21)
load("data/Genotype_data/AlaT_markers.Rdata")
AlaAT<-AlaT_markers%>%as.data.frame()%>%select(Entry,AlaT_Allele)%>%rename(gid=Entry)
MQ_WMB21<-MQ_WMB21%>%left_join(AlaAT,by="gid")
MQ_WMB21$AlaT_Allele<-as.factor(MQ_WMB21$AlaT_Allele)
#using TMC to check variance

#checking significance for traits
#Beta glucan
#View(MQ_WMB21)
summary(lm(ST~AlaT_Allele,data=MQ_WMB21))
summary(lm(BG~gid+Timepoint+Timepoint:Timepoint_group,data = MQ_WMB21))
#significant effects: Timepoint 1:Group 3 & 4 and Timepoint
summary(lm(DP~Check.x+Timepoint+gid,data=MQ_WMB21))
#Timepoint and Check are significant
summary(lm(AA~Check.x+Timepoint+Timepoint:Timepoint_group+gid,data=MQ_WMB21))
#Timepoint and Check are significant, Timepoint1:Group 3 4 and Timepoint 2:4 are significant
summary(lm(ST~gid+Check.x+Timepoint+Timepoint:Timepoint_group,data=MQ_WMB21))
#only timepoint2 group 3 and 4 are significant
summary(lm(FAN~gid+Check.y+Timepoint+Timepoint:Timepoint_group,data=MQ_WMB21))
#Timepoint is significant, so is check and TP1:2, TP1:4, TP2:4
summary(lm(ME~gid+Check.y+Timepoint+Timepoint:Timepoint_group,data=MQ_WMB21))
#Timepoint is significant, TP2:3 and TP2:4 is signficant


library(asreml)
MQ_WMB21$gid<-as.factor(MQ_WMB21$gid)
mod<-asreml(fixed=BG~Timepoint+Timepoint:Timepoint_group,random=~gid,residual=~units,data=MQ_WMB21)
summary(mod)$varcomp[1,1]/(summary(mod)$varcomp["units!R",1]/2+summary(mod)$varcomp[1,1])

mod<-asreml(fixed=AA~Timepoint+Timepoint:Timepoint_group,random=~gid,residual=~units,data=MQ_WMB21)
summary(mod)$varcomp[1,1]/(summary(mod)$varcomp["units!R",1]/2+summary(mod)$varcomp[1,1])

mod<-asreml(fixed=DP~Timepoint+Timepoint:Timepoint_group,random=~gid,residual=~units,data=MQ_WMB21)
summary(mod)$varcomp[1,1]/(summary(mod)$varcomp["units!R",1]/2+summary(mod)$varcomp[1,1])

mod<-asreml(fixed=FAN~Timepoint+Timepoint:Timepoint_group,random=~gid,residual=~units,data=MQ_WMB21)
summary(mod)$varcomp[1,1]/(summary(mod)$varcomp["units!R",1]/2+summary(mod)$varcomp[1,1])

mod<-asreml(fixed=ME~1,random=~gid,residual=~units,data=MQ_WMB21)
summary(mod)$varcomp[1,1]/(summary(mod)$varcomp["units!R",1]/2+summary(mod)$varcomp[1,1])
mod<-asreml(fixed=ST~1,random=~gid,residual=~units,data=MQ_WMB21)
summary(mod)$varcomp[1,1]/(summary(mod)$varcomp["units!R",1]/2+summary(mod)$varcomp[1,1])

save(MQ_WMB21,file="data/MQ/Maltinq_quality_proccesed.Rdata")
