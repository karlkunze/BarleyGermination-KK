load("data/MQ/WMB21_Master_Cornell.Rdata")#MQ_WMB21
getwd()
patht<-"/home/karl/git/NY-winter-barley-analysis/"
load(paste0(patht,"data/phenotypes/Field_data_2020_2022.Rdata"))#all
load(paste0(patht,"data/genotypes/wmb_GD_rrblup.Rdata"))

load(paste0(patht,"data/genotypes/GAPIT_wmb.Rdata"))
load(paste0(patht,"data/genotypes/wmb_pedigree.Rdata"))
library(readxl)
AlaT_markers<-read_excel(paste0(patht,"data/genotypes/WMB_DH_AlaAT_KASP_rawdata.xlsx"),sheet = 'Results')%>%
  select("Sample Name",Allele_bp)%>%rename(GID=1,AlaAT=2)%>%mutate(GID=toupper(GID))%>%filter(!GID%in%c("UNKNOWN-OMIT"))%>%tidyr::drop_na(AlaAT)

colnames(im)
im$Qs
library(dplyr)
MQ_WMB21$Entry

MQ_WMB21<-MQ_WMB21%>%rename(gid=Entry,N_percent=NitrogenPercent)%>%mutate(MP=N_percent*6.25/(1-Moisture_percent/100),sp_mean=abs(sp_mean),ST=sp_mean/MP)

MQ_WMB21$Timepoint<-substr(MQ_WMB21$Treatment, start = 9, stop =9)
MQ_WMB21$Timepoint_group<-substr(MQ_WMB21$Treatment, start = 11, stop =11)

MQ_WMB21[MQ_WMB21$Timepoint%in%c("1"),]$Timepoint<-"TP1"
MQ_WMB21[MQ_WMB21$Timepoint%in%c("2"),]$Timepoint<-"TP2"

MQ_WMB21[MQ_WMB21$gid=="TMC",]$PLOT<-paste0("TMC-",MQ_WMB21[MQ_WMB21$gid=="TMC",]$Timepoint_group,"-",MQ_WMB21[MQ_WMB21$gid=="TMC",]$TB)
MQ_WMB21<-MQ_WMB21%>%mutate_at(vars(Timepoint_group,Timepoint,gid),  list(factor))
library(lme4)
colnames(MQ_WMB21)
load("data/Genotype_data/AlaT_markers.Rdata")


#MQ_WMB21<-MQ_WMB21%>%left_join(AlaAT,by="gid")
#MQ_WMB21$AlaT_Allele<-as.factor(MQ_WMB21$AlaT_Allele)


MQ_WMB21$gid<-toupper(MQ_WMB21$gid)
MQ_WMB21[MQ_WMB21$gid=="KWS SCALA",]$gid<-"SCALA"
MQ_WMB21[MQ_WMB21$gid=="ENDEAVOR",]$gid<-"ENDEAVOR"
MQ_WMB21$TotalProtein<-MQ_WMB21$N_percent*6.25
#Split up into timepoints

MQ_TP<-MQ_WMB21%>%as.data.frame()%>%
  mutate_at(vars(gid,Treatment,Timepoint,PLOT,Timepoint),  list(factor))
save(MQ_TP,file="data/MQ/Maltinq_quality_proccesed.Rdata")
library(asreml)
library(asremlPlus)
MQ_TP1<-MQ_TP%>%filter(Timepoint%in%c("TP1"))

#BG
tp1_BG<-asreml(fixed=BG~Treatment,random = ~gid,residual = ~units,data=MQ_TP1,na.action = na.method(x = "include"))
#fixed_BG2<-asreml(fixed=BG~Treatment,random=~gid,residual = ~units,data=MQ_TP1,na.action = na.method(x = "include"))
BG_TP1<-predict(tp1_BG,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP1",trait="BG")
BG_TP1$pvals.predicted.value<-BG_TP1$pvals.predicted.value/(1-(BG_TP1$pvals.std.error/summary(tp1_BG)$varcomp["gid","component"]))
BG_TP1$pvals.predicted.value
#AA
model<-asreml(fixed=AA~Treatment,random = ~gid,residual = ~units,data=MQ_TP1,na.action = na.method(x = "include"))
model_TP1<-asreml(fixed=AA~Treatment,random = ~gid,residual = ~units,data=MQ_TP1,na.action = na.method(x = "include"))
wald.asreml(model,model_TP1)
modelAA_TP1<-predict(model_TP1,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP1",trait="AA")

#ST
modelTP1.T<-asreml(fixed=ST~1,random = ~gid,residual = ~units,data=MQ_TP1,na.action = na.method(x = "include"))
modelST_TP1<-predict(modelTP1.T,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP1",trait="ST")
modelST_TP1$pvals.predicted.value
#MP

modelTP1.P<-asreml(fixed=TotalProtein~1,random = ~gid,residual = ~units,data=MQ_TP1,na.action = na.method(x = "include"))
modelP_TP1<-predict(modelTP1.P,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP1",trait="P")
modelP_TP1$pvals.predicted.value
#DP
modelTP1.DP<-asreml(fixed=DP~Treatment,random = ~gid,residual = ~units,data=MQ_TP1,na.action = na.method(x = "include"))
summary(modelTP1.DP)$varcomp
modelDP_TP1<-predict(modelTP1.DP,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP1",trait="DP")
modelDP_TP1$pvals.predicted.value
#FAN
modelTP1.FAN<-asreml(fixed=FAN~Treatment,random = ~gid,residual = ~units,data=MQ_TP1,na.action = na.method(x = "include"))
summary(modelTP1.FAN)$varcomp
modelFAN_TP1<-predict(modelTP1.FAN,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP1",trait="FAN")
#ME
modelTP1.ME.0<-asreml(fixed=ME~1,random = ~gid,residual = ~units,data=MQ_TP1,na.action = na.method(x = "include"))
MQ_TP1
modelTP1.ME<-asreml(fixed=ME~Treatment,random = ~gid,residual = ~units,data=MQ_TP1,na.action = na.method(x = "include"))
wald.asreml(modelTP1.ME.0,modelTP1.ME)
modelME_TP1<-predict(modelTP1.ME,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP1",trait="ME")

TP1_MQ_values<-rbind(BG_TP1,modelAA_TP1,modelST_TP1,modelP_TP1,modelME_TP1,modelDP_TP1,modelFAN_TP1)
#Timepoint 2

MQ_TP2<-MQ_TP%>%filter(Timepoint%in%c("TP2"))

#BG
TP2_BG<-asreml(fixed=BG~Timepoint_group,random = ~gid,residual = ~units,data=MQ_TP2,na.action = na.method(x = "include"))
summary(TP2_BG)$varcomp
BG_TP2<-predict(TP2_BG,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP2",trait="BG")
BG_TP2$pvals.predicted.value<-BG_TP2$pvals.predicted.value/(1-(BG_TP2$pvals.std.error/summary(TP2_BG)$varcomp["gid","component"]))
BG_TP2$pvals.predicted.value
#AA
model<-asreml(fixed=AA~Treatment,random = ~gid,residual = ~units,data=MQ_TP2,na.action = na.method(x = "include"))
model_TP2<-asreml(fixed=AA~Treatment,random = ~gid,residual = ~units,data=MQ_TP2,na.action = na.method(x = "include"))
wald.asreml(model,model_TP2)
modelAA_TP2<-predict(model_TP2,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP2",trait="AA")
modelAA_TP2
#ST
modelTP2.T<-asreml(fixed=ST~1,random = ~gid,residual = ~units,data=MQ_TP2,na.action = na.method(x = "include"))
modelST_TP2<-predict(modelTP2.T,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP2",trait="ST")
modelST_TP2$pvals.predicted.value
#MP

modelTP2.P<-asreml(fixed=TotalProtein~1,random = ~gid,residual = ~units,data=MQ_TP2,na.action = na.method(x = "include"))
modelP_TP2<-predict(modelTP2.P,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP2",trait="P")
modelP_TP2$pvals.predicted.value
#DP
modelTP2.DP<-asreml(fixed=DP~Treatment,random = ~gid,residual = ~units,data=MQ_TP2,na.action = na.method(x = "include"))
summary(modelTP2.DP)$varcomp
modelDP_TP2<-predict(modelTP2.DP,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP2",trait="DP")
modelDP_TP2$pvals.predicted.value
#FAN
modelTP2.FAN<-asreml(fixed=FAN~Treatment,random = ~gid,residual = ~units,data=MQ_TP2,na.action = na.method(x = "include"))
summary(modelTP2.FAN)$varcomp
modelFAN_TP2<-predict(modelTP2.FAN,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP2",trait="FAN")
#ME
modelTP2.ME.0<-asreml(fixed=ME~1,random = ~gid,residual = ~units,data=MQ_TP2,na.action = na.method(x = "include"))
MQ_TP2
modelTP2.ME<-asreml(fixed=ME~Treatment,random = ~gid,residual = ~units,data=MQ_TP2,na.action = na.method(x = "include"))
wald.asreml(modelTP2.ME.0,modelTP2.ME)
modelME_TP2<-predict(modelTP2.ME,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP2",trait="ME")

TP2_MQ_values<-rbind(BG_TP2,modelAA_TP2,modelST_TP2,modelP_TP2,modelME_TP2,modelDP_TP2,modelFAN_TP2)

MQ_values<-rbind(TP1_MQ_values,TP2_MQ_values)%>%rename(GID=1,value=2)%>%mutate_at(vars(Timepoint,trait,GID),  list(factor))


#now analyze across both
#BG
modelBG<-asreml(fixed=value~Timepoint,random = ~GID,residual = ~units,data=MQ_values[MQ_values$trait=="BG",],na.action = na.method(x = "include"))
summary(modelBG)$varcomp
modelBG<-predict(modelBG,classify = "GID")
modelBG
#AA
modelAA<-asreml(fixed=value~Timepoint,random = ~GID,residual = ~units,data=MQ_values[MQ_values$trait=="AA",],na.action = na.method(x = "include"))
summary(modelAA)$varcomp
modelAA<-predict(modelAA,classify = "GID")
table(MQ_values$trait)
#P
modelP<-asreml(fixed=value~Timepoint,random = ~GID,residual = ~units,data=MQ_values[MQ_values$trait=="P",],na.action = na.method(x = "include"))
summary(modelP)$varcomp
modelP<-predict(modelP,classify = "GID")
#ST
modelST.0<-asreml(fixed=value~1,random = ~GID,residual = ~units,data=MQ_values[MQ_values$trait=="ST",],na.action = na.method(x = "include"))
modelST<-asreml(fixed=value~Timepoint,random = ~GID,residual = ~units,data=MQ_values[MQ_values$trait=="ST",],na.action = na.method(x = "include"))
summary(modelST)$varcomp
modelST<-predict(modelST,classify = "GID")
#ME
modelME.0<-asreml(fixed=value~1,random = ~GID,residual = ~units,data=MQ_values[MQ_values$trait=="ME",],na.action = na.method(x = "include"))
modelME<-asreml(fixed=value~Timepoint,random = ~GID,residual = ~units,data=MQ_values[MQ_values$trait=="ME",],na.action = na.method(x = "include"))
summary(modelME)$varcomp
modelME<-predict(modelME,classify = "GID")
modelME
#FAN
modelFAN.0<-asreml(fixed=value~1,random = ~GID,residual = ~units,data=MQ_values[MQ_values$trait=="FAN",],na.action = na.method(x = "include"))
modelFAN<-asreml(fixed=value~Timepoint,random = ~GID,residual = ~units,data=MQ_values[MQ_values$trait=="FAN",],na.action = na.method(x = "include"))
summary(modelFAN)$varcomp
modelME<-predict(modelFAN,classify = "GID")
modelME

#We need to graph this a bit
#per timepoint
library(ggplot2)
table(MQ_values$GID)
all_pheno
AlaT_markers$GID
AlaT_markers$AlaAT
AlaT_markers$AlaAT
AlaT_markers
AlaT_markers%>%filter(AlaAT%in%c("C","G"))
MQ_values$ratio<-MQ_values$pvals.std.error/MQ_values$value
df3 <- data.frame(trait= c('AA','BG','DP','FAN','ME','P','ST'),threshold=c(0.3,0.7,0.17,0.16,0.009,0.07,0.07),stringsAsFactors = FALSE)

t<-MQ_values%>%group_by(Timepoint,trait)%>%inner_join(df3,by="trait") %>% filter(ratio < threshold)%>%ungroup()



t %>%filter(!GID%in%c("NACL","TMC"))%>%left_join(AlaT_markers,by="GID")%>%filter(AlaAT%in%c("C","G"))%>%

  mutate( trait = plyr::mapvalues(trait, from = c('AA','BG','DP','FAN','ME','P','SP','ST'),
                            to = c('Alpha-amylase','Beta-glucan','Diastatic power','Free amino nitrogen','Malt extract','Malt protein',
                                   'Soluble protein','Soluble/total protein')),AlaAT=plyr::mapvalues(AlaAT,from=c("C","G"),to=c("Dormant","Nondormant")),
                                                                                                     Timepoint=plyr::mapvalues(Timepoint,from=c("TP1","TP2"),to=c("2 months","5 months"))) %>%
  ggplot(aes(x = as.factor(Timepoint), y = value, fill = AlaAT)) +
  geom_boxplot()+facet_wrap(~trait, scales = 'free') +labs(x="Timepoint",y="MQ value")+ggtitle("Winter Malting Barley 2021 Nano malts")+
  theme_bw() 




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



