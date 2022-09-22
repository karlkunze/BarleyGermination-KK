


Sys.info()['sysname']
if( .Platform$OS.type == "unix" )
  patht <- "/home/karl/git/"
if( .Platform$OS.type == "windows" )
  patht <- "C:/Users/kars8/git/"
#ent

getwd()

patht<-paste0(patht,"NY-winter-barley-analysis/")
load("data/MQ/WMB21_Master_Cornell.Rdata")
patht
load(paste0(patht,"data/phenotypes/Field_data_2020_2022.Rdata"))#all
load(paste0(patht,"data/genotypes/wmb_GD_rrblup.Rdata"))

load(paste0(patht,"data/genotypes/GAPIT_wmb.Rdata"))
load(paste0(patht,"data/genotypes/wmb_pedigree.Rdata"))
library(readxl)
library(dplyr)

library(asreml)
library(asremlPlus)
MQ_WMB21$Env
field_data<-all_pheno%>%filter(trial=="PYT",Year=="2021")%>%dplyr::select(SourcePLOT,GID,Env,Row,Column,PHS,Scald,winter_survival,yield_kgha)%>%rename(PLOT=SourcePLOT,taxa=GID)
colnames(field_data)
str(field_data$PLOT)
field_data$PLOT<-as.character(field_data$PLOT)
#View(MQ_WMB21)

MQ_TP<-MQ_WMB21%>%dplyr::select(PLOT,Treatment,gid,Timepoint_group,Timepoint,Check.x,TB,AA,BG,ME,FAN,ST,DP,TotalProtein,AlaAT)%>%rename(taxa=gid)%>%left_join(field_data,by=c("PLOT","taxa"))%>%
  tidyr::pivot_longer(Timepoint,names_to = "name",values_to = "Timepoint")%>%tidyr::pivot_longer(cols=c("AA","BG","ME","FAN","ST","DP","TotalProtein"),names_to="trait")%>%
  arrange(Timepoint,trait,PLOT)%>%mutate_at(vars(Timepoint_group,Timepoint,Treatment,taxa,Check.x,Env,Row,Column),  list(factor))%>%rename(gid=taxa)
table(MQ_TP[MQ_TP$trait=="ME",]$Timepoint)

MQ_TP1<-MQ_TP%>%filter(Timepoint%in%c("TP1"))
MQ_TP1$trait
#BG
tp1_BG<-asreml(fixed=value~Treatment,random = ~gid,residual = ~units,data=MQ_TP1%>%filter(trait=="BG"),na.action = na.method(x = "include"))
#fixed_BG2<-asreml(fixed=BG~Treatment,random=~gid,residual = ~units,data=MQ_TP1,na.action = na.method(x = "include"))
BG_TP1<-predict(tp1_BG,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP1",trait="BG")
BG_TP1$pvals.predicted.value<-BG_TP1$pvals.predicted.value/(1-(BG_TP1$pvals.std.error/summary(tp1_BG)$varcomp["gid","component"]))
BG_TP1$pvals.predicted.value
#AA
model<-asreml(fixed=value~Treatment,random = ~gid,residual = ~units,data=MQ_TP1%>%filter(trait=="AA"),na.action = na.method(x = "include"))
model_TP1<-asreml(fixed=value~Treatment,random = ~gid,residual = ~units,data=MQ_TP1%>%filter(trait=="AA"),na.action = na.method(x = "include"))
wald.asreml(model,model_TP1)
modelAA_TP1<-predict(model_TP1,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP1",trait="AA")

#ST
modelTP1.T<-asreml(fixed=value~1,random = ~gid,residual = ~units,data=MQ_TP1%>%filter(trait=="ST"),na.action = na.method(x = "include"))
modelST_TP1<-predict(modelTP1.T,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP1",trait="ST")
modelST_TP1$pvals.predicted.value
#MP

modelTP1.P<-asreml(fixed=value~1,random = ~gid,residual = ~units,data=MQ_TP1%>%filter(trait=="TotalProtein"),na.action = na.method(x = "include"))
modelP_TP1<-predict(modelTP1.P,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP1",trait="P")
modelP_TP1$pvals.predicted.value
#DP
modelTP1.DP<-asreml(fixed=value~Treatment,random = ~gid,residual = ~units,data=MQ_TP1%>%filter(trait=="DP"),na.action = na.method(x = "include"))
summary(modelTP1.DP)$varcomp
modelDP_TP1<-predict(modelTP1.DP,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP1",trait="DP")
modelDP_TP1$pvals.predicted.value
#FAN
modelTP1.FAN<-asreml(fixed=value~Treatment,random = ~gid,residual = ~units,data=MQ_TP1%>%filter(trait=="FAN"),na.action = na.method(x = "include"))
summary(modelTP1.FAN)$varcomp
modelFAN_TP1<-predict(modelTP1.FAN,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP1",trait="FAN")
modelFAN_TP1
#ME
modelTP1.ME.0<-asreml(fixed=value~1,random = ~gid,residual = ~units,data=MQ_TP1%>%filter(trait=="ME"),na.action = na.method(x = "include"))
MQ_TP1
modelTP1.ME<-asreml(fixed=value~Treatment,random = ~gid,residual = ~units,data=MQ_TP1%>%filter(trait=="ME"),na.action = na.method(x = "include"))
wald.asreml(modelTP1.ME.0,modelTP1.ME)
modelME_TP1<-predict(modelTP1.ME,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP1",trait="ME")

TP1_MQ_values<-rbind(BG_TP1,modelAA_TP1,modelST_TP1,modelP_TP1,modelME_TP1,modelDP_TP1,modelFAN_TP1)
View(TP1_MQ_values)
#Timepoint 2

MQ_TP2<-MQ_TP%>%filter(Timepoint%in%c("TP2"))

#BG
TP2_BG<-asreml(fixed=value~Timepoint_group,random = ~gid,residual = ~units,data=MQ_TP2%>%filter(trait=="BG"),na.action = na.method(x = "include"))
summary(TP2_BG)$varcomp
BG_TP2<-predict(TP2_BG,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP2",trait="BG")
BG_TP2$pvals.predicted.value<-BG_TP2$pvals.predicted.value/(1-(BG_TP2$pvals.std.error/summary(TP2_BG)$varcomp["gid","component"]))
BG_TP2$pvals.predicted.value
#AA
model<-asreml(fixed=value~Treatment,random = ~gid,residual = ~units,data=MQ_TP2%>%filter(trait=="AA"),na.action = na.method(x = "include"))
model_TP2<-asreml(fixed=value~Treatment,random = ~gid,residual = ~units,data=MQ_TP2%>%filter(trait=="AA"),na.action = na.method(x = "include"))
wald.asreml(model,model_TP2)
modelAA_TP2<-predict(model_TP2,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP2",trait="AA")
modelAA_TP2
#ST
modelTP2.T<-asreml(fixed=value~1,random = ~gid,residual = ~units,data=MQ_TP2%>%filter(trait=="ST"),na.action = na.method(x = "include"))
modelST_TP2<-predict(modelTP2.T,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP2",trait="ST")
modelST_TP2$pvals.predicted.value
#MP

modelTP2.P<-asreml(fixed=value~1,random = ~gid,residual = ~units,data=MQ_TP2%>%filter(trait=="TotalProtein"),na.action = na.method(x = "include"))
modelP_TP2<-predict(modelTP2.P,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP2",trait="P")
modelP_TP2$pvals.predicted.value
#DP
modelTP2.DP<-asreml(fixed=value~Treatment,random = ~gid,residual = ~units,data=MQ_TP2%>%filter(trait=="DP"),na.action = na.method(x = "include"))
summary(modelTP2.DP)$varcomp
modelDP_TP2<-predict(modelTP2.DP,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP2",trait="DP")
modelDP_TP2$pvals.predicted.value
#FAN
modelTP2.FAN<-asreml(fixed=value~1,random = ~gid,residual = ~units,data=MQ_TP2%>%filter(trait=="FAN"),na.action = na.method(x = "include"))
summary(modelTP2.FAN)$varcomp
modelFAN_TP2<-predict(modelTP2.FAN,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP2",trait="FAN")
modelFAN_TP2
#ME
modelTP2.ME.0<-asreml(fixed=value~1,random = ~gid,residual = ~units,data=MQ_TP2%>%filter(trait=="ME"),na.action = na.method(x = "include"))
MQ_TP2
modelTP2.ME<-asreml(fixed=value~Treatment,random = ~gid,residual = ~units,data=MQ_TP2%>%filter(trait=="ME"),na.action = na.method(x = "include"))
wald.asreml(modelTP2.ME.0,modelTP2.ME)
modelME_TP2<-predict(modelTP2.ME.0,classify ="gid")%>%as.data.frame()%>%select(1,2,3)%>%mutate(Timepoint="TP2",trait="ME")
modelME_TP2
TP2_MQ_values<-rbind(BG_TP2,modelAA_TP2,modelST_TP2,modelP_TP2,modelME_TP2,modelDP_TP2,modelFAN_TP2)
View(TP2_MQ_values)
MQ_values<-rbind(TP1_MQ_values,TP2_MQ_values)%>%rename(GID=1,value=2)%>%mutate_at(vars(Timepoint,trait,GID),  list(factor))
View(MQ_values)
MQ_values%>%filter(trait=="DP")
save(MQ_values,file="data/MQ/MQ_values_step1.Rdata")

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
#DP
modelDP<-asreml(fixed=value~Timepoint,random = ~GID,residual = ~units,data=MQ_values[MQ_values$trait=="DP",],na.action = na.method(x = "include"))
summary(modelDP)$varcomp
modelDP<-predict(modelDP,classify = "GID")
#ST
modelST.0<-asreml(fixed=value~1,random = ~GID,residual = ~units,data=MQ_values[MQ_values$trait=="ST",],na.action = na.method(x = "include"))
modelST<-asreml(fixed=value~Timepoint,random = ~GID,residual = ~units,data=MQ_values[MQ_values$trait=="ST",],na.action = na.method(x = "include"))
summary(modelST)$varcomp
modelST<-predict(modelST,classify = "GID")
#ME
table(MQ_values$Timepoint)
MQ_values[MQ_values$trait=="ME",]
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
MQ_values[MQ_values$trait=="DP",]
MQ_values$ratio<-MQ_values$pvals.std.error/MQ_values$value
df3 <- data.frame(trait= c('AA','BG','DP','FAN','ME','P','ST'),threshold=c(0.3,0.7,0.17,0.16,0.009,0.07,0.07),stringsAsFactors = FALSE)
df3
t<-MQ_values%>%group_by(Timepoint,trait)%>%inner_join(df3,by="trait") %>% filter(ratio < threshold)%>%ungroup()

load(paste0(patht,"data/genotypes/GAPIT_wmb.Rdata"))
GD_Q<-GAPIT_wmb$GD[,c("Qsd1","taxa")]
GD_Q$taxa<-rownames(GD_Q)
colnames(GD_Q)[2]<-'GID'
GD_Q[GD_Q$GID=="DH130910",]
t %>%filter(!GID%in%c("NACL","TMC"))%>%left_join(GD_Q,by="GID")%>%filter(!Qsd1==1)%>%
  
  mutate( trait = plyr::mapvalues(trait, from = c('AA','BG','DP','FAN','ME','P','SP','ST'),
                                  to = c('Alpha-amylase','Beta-glucan','Diastatic power','Free amino nitrogen','Malt extract','Malt protein',
                                         'Soluble protein','Soluble/total protein')),Qsd1=plyr::mapvalues(Qsd1,from=c(0,2),to=c("Dormant","Nondormant")),
          Timepoint=plyr::mapvalues(Timepoint,from=c("TP1","TP2"),to=c("67 days \n after Maturity","140 days \n after Maturity"))) %>%
  ggplot(aes(x = as.factor(Timepoint), y = value,fill=Qsd1)) +
  geom_boxplot(show.legend = "none")+facet_wrap(~trait, scales = 'free')+
  
labs(x="Timepoint",y="MQ value",face="bold")+ggtitle("WMB 21 Nano Malting Distributions")+
  theme_bw() 
t %>%filter(!GID%in%c("NACL","TMC"))%>%left_join(GD_Q,by="GID")%>%filter(!Qsd1==1)%>%

  mutate( trait = plyr::mapvalues(trait, from = c('AA','BG','DP','FAN','ME','P','SP','ST'),
                            to = c('Alpha-amylase','Beta-glucan','Diastatic power','Free amino nitrogen','Malt extract','Malt protein',
                                   'Soluble protein','Soluble/total protein')),Qsd1=plyr::mapvalues(Qsd1,from=c(0,2),to=c("Dormant","Nondormant")),
                                                                                                     Timepoint=plyr::mapvalues(Timepoint,from=c("TP1","TP2"),to=c("67 days post PM \n(2 months)","140 days post PM \n(5 months)"))) %>%
filter(!trait=="Malt protein")%>%  ggplot(aes(x = as.factor(Timepoint), y = value,fill=Qsd1)) +
  geom_boxplot(show.legend = "none")+facet_wrap(~trait, scales = 'free') +labs(x="Timepoint",y="MQ value")+ggtitle("WMB 21 Nano Malting Distributions")+
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



