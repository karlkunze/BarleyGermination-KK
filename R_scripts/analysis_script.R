
library(here)
library(lme4)
library(reshape2)
getwd()

load(here( "WMB DH 2020 all data for analysis.RData"))
WDH20_pheno<-All_pheno_TPall

str(WDH20_pheno)

WDHphs.lm=lm(phs~ PM + Entry, data=WDH20_pheno)
anova(WDHphs.lm) #all factors signif

WDHge1.lm=lm(GE~ Location + Entry + PM, data=WDH20_pheno[WDH20_pheno$Timepoint=="1",])
anova(WDHge1.lm) #location, #PM significant

WDHgil.lm=lm(GIscale~ Location+ Entry + PM, data=WDH20_pheno[WDH20_pheno$Timepoint=="1",])
anova(WDHgil.lm)#location sign, PM date ns

WDHge2.lm=lm(GE~ Location+ Entry + PM + PM:replication, data=WDH20_pheno[WDH20_pheno$Timepoint=="2",])
anova(WDHge2.lm) #PM:replication ns

WDHgI2.lm=lm(GI~ Location+ Entry + PM + PM:replication, data=WDH20_pheno[WDH20_pheno$Timepoint=="2",])
anova(WDHgI2.lm) #PM:replication ns
GGge3.lm=lm(GE~ Location + Entry + PM + PM:replication, data=WDH20_pheno[WDH20_pheno$Timepoint=="3",])
anova(GGge3.lm) #PMe and PMe:rep, loc:block ns
GGgi3.lm=lm(GIscale~ Location+ Location:Block + Entry + PMe + PMe:rep, data=TP3all_30k)
anova(GGgi3.lm) #all factors signif

# random effects

WDHphs.lmer=lmer(phs~ PM + (1|Entry), data=WDH20_pheno)
WDH20_pheno
WDHGE_TP1.lmer=lmer(GE~ Location + PM + (1|Entry), data=WDH20_pheno[WDH20_pheno$Timepoint=="1",])
WDHGE_TP2.lmer=lmer(GE~  Location + PM+ (1|Entry) , data=WDH20_pheno[WDH20_pheno$Timepoint=="2",])
WDHGI_TP2.lmer=lmer(GIscale~  Location +PM + (1|Entry) + (1|PM:replication) , data=WDH20_pheno[WDH20_pheno$Timepoint=="2",])
WDHGE_TP3.lmer=lmer(GE~  Location + PM+ (1|Entry) , data=WDH20_pheno[WDH20_pheno$Timepoint=="3",])
WDHGI_TP3.lmer=lmer(GIscale~  Location + PM + (1|Entry) + (1|PM:replication) , data=WDH20_pheno[WDH20_pheno$Timepoint=="3",])
#plotting
str(WDH20_pheno)
WDH20_sub<-WDH20_pheno[,c("Entry","Pedigree","Location","Timepoint","phs","GE","GI","GIscale","GE_5D")]
GEall_DH<-WDH20_sub[,c("Entry","Pedigree","Location","Timepoint","GE")]
GIall_DH<-WDH20_sub[,c("Entry","Pedigree","Location","Timepoint","GIscale")]
GE_5Dall_DH<-WDH20_sub[,c("Entry","Pedigree","Location","Timepoint","GE_5D")]



WDHmelt=rbind(melt(GEall_DH, id.vars =c("Entry","Timepoint"), measure.vars = c("GE")),
              melt(GIall_DH, id.vars =c("Entry","Timepoint"), measure.vars = c("GIscale")),
              melt(GE_5Dall_DH, id.vars =c("Entry","Timepoint"), measure.vars = c("GE_5D")))
WDHmelt=WDHmelt[-c(which(is.na(WDHmelt$Timepoint))),]

WDHmelt_Sny=rbind(melt(GEall_DH[GEall_DH$Location=="Snyder",], id.vars =c("Entry","Timepoint"), measure.vars = c("GE")),
              melt(GIall_DH[GEall_DH$Location=="Snyder",], id.vars =c("Entry","Timepoint"), measure.vars = c("GIscale")),
              melt(GE_5Dall_DH[GEall_DH$Location=="Snyder",], id.vars =c("Entry","Timepoint"), measure.vars = c("GE_5D")))

WDHmelt_Ket=rbind(melt(GEall_DH[GEall_DH$Location=="Ketola 5",], id.vars =c("Entry","Timepoint"), measure.vars = c("GE")),
                  melt(GIall_DH[GEall_DH$Location=="Ketola 5",], id.vars =c("Entry","Timepoint"), measure.vars = c("GIscale")),
                  melt(GE_5Dall_DH[GEall_DH$Location=="Ketola 5",], id.vars =c("Entry","Timepoint"), measure.vars = c("GE_5D")))
WDHmelt_Sny$dataset<-"Snyder"
WDHmelt_Ket$dataset<-"Ketola"
WDH_melt_L=rbind(WDHmelt_Sny, WDHmelt_Ket)
levels(WDHmelt$variable)= c("GE", "GI","5 Day GE")
levels(WDH_melt_L$variable)= c("GE", "GI","5 Day GE")

library(ggplot2)
str(WDHmelt$Timepoint)

ggplot(data=WDHmelt, aes(x=Timepoint, y=value)) +
  geom_boxplot()+
  facet_grid( scales="free", rows=vars(variable)) +
  theme_bw() + 
  ggtitle("Raw Germination pheno distributions")+
  scale_fill_manual(values=c("#636363")) +
  theme(legend.position = "none")+
  xlab("Time point") + ylab("")
  
ggplot(data=WDHmelt, aes(x=Timepoint, y=value)) +
  geom_boxplot()+
  facet_grid( scales="free", rows=vars(variable)) +
  theme_bw() + 
  ggtitle("Raw Germination pheno distributions")+
  scale_fill_manual(values=c("#636363")) +
  theme(legend.position = "none")+
  xlab("Time point") + ylab("")


ggplot(data=WDH_melt_L, aes(x=Timepoint, y=value, fill=dataset)) +
  geom_boxplot()+
  facet_grid(cols=vars(dataset), scales="free", rows=vars(variable)) +
  theme_bw() + 
  ggtitle("Germination phenotype distributions by location")+
  scale_fill_manual(values=c("#cccccc", "#636363")) +
  theme(legend.position = "none")+
  xlab("Time point") + ylab("")
str(WDHmelt)

WDHmelt_P=rbind(melt(GEall_DH, id.vars =c("Entry","Timepoint","Pedigree"), measure.vars = c("GE")),
              melt(GIall_DH, id.vars =c("Entry","Timepoint","Pedigree"), measure.vars = c("GIscale")),
              melt(GE_5Dall_DH, id.vars =c("Entry","Timepoint","Pedigree"), measure.vars = c("GE_5D")))
WDHmelt_P=WDHmelt_P[-c(which(is.na(WDHmelt_P$Timepoint))),]

WDHmelt_P[WDHmelt_P$Pedigree=="DH130910/Wintmalt",]$Pedigree<-"Wintmalt/DH130910"
ggplot(data=WDHmelt_P, aes(x=Timepoint, y=value, fill=Pedigree)) +
  geom_boxplot()+
  facet_grid(cols=vars(Pedigree), scales="free", rows=vars(variable)) +
  theme_bw() + 
  ggtitle("Germination phenotype distributions by Pedigree")+

  theme(legend.position = "none")+
  xlab("Time point") + ylab("")

