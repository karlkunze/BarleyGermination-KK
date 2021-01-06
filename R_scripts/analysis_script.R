
library(here)
library(lme4)
library(reshape2)
getwd()
library(here)
load(here("WMB DH 2020 all data for analysis.RData"))
WDH20_pheno<-All_pheno_TPall
table(WDH20_pheno$Entry)
All_pheno_TPall$GIscale
All_pheno_TPall$GIscale_5D
View(cbind(All_pheno_TPall$Plus_PM,All_pheno_TPall$Entry,All_pheno_TPall$GIscale,All_pheno_TPall$GIscale_5D))
str(WDH20_pheno)
#comments about the data table
#each entry is replicated 8 times to cover 4 TPs for each replication
#This means that phenotype data for each entry is replicated that many times
#This shouldn't be a problem if comparing phenotype data within a timepoint,
# but when comparing between ensure that total number of observations are correct

#####



WDH20_pheno$Plus_PM<-as.factor(WDH20_pheno$Plus_PM)
WDH20_pheno$PM<-as.factor(WDH20_pheno$PM)
WDHphs.lm=lm(phs~ Plus_PM + Entry, data=WDH20_pheno)
anova(WDHphs.lm) #all factors signif


#test for time point 1
# *note we only did half the samples for this timepoint
#all fac were retained
#Entries were scored based on ranges, not reps so Travis(or Dan) 
#scored plots 11001-11250 and I scored 11252-11543 for example
c<-WDH20_pheno[WDH20_pheno$Timepoint=="1",]

WDHge_1.lm=lm(GE~ Location + Entry + PM + UserDay1, data=WDH20_pheno[WDH20_pheno$Timepoint=="1",])
anova(WDHge_1.lm) #location, #PM significant
WDHgi_1.lm=lm(GI~ Location + Entry + PM + UserDay1, data=WDH20_pheno[WDH20_pheno$Timepoint=="1",])
anova(WDHgi_1.lm) #all sign
WDHgiScale_1.lm=lm(GIscale~ Location+ Entry + PM +UserDay1, data=WDH20_pheno[WDH20_pheno$Timepoint=="1",])
anova(WDHgiScale_1.lm) # all sing
# 5 day
WDHge5D_1.lm=lm(GE_5D~ Location + Entry + PM +UserDay1, data=WDH20_pheno[WDH20_pheno$Timepoint=="1",])
anova(WDHge5D_1.lm) #location, #PM significant

WDHgi5D_1.lm=lm(GI_5D~ Location + Entry + PM +UserDay1, data=WDH20_pheno[WDH20_pheno$Timepoint=="1",])
anova(WDHgi5D_1.lm)
WDHgiScale5D_1.lm=lm(GIscale_5D~ Location+ Entry + PM +UserDay1, data=WDH20_pheno[WDH20_pheno$Timepoint=="1",])
anova(WDHgiScale5D_1.lm)#all sign



#Time point 2
#for this and following time points the full sample amounts were used
WDHge_2.lm=lm(GE~ Location+ Entry + PM + UserDay1, data=WDH20_pheno[WDH20_pheno$Timepoint=="2",])
anova(WDHge_2.lm) #all sing

WDHgI_2.lm=lm(GI~ Location+ Entry + PM + UserDay1, data=WDH20_pheno[WDH20_pheno$Timepoint=="2",])
anova(WDHgI_2.lm) #PM:replication ns

WDHgScale_2.lm=lm(GI~ Location+ Entry + PM + UserDay1, data=WDH20_pheno[WDH20_pheno$Timepoint=="2",])
anova(WDHgScale_2.lm) #PM:replication ns
#5 Day
WDHge5D_2.lm=lm(GE_5D~ Location+ Entry + PM + UserDay4, data=WDH20_pheno[WDH20_pheno$Timepoint=="2",])
anova(WDHge5D_2.lm)  #rep ns

WDHgI5D_2.lm=lm(GI_5D~ Location+ Entry + PM + UserDay4, data=WDH20_pheno[WDH20_pheno$Timepoint=="2",])
anova(WDHgI5D_2.lm)  #rep ns

WDHgScale5D_2.lm=lm(GIscale_5D~ Location+ Entry + PM + UserDay4, data=WDH20_pheno[WDH20_pheno$Timepoint=="2",])
anova(WDHgScale_2.lm)  #rep ns


# Time point 3
#scoring was done based on reps
WDHge_3.lm=lm(GE~ Location + Entry + PM + replication, data=WDH20_pheno[WDH20_pheno$Timepoint=="3",])
anova(WDHge_3.lm) #all significant, rep is significant, could be a problem

WDHgI_3.lm=lm(GI~ Location+ Entry + PM + replication, data=WDH20_pheno[WDH20_pheno$Timepoint=="3",])
anova(WDHgI_3.lm) #all factors signif, rep highly significant

WDHgScale_3.lm=lm(GIscale~ Location+ Entry + PM + replication, data=WDH20_pheno[WDH20_pheno$Timepoint=="3",])
anova(WDHgScale_3.lm) #all factors signif
# 5 Day

WDHge5D_3.lm=lm(GE_5D~ Location + Entry + PM + replication, data=WDH20_pheno[WDH20_pheno$Timepoint=="3",])
anova(WDHge5D_3.lm) # PM ns, rep slightly significant

WDHgI5D_3.lm=lm(GI_5D~ Location+ Entry + PM + replication, data=WDH20_pheno[WDH20_pheno$Timepoint=="3",])
anova(WDHgI5D_3.lm) #rep slightly significnat

WDHgScale5D_3.lm=lm(GIscale_5D~ Location+ Entry + PM + replication, data=WDH20_pheno[WDH20_pheno$Timepoint=="3",])
anova(WDHgScale5D_3.lm) #all factors signif

#Time point 4 
WDHge_4.lm=lm(GE~ Location+ Entry + PM + replicaton, data=WDH20_pheno[WDH20_pheno$Timepoint=="4",])
anova(WDHge_4.lm) #PM:replication ns

WDHge_4.lm=lm(GI~ Location+ Entry + PM + replication, data=WDH20_pheno[WDH20_pheno$Timepoint=="4",])
anova(WDHge_4.lm)

WDHgScale_4.lm=lm(GIscale~ Location+ Entry + PM + replication, data=WDH20_pheno[WDH20_pheno$Timepoint=="4",])
anova(WDHgScale_4.lm)
# random effects

WDHphs.lmer=lmer(phs~ PM + (1|Entry), data=WDH20_pheno)
WDHphs.lmer
WDHGE_TP1.lmer=lmer(GE~ Location + PM + (1|Entry), data=WDH20_pheno[WDH20_pheno$Timepoint=="1",])
WDHGE_TP2.lmer=lmer(GE~  Location + PM+ (1|Entry) , data=WDH20_pheno[WDH20_pheno$Timepoint=="2",])
WDHGI_TP2.lmer=lmer(GIscale~  Location +PM + (1|Entry) + (1|PM:replication) , data=WDH20_pheno[WDH20_pheno$Timepoint=="2",])
WDHGE_TP3.lmer=lmer(GE~  Location + PM+ (1|Entry) , data=WDH20_pheno[WDH20_pheno$Timepoint=="3",])
WDHGI_TP3.lmer=lmer(GIscale~  Location + PM + (1|Entry) + (1|PM:replication) , data=WDH20_pheno[WDH20_pheno$Timepoint=="3",])

WDHGE_TP4.lmer=lmer(GE~  Location + PM+ (1|Entry) , data=WDH20_pheno[WDH20_pheno$Timepoint=="4",])
WDHGI_TP4.lmer=lmer(GIscale~  Location + PM + (1|Entry) + (1|PM:replication) , data=WDH20_pheno[WDH20_pheno$Timepoint=="4",])


#plotting
str(WDH20_pheno)
WDH20_sub<-WDH20_pheno[,c("Entry","Pedigree","Location","Timepoint","phs","GE","GI","GIscale","GE_5D","GI_5D","GIscale_5D")]
GEall_DH<-WDH20_sub[,c("Entry","Pedigree","Location","Timepoint","GE")]
GIall_DH<-WDH20_sub[,c("Entry","Pedigree","Location","Timepoint","GIscale")]
GE_5Dall_DH<-WDH20_sub[,c("Entry","Pedigree","Location","Timepoint","GE_5D")]
GI_5Dall_DH<-WDH20_sub[,c("Entry","Pedigree","Location","Timepoint","GIscale_5D")]
cor(GEall_DH$GE,GE_5Dall_DH$GE_5D,use = "complete.obs")
cor(GIall_DH$GE,GI_5Dall_DH$GE_5D,use = "complete.obs")
WDHmelt=rbind(melt(GEall_DH, id.vars =c("Entry","Timepoint"), measure.vars = c("GE")),
              melt(GIall_DH, id.vars =c("Entry","Timepoint"), measure.vars = c("GIscale")),
              melt(GE_5Dall_DH, id.vars =c("Entry","Timepoint"), measure.vars = c("GE_5D")),
              melt(GI_5Dall_DH, id.vars =c("Entry","Timepoint"), measure.vars = c("GIscale_5D")))
WDHmelt=WDHmelt[-c(which(is.na(WDHmelt$Timepoint))),]

WDHmelt_Sny=rbind(melt(GEall_DH[GEall_DH$Location=="Snyder",], id.vars =c("Entry","Timepoint"), measure.vars = c("GE")),
              melt(GIall_DH[GEall_DH$Location=="Snyder",], id.vars =c("Entry","Timepoint"), measure.vars = c("GIscale")),
              melt(GE_5Dall_DH[GEall_DH$Location=="Snyder",], id.vars =c("Entry","Timepoint"), measure.vars = c("GE_5D")),
              melt(GI_5Dall_DH[GEall_DH$Location=="Snyder",], id.vars =c("Entry","Timepoint"), measure.vars = c("GIscale_5D")))

WDHmelt_Ket=rbind(melt(GEall_DH[GEall_DH$Location=="Ketola 5",], id.vars =c("Entry","Timepoint"), measure.vars = c("GE")),
                  melt(GIall_DH[GIall_DH$Location=="Ketola 5",], id.vars =c("Entry","Timepoint"), measure.vars = c("GIscale")),
                  melt(GE_5Dall_DH[GE_5Dall_DH$Location=="Ketola 5",], id.vars =c("Entry","Timepoint"), measure.vars = c("GE_5D")),
                  melt(GI_5Dall_DH[GI_5Dall_DH$Location=="Ketola 5",], id.vars =c("Entry","Timepoint"), measure.vars = c("GIscale_5D")))
WDHmelt_Sny$dataset<-"Snyder"
WDHmelt_Ket$dataset<-"Ketola"
WDH_melt_L=rbind(WDHmelt_Sny, WDHmelt_Ket)
levels(WDHmelt$variable)= c("GE", "GI","5 Day GE", " 5 Day GI")
levels(WDH_melt_L$variable)= c("GE", "GI","5 Day GE","5 Day GI")

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
              melt(GE_5Dall_DH, id.vars =c("Entry","Timepoint","Pedigree"), measure.vars = c("GE_5D")),
              melt(GI_5Dall_DH, id.vars =c("Entry","Timepoint","Pedigree"), measure.vars = c("GIscale_5D")))
WDHmelt_P=WDHmelt_P[-c(which(is.na(WDHmelt_P$Timepoint))),]

WDHmelt_P[WDHmelt_P$Pedigree=="DH130910/Wintmalt",]$Pedigree<-"Wintmalt/DH130910"
ggplot(data=WDHmelt_P, aes(x=Timepoint, y=value, fill=Pedigree)) +
  geom_boxplot()+
  facet_grid(cols=vars(Pedigree), scales="free", rows=vars(variable)) +
  theme_bw() + 
  ggtitle("Germination phenotype distributions by Pedigree")+

  theme(legend.position = "none")+
  xlab("Time point") + ylab("")

