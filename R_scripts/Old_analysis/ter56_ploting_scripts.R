
getwd()
library(here)
library(readxl)
load(here("Data/GGS", "GGS19_PHS_phenotypes.RData"))
#S2MET_PHS <- read_csv("Data/S2MET_data/S2MET PHS.csv")
#JIC19 <- read_excel("Data/Spring Barley Sprout 19_DS.xlsx", sheet = "Spring Obs19-JIC")
#BBBreg19 <- read_excel("Data/Spring Barley Sprout 19_DS.xlsx", sheet = "BBB19-Snyder")
str(GGS19_PHS_phenotypes)
five_avg=function(df,pheno,new_pheno,missing_phenos=F){
  df[[pheno]]=as.numeric(df[[pheno]])
  df <- cbind(df, "ID" = outer(df[[pheno]], 10^c(4:0), function(a, b) a %/% b %% 10))
  if(missing_phenos==T){
    df$ID.1[which(df$note=="*")]=NA
  }
  df$new_pheno=apply(df[,c((ncol(df)-4):ncol(df))],1,function(x) mean(x, na.rm=T))
  df$new_pheno=as.numeric(df$new_pheno)
}

#JIC19$PHS=five_avg(JIC19, 'Sprout Score', "PHS")
#BBBreg19$PHS=five_avg(BBBreg19, 'Score',"PHS")

#BBBreg19=cbind(BBBreg19[,c(4,12)], as.factor("CU2_19")); colnames(BBBreg19)=c("Entry","PHS","Trial")
#JIC19=cbind(JIC19[,c(2,11)], as.factor("JIC1_19")); colnames(JIC19)=c("Entry","PHS","Trial")
#S2MET_PHS15=cbind(S2MET_PHS[which(S2MET_PHS$year=="2015"),c(5,10)], as.factor("S2MET_15")); colnames(S2MET_PHS15)=c("Entry","PHS","Trial")
#S2MET_PHS16=cbind(S2MET_PHS[which(S2MET_PHS$year=="2016"),c(5,10)], as.factor("S2MET_16")); colnames(S2MET_PHS16)=c("Entry","PHS","Trial")
GGS19_PHS_phenotypes=cbind(GGS19_PHS_phenotypes[,c(3,9)], as.factor("CU1_19")); colnames(GGS19_PHS_phenotypes)=c("Entry","PHS","Trial")
allPHS=rbind(S2MET_PHS15,S2MET_PHS16, GGS19_PHS_phenotypes, BBBreg19,JIC19)


#plot phenotypic distributions

#CU1/GGS data
load(here::here("Data/GGS", "TP1 all data for analysis.RData"))
load(here::here("Data/GGS", "TP2 all data for analysis.RData"))
load(here::here("Data/GGS", "TP3 30k all data for analysis.RData"))

#BBBReg/CU2 data
#BBBReg_allTP <- read_excel("Data/BBBRegional/BBBReg_allTP.xlsx")

BBBReg_allTP$TimePoint=as.factor(BBBReg_allTP$TimePoint)
BBBReg_allTP$GE=apply(BBBReg_allTP[,6:10],1,function(x) {sum(as.numeric(x[1:3]))/sum(as.numeric(x[1:5]))}) 
BBBReg_allTP$GI=apply(BBBReg_allTP[,6:8],1, function(x) {(sum(as.numeric(x[1:3]))/(as.numeric(x[1])+2*as.numeric(x[2])+ 3*as.numeric(x[3])))*10})
BBBReg_allTP$GIscale=apply(BBBReg_allTP[,c(16,17)],1, function(x) {as.numeric(x[1])*as.numeric(x[2])})
colnames(BBBReg_allTP)[1]="taxa"  
BBBReg_allTP$Location=as.factor(BBBReg_allTP$Location)
BBBReg_allTP$TimePoint[162]="TP3"

#GE
str(TP1all)
tp1_ge=TP1all[,c(2,11,13)]
str(TP1all)
tp2_ge=TP2all[,c(2,11,13)]
tp3_ge=TP3all_30k[,c(2,11,13)]
#breg_ge=BBBReg_allTP[,c(1,5,16)]; colnames(breg_ge)=c("Entry", "TP","GE")

GEall=rbind(tp1_ge, tp2_ge, tp3_ge)
# GEall$Dataset=c(rep("CU1", 3*1302), rep("CU2", 640))
# GEall=GEall[-c(which(is.na(GEall$TP))),]
str(GEall)
#GI
tp1_gi=TP1all[,c(2,11,16)]
tp2_gi=TP2all[,c(2,11,16)]
tp3_gi=TP3all_30k[,c(2,11,16)]
#breg_gi=BBBReg_allTP[,c(1,5,18)]; colnames(breg_gi)=c("Entry", "TP","GIscale")

GI_all=rbind(tp1_gi, tp2_gi, tp3_gi)
str(GI_all)
# GI_all$Dataset=c(rep("CU1", 3*1302), rep("CU2", 640))
#cu2melt=melt(BBBReg_allTP, id.vars =c("taxa","TimePoint"), measure.vars = c("GE", "GIscale"))
cu1melt=rbind(melt(GEall, id.vars =c("Entry","TP"), measure.vars = c("GE")),
              melt(GI_all, id.vars =c("Entry","TP"), measure.vars = c("GIscale")))
cu1melt=cu1melt[-c(which(is.na(cu1melt$TP))),]
cu2melt$dataset="CU2"; cu1melt$dataset="CU1"; colnames(cu1melt)[1:2]=c("taxa", "TimePoint")
cumelt=rbind(cu2melt, cu1melt)
levels(cumelt$variable)= c("Germination energy", "Germination index")

ggplot(data=allPHS, aes(x=Trial, y=PHS, fill=Trial)) +
  geom_boxplot()+
  theme_bw() + #ylim(1,6.8)+
  ggtitle("Raw PHS distributions")+ theme(legend.position = "none")+
  scale_fill_manual(values=c("#f7f7f7", "#cccccc","#969696","#636363", "#252525")) +
  xlab("Dataset") + ylab("Pre-harvest sprouting")


ggplot(data=cumelt, aes(x=TimePoint, y=value, fill=dataset)) +
  geom_boxplot()+
  facet_grid(cols=vars(dataset), scales="free", rows=vars(variable)) +
  theme_bw() + 
  #ggtitle("Raw GE distributions")+
  scale_fill_manual(values=c("#cccccc", "#636363")) +
  theme(legend.position = "none")+
  xlab("Time point") + ylab("")

#make PHS, GI, GE plots by haplo here. Use consistent color scheme and background theme

#line plot for CU2/BReg19 GI
load(here("Data","GG19 allGI haplotype data.RData"))
BBBReg_allTP=merge(BBBReg_allTP, allGI[,c(1,11)], by="taxa")
BBBregmelt=melt(BBBReg_allTP, id.vars =c("taxa","haplo","TimePoint", "Location"), measure.vars = c("GIscale"))
BBBregmelt$haplo=factor(BBBregmelt$haplo, levels = c( "202", "220","222","020",  "022"))

bbbreg2=aggregate(value~ haplo +TimePoint, data=BBBregmelt, mean)
ggplot(data=bbbreg2, aes(x=TimePoint, y=value, group=haplo, color=haplo)) +
  geom_point()+
  geom_line(size=1.5) +
  theme_bw() + 
  geom_hline(yintercept =5.5, color="red", linetype=3, size=1)+
  scale_color_manual(values=c("#b2182b","#fddbc7","#d1e5f0", "#67a9cf", "#2166ac"), labels=c("DNN", "DDD","DDN","NDD", "NDN"), name="Haplotype") + 
  ylab("CU2 germination index") + xlab("Time point")

#boxplots for CU2/BBBReg19 GI averaged across both locations
ggplot(data=BBBregmelt, aes(x=haplo, y=value, fill=haplo)) +
  geom_boxplot(show.legend = F)+
  facet_wrap(~TimePoint, nrow=2) + 
  ylim(0,7) +
  theme_bw() +
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "202") +
  scale_x_discrete(breaks=c("202","220", "222","020","022"), labels=c("DNN", "DDD","DDN","NDD", "NDN")) +
  scale_fill_manual(values=c("#b2182b","#fddbc7","#d1e5f0", "#67a9cf", "#2166ac")) + 
  ylab("CU2 germination index") + xlab("Seed dormancy haplotype (AlaAt/GA20ox1/MKK3)")

#separated by locations
ggplot(data=BBBregmelt, aes(x=haplo, y=value, fill=haplo)) +
  geom_boxplot(show.legend = F)+
  facet_wrap(Location~TimePoint, nrow=2) + 
  ylim(0,7) +
  theme_bw() +
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "202") +
  scale_x_discrete(breaks=c("202","220", "222","020","022"), labels=c("DNN", "DDD","DDN","NDD", "NDN")) +
  scale_fill_manual(values=c("#b2182b","#fddbc7","#d1e5f0", "#67a9cf", "#2166ac")) + 
  ylab("CU2 germination index") + xlab("Seed dormancy haplotype (AlaAt/GA20ox1/MKK3)")

## GE for CU2/BReg
BBBregmelt2=melt(BBBReg_allTP, id.vars =c("taxa","haplo","TimePoint", "Location"), measure.vars = c("GE"))
BBBregmelt2$haplo=factor(BBBregmelt2$haplo, levels = c( "202", "220","222","020",  "022"))

ggplot(data=BBBregmelt2, aes(x=haplo, y=value, fill=haplo)) +
  geom_boxplot(show.legend = F)+
  facet_wrap(~TimePoint, nrow=2) + 
  ylim(0,1.05)+
  theme_bw() +
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "202") +
  scale_x_discrete(breaks=c("202","220", "222","020","022"), labels=c("DNN", "DDD","DDN","NDD", "NDN")) +
  scale_fill_manual(values=c("#b2182b","#fddbc7","#d1e5f0", "#67a9cf", "#2166ac")) + 
  ylab("CU2 germination energy") + xlab("Seed dormancy haplotype (AlaAt/GA20ox1/MKK3)")

## GE and GI at the same time for CU2/BReg
BBBregmelt3=melt(BBBReg_allTP, id.vars =c("taxa","haplo","TimePoint", "Location"), measure.vars = c("GIscale", "GE"))
BBBregmelt3$haplo=factor(BBBregmelt3$haplo, levels = c( "202", "220","020","222",  "022"))

ggplot(data=BBBregmelt3, aes(x=haplo, y=value, fill=haplo)) +
  geom_boxplot(show.legend = F)+
  facet_grid(variable~TimePoint, scales = "free_y") + 
  theme_bw() +
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "202") +
  scale_x_discrete(breaks=c("202","220","020", "222","022"), labels=c("DNN", "DDD","NDD","DDN", "NDN")) +
  scale_fill_manual(values=c("#b2182b","#fddbc7","#d1e5f0", "#67a9cf", "#2166ac")) + 
  xlab("Seed dormancy haplotype (AlaAt/GA20ox1/MKK3)") + ylab("")

