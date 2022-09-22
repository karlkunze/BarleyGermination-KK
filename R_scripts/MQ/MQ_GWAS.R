

if( .Platform$OS.type == "unix" )
  patht <- "/home/karl/git/"
if( .Platform$OS.type == "windows" )
  patht <- "C:/Users/kars8/git/"
#ent
patht
getwd()
pathT
Sys.info()['sysname']
library(dplyr)
library(GAPIT3)
library(SNPRelate)
library(readxl)
library(dplyr)
load("data/MQ/WMB21_Master_Cornell.Rdata")#MQ_WMB21#load("data/MQ/Maltinq_quality_proccesed.Rdata")

patht<-paste0(patht,"NY-winter-barley-analysis/")
patht
load(paste0(patht,"data/phenotypes/Field_data_2020_2022.Rdata"))#all
load(paste0(patht,"data/genotypes/wmb_GD_rrblup.Rdata"))

load(paste0(patht,"data/genotypes/GAPIT_wmb.Rdata"))
load(paste0(patht,"data/genotypes/wmb_pedigree.Rdata"))




GD<-wmb_GD_rr$imputed;GD<-GD+1;GD<-as.data.frame(GD)
Am<-wmb_GD_rr$A
GM<-GAPIT_wmb$GM

snpgdsCreateGeno("test4.gds", genmat = as.matrix(GD), sample.id=rownames(GD), snp.id = colnames(GD), snp.chromosome = GM$Chromosome,
                 snp.position = as.integer(GM$Position), snpfirstdim = F)
(gdsfile <- snpgdsOpen("test4.gds"))
set.seed(1000)
snpset <- snpgdsLDpruning(gdsfile, method="corr", ld.threshold=0.95, slide.max.bp = 3000)
snpset_names=unlist(snpset)
GD_prune=GD[,c(1,which(colnames(GD) %in% snpset_names))]
GD_prune_gapit<-GD_prune%>%tibble::add_column(.before = 1,taxa=rownames(GD_prune))
?add_column()
GM_prune=GM[which(GM$SNP %in% colnames(GD_prune)),]
GM_prune<-GM_prune%>%arrange(Chromosome,Position)

DF_OperationsV3 = function(df){
df = df %>% arrange(P.value)  %>% mutate(logPercentileQQplot = -log10(c(1:length(df$SNP))/length(df$SNP)),
                                         rank = c(1:length(df$SNP))) %>% arrange(Chr, as.numeric(Pos)) %>%
  mutate(log10PVal = -log10(P.value),ordinal = c(1:length(df$SNP)))
return(df)
}
# #GWA_MLMM_fortidyR = function(df, groupvars) {
#   GWAS = GAPIT(Y = df %>% dplyr::select(taxa, value,Location,Year,Env) %>% as.data.frame(),GD=GD, GM=GM,PCA.total = 2,
#                Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05,h2=TRUE)
#   Out = DF_OperationsV3(GWAS$GWAS) %>% arrange(P.value) %>% slice_head(n=1000)
#   return(Out)
# }
MQ_WMB21<-MQ_TP

field_data<-all_pheno%>%filter(trial=="PYT",Year=="2021")%>%dplyr::select(SourcePLOT,GID,Env,Row,Column,PHS,Scald,winter_survival,yield_kgha)%>%rename(PLOT=SourcePLOT,taxa=GID)
colnames(field_data)
str(field_data$PLOT)
field_data$PLOT<-as.character(field_data$PLOT)
#View(MQ_WMB21)

data<-MQ_WMB21%>%dplyr::select(PLOT,Treatment,gid,Timepoint_group,Timepoint,Check.x,TB,AA,BG,ME,FAN,ST,DP,TotalProtein,AlaAT)%>%rename(taxa=gid)%>%left_join(field_data,by=c("PLOT","taxa"))%>%
  tidyr::pivot_longer(Timepoint,names_to = "name",values_to = "Timepoint")%>%tidyr::pivot_longer(cols=c("AA","BG","ME","FAN","ST","DP","TotalProtein"),names_to="trait")%>%
  arrange(Timepoint,trait,PLOT)%>%mutate_at(vars(Timepoint_group,Timepoint,taxa,Check.x,Env,Row,Column),  list(factor))
data%>%select(taxa)%>%unique()
library(asreml)
table(data$Timepoint)
#View(data)
rr

colnames(Am)
asreml_step1 <- function(d2, groupvars) {
#d2<-data%>%filter(trait=="AA",Timepoint=="TP2")
model0<-asreml(fixed=value~Timepoint_group,random=~taxa,residual=~units,data=d2,na.action = na.method(x = "include"))

model_narrow<-asreml(fixed=value~Timepoint_group,random=~vm(taxa,Am, singG = "PSD"),residual=~units,data=d2,na.action = na.method(x = "include"))
var_comp_narrow<-as.data.frame(summary(model_narrow)$varcomp)
var_comp_narrow
var_Comp<-as.data.frame(summary(model0)$varcomp)

t<-as.data.frame(table(d2$taxa));t<-t[!t$Var1%in%c("R","NACL"),];r<-mean(t$Freq)
r
var_Comp
H2<-round(var_Comp["taxa",]$component/((var_Comp["taxa",]$component+var_Comp["units!R",]$component)/r),4)
h2<-round(var_comp_narrow[1,]$component/((var_comp_narrow[1,]$component+var_comp_narrow["units!R",]$component)/r),4)
h2<-as.data.frame(h2)
H2<-as.data.frame(H2)
H2
predict
predictV<-predict(model0,classify = "taxa")

object<-as.data.frame(list(predictV,H2,h2))
object
return(object)
}


#pheno_MQ<-all_pheno%>%tidyr::pivot_longer(cols=c("yield_kgha","winter_survival","Scald","SpotBlotch","Ht","Lodging","PHS","TP_JD","HD_JD","Maturity_JD","TWT"),names_to = "trait",values_to = "value")%>%
 # rename(taxa=GID)
data$trait<-as.factor(data$trait)
data$Timepoint
step1_MQ.0<-data%>%group_by(Timepoint,trait)%>%group_modify(asreml_step1)%>%ungroup()

step1_MQ1<-step1_MQ.0%>%rename(taxa=pvals.taxa)#%>%left_join(field_data,by=c("taxa"))
asreml_step1_AlaAT <- function(d2, groupvars) {
  #d2<-data%>%filter(trait=="AA",Timepoint=="TP2")
  model0<-asreml(fixed=value~Qsd1+Timepoint:Qsd1,random=~taxa,residual=~units,data=d2,na.action = na.method(x = "include"))
  
  model_narrow<-asreml(fixed=value~Timepoint+Timepoint_group,random=~vm(taxa,Am, singG = "PSD"),residual=~units,data=d2,na.action = na.method(x = "include"))
  var_comp_narrow<-as.data.frame(summary(model_narrow)$varcomp)
  var_comp_narrow
  var_Comp<-as.data.frame(summary(model0)$varcomp)
  
  t<-as.data.frame(table(d2$taxa));t<-t[!t$Var1%in%c("R","NACL"),];r<-mean(t$Freq)
  r
  var_Comp
  H2<-round(var_Comp["taxa",]$component/((var_Comp["taxa",]$component+var_Comp["units!R",]$component)/r),4)
  h2<-round(var_comp_narrow[1,]$component/((var_comp_narrow[1,]$component+var_comp_narrow["units!R",]$component)/r),4)
  h2<-as.data.frame(h2)
  H2<-as.data.frame(H2)
  H2
  predict
  predictV<-predict(model0,classify = "taxa")
  
  object<-as.data.frame(list(predictV,H2,h2))
  object
  return(object)
}
library(asreml)
step1_MQ.AlaAT<-data%>%left_join(GD_prune_gapit[,c("taxa","Qsd1")])%>%mutate_at(vars(taxa,Qsd1),  list(factor))%>%group_by(trait)%>%group_modify(asreml_step1_AlaAT)%>%ungroup()

step1_MQ.AlaAT<-step1_MQ.AlaAT%>%rename(taxa=pvals.taxa)
table(step1_MQ.AlaAT$h2,step1_MQ.AlaAT$trait)
#cant really adjust here
#View(step1_MQ1)
#Broad_H2<-function(d2, groupvars)
#GWAS
chrTable=GM_prune%>%arrange(Chromosome,Position)
chrTable = c(as.data.frame(table(GM_prune$Chromosome))$Freq);chrLabel = c(1:7, 'UN')
winterOrdinalBreaks = c(chrTable[1]/2)
WinterChrLines = c(as.data.frame(table(GM_prune$Chromosome))$Freq[1])

for (i in 2:8){
  winterOrdinalBreaks[i] = sum(chrTable[1:i-1])+chrTable[i]/2
  WinterChrLines[i] = sum(chrTable[1:i])
}
#Sommer 

GD<-GD_prune
GM<-GM_prune
GD_gapit<-GD_prune%>%tibble::add_column(taxa=rownames(GD_prune),.before = 1)

rownames(GD)
data
pheno<-data%>%filter(taxa%in%rownames(GD))
DT <- pheno
#%>%filter(!Year=="2020",trial=="PYT")%>%mutate(Fac = ifelse(Fac_type == "Fac",1,0))%>%filter(trait%in%c("winter_survival"))
#filter(trait=="Fac_type",trial=="PYT",Env%in%c("KET-2021","WMB22 Master - Snyder"))

GT <- GD
MP <- GM
numo<-GD-1
common <- intersect(DT$taxa,rownames(numo))
table(common)
marks <- numo[common,]; marks[1:5,1:5]
dim(marks)
DT<- DT%>%filter(taxa%in%common);DT <- as.data.frame(DT)
marks<-as.matrix(marks)
A <- sommer::A.mat(marks)
dt3<-DT%>%filter(trait=="AA",Timepoint=="TP1")%>%filter(taxa%in%colnames(A))%>%droplevels()
ans2 <- GWAS(value~1,
             random=~vsr(taxa,Gu=A),
             rcov=~units,
             gTerm = "u:taxa",
             M=marks, data=dt3,getPEV=TRUE)
DF_Operations_sommer = function(df){
  df<-as.data.frame(df$scores[,1])
  df$SNP=rownames(df)
  colnames(df)[1]<-"P.value"
  df = df%>%left_join(GM,by="SNP")%>%arrange(P.value)  %>% mutate(logPercentileQQplot = -log10(c(1:length(df$SNP))/length(df$SNP)),
                                                                  rank = c(1:length(df$SNP))) %>% arrange(Chromosome, as.numeric(Position)) %>%
    mutate(log10PVal = -log10(P.value),ordinal = c(1:length(df$SNP)))
  return(df)
}
Out1 = DF_Operations_sommer(ans2)
Out1 %>%
  ggplot(aes(x = ordinal, y = log10PVal)) +
  geom_point()+  geom_vline(xintercept = 4780, color = 'red')+
  
  geom_vline(xintercept = WinterChrLines)+
  ylim(0,max(Out1$log10PVal)*2.2)+
  
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = winterOrdinalBreaks)+
  ylab('-log(p-value)')+xlab('Chromosome')+ geom_hline(yintercept = -log10(5e-5))+theme_bw()

Out1%>%filter(P.value<5e-5)
#GAPIT
library(GAPIT3)
GM_prune
DF_OperationsV3 = function(df){
  df = df %>% arrange(P.value)  %>% mutate(logPercentileQQplot = -log10(c(1:length(df$SNP))/length(df$SNP)),
                                           rank = c(1:length(df$SNP))) %>% arrange(Chromosome, Position) %>%
    mutate(log10PVal = -log10(P.value),ordinal = c(1:length(df$SNP)))
  return(df)
}
modelT<-c("MLMM","Blink", "FarmCPU")
GWA_MLMM_fortidyR = function(df, groupvars) {
  
 # df<-data0%>%filter(trait=="BG",Timepoint=="TP1",model=="FarmCPU")
   
    GWAS = GAPIT(Y = df %>% dplyr::select(taxa, value) %>% as.data.frame(),GD=GD_prune_gapit, GM=GM_prune,PCA.total = 2,
               Geno.View.output=F, model="FarmCPU", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
    DF_OperationsV3 
  str(GWAS)[4]
  GWAS=as.data.frame(GWAS$GWAS)
  return(GWAS)
}
data$model<-NA
data1=data;data2=data;data3=data
data1$model<-"MLMM"
data2$model<-"Blink"
data3$model<-"FarmCPU"
data0=rbind(data1,data2,data3)
rm(data1,data2,data3)
Gapit_results_base<-data0%>%group_by(Timepoint,trait,model)%>%group_modify(GWA_MLMM_fortidyR)

df<-Gapit_results %>%filter(Timepoint=="TP2",trait=="BG",model=="Blink")

Model_pvalues = function(df, groupvars) {
  
  df=DF_OperationsV3(df)
  return(df)
}
Gapit_results_base
object<-Gapit_results_base%>%group_by(model,Timepoint,trait)%>%group_modify(Model_pvalues)

save(object,file="data/GWA_results/MQ_results.Rdata")
ggtitle()
library(ggplot2)
object%>%filter(model=="MLMM")%>%ggplot(aes(ordinal, log10PVal, color = trait,shape=Timepoint))+geom_point(size=2)+
  geom_vline(xintercept = WinterChrLines, color = 'black')+
 geom_vline(xintercept = 4975, color = 'red')+
 annotate(geom= 'text', x = 4975, y = 15, label = 'Qsd1')+
geom_vline(xintercept = WinterChrLines)+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = winterOrdinalBreaks)+
  ggtitle(lab="Genome Wide Assocation of Malting Quality Traits")+
  ylab('-log(p-value)')+xlab('Chromosome')+ geom_hline(yintercept = -log10(5e-5)) +
  theme_bw()
table(object$model)
asreml.license.status(quiet = FALSE, task = "checkout", json = "")
#second step results
Gapit_results$model
second_step<-step1_MQ1%>%rename(value=4)%>%group_by(Timepoint,trait)%>%group_modify(GWA_MLMM_fortidyR)
second_step%>% ggplot(aes(ordinal, log10PVal, color = trait,shape=Timepoint))+geom_point(size=3)+
  geom_vline(xintercept = WinterChrLines, color = 'black')+
  geom_vline(xintercept = 4975, color = 'red')+
  annotate(geom= 'text', x = 4975, y = 15, label = 'AlaAT1')+
  geom_vline(xintercept = WinterChrLines)+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = winterOrdinalBreaks)+
  ylab('-log(p-value)')+xlab('Chromosome')+ geom_hline(yintercept = -log10(5e-5)) +
  theme_bw()
step1_MQ.AlaAT
AlAT_second_step<-step1_MQ.AlaAT%>%rename(value=3)%>%group_by(trait)%>%group_modify(GWA_MLMM_fortidyR)

AlAT_second_step%>% ggplot(aes(ordinal, log10PVal, color = trait))+geom_point(size=1)+
  geom_vline(xintercept = WinterChrLines, color = 'black')+
  geom_vline(xintercept = 4975, color = 'red')+
  annotate(geom= 'text', x = 4975, y = 15, label = 'AlaAT1')+
  geom_vline(xintercept = WinterChrLines)+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = winterOrdinalBreaks)+
  ylab('-log(p-value)')+xlab('Chromosome')+ geom_hline(yintercept = -log10(5e-5)) +
  theme_bw()
