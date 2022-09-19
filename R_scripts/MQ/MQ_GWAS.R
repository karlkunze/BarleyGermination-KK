

if( .Platform$OS.type == "unix" )
  pathT <- "/home/karl/git/"
if( .Platform$OS.type == "windows" )
  pathT <- "windows"
#end
pathT

pathT
Sys.info()['sysname']
library(dplyr)
library(GAPIT3)
library(SNPRelate)
library(readxl)
library(dplyr)
load("data/MQ/WMB21_Master_Cornell.Rdata")#MQ_WMB21
getwd()

patht<-paste0(pathT,"NY-winter-barley-analysis/")
patht
load(paste0(patht,"data/phenotypes/Field_data_2020_2022.Rdata"))#all
load(paste0(patht,"data/genotypes/wmb_GD_rrblup.Rdata"))

load(paste0(patht,"data/genotypes/GAPIT_wmb.Rdata"))
load(paste0(patht,"data/genotypes/wmb_pedigree.Rdata"))

AlaT_markers<-read_excel(paste0(patht,"data/genotypes/WMB_DH_AlaAT_KASP_rawdata.xlsx"),sheet = 'Results')%>%
  select("Sample Name",Allele_bp)%>%rename(GID=1,AlaAT=2)%>%mutate(GID=toupper(GID))%>%filter(!GID%in%c("UNKNOWN-OMIT"))%>%tidyr::drop_na(AlaAT)



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

GD<-wmb_GD_rr$imputed;GD<-GD+1;GD<-as.data.frame(GD)
GM<-GAPIT_wmb$GM

snpgdsCreateGeno("test4.gds", genmat = as.matrix(GD), sample.id=rownames(GD), snp.id = colnames(GD), snp.chromosome = GM$Chromosome,
                 snp.position = as.integer(GM$Position), snpfirstdim = F)
(gdsfile <- snpgdsOpen("test4.gds"))
set.seed(1000)
snpset <- snpgdsLDpruning(gdsfile, method="corr", ld.threshold=0.95, slide.max.bp = 3000)
snpset_names=unlist(snpset)
GD_prune=GD[,c(1,which(colnames(GD) %in% snpset_names))]
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

all_pheno$winter_survival<-as.numeric(all_pheno$winter_survival)
all_pheno$Ht<-as.numeric(all_pheno$Ht)
colnames(MQ_WMB21$sp_mean)
MQ_WMB21$ST
MQ_WMB21$Timepoint
all_pheno$Scald
field_data<-all_pheno%>%filter(trial=="PYT",Year=="2021")%>%dplyr::select(SourcePLOT,GID,Env,Row,Column,PHS,Scald,winter_survival,yield_kgha)%>%rename(PLOT=SourcePLOT,taxa=GID)
colnames(field_data)
str(field_data$PLOT)
field_data$PLOT<-as.character(field_data$PLOT)
data<-MQ_WMB21%>%dplyr::select(PLOT,Treatment,gid,Timepoint_group,Timepoint,Check.x,TB,AA,BG,ME,FAN,ST,DP,TotalProtein)%>%rename(taxa=gid)%>%left_join(field_data,by=c("PLOT","taxa"))%>%
  tidyr::pivot_longer(Timepoint,names_to = "name",values_to = "Timepoint")%>%tidyr::pivot_longer(cols=c("AA","BG","ME","FAN","ST","DP","TotalProtein"),names_to="trait")%>%arrange(Timepoint,trait,PLOT)
data
View(data)
#pheno_MQ<-all_pheno%>%tidyr::pivot_longer(cols=c("yield_kgha","winter_survival","Scald","SpotBlotch","Ht","Lodging","PHS","TP_JD","HD_JD","Maturity_JD","TWT"),names_to = "trait",values_to = "value")%>%
 # rename(taxa=GID)

GD<-GD_prune
GM<-GM_prune
GD_gapit<-GD_prune%>%tibble::add_column(taxa=rownames(GD_prune),.before = 1)


length(table(pheno$taxa))
rownames(GD)
pheno<-data%>%filter(taxa%in%rownames(GD))
pheno
library(dplyr)
t
df = pheno#%>%filter(!Year=="2020",trial=="PYT")%>%mutate(Fac = ifelse(Fac_type == "Fac",1,0))%>%filter(trait%in%c("winter_survival"))
#filter(trait=="Fac_type",trial=="PYT",Env%in%c("KET-2021","WMB22 Master - Snyder"))


library(sommer)
?GWAS()
#install_github('covaruber/sommer')



data(DT_polyploid)
table(df$Env)
DT <- df
DT
GT <- GD
MP <- GM

####=========================================####
####### convert markers to numeric format
####=========================================####
GT[1:5,1:5]
#numo <- atcg1234(data=GT, ploidy=4)

numo<-GD-1
numo[1:5,1:5]
?atcg1234()
###=========================================####
###### plants with both genotypes and phenotypes
###=========================================####
common <- intersect(DT$taxa,rownames(numo))
table(common)
###=========================================####
### get the markers and phenotypes for such inds
###=========================================####
common
marks <- numo[common,]; marks[1:5,1:5]
dim(marks)
DT2 <- DT%>%filter(taxa%in%common)

DT2 <- as.data.frame(DT2)
table(DT2)
DT2[1:5,]

###=========================================####
###### Additive relationship matrix, specify ploidy
###=========================================####
marks<-as.matrix(marks)
A <- sommer::A.mat(marks)
###=========================================####
### run it as GWAS model
###=========================================####
marks

table(dt3$Env)
dt3<-DT2%>%filter(trait=="FAN")
dt3$Timepoint
ans2 <- GWAS(value~Timepoint+Timepoint_group,
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

library(ggplot2)
chrTable=GM%>%arrange(Chromosome,Position)
chrTable = c(as.data.frame(table(GM$Chromosome))$Freq)
chrTable
chrTable
chrLabel = c(1:7, 'UN')
winterOrdinalBreaks = c(chrTable[1]/2)
WinterChrLines = c(as.data.frame(table(GM_prune$Chromosome))$Freq[1])
for (i in 2:8){
  winterOrdinalBreaks[i] = sum(chrTable[1:i-1])+chrTable[i]/2
  WinterChrLines[i] = sum(chrTable[1:i])
}

##

Out1 %>%
  ggplot(aes(x = ordinal, y = log10PVal)) +
  geom_point()+  geom_vline(xintercept = 4780, color = 'red')+
  
  geom_vline(xintercept = WinterChrLines)+
  ylim(0,max(Out1$log10PVal)*2.2)+
  
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = winterOrdinalBreaks)+
  ylab('-log(p-value)')+xlab('Chromosome')+ geom_hline(yintercept = -log10(5e-5))+theme_bw()

View(Out%>%filter(P.value<5e-5))
Out1$P.value
