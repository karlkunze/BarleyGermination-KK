# Script for the analysis, time series, genomic prediction, etc... Of the
# Cornell Winter Maling Barley DH families. 

library(sommer);library(arm);library(lme4);library(Hmisc);library(plyr);library(readxl);
library(tibble);library(patchwork);library(ggplot2);library(fda) ; library(magic); 
library(drc);library(rrBLUP);library(tidyr);library(ggh4x);library(dplyr);
setwd(rprojroot::find_rstudio_root_file())
library(dplyr)
library(tidyr)
getwd()
path_n<-"/home/karl/git/NY-winter-barley-analysis/"# change this based on the path for the neighboring git repo
path_n<-"C:/Users/kars8/git/NY-winter-barley-analysis/"
getwd()
load(paste0(path_n,"data/genotypes/wmb_GD_rrblup.Rdata"))
load(paste0(path_n,"data/genotypes/GAPIT_wmb.Rdata"))
load(paste0(path_n,"data/genotypes/wmb_pedigree.Rdata"))
load(paste0(path_n,"data/phenotypes/Field_data_2020_2022.Rdata"))

dh_list<-wmb_ped%>%filter(!Population%in%c("Recombinant_inbred_lines"))
wmb_dh_imputed<-wmb_GD_rr$imputed[rownames(wmb_GD_rr$imputed)%in%dh_list$Ind,]
WinterGD_1=wmb_dh_imputed%>%as.data.frame()

load('data/Genotype_data/WinterGD_ter.RData')

WinterGM=GAPIT_wmb$GM%>%filter(SNP%in%colnames(WinterGD))
#load('data/Genotype_data/WinterGM_ter.RData')
WinterGD[1:5,1:5]
WinterGD_1$taxa<-rownames(WinterGD_1)

WinterGM = WinterGM %>% arrange(Chromosome, Position)

#WinterGD<-WinterGD_1%>%relocate("taxa",.before=1)
#taxa = gsub(pattern = '-',replacement = '_',taxa),
WinterGD = WinterGD %>%
  mutate(
         taxa = gsub(pattern = ' ', replacement = '_',taxa),
         taxa1 = taxa) %>% remove_rownames()%>% 
  column_to_rownames('taxa1')
sum(WinterGM)
sum(colnames(WinterGD[,-1]))
sum(WinterGM$SNP == colnames(WinterGD[,-1]))
# Make sure things are in the right order
# Sum should = 8384
library(ggplot2)
# PCA plot of pop structure #####
WinterRelationship = rrBLUP::A.mat(WinterGD[,-1]-1, impute.method = 'EM', return.imputed = F)
WinterPCA = eigen(WinterRelationship)
WinterPVEPCA = WinterPCA$values/sum(WinterPCA$values)
data.frame(ordinal = 1:10, PVE = WinterPVEPCA[1:10]) %>%plot(., xlab = 'PC', col = 'red') 
winterlinePCAvalues = WinterPCA$vectors %>% data.frame()%>% 
  mutate(family = mapvalues(substr(WinterGD$taxa,1,3), from = c('BS6','BS7','BS8','BS9','DH1','Fla','SY_','Sca','Win'), 
                           to = c('Flavia/DH130910','Scala/DH130910','SY_Tepee/DH130910','Wintmalt/DH130910',
                                  'DH130910','Flavia/DH130910','SY_Tepee/DH130910','Scala/DH130910','Wintmalt/DH130910')),
         taxa = WinterGD$taxa,
         shapes = ifelse(taxa %in% c('DH130910', 'FLAVIA','TEPEE','WINTMALT','SCALA'), taxa, 'Lines'),
         size = ifelse(taxa %in% c('DH130910', 'FLAVIA','TEPEE','WINTMALT','SCALA'), 3, 2))
 
winterlinePCAvalues %>% ggplot(aes(x = X1, y = X2, color = family)) + geom_point()+
  winterlinePCAvalues%>% ggplot(aes(x = X1, y = X3,color = family)) + geom_point()+
  winterlinePCAvalues%>% ggplot(aes(x = X2, y = X3,color = family)) + geom_point()

winterlinePCAvalues %>% filter(family != 'Cha') %>%
  ggplot(aes(x = X1, y = X2, color = family, shape = shapes)) + geom_point(aes(size = size))+theme_bw() +guides(size = "none")+
  xlab('PC1')+ylab('PC2')
  
# max 2 PCs accounts for population structure, though 1 would likely work
# Functions H2 get_blues #####
BLUPH2 = function(trait.lm) {
  ses<- se.ranef(trait.lm)$'taxa' #where 'm' is your model object from 'lmer' (replace 'genotypes' with whatever you call your individuals in the data)
  v_BLUP<- ses^2
  sigma2_g=VarCorr(trait.lm, comp="Variance")$'taxa'[1]
  Reliability<- 1- v_BLUP/ (2*sigma2_g)  #where sigma2_g is the genetic variance estimated
  H2<- round(mean(Reliability),3)
  return(H2)
}
# Load and modify input data, calcuate the GI and GE for all assays ####

setwd(rprojroot::find_rstudio_root_file())
getwd()
DHs2020 = rbind(read_excel("data/Phenotype_Data/2020/DHtp1_all.xlsx", guess_max = 1723)%>%mutate(TP = 'TP1',PM_date =5),
            read_excel("data/Phenotype_Data/2020/DHtp2_all.xlsx", guess_max = 1723)%>%mutate(TP = 'TP2',PM_date =19),
            read_excel("data/Phenotype_Data/2020/DHtp3_all.xlsx", guess_max = 1723)%>%mutate(TP = 'TP3',PM_date =47),
            read_excel("data/Phenotype_Data/2020/DHtp4_all.xlsx",  guess_max = 1723)%>%mutate(TP = 'TP4',PM_date =96),
            read_excel("data/Phenotype_Data/2020/DHtp5_all.xlsx", guess_max = 1723)%>%mutate(TP = 'TP5',PM_date =152)) %>%
  dplyr::mutate(GE =(Day1Germ+Day2Germ+Day3Germ)/(Day1Germ+Day2Germ+Day3Germ+Day4Germ+Day5Germ+KernelsLeft),
         GI = 10*GE*(Day1Germ+Day2Germ+Day3Germ)/(Day1Germ+2*Day2Germ+3*Day3Germ),
         GI = ifelse(is.nan(GI),0,GI),
         GE5 = (Day1Germ+Day2Germ+Day3Germ+Day4Germ+Day5Germ)/(Day1Germ+Day2Germ+Day3Germ+Day4Germ+Day5Germ+KernelsLeft),
         GI5 = 10*GE5*(Day1Germ+Day2Germ+Day3Germ+Day4Germ+Day5Germ)/(Day1Germ+Day2Germ*2+Day3Germ*3+Day4Germ*4+Day5Germ*5),
         Location = ifelse(substr(PLOT,1,2)=='11' | substr(PLOT,1,2)=='a1','2020Snyder','2020Caldwell'),
         rep = as.factor(replication),
         taxa = mapvalues(Entry, from = c('Check 1','Check 2','Check 3','Check 4','Check 5','Check 6'),
                          to = c('Flavia', 'Scala','DH130910','SY_Tepee','Wintmalt','Charles')),
         taxa =  mapvalues(taxa, from = c("Check 1-Flavia","Check 2-Scala","Check 3-DH130910","Check 3-SY Tepee",
                                          "Check 4-SY Tepee","Check 5-Wintmalt","Check 6-Charles"),
                           to = c('Flavia', 'Scala','DH130910','DH130910','SY_Tepee','Wintmalt','Charles')),
        
         taxa = gsub(pattern = ' ', replacement = '_',taxa),
         year = '2020',
         Family = mapvalues(substr(taxa,1,3), from = c('BS6','BS7','BS8','BS9','DH1','Fla','SY_','Sca','Win','KWS'), 
                            to = c('Flavia/DH130910','Scala/DH130910','SY_Tepee/DH130910','Wintmalt/DH130910',
                                   'DH130910','Flavia/DH130910','SY_Tepee/DH130910','Scala/DH130910','Wintmalt/DH130910','Scala/DH130910')),
         taxa = toupper(taxa),
         taxa = ifelse(taxa == 'SY_TEPEE','TEPEE',taxa)) %>%
  filter(Entry %nin% c('BBBdormSNLine', 'BBBdormSRLine'))

DHs2021 = rbind(read_excel('data/Phenotype_Data/2021/DHs_GGS_TP1_PM5_full.xlsx') %>% mutate(notes = NA, TP = 'TP1'),
                read_excel('data/Phenotype_Data/2021/DHs_PM_12_full.xlsx')%>% mutate(notes = NA,TP = 'TP1.5'),
                read_excel('data/Phenotype_Data/2021/DHs_PM_19_GGS_TP2_full.xlsx') %>% mutate(notes = NA, TP = 'TP2'),
                read_excel('data/Phenotype_Data/2021/DHs_PM_33_GGS_TP3_full.xlsx') %>% mutate(TP = 'TP2.5'),
                read_excel('data/Phenotype_Data/2021/DHs_PM_47_GGS_TP4_full.xlsx') %>% mutate(TP = 'TP3'),
                read_excel('data/Phenotype_Data/2021/DHs_PM_68_GGS_TP5_full.xlsx')%>% mutate(TP='TP3.5'),
                read_excel('data/Phenotype_Data/2021/DHs_PM_96_full.xlsx') %>% mutate(TP = 'TP4'),
                read_excel('data/Phenotype_Data/2021/DHs_PM_150_GGS_TP7_full.xlsx') %>% mutate(Pmdate = 152,notes = NA,TP = 'TP5')) %>%
  filter(!is.na(Day2Germ)) %>% filter(AssayNumber != 10 & AssayNumber != 318 & AssayNumber != 190) %>% #These look to be contamination
  filter(Location == 'McGowan' | Location == 'Ketola') %>%
  dplyr::select(!c('seed_mass', 'GerminationAssay', 'MaltingQualitySampling')) %>%
  dplyr::mutate(GE =(Day1Germ+Day2Germ+Day3Germ)/(Day1Germ+Day2Germ+Day3Germ+Day4Germ+KernelsLeft),
         GI = 10*GE*(Day1Germ+Day2Germ+Day3Germ)/(Day1Germ+2*Day2Germ+3*Day3Germ), Rep = as.factor(Rep)) %>%
  dplyr::rename('taxa' = Entry, 'PM_date' = Pmdate, 'rep' = Rep)  %>%
  mutate(Family = mapvalues(substr(taxa,1,3), from = c('BS6','BS7','BS8','BS9','DH1','Fla','SY_','KWS','Win'), 
                            to = c('Flavia/DH130910','Scala/DH130910','SY_Tepee/DH130910','Wintmalt/DH130910',
                                   'DH130910','Flavia/DH130910','SY_Tepee/DH130910','Scala/DH130910','Wintmalt/DH130910')),
         GI = ifelse(is.nan(GI),0,GI), 
         taxa = ifelse(AssayNumber == 408, 'Charles',taxa), #AssayNumber 408 is 'Charles' not the entry it claims to be!
         taxa = toupper(taxa), 
         taxa = gsub(pattern = ' ', replacement = '_',taxa),
         taxa = ifelse(taxa == 'KWS_SCALA','SCALA',taxa), year = '2021')

DHs2021 %>% filter(substr(taxa,1,2) !='BS') %>% select(taxa) %>% unique()
#
##

#View(DHs2021)

#PHS 
#DH_Sprout = read_excel('data/Phenotype_Data/2021/Data/DHs_PHS_2021.xlsx')
#View(DH_Sprout)
#DH_Sprout = read_excel('data/Phenotype_Data/2021/Data/DHs_PHS_2021.xlsx') %>%
 # mutate(location = ifelse(PLOT>6999, 'Ketola','McGowan'), year = '2021') %>% 
  # rbind(read_excel('WinterBarley/PhenotypeData/2020Harvest/2020WinterDH_PHS.xlsx')%>%
  # #         mutate(location = 'Ketola2020'), year = '2020')  #Not sure if this should be included as these were planed as facultatives. 
  # select(!Comment) %>% mutate(taxa = gsub(pattern = '-', replacement = '_',Entry),taxa = gsub(pattern = ' ', replacement = '_',taxa)) %>%
  # rename(score = `Sprout Score`) %>% 
  # separate(score, into =c('p0','p1','p2','p3','p4','p5'), sep = '') %>% select(!p0) %>% pivot_longer(cols = c(p1,p2,p3,p4,p5)) %>%
  # mutate(value = as.numeric(value)) %>%
  # group_by(taxa, location, Harv, year,PLOT) %>% summarise(PHS = mean(value,na.rm = T))
#source from all df
#phs<-all_pheno%>%filter(GID%in%dh_list$Ind)%>%filter(!trial%in%"Screening")%>%dplyr::select(GID,SourcePLOT,Location,Year,PHS,PHS_in,PHS_var)
DH_Sprout
listt<-c(DHs2021$taxa,DHs2020$taxa)%>%unique()
rm(list)
all_pheno
DH_Sprout<-all_pheno%>%rename(taxa=GID,year=Year)%>%filter(!year=="2020",trial=="PYT",taxa%in%listt)%>%
  dplyr::select(taxa,Location,SourcePLOT,year,PHS,PHS_in,PHS_var)%>%mutate_at(vars(taxa,year,Location),  list(factor))
DH_Sprout$SourcePLOT<-as.character(DH_Sprout$SourcePLOT)
#we added 2022 here
library(dplyr)
#colnames(DH_Sprout)[5]<-"SourcePLOT"
#DHs2021$SourcePLOT<-as.numeric(DHs2021$SourcePLOT)
DH_Sprout
#going to run a model here quuck for the phs 
library(asreml)
library
str(DH_Sprout)
str(DH_Sprout)
df<-DH_Sprout%>%as.data.frame()
df
library(lme4)
lmer(PHS~(1|taxa)+year,data=df)
asreml()
asreml(fixed=GE~taxa,residual =~units,data=df,asr_gaussian(link = "identity"))

mod_MET1<-asreml(fixed=ME~GID,residual=~units,data=df)


DHs2021_PHS=full_join(DHs2021,DH_Sprout[,c("SourcePLOT","PHS","year")],by=c("SourcePLOT","year"))
table(DHs2021_PHS$taxa)
checks=c("DH130910","Endeavor","Scala","BS908_25")


#DHs2021_PHS=DHs2021_PHS[!DHs2021_PHS$taxa%in%checks,]
DHs2021_PHS
test1<-as.data.frame(DHs2021_PHS)%>%filter(Year=="2021")
test1$taxa
test1$c<-"Test Entry"
test1[test1$taxa=="ENDEAVOR",]$c="Endeavor"
test1[test1$taxa=="SCALA",]$c="KWS Scala"
test1[test1$taxa=="DH130910",]$c="DH130910(Lightning)"
test1$year<-as.factor(test1$year)
ggplot(data = test1, aes(x =GE, y =PHS,color=c,shape=year)) + 
  ggtitle("Correlation of GE TP1 and PHS")+
  xlab("Germination Energy(TP 1)") + ylab("Pre-harvest-sprouting")+
 geom_point(show.legend = FALSE)+
  
 # geom_smooth(orientation = "x",method = "lm",se=TRUE)+
  scale_color_manual(values = c("green","#E69F00","blue","black"))+
  annotate("rect", xmin =0, xmax = 1, ymin = 0, ymax = 3,
           alpha = .1,fill="green")+
  annotate("rect", xmin =0, xmax = 1, ymin = 3, ymax = 5,
           alpha = .1,fill="yellow")+
  annotate("rect", xmin =0, xmax = 1, ymin = 5, ymax = 9.1,
           alpha = .1,fill="red")+
  theme_bw()
#GI
ggplot(data = test1, aes(x =GI, y =PHS,color=c)) + 
  ggtitle("Correlation of GI for TP1 and PHS")+
  xlab("Germination Index(TP 1)") + ylab("Pre-harvest-sprouting")+
  geom_point(show.legend = FALSE)+
  
  # geom_smooth(orientation = "x",method = "lm",se=TRUE)+
  scale_color_manual(values = c("green","#E69F00","blue","black"))+
  annotate("rect", xmin =0, xmax = 8, ymin = 0, ymax = 3,
           alpha = .1,fill="green")+
  annotate("rect", xmin =0, xmax = 8, ymin = 3, ymax = 5,
           alpha = .1,fill="yellow")+
  annotate("rect", xmin =0, xmax = 8, ymin = 5, ymax = 9.1,
           alpha = .1,fill="red")+
  theme_bw()
##
cor(test1$GI,test1$PHS,use = "pairwise.complete.obs")
cor(test1$GE,test1$PHS,use = "pairwise.complete.obs")
cor(test1$Day2Germ,test1$PHS,use = "pairwise.complete.obs")


test2<-as.data.frame(DHs2021_PHS[DHs2021_PHS$TP=="TP2",])
cor(test2$GI,test2$PHS,use = "pairwise.complete.obs")
cor(test2$GE,test2$PHS,use = "pairwise.complete.obs")
test3<-as.data.frame(DHs2021_PHS[DHs2021_PHS$TP=="TP3",])
cor(test3$GI,test3$PHS,use = "pairwise.complete.obs")
cor(test3$GE,test3$PHS,use = "pairwise.complete.obs")
test4<-as.data.frame(DHs2021_PHS[DHs2021_PHS$TP=="TP4",])
cor(test4$GI,test4$PHS,use = "pairwise.complete.obs")
cor(test4$GE,test4$PHS,use = "pairwise.complete.obs")

DH.PHS = (ranef(PHS.lmer)$taxa + fixef(PHS.lmer)[1]) %>% as.data.frame() %>% rownames_to_column('taxa') %>% rename(PHS =`(Intercept)`)
length(table(DH.PHS[DH.PHS$PHS<1,]$PHS))
View(DH.PHS)


# 2020 GE and GI anova #####
anova(lm(GI~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==5))) #rep insignificant
anova(lm(GI ~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==19)))
anova(lm(GI ~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==47)))
anova(lm(GI ~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==96)))
anova(lm(GI ~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==152)))

anova(lm(GE ~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==5)))
anova(lm(GE ~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==19)))
anova(lm(GE ~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==47)))
anova(lm(GE ~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==96)))
anova(lm(GE ~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==152)))

# 2021 GE and  GI anova #####
anova(lm(GE ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==5)))
anova(lm(GE ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==12)))
anova(lm(GE ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==19)))
anova(lm(GE ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==33)))
anova(lm(GE ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==47)))
anova(lm(GE ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==68)))
anova(lm(GE ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==96)))
anova(lm(GE ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==152)))

anova(lm(GI ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==5)))
anova(lm(GI ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==12)))
anova(lm(GI ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==19)))
anova(lm(GI ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==33)))
anova(lm(GI ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==47)))
anova(lm(GI ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==68)))
anova(lm(GI ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==96)))
anova(lm(GI ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==152)))

# Location is pretty much always significant, rep bounces around. 
# 2020,2021 combined, only for TPs that are phenotypes at both years #####
DHs20202021 = DHs2021 %>% select(taxa, rep, Location, TP, GE,GI,PM_date, year) %>%
  rbind(.,DHs2020 %>% select(taxa, rep, Location,TP, GE, GI, PM_date, year)) %>%
  mutate(year = factor(year,levels = c('2020','2021')))

anova(lm(GE ~ taxa +year+year %in% Location+rep, DHs20202021 %>% filter(PM_date ==5)))
anova(lm(GE ~ taxa +Location+rep+year, DHs20202021 %>% filter(PM_date ==19)))
anova(lm(GE ~ taxa +Location+rep+year, DHs20202021 %>% filter(PM_date ==47)))
anova(lm(GE ~ taxa +Location+rep+year, DHs20202021 %>% filter(PM_date ==96)))
anova(lm(GE ~ taxa +Location+rep+year, DHs20202021 %>% filter(PM_date ==152)))

anova(lm(GI ~ taxa +Location+rep, DHs20202021 %>% filter(PM_date ==5)))
anova(lm(GI ~ taxa +Location+rep, DHs20202021 %>% filter(PM_date ==19)))
anova(lm(GI ~ taxa +Location+rep, DHs20202021 %>% filter(PM_date ==47)))
anova(lm(GI ~ taxa +Location+rep, DHs20202021 %>% filter(PM_date ==96)))
anova(lm(GI ~ taxa +Location+rep, DHs20202021 %>% filter(PM_date ==152)))
# Location very significant!
# Lets plot at the raw phentypes for the sake of the story: ######
DHs20202021 %>% pivot_longer(cols = c(GE,GI), names_to = 'trait') %>% filter(trait =='GE') %>%
  ggplot(aes(x = value, fill = TP))+geom_density()+facet_nested(TP ~trait+year, scales = 'free')+
  theme_bw()+guides(fill = "none")+
    DHs20202021 %>% pivot_longer(cols = c(GE,GI), names_to = 'trait') %>% filter(trait =='GI') %>%
    ggplot(aes(x = value, fill = TP))+geom_density()+facet_nested(TP ~trait+year, scales = 'free')+
    theme_bw()+ guides(fill = "none")+
    plot_layout(ncol = 2)

DHs20202021 %>%pivot_longer(cols = c(GE,GI), names_to = 'trait') %>% 
  ggplot(aes(x = PM_date, y = value, group = taxa))+geom_line()+facet_grid(trait~year, scales = 'free')

#Malting Quality

#View(MQ_WMB21)
TechnicalReps%>%filter(year=="2021")
# Lets extract out 2020, 2021, and 2020/2021 values for GWA: Always including rep and location #####

BlueBlupsH2_Location_rep_taxa <- function(d, groupvars) {
  trait.lmer <- lmer(formula = value ~(1|taxa)+Location+rep, 
                     data = d)
  lineEf = (ranef(trait.lmer)$taxa + fixef(trait.lmer)[1]) %>% as.data.frame() %>% rownames_to_column('taxa') %>% 
    mutate(type = 'BLUP') %>% rename(value = '(Intercept)')
  trait.lm = broom::tidy(lm(value~ taxa+Location+rep, data=d))
  
  first_taxa = d %>% arrange(taxa) %>% slice_head( n = 1) %>% select(taxa) %>% as.character()
  Intercept = trait.lm %>% filter(term == '(Intercept)') %>% select(estimate) %>% as.numeric()
  lineBLUE = trait.lm %>% filter(substr(term,1,4)=='taxa') %>% 
    add_row(term = paste0('taxa',first_taxa),
            estimate = 0) %>% mutate(BLUE = estimate + Intercept) %>%
    transmute(taxa = gsub(pattern = 'taxa',replacement = '',x = term),
              value = BLUE, type = 'BLUE')
 H2 = BLUPH2(trait.lmer)
  return(rbind(lineEf, lineBLUE) %>% add_row(value = H2, type = 'H2') %>% arrange(type, taxa) %>%
           join(d %>% select(taxa, Family) %>% unique()))
}

d2 = DHs2020 %>% select(taxa, rep, Location,TP, GE, GI,PM_date,year) %>%
  rbind(., DHs2021 %>% select(taxa, rep, Location, TP, GE,GI,PM_date,year)) %>% 
  mutate(year = factor(year, levels = c('2021','2020'))) %>%
  pivot_longer(cols = c(GE, GI), names_to = 'trait') %>% filter(TP == 'TP2', trait == 'GE')

BlueBlupsH2_Year_rep_taxa <- function(d2, groupvars) {
  if (length(unique(d2$year))==2) {
  trait.lmer <- lmer(formula = value ~(1|taxa)+Location + rep, data = d2)
  # fixef(trait.lmer)[2]/need to think about this more. 
  
  Cept = (fixef(trait.lmer)[1]*4+sum(fixef(trait.lmer)[2:4]))/4
  lineEf = (ranef(trait.lmer)$taxa + Cept) %>% as.data.frame() %>% rownames_to_column('taxa') %>% 
    mutate(type = 'BLUP') %>% rename(value = '(Intercept)')
  
  trait.lm = broom::tidy(lm(value~ taxa+Location+ rep, data=d2))
  first_taxa = d2 %>% arrange(taxa) %>% slice_head( n = 1) %>% select(taxa) %>% as.character()
  Intercept_list = trait.lm %>% filter(term == '(Intercept)'|substr(term,1,8)=='Location')
  Intercept = (Intercept_list$estimate[1]*4+sum(Intercept_list$estimate[2:4]))/4

  lineBLUE = trait.lm %>% filter(substr(term,1,4)=='taxa') %>% 
    add_row(term = paste0('taxa',first_taxa),
            estimate = 0) %>% mutate(BLUE = estimate + Intercept) %>%
    transmute(taxa = gsub(pattern = 'taxa',replacement = '',x = term),
              value = BLUE, type = 'BLUE')
  H2 = BLUPH2(trait.lmer)
  return(rbind(lineEf, lineBLUE) %>% add_row(value = H2, type = 'H2')%>% join(d2 %>% select(taxa, Family)%>% unique()))
  }
  else {
    trait.lmer <- lmer(formula = value ~(1|taxa)+Location+rep, 
                       data = d2)
    lineEf = (ranef(trait.lmer)$taxa + fixef(trait.lmer)[1]) %>% as.data.frame() %>% rownames_to_column('taxa') %>% 
      mutate(type = 'BLUP') %>% rename(value = '(Intercept)')
    trait.lm = broom::tidy(lm(value~ taxa+Location+rep, data=d2))
    
    first_taxa = d2 %>% arrange(taxa) %>% slice_head( n = 1) %>% select(taxa) %>% as.character()
    Intercept = trait.lm %>% filter(term == '(Intercept)') %>% select(estimate) %>% as.numeric()
    lineBLUE = trait.lm %>% filter(substr(term,1,4)=='taxa') %>% 
      add_row(term = paste0('taxa',first_taxa),
              estimate = 0) %>% mutate(BLUE = estimate + Intercept) %>%
      transmute(taxa = gsub(pattern = 'taxa',replacement = '',x = term),
                value = BLUE,type = 'BLUE')
    H2 = BLUPH2(trait.lmer)
    return(rbind(lineEf, lineBLUE) %>% add_row(value = H2, type = 'H2') %>% arrange(type, taxa) %>%
             join(d2 %>% select(taxa, Family) %>% unique()))
  }
}

DH2020Estimates = DHs2020 %>% select(taxa, rep, Location,TP, GE, GI, PM_date,year, Family) %>% 
  pivot_longer(cols = c(GE, GI), names_to = 'trait') %>%
  group_by(TP, PM_date, trait, year) %>% group_modify(BlueBlupsH2_Location_rep_taxa) %>% ungroup()
 
DH2021Estimates = DHs2021 %>% select(taxa, rep, Location, TP, GE,GI,PM_date,year,Family) %>% pivot_longer(cols = c(GE, GI), names_to = 'trait') %>%
  group_by(TP,PM_date, trait, year) %>% group_modify(BlueBlupsH2_Location_rep_taxa) %>% ungroup()
 
DHCombined = DHs2020 %>% select(taxa, rep, Location,TP, GE, GI,PM_date,year,Family) %>%
  rbind(., DHs2021 %>% select(taxa, rep, Location, TP, GE,GI,PM_date,year,Family)) %>% 
  mutate(year = factor(year, levels = c('2021','2020'))) %>%  pivot_longer(cols = c(GE, GI), names_to = 'trait') %>%
  group_by(TP,PM_date, trait) %>% group_modify(BlueBlupsH2_Year_rep_taxa)  %>% mutate(year = '2020/2021')

AllDHBluesPerYear = rbind(DH2020Estimates, DH2021Estimates,DHCombined) %>% filter(type =='BLUE') %>% ungroup()
save(AllDHBluesPerYear, file = 'data/Analysis/AllDHBluesPerYear.RData')
save(AllDHBluesPerYear, file = 'C:/Users/kars8/git/Cornell-WMB21-selections/data/phenotypes/Germination/AllDHBluesGerminationPerYear.RData')
View(AllDHBluesPerYear)


# Plotting of 2020 vs 2021 to find outliers and removed contamination #####
DH2020Estimates %>% join(DH2021Estimates  %>% select(!year)%>% rename(value2021 = value)) %>%
  filter(!is.na(value2021)) %>% filter(trait == 'GE') %>% filter(type !='H2') %>%
  ggplot(aes(x = value, y= value2021)) +geom_point()+facet_grid(rows = vars(TP), cols = vars(type))+
  xlab('2020 BLUE or BLUP')+ylab('2021 BLUP or BLUP')

DH2020Estimates %>% join(DH2021Estimates  %>% select(!year)%>% rename(value2021 = value))  %>%
  filter(!is.na(value2021)) %>% filter(trait == 'GI') %>% filter(type !='H2') %>%
  ggplot(aes(x = value, y= value2021)) +geom_point()+facet_grid(rows = vars(TP), cols = vars(type))+
  xlab('2020 BLUE or BLUP')+ylab('2021 BLUP or BLUP')

# What is the correlation between the 2020 values and the 2021 values?
DH2020Estimates %>% join(DH2021Estimates  %>% select(!year)%>% rename(value2021 = value))  %>%
  filter(!is.na(value2021), type =='BLUE') %>% filter(type !='H2') %>% group_by(type, TP, trait) %>%
  summarise(correlation = cor(value, value2021)) 
#could add phs data here
# heritability plots and over time points #####
DH2020Estimates %>% rbind(DH2021Estimates, DHCombined) %>%
  filter(type == 'H2') %>% ggplot(aes(x = TP, y = value, fill = trait)) +geom_bar(stat = 'identity', position = 'dodge')+
  facet_wrap(vars(year), ncol = 1)+theme_bw()+labs(title= 'Broard sense heritability\nover time and datasets')



obj<-rbind(DH2021Estimates, DH2020Estimates) %>%filter(type == 'BLUE')%>% mutate(headerFacet = 'Data source')%>%

 group_by(TP,trait,year,Family) %>%
 summarise_at(vars(value), list(mean = mean))
obj
#plots over time - reveals some issues with combining data from a time series perspective
c<-rbind(DH2021Estimates, DH2020Estimates) %>%  filter(type == 'BLUE')%>%join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1!= 1) %>% mutate(Qsd1= as.factor(Qsd1)) %>%
  filter(Family %nin% c('Cha','End','DH130910')) %>% 
  mutate(Qsd1= ifelse(Qsd1==2,'Dormant','Nondormant')) %>% mutate(headerFacet = 'Data source')
  
  
dfr<-rbind(DH2021Estimates, DH2020Estimates) %>%  filter(type == 'BLUE')%>%join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1!= 1) %>% mutate(Qsd1= as.factor(Qsd1)) %>%
  filter(Family %nin% c('Cha','End','DH130910')) %>% 
  mutate(Qsd1= ifelse(Qsd1==2,'Dormant','Nondormant')) %>% mutate(headerFacet = 'Year')
 dfr%>%ggplot(aes(x = PM_date, y = value,fill=Qsd1, group = taxa))+
  
  geom_line()+
facet_nested(trait~headerFacet+year, scales = 'free')+
  #geom_vline(xintercept = c(12,33,68), color = 'red',linetype="longdash",alpha=0.5)+
  geom_vline(xintercept = c(68,140),color = 'red',linetype="longdash",alpha=1,size=0.8)+
 # annotate("rect", xmin =0, xmax = 150, ymin = 0.5, ymax = 1,
 # color="black")+
  labs(title = 'GE and GI of taxa over time across two years')+theme_bw()
#box plot this on TP

rbind(DH2021Estimates, DH2020Estimates) %>%  filter(type == 'BLUE')%>%join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1!= 1) %>% mutate(Qsd1= as.factor(Qsd1)) %>%
  filter(Family %nin% c('Cha','End','DH130910')) %>% 
  mutate(Qsd1= ifelse(Qsd1==2,'Dormant','Nondormant')) %>% mutate(headerFacet = 'Year') %>%
  ggplot(aes(x = PM_date, y = value,group=taxa,color=Qsd1))+
  
  geom_line(show.legend = FALSE)+
  facet_nested(trait~headerFacet+year, scales = 'free')+
  geom_vline(xintercept = c(68,140),color = 'red',linetype="longdash",alpha=1,size=0.8)+
  # annotate("rect", xmin =0, xmax = 150, ymin = 0.5, ymax = 1,
  # color="black")+
  ylab(label = "GE/GI value")+
  labs(title = 'GE and GI of taxa over time across two years based on Qsd1')+theme_bw()



##

DHCombined %>% filter(type == 'BLUE' & trait == 'GI') %>% ggplot(aes(x = PM_date, y = value, group = taxa))+geom_line()+
  geom_vline(xintercept = c(12,33,68), color = 'red')
dev.off()
png('WinterBarley/WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/Sup_BluesByFamilyQsd1AllYears.png', 2400, 1500, res =120)
#WinterGD[-1]=WinterGD[-1]+1
WinterGD[1:5,1:5]

AllDHBluesPerYear
Qsd1<-WinterGD[,c('taxa','Qsd1')]
Qsd1[Qsd1$Qsd1==1,]$Qsd1<-2
Qsd1[Qsd1$taxa=="FLAVIA",]$Qsd1<-0
Qsd1
table(AllDHBluesPerYear$Family)
AllDHBluesPerYear %>%  filter(type == 'BLUE') %>% mutate(year = factor(year, levels = c('2020','2021','2020/2021')))%>%
  left_join(Qsd1,by="taxa") %>% filter(Qsd1!= 1) %>% mutate(Qsd1= as.factor(Qsd1)) %>%
  filter(Family %nin% c('Cha','End','DH130910')) %>% 
  mutate(Qsd1= ifelse(Qsd1==2,'Dormant','Nondormant')) %>%
  filter(!(trait =='GE' &value>1.05)) %>%
  filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  ggplot(aes(x = TP, y = value, fill = Qsd1))  +
  geom_boxplot()+facet_nested(trait~year+Family, scales = 'free', space = 'free_x')

dev.off()
#2020
png('plots/Presentation2022_01/BluesByFamilyQsd1_2020.png', 1400, 800, res =120)
AllDHBluesPerYear$PM_date<-as.factor(AllDHBluesPerYear$PM_date)
AllDHBluesPerYear %>%  filter(type == 'BLUE') 
AllDHBluesPerYear %>%  filter(type == 'BLUE') %>% mutate(year = factor(year, levels = c('2020','2021','2020/2021')))%>%
  left_join(Qsd1,by="taxa") %>% filter(Qsd1!= 1) %>% mutate(Qsd1= ifelse(Qsd1==2,'Dormant','Nondormant')) %>%
  filter(Family %nin% c('Cha','End','DH130910')) %>% filter(year %in% c('2020','2021'))  %>%
  filter(!(year == '2020' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  filter(!(trait =='GE' &value>1.05)) %>%
  ggplot(aes(x = PM_date, y = value, fill = Qsd1))+  
  geom_vline(xintercept = c(1,5,color = 'red',linetype="longdash",alpha=1,size=0.8)+
  geom_boxplot(outlier.shape = NA,show.legend = NA)+facet_nested(trait~year, scales = 'free')+
 
  
  geom_vline(xintercept=68)+
  ggtitle(label="Germination over Timepoints by Qsd1")+
  ylab(label = "GE/GI value")+
  xlab(label="Post Days after Maturity")+

 
 theme_bw()+

  theme(legend.position = "none")
dev.off()
png('plots/Presentation2022_01/BluesByFamilyQsd1_2021.png', 1400, 800, res =120)
AllDHBluesPerYear %>%  filter(type == 'BLUE') %>% mutate(year = factor(year, levels = c('2020','2021','2020/2021')))%>%
  join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1!= 1) %>% mutate(Qsd1= ifelse(Qsd1==2,'Dormant','Nondormant')) %>%
  filter(Family %nin% c('Cha','End','DH130910')) %>% filter(year %in% c('2021'))  %>%
  filter(!(year == '2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  filter(!(trait =='GE' &value>1.05)) %>%
  ggplot(aes(x = TP, y = value, fill = Qsd1))+  
  geom_boxplot()+facet_nested(trait~year, scales = 'free')

AllDHBluesPerYear %>%  filter(type == 'BLUE') %>% mutate(year = factor(year, levels = c('2020','2021','2020/2021')))%>%
  join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1!= 1) %>% mutate(Qsd1= ifelse(Qsd1==2,'Dormant','Nondormant')) %>%
  filter(Family %nin% c('Cha','End','DH130910')) %>% filter(year %in% c('2020'))  %>%
  filter(!(year == '2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  filter(!(trait =='GE' &value>1.05)) %>%
  ggplot(aes(x = TP, y = value, fill = Qsd1))+  
  geom_boxplot(show.legend = NULL)+facet_nested(trait~year, scales = 'free')+theme_bw()

dev.off()
png('plots/Presentation2022_01/BluesByFamilyQsd1_2020_2021.png', 1400, 800, res =120)
AllDHBluesPerYear %>%  filter(type == 'BLUE') %>% mutate(year = factor(year, levels = c('2020','2021','2020/2021')))%>%
  join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1!= 1) %>% mutate(Qsd1= ifelse(Qsd1==2,'Dormant','Nondormant')) %>%
  filter(Family %nin% c('Cha','End','DH130910')) %>% filter(year %in% c('2020/2021'))  %>%
  filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  ggplot(aes(x = TP, y = value, fill = Qsd1))+  
  
  geom_boxplot(show.legend = NULL)+facet_nested(trait~year, scales = 'free')+theme_bw()
# what is the heritability if Qsd1 is accounted for in the model? 
dev.off()
DHs2020 %>% select(taxa, rep, Location,TP, GE, GI,PM_date,year,Family) %>%
  rbind(., DHs2021 %>% select(taxa, rep, Location, TP, GE,GI,PM_date,year,Family)) %>% 
  mutate(year = factor(year, levels = c('2021','2020'))) %>%  pivot_longer(cols = c(GE, GI), names_to = 'trait') %>% 
  filter(Family %nin%  c('Cha','End','DH130910')) %>%
  group_by(TP,PM_date, trait, Family) %>% join(WinterGD[,c('taxa', 'Qsd1')]) %>% filter(Qsd1 != 1) %>%
  group_modify(~{data.frame(H2 = BLUPH2(lmer(value~Location+rep+ Qsd1+(1|taxa),data = .x)))}) %>% ungroup() %>% select(TP,trait,Family,H2) %>%
  pivot_wider(values_from = H2, names_from = c(TP,trait))


# Now that we have the per TP values we can run GWA -  
# Single TP, TP/Family, TP/Family/year, TP/year GWAS #####
library(GAPIT3)
DF_OperationsV3 = function(df){
  df = df %>% arrange(P.value)  %>% mutate(logPercentileQQplot = -log10(c(1:length(df$SNP))/length(df$SNP)),
           rank = c(1:length(df$SNP))) %>% arrange(Chromosome, Position) %>%
    mutate(log10PVal = -log10(P.value),ordinal = c(1:length(df$SNP)))
  return(df)
}
GWA_MLMM_fortidyR = function(df, groupvars) {
   GWAS = GAPIT(Y = df %>% select(taxa, value) %>% as.data.frame(),GD=WinterGD, GM=WinterGM,PCA.total = 2,
               Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
  Out = DF_OperationsV3(GWAS$GWAS) %>% arrange(P.value) %>% slice_head(n=1000)
  return(Out)
}
GWA_MLMM_fortidyR.perFamily = function(df, groupvars) {
  WinterGD_sub = WinterGD %>% filter(taxa %in% df$taxa)
  Out = tryCatch(
    {GWAS = suppressWarnings(GAPIT(Y = df %>%ungroup() %>% select(taxa, value) %>% as.data.frame(),GD=WinterGD_sub, GM=WinterGM,PCA.total = 0,
               Geno.View.output=F, model="MLMM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05))
      DF_OperationsV3(GWAS$GWAS) %>%
        arrange(P.value) %>% slice_head(n = 100)
    }, error = function(e){data.frame() } )
   return(Out)
}
#GWA_MLM_fortidyR = function(df, groupvars) {
  GWAS = GAPIT(Y = df %>% ungroup() %>% select(taxa, value) %>% as.data.frame(),GD=WinterGD, GM=WinterGM,PCA.total = 2,
               Geno.View.output=F, model="MLM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)
  Out = (GWAS$GWAS) %>% arrange(P.value)  %>%
    mutate(logPercentileQQplot = -log10(c(1:length(GWAS$GWAS$SNP))/length(GWAS$GWAS$SNP)),
           rank = c(1:length(GWAS$GWAS$SNP))) %>% arrange(Chromosome, Position) %>%
    mutate(log10PVal = -log10(P.value), ordinal = c(1:length(GWAS$GWAS$SNP)))%>% 
    arrange(P.value) %>% slice_head(n = 1000)
  return(Out)
}
#GWA_MLM_fortidyR.perFamily = function(df, groupvars) {
  WinterGD_sub = WinterGD %>% filter(taxa %in% df$taxa)
  Out = tryCatch(
    {GWAS = suppressWarnings(GAPIT(Y = df %>% ungroup() %>% select(taxa, value) %>% as.data.frame(),GD=WinterGD_sub, GM=WinterGM,PCA.total = 0,
                                   Geno.View.output=F, model="MLM", Major.allele.zero = F, file.output=F,SNP.MAF = 0.05))
    (GWAS$GWAS) %>% arrange(P.value)  %>%
      mutate(logPercentileQQplot = -log10(c(1:length(GWAS$GWAS$SNP))/length(GWAS$GWAS$SNP)),
             rank = c(1:length(GWAS$GWAS$SNP))) %>%  arrange(Chromosome, 'Position') %>%
      mutate(log10PVal = -log10(P.value), ordinal = c(1:length(GWAS$GWAS$SNP)))%>% 
      arrange(P.value) %>% filter(P.value < 0.01)
    }, error = function(e){data.frame() } )
  return(Out)
}
Add_DH130910_familyStructure = function(df3, groupvar){
  DH130910_rows = df3 %>% filter(taxa == 'DH130910') %>% select(!Family)
  df3= df3 %>% filter(taxa != 'DH130910') %>%
    rbind(DH130910_rows %>% mutate(Family='Flavia/DH130910'),
          DH130910_rows %>% mutate(Family='SY_Tepee/DH130910'),
          DH130910_rows %>% mutate(Family='Scala/DH130910'),
          DH130910_rows %>% mutate(Family='Wintmalt/DH130910')
          ) %>%
    filter(Family %in% c('Flavia/DH130910','SY_Tepee/DH130910','Scala/DH130910','Wintmalt/DH130910'))
  return(df3)
}

chrTable = c(785, 1534, 1187,  760, 1348,  823, 1334,16 )
chrLabel = c(1:7, 'UN')
winterOrdinalBreaks = c(chrTable[1]/2)
WinterChrLines = c(785)
for (i in 2:8){
  winterOrdinalBreaks[i] = sum(chrTable[1:i-1])+chrTable[i]/2
  WinterChrLines[i] = sum(chrTable[1:i])
}

# per time point GWAS with MLMM
WinterPerTPGWAS = DH2020Estimates %>% rbind(DH2021Estimates, DHCombined) %>%  filter(type == 'BLUE') %>%
  ungroup() %>% group_by(year, TP, trait) %>% group_modify(GWA_MLMM_fortidyR)

WinterPerTPGWAS %>% arrange(P.value) %>% ungroup()%>% select(SNP, Chromosome, Position) %>% unique()
WinterPerTPGWAS %>% arrange(P.value) %>% view()
save(WinterPerTPGWAS,file = "data/WinterPerTPGWAS.Rdata")
WinterPerTPGWAS %>% ggplot(aes(ordinal, log10PVal, color = TP, shape = year))+geom_point()+
  geom_vline(xintercept = WinterChrLines, color = 'black')+
  geom_vline(xintercept = 4780, color = 'red')+
  annotate(geom= 'text', x = 4780, y = 30, label = 'AlaAT1')+
  geom_vline(xintercept = WinterChrLines)+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = winterOrdinalBreaks)+
  ylab('-log(p-value)')+xlab('Chromosome')+ geom_hline(yintercept = -log10(5e-5)) +
  facet_grid(rows = vars(trait), scales = 'free_y')+theme_bw()

# Per TP MLM 
#WinterPerTPGWAS.MLM = DH2020Estimates %>% rbind(DH2021Estimates, DHCombined) %>% 
 # filter(type == 'BLUE') %>% ungroup() %>% group_by(year, TP, trait) %>% group_modify(GWA_MLM_fortidyR)
# Plot over time the MLM values as MLMM may swap from one loci to another
#WinterPerTPGWAS.MLM %>% filter(SNP %in% c('Qsd1','JHI-Hv50k-2016-311908','JHI-Hv50k-2016-276836')) %>%
 # mutate(label = paste0(SNP,'\n',Chromosome,':',`Position `)) %>%
  #ggplot(aes(TP, log10PVal, group = SNP, color = label))+geom_line()+
#  ylab('-log(p-value)')+xlab('TP')+ geom_hline(yintercept = -log10(5e-5), color = 'green')+
 # facet_grid(rows = vars(trait), cols = vars(year), scales = 'free_y')+theme_bw()

# per TP per Family with MLMM 
WinterPerTPGWAS.PerFamily = DH2020Estimates %>% rbind(DH2021Estimates, DHCombined) %>% ungroup() %>%  Add_DH130910_familyStructure() %>%
  filter(type == 'BLUE') %>%  ungroup() %>% group_by(Family, year, TP, trait) %>% group_modify(GWA_MLMM_fortidyR.perFamily)
WinterPerTPGWAS.PerFamily %>% arrange(P.value) %>% view()

WinterPerTPGWAS.PerFamily%>% filter(SNP %in% c('Qsd1','JHI-Hv50k-2016-311908','JHI-Hv50k-2016-276836')) %>%
  mutate(label = paste0(SNP,'\n',Chromosome,':',Position)) %>%
  ggplot(aes(TP, log10PVal, group = SNP, color = label))+geom_line()+
  ylab('-log(p-value)')+xlab('TP')+ geom_hline(yintercept = -log10(5e-5), color = 'green')+
  facet_nested(year~Family+trait, scales = 'free_y')+theme_bw()

# Per TP per family MLM 
#WinterPerTPGWAS.PerFamily.mlm = DH2020Estimates %>% rbind(DH2021Estimates, DHCombined)%>% ungroup() %>% 
 # Add_DH130910_familyStructure() %>% filter(type == 'BLUE') %>% ungroup() %>% group_by(Family, year, TP, trait) %>%
  #group_modify(GWA_MLM_fortidyR.perFamily)

#WinterPerTPGWAS.PerFamily.mlm %>% filter(SNP %in% c('Qsd1','JHI-Hv50k-2016-311908','JHI-Hv50k-2016-276836','JHI-Hv50k-2016-308652')) %>%
 # mutate(label = paste0(SNP,'\n',Chromosome,':',`Position `)) %>% filter(year =='2021') %>%
  #ggplot(aes(TP, log10PVal, group = SNP, color = label))+geom_line()+
#  ylab('-log(p-value)')+xlab('TP')+ geom_hline(yintercept = -log10(5e-5), color = 'green')+
 # facet_nested_wrap(~Family+trait, nrow = 2, scales = 'free_y')+theme_bw()

#WinterPerTPGWAS.PerFamily.mlm %>% filter(Family == 'Flavia/DH130910') %>% arrange(P.value) %>% View()
#WinterPerTPGWAS.PerFamily.mlm %>% filter(Family == 'Scala/DH130910') %>% arrange(P.value) %>% View()
#WinterPerTPGWAS.PerFamily.mlm %>% arrange(P.value) %>% View()


Dprime_LD = function(SNPs) { #SNPs is a n rowed by 2 column dataframe with the 2 snps for LD calc
  AlleleFreq1 = sum(SNPs[,1])/length(SNPs[,1])/2
  AlleleFreq2 = sum(SNPs[,2])/length(SNPs[,2])/2
  P_11 = table(SNPs)[1,1]/length(SNPs[,1])
  P_22 = table(SNPs)[2,2]/length(SNPs[,1])
  P_12 = table(SNPs)[1,2]/length(SNPs[,1])
  P_21 = table(SNPs)[2,1]/length(SNPs[,1])
  
  D = P_11*P_22 - P_21*P_12
  Dmax = ifelse(D>0, min(AlleleFreq1*(1-AlleleFreq2),(1-AlleleFreq1)*AlleleFreq2),
                min(AlleleFreq1*AlleleFreq2,(1-AlleleFreq1)*(1-AlleleFreq2)))
  D_prime = D/Dmax
  return(D_prime)
  }

#WinterPerTPGWAS.PerFamily.mlm %>% rbind(WinterPerTPGWAS.MLM%>% mutate(Family = 'Overall', .before = year)) %>%
 # mutate(ImportantSNPs = ifelse(SNP %in% c('Qsd1','JHI-Hv50k-2016-311497',
#                                           'JHI-Hv50k-2016-311574','JHI-Hv50k-2016-337656',
 #                                          'JHI-Hv50k-2016-273541'), SNP, 'others'),
  #       ImportantSNPs = factor(ImportantSNPs, levels =c('Qsd1','JHI-Hv50k-2016-311497','JHI-Hv50k-2016-311574', 'others') ))%>%
 # filter(SNP %in% c('Qsd1','JHI-Hv50k-2016-311497','JHI-Hv50k-2016-311574','JHI-Hv50k-2016-337656','JHI-Hv50k-2016-273541'))%>%
  ggplot(aes(x = TP, y = log10PVal, group = SNP, color = SNP)) +
 # geom_line()+facet_nested(trait+year~Family, scales = 'free')+theme_bw()

# Plots etc, not the most useful but we will see here #####

WinterPerTPGWAS %>% filter(P.value<1E-4) %>% filter(year == '2020/2021') %>%
  ggplot(aes(x = ordinal, y = log10PVal, color = TP, shape = trait)) +
    geom_point()+  geom_vline(xintercept = 4780, color = 'red')+
    annotate(geom= 'text', x = 4780, y = 30, label = 'AlaAT1')+
    geom_vline(xintercept = WinterChrLines)+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = winterOrdinalBreaks)+
  ylab('-log(p-value)')+xlab('Chromosome')+ geom_hline(yintercept = -log10(5e-5))

WinterPerTPGWAS %>% filter(P.value<1E-4) %>% filter(year == '2020') %>%
  ggplot(aes(x = ordinal, y = log10PVal, color = TP, shape = trait)) +
  geom_point()+ geom_vline(xintercept = WinterChrLines)+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = winterOrdinalBreaks)+
  ylab('-log(p-value)')+xlab('Chromosome')+
  geom_hline(yintercept = -log10(5e-5))+geom_vline(xintercept = 4772, color = 'red')+
  annotate(geom= 'text', x = 4772, y = 30, label = 'AlaAT1')

WinterPerTPGWAS %>% filter(P.value<1E-4) %>% filter(year == '2021') %>%
  ggplot(aes(x = ordinal, y = log10PVal, color = TP, shape = trait)) +
  geom_point()+  geom_vline(xintercept = WinterChrLines)+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = winterOrdinalBreaks)+
  ylab('-log(p-value)')+xlab('Chromosome')+
  geom_hline(yintercept = -log10(5e-5))+geom_vline(xintercept = 4772, color = 'red')+
  annotate(geom= 'text', x = 4772, y = 30, label = 'AlaAT1')

# PResent this #####
Winter_traits_withsigMarkers = WinterGD %>% select(taxa,'Qsd1','JHI-Hv50k-2016-308652' ) %>%
  filter(Qsd1!= 1) %>% filter(`JHI-Hv50k-2016-308652` != 1) %>%
  mutate(haplo = paste0(Qsd1,`JHI-Hv50k-2016-308652`))%>% select(taxa, haplo, Qsd1) %>% 
  merge(.,rbind(DH2020Estimates, DH2021Estimates, DHCombined %>% mutate(year = '2020/2021')),by = 'taxa') %>% 
  filter(trait %in% c('GE','GI') & type == 'BLUE') %>%
  filter(Family %nin% c('DH130910','Flavia','SY_Tepee','Scala','Wintmalt', 'Cha'))

Winter_traits_withsigMarkers %>% ggplot(aes(x = TP, y= value, fill = as.factor(Qsd1)))+
  geom_boxplot() +facet_nested(trait+year~Family, scales = 'free')+labs(title = 'by Qds1 at 442 Mbp', fill ='Qsd1')

Winter_traits_withsigMarkers %>%filter(year =='2021') %>%ggplot(aes(x = as.factor(PM_date), y= value, fill = as.factor(Qsd1)))+
  geom_boxplot() +facet_nested(trait+year~Family, scales = 'free')+labs(title = 'by Qds1 at 442 Mbp', fill ='Qsd1')+xlab('PM_date')
# Time series startings #########
# Lets predict the first TP 2020 unobserved lines using rrblup
FoldCVGPv2 = function(df, myGD=WinterGD, numfolds,datasplit, trait_name){
  # Phenotype is full vector of phenotypes
  # myGD is a n taxa x N markers dataframe with -1,0,1 coding and no 'taxa' column in it
  # numfolds is the number of folds to cv (ie 5 for five-fold cv)
  # datasplit is the percentage as a decimal for the training/testing split.
  
  df = df %>% filter(taxa %in% myGD$taxa) %>% arrange(taxa) 
  myGD = myGD %>% filter(taxa %in% df$taxa) %>% arrange(taxa)
  dim(df)
  dim(myGD)
  if(!length(unique(df$taxa))==sum(df$taxa == myGD$taxa)){
    stop('taxa lists are not correctly aligning')
  }
  
  TrueandPredicted = data.frame()
  
  Training_size = round(datasplit*num_entries,0)
  start_time <- Sys.time() 
  set.seed(1)
  Phenotype = as.vector(df$value)
  num_entries = length(Phenotype)
  myGD = myGD[,-1]-1
  for (i in 1:numfolds){
    
    trainingSet = sort(sample(1:num_entries,Training_size))
    testingSet = setdiff(1:num_entries,trainingSet)
    
    y_train = Phenotype[trainingSet] 
    y_test = Phenotype[testingSet]
    
    marker_train = myGD[trainingSet,]
    marker_test = myGD[testingSet,]
    
    trained.Model = mixed.solve(y_train, Z = marker_train, K = NULL, SE =FALSE)
    
    PredictedPheno = as.matrix(marker_test) %*% as.matrix(trained.Model$u)+ as.numeric(trained.Model$beta)
    print(cor(y_test,PredictedPheno))
    
    TrueandPredicted = rbind(TrueandPredicted, 
                             data.frame(TruePheno = y_test,
                                        PredPheno = PredictedPheno,
                                        taxa = df$taxa[testingSet],
                                        fold = i,
                                        trait = trait_name,
                                        correlation = cor(y_test,PredictedPheno)))
  }
  end_time <- Sys.time()
  print(end_time-start_time)
  return(TrueandPredicted)
}

GETP22020= FoldCVGPv2(df = DH2020Estimates %>% filter(TP=='TP2', trait == 'GE',type =='BLUE'),myGD = WinterGD, numfolds = 10,datasplit = .5,trait_name ='GE' )
# We are going to use our DH2020Estimates to fit time series. We need to GP the first TP of GE and GI values on the 
# half of the DHs that were not phenotyped. 
Taxacounts2020 = DH2020Estimates %>% filter(type == 'BLUE', trait == 'GE') %>% ungroup() %>% group_by(taxa) %>% 
  summarise(count = n()) %>% mutate(DataC = ifelse(count == 4,'Pred','Full')) %>%
  filter(taxa %in% WinterGD$taxa)
#everything has either 4 or 5 observations the ones with 4 need GP applied on them to get TP1 values. ######
GE_TP1_2020 = DH2020Estimates %>% filter(TP == 'TP1', trait == 'GE', type == 'BLUE') 
GE_TP1_2020_model_Train = mixed.solve(y = GE_TP1_2020 %>% arrange(taxa) %>% ungroup() %>%
                                        select(value) %>% as.matrix(),
                                      Z = WinterGD %>% filter(taxa%in%GE_TP1_2020$taxa)%>% arrange(taxa) %>%
                                        select(!taxa) %>% as.matrix()-1,
                                      K=NULL,
                                      SE=FALSE)
GE_TP1_2020_predValues = (as.matrix(WinterGD %>% filter(taxa %nin% GE_TP1_2020$taxa) %>% arrange(taxa) %>%
                                     select(!taxa) %>% as.matrix()-1) %*% as.matrix(GE_TP1_2020_model_Train$u)+
  as.numeric(GE_TP1_2020_model_Train$beta)) %>% as.data.frame() %>%
  rownames_to_column(var = 'taxa') %>% transmute(TP = 'TP1',PM_date = 5, trait = 'GE', taxa = taxa,
                                                 value = V1,type = 'BLUE',
                                                 Family = mapvalues(substr(taxa,1,3), 
                                                                    from = c('BS6','BS7','BS8','BS9','DH1','Fla','SY_','Sca','Win','KWS'), 
                                                                    to = c('Flavia/DH130910','Scala/DH130910',
                                                                           'SY_Tepee/DH130910','Wintmalt/DH130910',
                                                                           'DH130910','Flavia','SY_Tepee','Scala','Wintmalt','Scala'))) %>%
  filter(taxa %in% Taxacounts2020$taxa)%>% mutate(year = '2020')

GI_TP1_2020 = DH2020Estimates %>% filter(TP == 'TP1', trait == 'GI', type == 'BLUE') 
GI_TP1_2020_model_Train = mixed.solve(y = GI_TP1_2020 %>% arrange(taxa) %>% ungroup() %>%
                                        select(value) %>% as.matrix(),
                                      Z = WinterGD %>% filter(taxa%in%GI_TP1_2020$taxa)%>% arrange(taxa) %>%
                                        select(!taxa) %>% as.matrix()-1,
                                      K=NULL,
                                      SE=FALSE)
GI_TP1_2020_predValues = (as.matrix(WinterGD %>% filter(taxa %nin% GI_TP1_2020$taxa) %>% arrange(taxa) %>%
                                      select(!taxa) %>% as.matrix()-1) %*% as.matrix(GI_TP1_2020_model_Train$u)+
                            as.numeric(GI_TP1_2020_model_Train$beta)) %>% as.data.frame() %>%
  rownames_to_column(var = 'taxa') %>% transmute(TP = 'TP1',PM_date = 5, trait = 'GI', taxa = taxa,
                                                 value = V1, type = 'BLUE',
                                                 Family = mapvalues(substr(taxa,1,3), 
                                                                    from = c('BS6','BS7','BS8','BS9','DH1','Fla','SY_','Sca','Win','KWS'), 
                                                                    to = c('Flavia/DH130910','Scala/DH130910',
                                                                           'SY_Tepee/DH130910','Wintmalt/DH130910',
                                                                           'DH130910','Flavia','SY_Tepee','Scala','Wintmalt','Scala'))) %>%
  filter(taxa %in% Taxacounts2020$taxa) %>% mutate(year = '2020')
# Now we have predictions for our 2020 dataset. Move to fitting ######
# GE logfits #####
G2020.TSF.GE = DH2020Estimates%>% ungroup() %>% rbind(GI_TP1_2020_predValues, GE_TP1_2020_predValues) %>% 
  filter(type == 'BLUE', trait =='GE') %>%
  join(Taxacounts2020) %>% mutate(value = ifelse(value<0, 0, value))%>% 
  mutate(PM_date = PM_date-5,  year = '2020')

GE_logFits_trycatch = function(df, groupvars){
  Out = tryCatch(
    {suppressWarnings(broom::tidy(drm(value~PM_date, data = df,
                      fct=LL.4(fixed = c(NA, NA, 1, NA),
                               names = c('Rate','Lower','Upper','Centering'))))) %>%
        add_row(.,term = 'TimeTo95',curve ='Derived',
                estimate = as.double(exp(log((1-.[2,3])/(0.95-.[2,3])-1)/.[1,3]+log(.[3,3])))) %>%
         add_row(.,term = 'rTimeTo95',curve ='Derived',
             estimate = as.double(.[3,3]*exp(log((1-0.95)/0.95)/.[1,3])))
    }, error = function(e){
      data.frame()
    }
  )
  return(Out)
}
G2020.TSF.GE %>% group_by(taxa) %>% summarise(count = n()) %>% filter(count==4)
G2020.TSF.GI %>% group_by(taxa) %>% summarise(count = n()) %>% filter(count==4)

DHGE.logfits = G2020.TSF.GE %>% ungroup()  %>%
  # rbind(G2021.TSF.GE %>% mutate(PM_date = PM_date-5,year = '2021'),
  #        G2020_2021.TSF.GE %>% mutate(PM_date = PM_date-5,year = '2020/2021'))
  arrange(taxa, TP) %>% group_by(year, taxa) %>% group_modify(GE_logFits_trycatch) %>% ungroup() 
# There are a few that fail of course, but overall all looks really good. 
TaxaToFilter2020GE = DHGE.logfits %>% filter(year == '2020') %>% 
  filter(term=='Centering' & estimate> 200 | term=='TimeTo95' & estimate>250)
DHGE.logfits %>% filter(taxa %nin% TaxaToFilter2020GE$taxa) %>% ggplot(aes(x = estimate)) +geom_density()+facet_wrap(vars(term), scales = 'free')
# figure out what to set filters on for the 2021 and combined 2020/2021 data.sets
DH.GElogfitGWA.mlm = DHGE.logfits %>% filter(!(year=='2020' & taxa %in% TaxaToFilter2020GE$taxa)) %>%
  rename(value = estimate) %>%
  group_by(year,term) %>% group_modify(GWA_MLM_fortidyR)

DH.GElogfitGWA.mlmm = DHGE.logfits %>% filter(!(year=='2020' & taxa %in% TaxaToFilter2020GE$taxa)) %>%
  rename(value = estimate) %>%
  group_by(year,term) %>%  group_modify(GWA_MLMM_fortidyR)

DH.GElogfitGWA.mlm %>% group_by(year, term) %>% slice_head(n = 5) %>% view()
DH.GElogfitGWA.mlmm %>% group_by(year, term) %>% slice_head(n = 5) %>% view()
#lets Plot the fitted 
t = seq(1,150, by = 3)
DHGE.logestimates = DHGE.logfits %>% filter(term %in% c('Centering','Lower','Rate')) %>% group_by(year,taxa) %>% 
  group_modify(~{data.frame(time = t,GE_est = .x$estimate[2]+(1-.x$estimate[2])/(1+exp(.x$estimate[1]*(log(t)-log(.x$estimate[3])))))})
DHGE.logestimates %>% filter(taxa %nin% TaxaToFilter2020GE$taxa) %>% ggplot(aes(x = time, y = GE_est, group = taxa)) +geom_line()

# GI log fits #####
G2020.TSF.GI = DH2020Estimates%>% ungroup() %>% rbind(GI_TP1_2020_predValues, GE_TP1_2020_predValues) %>% 
  filter(type == 'BLUE', trait =='GI') %>% join(Taxacounts2020) %>% mutate(value = ifelse(value<0, 0, value))

GI_logfits_trycatch = function(df, groupvars) {
  Out = tryCatch(
    {suppressWarnings(broom::tidy(drm(value~PM_date, data = df,
                                       fct =LL.4(names = c('Rate','Lower','Upper','Centering'))))
    ) %>%
        add_row(.,term = 'TimeTo5.0',curve ='Derived',
                estimate = ifelse(.[2,3] > 5.0, 0,as.double((((.[3,3]-.[2,3])/(5.0-.[2,3])-1)^(1/.[1,3]))*.[4,3])))
    }, error = function(e){
      data.frame()
    }
  )
}

DHGIlogfits = G2020.TSF.GI %>% mutate(PM_date = PM_date-5, year = '2020') %>%
  # rbind(G2021.TSF.GI %>% mutate(PM_date = PM_date-5,year = '2021'),
  #        G2020_2021.TSF.GI %>% mutate(PM_date = PM_date-5,year = '2020/2021'))%>%
  arrange(taxa, TP) %>% group_by(year, taxa) %>%
  group_modify(GI_logfits_trycatch) %>% ungroup() 

taxaToFilter2020GI = DHGIlogfits %>% filter(year=='2020') %>% 
  filter(term == 'Centering' & estimate > 177 |
           term == 'Lower' & estimate < 0 |
           term == 'TimeTo5.0' & estimate > 250)

 DH.GIlogfitGWA.mlm = DHGIlogfits %>% filter(!(year == '2020' & taxa %in% taxaToFilter2020GI$taxa)) %>%
   # Add filtering steps here
   group_by(year,term) %>% rename(value = estimate) %>%
   group_modify(GWA_MLM_fortidyR)
 DH.GIlogfitGWA.mlm %>% group_by(term) %>% slice_head(n=10) %>% View()

 DH.GIlogfitGWA.mlmm = DHGIlogfits %>% 
   filter(!(year == '2020' & taxa %in% taxaToFilter2020GI$taxa)) %>%
   # Add filtering steps here
   group_by(year,term) %>%  rename(value = estimate) %>%
   group_modify(GWA_MLMM_fortidyR)
 
DHGI.logestimates = DHGIlogfits %>% filter(term %in% c('Centering','Rate','Upper','Lower')) %>% group_by(year, taxa) %>% 
   group_modify(~{data.frame(time = t,GI_est= .x$estimate[2]+(.x$estimate[3]-.x$estimate[2])/(1+exp(.x$estimate[1]*(log(t)-log(.x$estimate[4])))))})

DHGI.logestimates %>% filter(taxa %nin% taxaToFilter2020GI$taxa) %>% ggplot(aes(x = time, y = GI_est, group = taxa)) +geom_line()

DH.GIlogfitGWA.mlmm %>% group_by(year, term) %>% slice_head(n = 5) %>% view()

# FPCA on the 2020 data. ######
source('SpringBarley/Analysis/All_taxa/FPCA/Functions/FPCA_function.R')
source('SpringBarley/Analysis/All_taxa/FPCA/Functions/pca_fun.R')
source('SpringBarley/Analysis/All_taxa/FPCA/Functions/pca_score.R')
source('SpringBarley/Analysis/All_taxa/FPCA/Functions/tuning_nointer.R')

GI2020FPCA = G2020.TSF.GI %>% group_by(taxa) %>% dplyr::add_tally() %>% filter(n ==5) %>% ungroup() %>% 
  transmute(taxa = taxa, time = PM_date, GI = value) %>%
  FPCA_function(dfTaxaTraitTime = .,
                             Trait = 'GI', #Trait name must be entered as a character 
                             NumKnots = 0, # NumKnots is the number of interior knots to be fitted
                             order = 3, # Order is the dergree of the polynomial to be fit to the data.
                             NumObsevationsPerLine = 5)

GE2020FPCA = G2020.TSF.GE %>% group_by(taxa) %>% dplyr::add_tally() %>% filter(n ==5) %>% ungroup() %>% 
  transmute(taxa = taxa, time = PM_date, GE = value) %>%
  FPCA_function(dfTaxaTraitTime = .,
                Trait = 'GE', #Trait name must be entered as a character
                NumKnots = 0, # NumKnots is the number of interior knots to be fitted
                order = 3, # Order is the dergree of the polynomial to be fit to the data.
                NumObsevationsPerLine = 5)
GE2020FPCA$v1
GI2020FPCA$v1
GEGI2020FPCATimeTo = data.frame(time = GI2020FPCA$EstimatedAndEmpiricalMu$time[-c(1:5)]) %>% cbind(GI2020FPCA$RecoveredCurves)  %>% pivot_longer(cols = !time, names_to='taxa') %>%
  group_by(taxa) %>% group_modify(~{.x %>% filter(value>5.0) %>% slice_head(n=1)}) %>% ungroup() %>% mutate(Param = 'fTimeTo5.0', year= '2020', trait = 'GI') %>%
  rbind(data.frame(time = GE2020FPCA$EstimatedAndEmpiricalMu$time[-c(1:5)]) %>% cbind(GE2020FPCA$RecoveredCurves)  %>% pivot_longer(cols = !time, names_to='taxa') %>%
      group_by(taxa) %>% group_modify(~{.x %>% filter(value>.95) %>% slice_head(n=1)}) %>% ungroup() %>% mutate(Param = 'fTimeTo95', year= '2020', trait= 'GE')) %>%
  select(taxa, trait, year, Param, time) %>% rename(value=time)

#GWA using the FPCs without interpolation and penalized spline fits to the data
GWA2020.FPCA.GEGI = rbind(GI2020FPCA$PCs_withTaxa %>% select(taxa, FPC1,FPC2) %>% mutate(trait='GI', year= '2020') %>%
        pivot_longer(cols = c(FPC1,FPC2), values_to = 'value',names_to='Param'),
      GE2020FPCA$PCs_withTaxa %>% select(taxa, FPC1,FPC2) %>% mutate(trait='GE', year= '2020') %>%
        pivot_longer(cols = c(FPC1,FPC2), values_to = 'value',names_to='Param')) %>% rbind(GEGI2020FPCATimeTo) %>%
  group_by(year,trait,Param) %>%
  group_modify(GWA_MLMM_fortidyR)

GWA2020.FPCA.GEGI %>% slice_head(n=5) %>%view()


# Lets try to fit a penalized regression spline to our obsevations to try and interpolate between out points both ##############
# to give FPCA more to chew on, but also so we have more balanced datasets to combine later in analysis when we have 
tuning_nointer = function(lower, upper, Omega, Xmat, Y.vec, xlen){
  lam.list=exp(seq(lower,upper,1))
  gcv=rep(0,length(lam.list))
  for(ii in 1:length(lam.list))
  {
    A <- solve(t(Xmat)%*%Xmat+adiag(Omega*lam.list[ii]))
    Y.vec.hat <- (Xmat%*%A) %*% (t(Xmat)%*%Y.vec)
    diag.mean <- sum(diag(t(Xmat)%*%Xmat%*%A))/(dim(Xmat)[1]) # the original mean(diag(Hmat))
    gcv[ii] <- mean((Y.vec-Y.vec.hat)^2)/(1-diag.mean)^2
  }
  ind=which(gcv==min(gcv))
  lam.list[ind]
}
PenalizedSplineFitting = function(df, groupvars, K.int = 2,order = 4 ) {
  # K.int = number of interior knots
    tt = df$PM_date	# evaluation points
    J = length(tt)		# number of evaluation points
    t.min = min(tt)
    t.max = max(tt)
    ### basis functions ###
    # cubic splines, for second derivative, using order 6 rather than 4
    knots = t.min + (t.max-t.min)*(1:K.int)/(1+K.int)
    # knots = 45
    K = length(knots) + order # number of basis functions
    basis = create.bspline.basis(c(t.min,t.max),K,norder=order)
    
    # evaluate original, the first, and second derivatives of basis functions
    BS = eval.basis(tt,basis,0)
    BS1 = eval.basis(tt,basis,1)
    BS2 = eval.basis(tt,basis,2)
    
    # penalty matrix
    Omega = inprod(basis,basis,2,2)	# for second derivative, using 4 rather than 2
    
    # function for tuning parameter selection using GCV
    
    Y.vec = as.vector(df$value)
    N = length(Y.vec)
    # design matrix
    Xmat = BS
    
    dim(t(Xmat)%*%Xmat)
      
    ### Penalized least squares estimates ###
    lam = tuning_nointer(-10,15,Omega,Xmat,Y.vec,xlen)	# tunning parameter
    # lam = 1000
    # lam= diag(c(10000,10000,10000,10000,10000,10000,10000,10000,10000,10000))
    bhat = solve(t(Xmat)%*%Xmat+adiag(Omega * lam))%*%t(Xmat)%*%Y.vec
    n.eval = 50 #range(tt)[2]-range(tt)[1]	
    t.eval = sort(c(seq(t.min,t.max,by=(t.max-t.min)/n.eval),7,14,28, 42,63, 91, 115, 144))
    mu.eval = eval.basis(t.eval,basis,0)%*%bhat
    y_lim = range(Y.vec)
    
 plot(t.eval,mu.eval)+points(x = df[,c('PM_date','value')], col= 'blue')

  return(data.frame(time = t.eval, value = mu.eval))
}
df = DH2020Estimates %>% ungroup() %>% group_by(taxa) %>%  rbind(GI_TP1_2020_predValues, GE_TP1_2020_predValues) %>% 
  filter(type == 'BLUE')%>% dplyr::add_tally() %>% ungroup() %>% filter(n == 10)  %>% group_by(year, taxa, trait)%>%
  mutate(PM_date=PM_date-5) %>% 
  filter(taxa == 'DH130910', trait== 'GI')
  
DH2020SplineFits = DH2020Estimates %>% ungroup() %>% group_by(taxa) %>%  rbind(GI_TP1_2020_predValues, GE_TP1_2020_predValues) %>% 
  filter(type == 'BLUE')%>% dplyr::add_tally() %>% ungroup() %>% filter(n ==10) %>% ungroup()  %>% group_by(year, taxa, trait)%>% mutate(PM_date=PM_date-5) %>%
  group_modify(PenalizedSplineFitting)
DH2020SplineFits %>% ggplot(aes( x= time, y= value, group = taxa)) +
  geom_line()+facet_wrap(~trait, scales = 'free') +geom_hline(yintercept = 1, color = 'red')

DH2020SplineFits%>%ungroup() %>%
  join(DH2020Estimates %>%ungroup() %>% filter(type == 'BLUE') %>% transmute(taxa = taxa, trait = trait, time = PM_date-5, Tvalue = value)) %>%
  filter(!(is.na(Tvalue))) %>% group_by(trait, time) %>% group_modify(~{data.frame(Cor = cor(.$value,.$Tvalue))})


# FPCA on the b-spline fit stuff ######
# each of these take about 25 minutes to run so be careful. I am going to comment them out to prevent accedental running
# SplineFitFPCA.GE.2020 = DH2020SplineFits %>% ungroup() %>% group_by(taxa) %>% filter(trait=='GE') %>% add_tally() %>% filter(n==59) %>%
#   ungroup() %>% transmute(taxa = taxa, time = time, GE=value) %>%
#   FPCA_function(dfTaxaTraitTime = .,
#                 Trait = 'GE', #Trait name must be entered as a character ie Trait = 'GE3'
#                 NumKnots = 2, # NumKnots is the number of interior knots to be fitted
#                 order = 3, # Order is the dergree of the polynomial to be fit to the data.
#                 NumObsevationsPerLine = 59)
# SplineFitFPCA.GE.2020$phi.fun.plot
# 
# start_time = Sys.time()
# SplineFitFPCA.GI.2020 = DH2020SplineFits %>% ungroup() %>% group_by(taxa) %>% filter(trait=='GI') %>% add_tally() %>% filter(n==59) %>%
#   ungroup() %>% transmute(taxa = taxa, time = time, GI=value) %>%
#   FPCA_function(dfTaxaTraitTime = .,
#                 Trait = 'GI', #Trait name must be entered as a character ie Trait = 'GE3'
#                 NumKnots = 2, # NumKnots is the number of interior knots to be fitted
#                 order = 3, # Order is the dergree of the polynomial to be fit to the data.
#                 NumObsevationsPerLine = 59)
# stoptime=Sys.time()
print(stoptime-start_time)
SplineFitFPCA.GE.2020$v1
SplineFitFPCA.GI.2020$v1
#first two PCs look like they explain the most. 
GEGI2020FPCATimeTo.SplineFit = data.frame(time = SplineFitFPCA.GI.2020$EstimatedAndEmpiricalMu$time[-c(1:59)]) %>% cbind(SplineFitFPCA.GI.2020$RecoveredCurves)  %>% 
  pivot_longer(cols = !time, names_to='taxa') %>%
  group_by(taxa) %>% group_modify(~{.x %>% filter(value>4.5) %>% slice_head(n=1)}) %>% ungroup() %>% mutate(Param = 'fTimeTo5.0', year= '2020', trait = 'GI') %>%
  rbind(data.frame(time = SplineFitFPCA.GE.2020$EstimatedAndEmpiricalMu$time[-c(1:59)]) %>% cbind(SplineFitFPCA.GE.2020$RecoveredCurves)  %>% pivot_longer(cols = !time, names_to='taxa') %>%
          group_by(taxa) %>% group_modify(~{.x %>% filter(value>.95) %>% slice_head(n=1)}) %>% ungroup() %>% mutate(Param = 'fTimeTo95', year= '2020', trait= 'GE')) %>%
  select(taxa, trait, year, Param, time) %>% rename(value=time)

GEGI2020FPCATimeTo.SplineFit %>% ggplot(aes(x = value)) +facet_wrap(trait~Param) +geom_density()
GEGI2020FPCATimeTo%>% ggplot(aes(x = value)) +facet_wrap(trait~Param) +geom_density()

GWA2020.FPCA.GEGI.SplineFit = 
  rbind(SplineFitFPCA.GI.2020$PCs_withTaxa %>% select(taxa, FPC1,FPC2) %>% mutate(trait='GI', year= '2020') %>%
          pivot_longer(cols = c(FPC1,FPC2), values_to = 'value',names_to='Param'),
        SplineFitFPCA.GE.2020$PCs_withTaxa %>% select(taxa, FPC1,FPC2) %>% mutate(trait='GE', year= '2020') %>%
          pivot_longer(cols = c(FPC1,FPC2), values_to = 'value',names_to='Param')) %>% rbind(GEGI2020FPCATimeTo.SplineFit) %>%
  group_by(year,trait,Param) %>%
  group_modify(GWA_MLMM_fortidyR)

GWA2020.FPCA.GEGI.SplineFit %>% slice_head(n=5) %>%view()
# Now we have both the niave FPCA with only the 5 TP that we collected in 2020, and also the more advanced with interpolation according to
# a penalized spline fit, followed by FPCA. We can check the GWAs for differences. etc. 
# The next steps really can only take place once we have the 2021 data collection finished, as the log fits to the 2021 data need all timepoints
# and the combined as needs all the timepoint values. 
# What happens when we fit logistic models to the b-spline fitted data? would give us more interpreability than the fpcs? #####
DHGE.logfits.Spline = DH2020SplineFits %>% ungroup()  %>% filter(trait=='GE')  %>% group_by(year, trait, taxa) %>%  rename(PM_date = time)%>%
  group_modify(GE_logFits_trycatch) %>% ungroup() 
DHGE.logfits.Spline.taxaToFilter = DHGE.logfits.Spline %>% filter(term =='Centering' &estimate > 200 |
                                                                    term == 'TimeTo95' & estimate >250)
DH.GElogfitGWA.mlm.Spline = DHGE.logfits.Spline %>% filter(taxa %nin% DHGE.logfits.Spline.taxaToFilter$taxa) %>% rename(value= estimate)%>% 
  group_by(year, term) %>% group_modify(GWA_MLM_fortidyR) 
DH.GElogfitGWA.mlmm.Spline = DHGE.logfits.Spline %>% filter(taxa %nin% DHGE.logfits.Spline.taxaToFilter$taxa)%>% rename(value= estimate) %>% 
  group_by(year, term) %>% group_modify(GWA_MLMM_fortidyR) 

DH.GElogfitGWA.mlmm.Spline %>% slice_head(n=5) %>%view()

  
DHGI.logfits.Spline = DH2020SplineFits %>% ungroup() %>% filter(trait == 'GI') %>% group_by(year, trait, taxa) %>% rename(PM_date = time) %>%
  group_modify(GI_logfits_trycatch)
DHGI.logfits.Spline.taxaToFilter = DHGI.logfits.Spline %>% filter(term =='Centering' & estimate> 250 | term =='TimeTo5.0' & estimate > 300)

DH.GIlogfitGWA.mlm.Spline = DHGI.logfits.Spline %>% filter(taxa %nin% DHGI.logfits.Spline.taxaToFilter$taxa) %>% rename(value= estimate) %>%
  group_by(year, term) %>% group_modify(GWA_MLM_fortidyR) 
DH.GIlogfitGWA.mlmm.Spline = DHGI.logfits.Spline %>% filter(taxa %nin% DHGI.logfits.Spline.taxaToFilter$taxa) %>% rename(value= estimate) %>% 
  group_by(year, term) %>% group_modify(GWA_MLMM_fortidyR) 

DH.GIlogfitGWA.mlmm.Spline %>% slice_head(n=5) %>%view()


# Spline fits per PLOT and associted analysis 
PenalizedSplineFitting.temp = function(df, groupvars, K.int = 1,order = 3 ) {
  # K.int = number of interior knots
  # print(groupvars)
  tt = df$PM_date	# evaluation points
  J = length(tt)		# number of evaluation points
  t.min = min(tt)
  t.max = max(tt)
  ### basis functions ###
  # cubic splines, for second derivative, using order 6 rather than 4
  knots = t.min + (t.max-t.min)*(1:K.int)/(1+K.int)
  # knots = 45
  K = length(knots) + order # number of basis functions
  basis = create.bspline.basis(c(t.min,t.max),K,norder=order)
  
  # evaluate original, the first, and second derivatives of basis functions
  BS = eval.basis(tt,basis,0)
  BS1 = eval.basis(tt,basis,1)
  BS2 = eval.basis(tt,basis,2)
  
  # penalty matrix
  Omega = inprod(basis,basis,2,2)	# for second derivative, using 4 rather than 2
  
  # function for tuning parameter selection using GCV
  
  Y.vec = as.vector(df$value)
  N = length(Y.vec)
  # design matrix
  Xmat = BS
  
  ### Penalized least squares estimates ###
  lam = 50000#
  # lam = tuning_nointer(-10,15,Omega,Xmat,Y.vec,xlen)	# tunning parameter
  bhat = solve(t(Xmat)%*%Xmat+adiag(Omega*lam))%*%t(Xmat)%*%Y.vec
  # n.eval = 50 #range(tt)[2]-range(tt)[1]	
  t.eval = sort(c(seq(t.min,t.max,by=1))) #,7,14,28, 42,63)) #, 91, 115, 144))
  mu.eval = eval.basis(t.eval,basis,0)%*%bhat
  y_lim = range(Y.vec)
  
  # plot(t.eval,mu.eval)+points(x = df[,c('PM_date','value')], col= 'blue')
  return(data.frame(time = t.eval, value = mu.eval))
}

test.taxa = sample(unique(DHs2020$taxa),size = 20)

# we are going to take all the GE and GI predictions for TP1 2020 and use them in our model.
PerPlot.SplineFit = DHs2020 %>% select(year, Family, taxa, Location, PLOT, PM_date, GE, GI) %>%
  rbind(DHs2021 %>% rename(PLOT = SourcePLOT) %>% select(year, Family, taxa, Location,PLOT, PM_date, GE, GI)) %>%
  pivot_longer(cols = c(GE,GI), names_to = 'trait', values_to = 'value') %>% 
  rbind(rbind(GE_TP1_2020_predValues,GI_TP1_2020_predValues) %>% select(year, Family, taxa, trait, value, PM_date) %>%
          join(DHs2020 %>% filter(PM_date==19) %>% mutate(PM_date=5)%>% select(taxa, Location, PLOT, PM_date)) %>%
          select(year, Family, taxa, Location, PLOT, PM_date,trait, value))  %>%
  mutate(PM_date = PM_date -5) %>%
  group_by(year, Family, Location,taxa, PLOT, trait) %>% filter(!is.na(value) )%>% add_tally() %>% filter(n>6) %>%
  group_modify(PenalizedSplineFitting.temp)

PerPlot.SplineFit %>%# mutate(value = ifelse(trait =='GE' & value > 1,1,value)) %>%# filter(time %% 5 == 0) %>%
  ggplot(aes( x= time, y= value, group = PLOT)) + geom_hline(yintercept = 1, color = 'red')+
  geom_line()+facet_wrap(~trait,scales = 'free')

PerPlot.SplineFitTimeTo = PerPlot.SplineFit %>% group_modify(~{.x %>% filter(value>4.5) %>% slice_head(n=1)}) %>% 
  rbind(PerPlot.SplineFit %>% filter(trait == 'GE') %>% group_modify(~{.x %>% filter(value>.95) %>% slice_head(n=1)}))
BLUPH2(lmer(time~(1|taxa)+Location, data = PerPlot.SplineFitTimeTo %>% filter(trait == 'GI')))
BLUPH2(lmer(time~(1|taxa)+Location, data = PerPlot.SplineFitTimeTo %>% filter(trait == 'GE')))

BLUE_BLUPs.Spline.tests <- function(d, groupvars) {
  trait.lmer <- lmer(formula = value ~(1|taxa)+Location, 
                     data = d)
  
  lineEf = (ranef(trait.lmer)$taxa + fixef(trait.lmer)[1]) %>% as.data.frame() %>% rownames_to_column('taxa') %>% 
    mutate(type = 'BLUP') %>% rename(value = '(Intercept)')
  trait.lm = broom::tidy(lm(value~ taxa+Location, data=d))
  
  first_taxa = d %>% arrange(taxa) %>% slice_head( n = 1) %>% select(taxa) %>% as.character()
  Intercept = trait.lm %>% filter(term == '(Intercept)') %>% select(estimate) %>% as.numeric()
  lineBLUE = trait.lm %>% filter(substr(term,1,4)=='taxa') %>% 
    add_row(term = paste0('taxa',first_taxa),
            estimate = 0) %>% mutate(BLUE = estimate + Intercept) %>%
    transmute(taxa = gsub(pattern = 'taxa',replacement = '',x = term),
              value = BLUE,
              type = 'BLUE')
  H2 = BLUPH2(trait.lmer)
  return(rbind(lineEf, lineBLUE) %>% add_row(value = H2, type = 'H2') %>% arrange(type, taxa))
}
DH2020SplineFits %>% ggplot(aes( x= time, y= value, group = taxa)) +
  geom_line()+facet_wrap(~trait, scales = 'free')

spline.taxaBLH2 = PerPlot.SplineFit %>% ungroup() %>% mutate(value = ifelse(trait =='GE' & value > 1,1,value))%>% filter(time %% 5 == 0) %>%
  group_by(time, trait) %>% group_modify(BLUE_BLUPs.Spline.tests )
spline.taxaBLH2 %>% filter(type == 'H2') %>% ggplot(aes(time, value))+geom_line()+facet_wrap(~trait, scales = 'free')
spline.taxaBLH2 %>% filter(type == 'BLUE') %>% ggplot(aes(time, value, group = taxa)) +geom_line() + facet_wrap(~trait, scales = 'free')

df = PerPlot.SplineFit %>% filter(PLOT == 'a11450', trait == 'GI')

df = DHs2020 %>% select(year, Family, taxa, Location, PLOT, PM_date, GE, GI) %>%
  # rbind(DHs2021 %>% rename(PLOT = SourcePLOT) %>% select(year, Family, taxa, Location,PLOT, PM_date, GE, GI)) %>%
  pivot_longer(cols = c(GE,GI), names_to = 'trait', values_to = 'value') %>% mutate(PM_date = PM_date -5) %>%
  group_by(year, Family, Location, PLOT, trait) %>% tally()

df$n %>% table()
filter(PLOT == '6011', trait =='GE')#%>% group_modify(PenalizedSplineFitting)
# Genomic Prediction ################################################################
#All will use the data from DHCombined for the non .5 TP. 
WinterGPdata = DHCombined %>% filter(TP %in% c('TP1','TP2','TP3','TP4','TP5')) %>% filter(type == 'BLUE',Family != 'Cha')
# Question can we reduce SNP numbers? #######
# Based on Dans paper: yes ~3k gives a prediction platue within 7 families. What about ours?
dim(WinterGD)
MarkerNumberRandonTrainingTestingrrBLUP = function(df, groupvars,GD=WinterGD, numfolds = 50){ 
  Phenotype = df %>% filter(taxa %in% GD$taxa) %>% arrange(taxa) 
  myGD = GD %>% filter(taxa %in% Phenotype$taxa) %>% arrange(taxa)
  if(!length(unique(Phenotype$taxa))==sum(Phenotype$taxa == myGD$taxa)){
    stop('taxa lists are not correctly aligning')
  }
  start_time <- Sys.time() 
  set.seed(1)
  Qsd1 = myGD$Qsd1
  myGDrr = myGD[,-1]-1 
  Y = as.vector(Phenotype$value)
  num_entries = nrow(Phenotype)
  Num_train = .8*num_entries
  TrueandPredicted = data.frame()
  for (ii in c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)){
    for (i in 1:numfolds){
      trainingSet = sort(sample(1:num_entries,Num_train))
      testingSet = setdiff(1:num_entries,trainingSet)
      
      fam_test = Phenotype$Family[testingSet]
      fam_train = Phenotype$Family[trainingSet]
      
      y_train = Y[trainingSet] 
      y_test = Y[testingSet]
      # Qsd1 is 5148
      MarkerSet = sort(c(sample(x = 1:8384,size =round(8384*ii), replace = F),5148))
      marker_train = myGDrr[trainingSet,MarkerSet]
      marker_test = myGDrr[testingSet,MarkerSet]
      
      trained.Model = mixed.solve(y_train, Z = marker_train, K = NULL, SE =FALSE)
      
      PredictedPheno = as.matrix(marker_test) %*% as.matrix(trained.Model$u)+ as.numeric(trained.Model$beta)
      print(cor(y_test,PredictedPheno))
      
      pred = data.frame(Family = fam_test, Predicted = PredictedPheno, TrueValue = y_test, Qsd1 = 'AllAlleles') %>% 
        rbind(data.frame(Family = fam_test, Predicted = PredictedPheno, TrueValue = y_test, Qsd1 = Qsd1[testingSet])) %>%
        rbind(data.frame(Family = 'Overall', Predicted = PredictedPheno, TrueValue = y_test, Qsd1 = Qsd1[testingSet])) %>%
        rbind(data.frame(Family = 'Overall', Predicted = PredictedPheno, TrueValue = y_test, Qsd1 = 'AllAlleles')) %>%
        filter(Qsd1 != '1', Family != 'DH130910') %>%
        group_by(Family, Qsd1) %>%
        summarize(correlation = cor(Predicted, TrueValue),n_test = n()) %>%
        ungroup() %>%
        join(as.data.frame(table(Phenotype$Family[trainingSet])) %>%
               rename(Family= Var1, n_train =Freq) %>% filter(Family!='DH130910') %>%
               add_row(Family = 'Overall',n_train = Num_train)) %>% mutate(fold = i, MarkerProportion = ii)
      TrueandPredicted = rbind(TrueandPredicted, pred)
      
    }
  }
  return(TrueandPredicted)
}
start_time <- Sys.time() 
MarkerReductionTests = WinterGPdata %>% group_modify(MarkerNumberRandonTrainingTestingrrBLUP)
Sys.time() - start_time
MarkerReductionTests %>% ungroup() %>% filter(Family == 'Overall') %>%
  ggplot(aes(x = MarkerProportion, y = correlation, fill = as.factor(MarkerProportion)))+geom_boxplot()+
  facet_grid(trait~TP)
MarkerReductionTests 
# save(MarkerReductionTests, file = 'WinterBarley/Analysis/MarkerReductionTests.RData')

# Question: Which model to use for our germination data? ######
# rrblup, rrblup+fixed effect for Qsd1, Bayesian LASSO, BayesA, BayesB, BayesC, RKHS.
library(BGLR)

variousModelPAdetermination = function(df, groupvars, GD = WinterGD, numfolds = 50){
  setwd(rprojroot::find_rstudio_root_file())
  Phenotype = df %>% filter(taxa %in% GD$taxa) %>% arrange(taxa) 
  myGD = GD %>% filter(taxa %in% Phenotype$taxa) %>% arrange(taxa)
  if(!length(unique(Phenotype$taxa))==sum(Phenotype$taxa == myGD$taxa)){
    stop('taxa lists are not correctly aligning')
  }
  start_time <- Sys.time() 
  set.seed(1)
  Qsd1 = myGD$Qsd1
  myGDrr = myGD[,-1]-1 
  Y = as.vector(Phenotype$value)
  num_entries = nrow(Phenotype)
  PredictionsA = data.frame()
  for (ii in c(40,70,100,150,200,250,300,350,0.5,0.75,0.9)) {
    if(ii<1){
      Num_train = round(num_entries*ii,0)
    } else{
      Num_train = ii
    }
    for (i in 1:numfolds){
      trainingSet = sort(sample(1:num_entries,Num_train, replace = FALSE))
      testingSet = setdiff(1:num_entries,trainingSet)
      yNA = Y
      yNA[testingSet] <- NA
      
      # BRR
      ETA_RR = list(list(X=myGDrr, model = 'BRR'))
      RR = BGLR(y = yNA,ETA = ETA_RR,saveAt ='WinterBarley/Analysis/BGLROutput/') 
      #BRR+fixed effect of Qsd1
      ETA_RR_qsd1Fixed = list(list(X=Qsd1, model = 'FIXED'),list(X=myGDrr, model = 'BRR'))
      RR_qsd1_fixed = BGLR(y = yNA,ETA = ETA_RR_qsd1Fixed,saveAt ='WinterBarley/Analysis/BGLROutput/') 
      #LASSO
      ETA_LASSO = list(list(X = myGDrr, model = 'BL'))
      BL = BGLR(y = yNA,ETA = ETA_LASSO,saveAt ='WinterBarley/Analysis/BGLROutput/') 
      # LASSE with Qsd1
      ETA_LASSO_qsd1Fixed = list(list(X = Qsd1, model = "FIXED"), list(X = myGDrr, model = 'BL'))
      BL_qsd1_fixed = BGLR(y = yNA,ETA = ETA_LASSO_qsd1Fixed,saveAt ='WinterBarley/Analysis/BGLROutput/') 
      #BayesB
      ETA_BayesB = list(list(X = myGDrr, model = 'BayesB'))
      BayesB = BGLR(y = yNA, ETA = ETA_BayesB,saveAt ='WinterBarley/Analysis/BGLROutput/')
      #BayesB with Qsd1 Fixed
      ETA_BayesB_Qsd1Fixed = list(list(X = Qsd1, model = "FIXED"),list(X = myGDrr, model = 'BayesB'))
      BayesB_Qsd1Fixed = BGLR(y = yNA, ETA = ETA_BayesB_Qsd1Fixed,saveAt ='WinterBarley/Analysis/BGLROutput/')
      #BayesC
      ETA_BayesC = list(list(X = myGDrr, model = 'BayesC'))
      BayesC = BGLR(y = yNA, ETA = ETA_BayesC,saveAt ='WinterBarley/Analysis/BGLROutput/')
      #BayesC with Qsd1 Fixed
      ETA_BayesC_Qsd1Fixed = list(list(X = Qsd1, model = "FIXED"),list(X = myGDrr, model = 'BayesC'))
      BayesC_Qsd1Fixed = BGLR(y = yNA, ETA = ETA_BayesC_Qsd1Fixed,saveAt ='WinterBarley/Analysis/BGLROutput/')
      
      results = data.frame(Family= 'Overall',True = Y[testingSet], Qsd1Status = 'AllAlleles',
                           BRR = RR$yHat[testingSet], BRR_Qsd1_Fixed = RR_qsd1_fixed$yHat[testingSet],
                           LASSO = BL$yHat[testingSet], LASSO_Qsd1_Fixed = BL_qsd1_fixed$yHat[testingSet],
                           BayesB = BayesB$yHat[testingSet], BayesB_Qsd1_Fixed = BayesB_Qsd1Fixed$yHat[testingSet],
                           BayesC = BayesC$yHat[testingSet], BayesC_Qsd1_Fixed = BayesC_Qsd1Fixed$yHat[testingSet]) %>%
        rbind(data.frame(Family= 'Overall',True = Y[testingSet], Qsd1Status = Qsd1[testingSet],
                         BRR = RR$yHat[testingSet], BRR_Qsd1_Fixed = RR_qsd1_fixed$yHat[testingSet],
                         LASSO = BL$yHat[testingSet], LASSO_Qsd1_Fixed = BL_qsd1_fixed$yHat[testingSet],
                         BayesB = BayesB$yHat[testingSet], BayesB_Qsd1_Fixed = BayesB_Qsd1Fixed$yHat[testingSet],
                         BayesC = BayesC$yHat[testingSet], BayesC_Qsd1_Fixed = BayesC_Qsd1Fixed$yHat[testingSet])) %>%
        rbind(data.frame(Family= Phenotype$Family[testingSet],True = Y[testingSet],Qsd1Status = 'AllAlleles',
                         BRR = RR$yHat[testingSet], BRR_Qsd1_Fixed = RR_qsd1_fixed$yHat[testingSet],
                         LASSO = BL$yHat[testingSet], LASSO_Qsd1_Fixed = BL_qsd1_fixed$yHat[testingSet],
                         BayesB = BayesB$yHat[testingSet], BayesB_Qsd1_Fixed = BayesB_Qsd1Fixed$yHat[testingSet],
                         BayesC = BayesC$yHat[testingSet], BayesC_Qsd1_Fixed = BayesC_Qsd1Fixed$yHat[testingSet]))%>%
        rbind(data.frame(Family= Phenotype$Family[testingSet],True = Y[testingSet],Qsd1Status = Qsd1[testingSet],
                         BRR = RR$yHat[testingSet], BRR_Qsd1_Fixed = RR_qsd1_fixed$yHat[testingSet],
                         LASSO = BL$yHat[testingSet], LASSO_Qsd1_Fixed = BL_qsd1_fixed$yHat[testingSet],
                         BayesB = BayesB$yHat[testingSet], BayesB_Qsd1_Fixed = BayesB_Qsd1Fixed$yHat[testingSet],
                         BayesC = BayesC$yHat[testingSet], BayesC_Qsd1_Fixed = BayesC_Qsd1Fixed$yHat[testingSet])) %>%
        pivot_longer(cols = !c('Family','True', Qsd1Status), names_to = 'Model', values_to = 'Predicted') %>%
        filter(Family!='DH130910') %>%
          group_by(Family, Model, Qsd1Status) %>%
        summarise(correlation = cor(True,Predicted), n_test = n()) %>% ungroup() %>%
        join(as.data.frame(table(Phenotype$Family[trainingSet])) %>%
               rename(Family= Var1, n_train =Freq) %>% filter(Family!='DH130910') %>%
             add_row(Family = 'Overall',n_train = Num_train)) %>% 
        mutate(fold = i, trainingProportion = ii)
     
      PredictionsA = rbind(PredictionsA, results)     
    }
  }
return(PredictionsA)
}
# run on mac or other compter with multiple cores
ModelFrameWorkTests = WinterGPdata %>% group_by(trait, TP)%>% group_modify(variousModelPAdetermination) %>% collect()


# GP: RADOM TRAINING AND TEST SETS and grouped by Qsd1 ############
df = DH2020Estimates %>% ungroup() %>% filter(TP=='TP2',trait=='GI',type=='BLUE')%>% filter(Family != 'Cha') 
df$Family %>% table()

FoldCVGP_tidy_randomTrainS = function(df, groupvars, myGD=WinterGD, numfolds = 5){
  Phenotype = df %>% filter(taxa %in% myGD$taxa) %>% arrange(taxa) 
  myGD = myGD %>% filter(taxa %in% Phenotype$taxa) %>% arrange(taxa)
  dim(Phenotype)
  dim(myGD)
  if(!length(unique(Phenotype$taxa))==sum(Phenotype$taxa == myGD$taxa)){
    stop('taxa lists are not correctly aligning')
  }
  
  TrueandPredicted = data.frame()
  Qsd1 = myGD$Qsd1
  start_time <- Sys.time() 
  set.seed(1)
  myGD = myGD[,-1]-1 
  phenotype = as.vector(Phenotype$value)
  num_entries = nrow(Phenotype)

  for (iii in c(40,70,100,150,200,250,300,350,0.5,0.75,0.9)){
    if(iii<1){
      Num_train = round(num_entries*iii,0)
    } else{
      Num_train = iii
    }
    
    if( num_entries < iii){
      next
    }
    
    
  for (i in 1:numfolds){
    
    trainingSet = sort(sample(1:num_entries,Num_train))
    testingSet = setdiff(1:num_entries,trainingSet)
    
    fam_test = Phenotype$Family[testingSet]
    fam_train = Phenotype$Family[trainingSet]
    
    y_train = phenotype[trainingSet] 
    y_test = phenotype[testingSet]
    
    marker_train = myGD[trainingSet,]
    marker_test = myGD[testingSet,]
    
    trained.Model = mixed.solve(y_train, Z = marker_train, K = NULL, SE =FALSE)
    
    PredictedPheno = as.matrix(marker_test) %*% as.matrix(trained.Model$u)+ as.numeric(trained.Model$beta)
    print(cor(y_test,PredictedPheno))
    
    pred = data.frame(Family = fam_test, Predicted = PredictedPheno, TrueValue = y_test, Qsd1 = 'AllAlleles') %>% 
      rbind(data.frame(Family = fam_test, Predicted = PredictedPheno, TrueValue = y_test, Qsd1 = Qsd1[testingSet])) %>%
      rbind(data.frame(Family = 'Overall', Predicted = PredictedPheno, TrueValue = y_test, Qsd1 = Qsd1[testingSet])) %>%
      rbind(data.frame(Family = 'Overall', Predicted = PredictedPheno, TrueValue = y_test, Qsd1 = 'AllAlleles')) %>%
      filter(Qsd1 != '1', Family != 'DH130910') %>%
      group_by(Family, Qsd1) %>%
      summarize(correlation = cor(Predicted, TrueValue),n_test = n()) %>%
       ungroup() %>%
      join(as.data.frame(table(Phenotype$Family[trainingSet])) %>%
             rename(Family= Var1, n_train =Freq) %>% filter(Family!='DH130910') %>%
             add_row(Family = 'Overall',n_train = Num_train)) %>% mutate(fold = i, trainingProportion = iii)
    
    
    TrueandPredicted = rbind(TrueandPredicted, pred)
  }}
  end_time <- Sys.time()
  print(end_time-start_time)
  return(TrueandPredicted)
}

DH2020_predictions = DH2020Estimates %>% filter(Family != 'Cha')%>% ungroup() %>% 
  filter(!(trait =='GE' & TP %in% c('TP4','TP5'))) %>% 
  filter(type =='BLUE') %>% group_by(year, TP, PM_date, trait)%>%
  group_modify(FoldCVGP_tidy_randomTrainS) 

DH2020_predictions_Qsd1Group = DH2020Estimates %>% filter(Family != 'Cha') %>% ungroup() %>% filter(!(trait =='GE' & TP %in% c('TP4','TP5'))) %>% 
  join(WinterGD[,c('taxa','Qsd1')]) %>%
  filter(Qsd1 != 1) %>%  filter(type =='BLUE') %>% group_by(year, TP, PM_date, trait, Qsd1)%>%
  group_modify(FoldCVGP_tidy_randomTrainS) 

# GP: Stratified random Sampling #####
FoldCVGP_tidy_StratRandSampling = function(df, groupvars, myGD=WinterGD, numFolds = 50){
  print(groupvars)
  Phenotype = df %>% filter(taxa %in% myGD$taxa) %>% arrange(taxa) 
  myGD = myGD %>% filter(taxa %in% Phenotype$taxa) %>% arrange(taxa)
  dim(Phenotype)
  dim(myGD)
  if(!length(unique(Phenotype$taxa))==sum(Phenotype$taxa == myGD$taxa)){
    stop('taxa lists are not correctly aligning')
  }
  start_time <- Sys.time() 
  set.seed(1)
  myGD = myGD[,-1]-1
  
  TrueandPredicted = data.frame()
    num_entries = nrow(Phenotype)

  Familycounts = Phenotype %>% group_by(Family) %>% select(Family) %>% tally() %>% filter(Family != 'DH130910') %>%
    mutate(prop = n/num_entries)
  

  for (iii in c(40,70,100,150,200,250,300,350)){
      Num_train = iii-1
      
      if(num_entries < iii){
        next
      }
      
      Familycounts = Familycounts %>% mutate(NperFam = round(Familycounts$prop*Num_train,0))
    
   
      for (i in 1:numFolds){
        F1_t = sample(which(Phenotype$Family == 'Flavia/DH130910'),size = Familycounts$NperFam[1])
        F2_t = sample(which(Phenotype$Family == 'Scala/DH130910'),size = Familycounts$NperFam[2])
        F3_t = sample(which(Phenotype$Family == 'SY_Tepee/DH130910'),size = Familycounts$NperFam[3])
        F4_t = sample(which(Phenotype$Family == 'Wintmalt/DH130910'),size = Familycounts$NperFam[4])
        DH130910 = which(Phenotype$taxa=='DH130910')
        training_set = sort(c(F1_t,F2_t,F3_t,F4_t,DH130910))
        testing_set =setdiff(c(1:nrow(Phenotype)),training_set)
    
        y_train = Phenotype$value[training_set]
        y_test = Phenotype$value[testing_set]
        
        fam_test = Phenotype$Family[testing_set]
        fam_train = Phenotype$Family[training_set]

        marker_train = myGD[training_set,]
        marker_test = myGD[testing_set,]
        
        trained.Model = mixed.solve(y_train, Z = marker_train, K = NULL, SE =FALSE)
        
        PredictedPheno = as.matrix(marker_test) %*% as.matrix(trained.Model$u)+ as.numeric(trained.Model$beta)
        
        print(cor(y_test,PredictedPheno))
        
        pred = data.frame(Family = fam_test, Predicted = PredictedPheno, TrueValue = y_test) %>% group_by(Family) %>%
          add_tally() %>%
          group_modify(~{data.frame(correlation =cor(.$Predicted,.$TrueValue),n_test = .$n[1])}) %>%ungroup() %>%
          join(as.data.frame(table(fam_train)) %>% rename(Family= fam_train, n_train =Freq)) %>%
          add_row(Family = 'Overall', correlation = cor(y_test,PredictedPheno),n_test = num_entries-Num_train, n_train = Num_train+1) %>%
          mutate(fold = i, trainingProportion = iii)
        
        TrueandPredicted = rbind(TrueandPredicted,pred)
     
    }
  }
  return(TrueandPredicted)
}

DH2020StratRandomSam = DH2020Estimates %>% filter(Family != 'Cha') %>% filter(type=='BLUE') %>%
  filter(!(trait =='GE' & TP %in% c('TP4','TP5'))) %>% ungroup() %>%
  group_by(year, TP, trait) %>% group_modify(FoldCVGP_tidy_StratRandSampling)

# GP: structured GBLUP #####

# GP: CDmean Optimal training population #####
size <- c(148, 197, 246, 295, 344, 394, 443, 492, 541, 590, 640, 689, 738, 787, 836, 886, 935, 983)


# GP: Train with 3 families predict 1 #####
FoldCVGP_tidy_ThreeFamilysPredictOne= function(df, groupvars, myGD=WinterGD, numFolds = 50){
  print(groupvars)
  Phenotype = df %>% filter(taxa %in% myGD$taxa) %>% arrange(taxa) 
  myGD = myGD %>% filter(taxa %in% Phenotype$taxa) %>% arrange(taxa)
  dim(Phenotype)
  dim(myGD)
  if(!length(unique(Phenotype$taxa))==sum(Phenotype$taxa == myGD$taxa)){
    stop('taxa lists are not correctly aligning')
  }
  start_time <- Sys.time() 
  set.seed(1)
  myGD = myGD[,-1]-1

  TrueandPredicted = data.frame()
  for (ii in c("Flavia/DH130910", "Scala/DH130910", "SY_Tepee/DH130910", "Wintmalt/DH130910")){
    training_set = which(Phenotype$Family != ii)
    testing_set = which(Phenotype$Family == ii)
    for (iii in c(40,70,100,150,200,250,300,0.5,0.75,0.9)){
      if(iii<1){
        Num_train = round(length(training_set)*iii,0)
      } else{
        Num_train = iii
      }
      
      if(iii > length(training_set)){
        next
      }
      
      num_entries = nrow(Phenotype)
    
      for (i in 1:numFolds){
        train_rows = sort(sample(training_set,Num_train, replace = F))
        test_rows =sort(c(testing_set, setdiff(training_set,train_rows)))
        
        y_train = Phenotype$value[train_rows]
        y_test = Phenotype$value[test_rows]
        
        fam_test = Phenotype$Family[test_rows]
        fam_train = Phenotype$Family[train_rows]
        
        marker_train = myGD[train_rows,]
        marker_test = myGD[test_rows,]
        
        trained.Model = mixed.solve(y_train, Z = marker_train, K = NULL, SE =FALSE)
        
        PredictedPheno = as.matrix(marker_test) %*% as.matrix(trained.Model$u)+ as.numeric(trained.Model$beta)
        
        print(cor(y_test,PredictedPheno))
        
        #DH130910 will always be excluded from the per family correlations, but included in the overall correlation. 
        pred = data.frame(Family = fam_test, Predicted = PredictedPheno, TrueValue = y_test) %>% group_by(Family) %>%
          add_tally() %>%
          group_modify(~{data.frame(correlation =cor(.$Predicted,.$TrueValue),n_test = .$n[1])}) %>%ungroup() %>%
          join(as.data.frame(table(fam_train)) %>% rename(Family= fam_train, n_train =Freq)) %>%
          add_row(Family = 'Overall', correlation = cor(y_test,PredictedPheno),n_test = num_entries-Num_train, n_train = Num_train) %>%
          mutate(fold = i, trainingProportion = iii, FamilyPredicted = ii)
        
        TrueandPredicted = rbind(TrueandPredicted,pred)
        }
    }
  }
    return(TrueandPredicted)
}

DH20203FamPred1 = DH2020Estimates %>% filter(Family != 'Cha') %>% filter(type=='BLUE') %>%
  filter(!(trait =='GE' & TP %in% c('TP4','TP5'))) %>% ungroup() %>%
  group_by(year, TP, trait) %>% group_modify(FoldCVGP_tidy_ThreeFamilysPredictOne)

DH20203FamPred1 %>% filter(n_test>4) %>% filter(trainingProportion==.75) %>%
  mutate(Family = factor(Family, levels = c('Overall',"Flavia/DH130910", "Scala/DH130910", "SY_Tepee/DH130910", "Wintmalt/DH130910")),
         yfacetlabel = 'Prediction accuracy for:') %>%
  ggplot(aes(x =TP, y = correlation, fill = FamilyPredicted)) +geom_boxplot()+facet_nested(yfacetlabel+Family~trait, scales = 'free',space = 'free_x')+
  geom_hline(yintercept = 0)+
  labs(title = 'Train with three families, predict outgroup',fill = 'Family Excluded \nfrom training set')+
  theme_bw()

DH20203FamPred1 %>% filter(n_test>4) %>% filter(trainingProportion==70) %>% filter(trait =='GI',TP == 'TP2') %>%
  mutate(Family = factor(Family, levels = c('Overall',"Flavia/DH130910", "Scala/DH130910", "SY_Tepee/DH130910", "Wintmalt/DH130910")),
         yfacetlabel = 'Prediction accuracy for:') %>% group_by(Family, FamilyPredicted) %>% summarize(correlation = mean(correlation)) %>%
  ggplot(aes(x =Family, y = FamilyPredicted, fill = correlation)) +geom_tile()+ ylab('Family Excluded from training') +
  labs(title = 'Train with three families, predict outgroup')+geom_text(aes(label = round(correlation,3)), color = 'white')+
  theme_bw()


# GP: Train with 2 families predict 2 #####
Family_combos = combn(c("Flavia/DH130910", "Scala/DH130910", "SY_Tepee/DH130910", "Wintmalt/DH130910"),2)
FoldCVGP_tidy_twoFamilysPredictTwo = function(df, groupvars, myGD=WinterGD, numFolds = 50){
  print(groupvars)
  Phenotype = df %>% filter(taxa %in% myGD$taxa) %>% arrange(taxa) 
  myGD = myGD %>% filter(taxa %in% Phenotype$taxa) %>% arrange(taxa)
  dim(Phenotype)
  dim(myGD)
  if(!length(unique(Phenotype$taxa))==sum(Phenotype$taxa == myGD$taxa)){
    stop('taxa lists are not correctly aligning')
  }
  start_time <- Sys.time() 
  set.seed(1)
  myGD = myGD[,-1]-1
  Family_combos = combn(c("Flavia/DH130910", "Scala/DH130910", "SY_Tepee/DH130910", "Wintmalt/DH130910"),2)
  
  TrueandPredicted = data.frame()
  for (ii in c(1:6)){
    testing_set = which(Phenotype$Family %nin% Family_combos[,ii])
    training_set = which(Phenotype$Family %in% Family_combos[,ii])
    
      for (iii in c(40,70,100,150,200,0.5,0.75,0.9)){
      if(iii<1){
        Num_train = round(length(training_set)*iii,0)
      } else{
        Num_train = iii
      }
      
      if(iii > length(training_set)){
        next
      }
      
      for (i in 1:numFolds){
        train_rows = sort(sample(training_set,Num_train, replace = F))
        test_rows = sort(c(testing_set, setdiff(training_set,train_rows)))
        
        y_train = Phenotype$value[train_rows]
        y_test = Phenotype$value[test_rows]
        
        marker_train = myGD[train_rows,]
        marker_test = myGD[test_rows,]
        
        fam_train = Phenotype$Family[train_rows]
        fam_test = Phenotype$Family[test_rows]
        
        trained.Model = mixed.solve(y_train, Z = marker_train, K = NULL, SE =FALSE)
        
        PredictedPheno = as.matrix(marker_test) %*% as.matrix(trained.Model$u)+ as.numeric(trained.Model$beta)
        
        #DH130910 will always be excluded from the per family correlations, but included in the overall correlation. 
        pred = data.frame(Family = fam_test, Predicted = PredictedPheno, TrueValue = y_test) %>% group_by(Family) %>%
          add_tally() %>%
          group_modify(~{data.frame(correlation =cor(.$Predicted,.$TrueValue), n_test = .$n[1])}) %>%ungroup() %>%
          join(as.data.frame(table(fam_train)) %>% rename(Family= fam_train, n_train =Freq)) %>%
          add_row(Family = 'Overall', correlation = cor(y_test,PredictedPheno), n_train = Num_train, n_test = length(training_set)-Num_train) %>%
          mutate(fold = i, trainingProportion = iii, FamilyTrain = ii)
        
        TrueandPredicted = rbind(TrueandPredicted,pred)
        print(cor(y_test,PredictedPheno))
      }
    }
  }
  return(TrueandPredicted)
}

DH20202FamPred2 = DH2020Estimates %>% filter(Family != 'Cha') %>% filter(type=='BLUE') %>%
  filter(!(trait =='GE' & TP %in% c('TP4','TP5'))) %>% ungroup() %>%
  group_by(year, TP, trait) %>% group_modify(FoldCVGP_tidy_twoFamilysPredictTwo)%>%
    mutate(trainingFamilies = mapvalues(FamilyTrain, from = c(1:6),
                                        to = c(paste0(Family_combos[1,1],'\n', Family_combos[2,1]),
                                               paste0(Family_combos[1,2],'\n', Family_combos[2,2]),
                                               paste0(Family_combos[1,3],'\n', Family_combos[2,3]),
                                               paste0(Family_combos[1,4],'\n', Family_combos[2,4]),
                                               paste0(Family_combos[1,5],'\n', Family_combos[2,5]),
                                               paste0(Family_combos[1,6],'\n', Family_combos[2,6]))))

DH20202FamPred2 %>% ggplot(aes(x = TP, y = correlation, color = trainingFamilies))+geom_boxplot()+
  facet_grid(trainingFamilies~trait, scales = 'free_x')+ylim(0,1)

DH20202FamPred2 %>% filter(n_test>4) %>% filter(trainingProportion==.75) %>%
  separate(trainingFamilies, into =c('TF1','TF2'), sep = '\n', remove = F) %>%
  mutate(Family = factor(Family, levels = c('Overall',"Flavia/DH130910", "Scala/DH130910", "SY_Tepee/DH130910", "Wintmalt/DH130910")),
         yfacetlabel = 'Prediction accuracy for:',
         IntrainingSet = (Family == TF1 |Family==TF2)) %>%
  ggplot(aes(x =TP, y = correlation, fill =  IntrainingSet)) +geom_boxplot()+
  facet_nested(yfacetlabel+Family~trait, scales = 'free',space = 'free_x')+
  geom_hline(yintercept = 0)+
  labs(title = 'Train with two families, predict other two', 
       subtitle = 'Prediction accuracy per family along side facets',
       fill = 'Family predicted\nin the training set?')+theme_bw()

DH20202FamPred2 %>% filter(n_test>4) %>% filter(trainingProportion==.75) %>%
  separate(trainingFamilies, into =c('TF1','TF2'), sep = '\n', remove = F) %>%
  mutate(Family = factor(Family, levels = c('Overall',"Flavia/DH130910", "Scala/DH130910", "SY_Tepee/DH130910", "Wintmalt/DH130910")),
         yfacetlabel = 'Prediction accuracy for:',
         IntrainingSet = (Family == TF1 |Family==TF2)) %>%
  ggplot(aes(x =TP, y = correlation, fill =  trainingFamilies)) +geom_boxplot()+
  facet_nested(yfacetlabel+Family~trait, scales = 'free',space = 'free_x')+
  geom_hline(yintercept = 0)+theme_bw()+
  labs(fill = 'Families in training set', title = 'Train on two families Predict other two')


# GP: Train with 1 family   predict 3 #######
FoldCVGP_tidy_OneFamilysPredictThree = function(df, groupvars, myGD=WinterGD, numFolds = 50){
  print(groupvars)
  Phenotype = df %>% filter(taxa %in% myGD$taxa) %>% arrange(taxa) 
  myGD = myGD %>% filter(taxa %in% Phenotype$taxa) %>% arrange(taxa)
  dim(Phenotype)
  dim(myGD)
  if(!length(unique(Phenotype$taxa))==sum(Phenotype$taxa == myGD$taxa)){
    stop('taxa lists are not correctly aligning')
  }
  start_time <- Sys.time() 
  set.seed(1)
  myGD = myGD[,-1]-1
  
  TrueandPredicted = data.frame()
  for (ii in c("Flavia/DH130910", "Scala/DH130910", "SY_Tepee/DH130910", "Wintmalt/DH130910")){
    testing_set = which(Phenotype$Family != ii)
    training_set = which(Phenotype$Family == ii)
    for (iii in c(40,70,100,0.5,0.75,0.9)){
      if(iii<1){
        Num_train = round(length(training_set)*iii,0)
      } else{
        Num_train = iii
      }
      
      if(iii > length(training_set)){
        next
      }
      
      for (i in 1:numFolds){
        train_rows = sort(sample(training_set,Num_train, replace = F))
        test_rows = sort(c(testing_set, setdiff(training_set,train_rows)))
        
        y_train = Phenotype$value[train_rows]
        y_test = Phenotype$value[test_rows]
        
        marker_train = myGD[train_rows,]
        marker_test = myGD[test_rows,]
        
        fam_train = Phenotype$Family[train_rows]
        fam_test = Phenotype$Family[test_rows]
        
        trained.Model = mixed.solve(y_train, Z = marker_train, K = NULL, SE =FALSE)
        
        PredictedPheno = as.matrix(marker_test) %*% as.matrix(trained.Model$u)+ as.numeric(trained.Model$beta)
        
        #DH130910 will always be excluded from the per family correlations, but included in the overall correlation. 
        pred = data.frame(Family = fam_test, Predicted = PredictedPheno, TrueValue = y_test) %>% group_by(Family) %>%
          add_tally() %>%
          group_modify(~{data.frame(correlation =cor(.$Predicted,.$TrueValue), n_test = .$n[1])}) %>%ungroup() %>%
          join(as.data.frame(table(fam_train)) %>% rename(Family= fam_train, n_train =Freq)) %>%
          add_row(Family = 'Overall', correlation = cor(y_test,PredictedPheno), n_train = Num_train, n_test = length(training_set)-Num_train) %>%
          mutate(fold = i, trainingProportion = iii, FamilyTrain = ii)
        
        TrueandPredicted = rbind(pred,TrueandPredicted)
                
        print(cor(y_test,PredictedPheno))
       }
    }
  }
  return(TrueandPredicted)
}

DH20201FamPred3 = DH2020Estimates %>% filter(Family != 'Cha') %>% filter(type=='BLUE') %>%
  filter(!(trait =='GE' & TP %in% c('TP4','TP5'))) %>% ungroup() %>%
  group_by(year, TP, trait) %>% group_modify(FoldCVGP_tidy_OneFamilysPredictThree)

DH20201FamPred3 %>% filter(n_test>5) %>%
  mutate(Family = factor(Family, levels = c('Overall',"Flavia/DH130910", "Scala/DH130910", "SY_Tepee/DH130910", "Wintmalt/DH130910")),
         yfacetLabel = 'Family used for training') %>% ggplot(aes(x = TP, y = correlation, fill = FamilyTrain))+
  geom_boxplot()+theme_bw()+ labs(title = 'Train with one famliy predict other three', fill = 'Training Family')+
  facet_nested(yfacetLabel + Family~trait, scales = 'free',space = 'free_x')+geom_hline(yintercept = 0)

DH20201FamPred3 %>% ungroup()%>% filter(trait == 'GI', TP == 'TP2') %>% filter(trainingProportion == 70.00) %>% filter(Family != 'DH130910') %>%
  group_by(year, TP, trait, Family, FamilyTrain) %>% summarise(correlation = mean(correlation)) %>%
  mutate(Family = factor(Family, levels = c("Flavia/DH130910", "Scala/DH130910", "SY_Tepee/DH130910", "Wintmalt/DH130910",'Overall'))) %>%
  ggplot(aes(x = Family, y = FamilyTrain, fill = correlation, label = round(correlation,3)))+geom_tile()+geom_text()

# GP Plotting of single timepoint and follow-up analysis ######
# Combined Predictions stratified random sample and random sampling  plot
RandSampleStratCombined = DH2020_predictions  %>% filter(Family %nin% c('Cha','DH130910')) %>% 
  group_by(Family, trait, TP,trainingProportion) %>% 
  summarise(correlation = mean(correlation,na.rm = T), stdev = sd(correlation, na.rm = T), n_train = mean(n_train)) %>%
  filter(!(TP=='TP1' &trainingProportion==200)) %>% mutate(type = 'Random Sampling') %>%
  rbind(DH2020StratRandomSam %>% group_by(Family, trait, TP, trainingProportion) %>% 
          summarise(stdev = sd(correlation, na.rm = T),correlation = mean(correlation,na.rm = T), n_train = mean(n_train)) %>%
          filter(!(TP=='TP1' &trainingProportion==200))%>% mutate(type = 'Stratified Sampling')) %>% 
  mutate(Family = factor(Family, levels = c('Overall',"Flavia/DH130910", "Scala/DH130910", "SY_Tepee/DH130910", "Wintmalt/DH130910")))

png('Presentations/Graphs/RandSample_vsStrat.png', 1400, 800, res =120)
RandSampleStratCombined %>% ggplot(aes(x = n_train, y = correlation, color = TP, linetype = type)) +geom_line(size = 1) + 
  geom_errorbar(aes(ymin = correlation - stdev, ymax = correlation+stdev))+theme_bw()+
  facet_grid(trait~Family, scales = 'free_x')+theme_bw() +ylim(0,1) +xlab('Number of training taxa within group') +
  scale_x_continuous(sec.axis = sec_axis(~./max(.)*(400*1.045), name = 'Total training population size') )+
  labs(linetype = 'Sampling type', color = 'Timepoint')+geom_hline(yintercept = 0)
dev.off()

DH2020_predictions %>% filter(Family %nin% c('Cha','DH130910')) %>% 
  group_by(trait, TP, Family, Qsd1, trainingProportion) %>%
  summarise(correlation = mean(correlation, na.rm = T), n_train = mean(n_train))%>% filter(trait=='GI') %>%
  ggplot(aes(x = n_train, y = correlation, color = Qsd1))+geom_line()+
  facet_grid(trait+TP~Family, scales = 'free_x')+theme_bw() +ylim(0,1) +xlab('Number of training taxa within group')


# Predictions within Qsd1 grouping. Dormant is more predictable group. 
png('Presentations/Graphs/Qsd1GroupingPredictions.png', 1400, 800, res =120)
DH2020_predictions_Qsd1Group %>% group_by(Family, trait, TP,trainingProportion, Qsd1) %>% filter(trainingProportion<149) %>%
  summarise(correlation = mean(correlation,na.rm = T), stdev = sd(correlation, na.rm = T), n_train = mean(n_train)) %>%
  filter(!(TP=='TP1' &trainingProportion==200)) %>% filter(Family != 'DH130910') %>%
  mutate(Family = factor(Family, levels = c('Overall',"Flavia/DH130910", "Scala/DH130910", "SY_Tepee/DH130910", "Wintmalt/DH130910")),
         Qsd1lab = ifelse(Qsd1 == 0, 'Non Dormant', 'Dormant'),
         training = 'Qsd1Specific') %>%
 ggplot(aes(x = n_train, y = correlation, color = TP)) +geom_line(aes(linetype = Qsd1lab), size = 1) + 
  facet_nested(trait~Family, scales = 'free_x')+theme_bw() +geom_hline(yintercept = 0)+theme_bw()+
  labs(linetype = 'Qsd1 status', color = 'Timepoint')+xlab('Number of training taxa within group')+
  scale_x_continuous(sec.axis = sec_axis(~./max(.)*(250*1.045), name = 'Total training population size'))+
  ylim(0,1)
dev.off()



png('Presentations/Graphs/RandSampleVsQsd1Strat.png', 1400, 800, res =120)
DH2020_predictions_Qsd1Group %>% group_by(Family, trait, TP,trainingProportion, Qsd1) %>% filter(trainingProportion<149) %>%
  summarise(correlation = mean(correlation,na.rm = T), stdev = sd(correlation, na.rm = T), n_train = mean(n_train)) %>%
  filter(!(TP=='TP1' &trainingProportion==200)) %>% filter(Family != 'DH130910') %>%
  mutate(Family = factor(Family, levels = c('Overall',"Flavia/DH130910", "Scala/DH130910", "SY_Tepee/DH130910", "Wintmalt/DH130910")),
         Qsd1lab = ifelse(Qsd1 == 0, 'Non Dormant', 'Dormant')) %>%select(!Qsd1) %>%
  rbind(RandSampleStratCombined %>% filter(type == 'Random Sampling') %>% rename(Qsd1lab = type)) %>%
  filter(TP %in% c('TP2','TP3')) %>%
  ggplot(aes(x = n_train, y = correlation, color = TP)) +geom_line(aes(linetype = Qsd1lab), size = 1) + 
  facet_nested(trait~Family, scales = 'free_x')+theme_bw() +geom_hline(yintercept = 0)+theme_bw()+
  labs(linetype = 'Qsd1 status', color = 'Timepoint')+xlab('Number of training taxa within group')+
  scale_x_continuous(sec.axis = sec_axis(~./max(.)*(250*1.045), name = 'Total training population size'))+
  ylim(0,1)
dev.off()

GITP22020= FoldCVGPv2(df = DH2020Estimates %>% filter(TP=='TP2', trait == 'GI',type =='BLUE'),myGD = WinterGD, numfolds = 10,datasplit = .5,trait_name ='GE' )
bysd1 = GITP22020 %>% join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1 != 1) %>% group_by(Qsd1) %>% summarise(cor = round(cor(TruePheno, PredPheno),3))
overall = GITP22020 %>% join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1 != 1) %>% summarise(cor = round(cor(TruePheno, PredPheno),3))
png('Presentations/Graphs/ExampleGroupPredictions.png', 1400, 800, res =120)
GITP22020 %>% join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1 != 1) %>% mutate(Qsd1 = ifelse(Qsd1==0,'NonDormant','Dormant')) %>%
  ggplot(aes(x = TruePheno, y = PredPheno, color = as.factor(Qsd1)))+geom_point()+theme_bw() +
  labs(color = 'Qsd1 Status', title = 'GI TP2, 2020 data predictions, training size = 200')+
  annotate(geom='text',x = 1, y = 6, label = paste('Overall:', overall,'\nQsd1 Dormant:',bysd1$cor[2], '\nQsd1 Non-Dormant:',bysd1$cor[1]))
dev.off()  

# three families predicting one plot
DH20203FamPred1 %>% filter(n_test>4) %>% filter(trainingProportion==.75) %>%
  mutate(Family = factor(Family, levels = c('Overall',"Flavia/DH130910", "Scala/DH130910", "SY_Tepee/DH130910", "Wintmalt/DH130910")),
         yfacetlabel = 'Prediction accuracy for:') %>%
  ggplot(aes(x =TP, y = correlation, fill = FamilyPredicted)) +geom_boxplot()+facet_nested(yfacetlabel+Family~trait, scales = 'free',space = 'free_x')+
  geom_hline(yintercept = 0)+
  labs(title = 'Train with three families, predict outgroup',fill = 'Family Excluded \nfrom training set')+
  theme_bw()

DH20203FamPred1 %>% filter(n_test>4) %>% filter(trainingProportion %in% c(40,70,100,150)) %>% filter(trait =='GI') %>% 
  filter(TP %in% c('TP1','TP3','TP5')) %>%
  mutate(Family = factor(Family, levels = c('Overall',"Flavia/DH130910", "Scala/DH130910", "SY_Tepee/DH130910", "Wintmalt/DH130910")),
         yfacetlabel = 'Training Size', xfacetlabel = 'Correlation within grouping for GI prediction') %>% 
  ggplot(aes(x = TP, y = correlation, fill = FamilyPredicted))+geom_boxplot()+facet_nested(yfacetlabel + trainingProportion~ xfacetlabel+Family)+
  labs(fill = 'Family excluded from\nthe training population')+ylim(0,1)+geom_hline(yintercept = c(0,0.5), color = 'black')

png('Presentations/Graphs/PAThreeFamiliesPredictOne.png', 1400, 800, res =120)
DH20203FamPred1 %>% filter(n_test>7) %>% filter(trainingProportion %in% c(70,100,150, 200)) %>% filter(trait =='GI') %>% 
  group_by(trait,trainingProportion, TP, Family, FamilyPredicted) %>% summarize(correlation = mean(correlation)) %>%
  mutate(yfacetlabel = 'Training Size', xfacetlabel = 'Correlation within grouping for GI prediction',
         Family = gsub(Family, pattern = '/DH130910', replacement = ''),
         Family = factor(Family, levels = c('Overall',"Flavia", "Scala", "SY_Tepee", "Wintmalt")),
         FamilyPredicted = gsub(FamilyPredicted, pattern = '/DH130910',replacement = ''),
         Box = ifelse(FamilyPredicted == Family, 'Y','N')) %>% 
  ggplot(aes(x = Family,  y= FamilyPredicted, fill = correlation))+
  geom_tile(aes(color = Box), height = .9, width = .9, size = 1)+
  scale_color_manual(values = c('white','black'), guide = 'none')+ 
  facet_nested(yfacetlabel+trainingProportion~TP)+theme_bw()+
  scale_fill_gradient2(limits=c(0,1),low = 'blue',high = 'red', midpoint =.5)+  
  geom_text(aes(label = round(correlation,2)), size = 2.5)+theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))+
  labs(y = 'Family excluded from model training', x = 'Family Predicted', title = 'Train on three families predict the other') 
dev.off()
# Two Families predict two 

DH20202FamPred2 %>% filter(n_test>4) %>% filter(trainingProportion %in% c(70,100,150)) %>% filter(trait == 'GI') %>%
  group_by(trait,trainingProportion, TP, Family, trainingFamilies) %>% summarize(correlation = mean(correlation)) %>%
  mutate(trainingFamilies= gsub(trainingFamilies, pattern = '/DH130910',replacement = '')) %>%
  separate(trainingFamilies, into =c('TF1','TF2'), sep = '\n', remove = F) %>%
  mutate(Family = gsub(Family, pattern = '/DH130910', replacement = ''),
         Family = factor(Family, levels = c('Overall',"Flavia", "Scala", "SY_Tepee", "Wintmalt")),
         yfacetlabel = 'Training population size',
         IntrainingSet = (Family == TF1 | Family==TF2),
         Box = ifelse(IntrainingSet==TRUE,'Y','N')) %>%
  ggplot(aes( x= Family, y =trainingFamilies, fill = correlation))+
  geom_tile(aes(color = Box), height = .9, width = .9, size = 1.5)+
  scale_color_manual(values = c('white','black'), guide = 'none')+ 
  facet_nested(yfacetlabel+trainingProportion~TP)+theme_bw()+
  scale_fill_gradient2(limits=c(0,1),low = 'blue',high = 'red', midpoint =.5)+  
  geom_text(aes(label = round(correlation,2)), size = 3)+theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
  

# one family train predict the other three
png('Presentations/Graphs/PAOneFamiliyPredictthree.png', 1400, 800, res =120)
DH20201FamPred3 %>% filter(n_test>4) %>% filter(trainingProportion %in% c(40,70,100)) %>% filter(trait =='GI') %>% 
  group_by(trait,trainingProportion, TP, Family, FamilyTrain) %>% summarize(correlation = mean(correlation)) %>%
  mutate(yfacetlabel = 'Training Size', xfacetlabel = 'Correlation within grouping for GI prediction',
         Family = gsub(Family, pattern = '/DH130910', replacement = ''),
         Family = factor(Family, levels = c('Overall',"Flavia", "Scala", "SY_Tepee", "Wintmalt")),
         FamilyTrain = gsub(FamilyTrain, pattern = '/DH130910',replacement = ''),
         Box = ifelse(FamilyTrain == Family, 'Y','N')) %>%
  ggplot(aes(x = Family,  y= FamilyTrain, fill = correlation))+
  geom_tile(aes(color = Box), height = .9, width = .9, size = 1)+
  scale_color_manual(values = c('white','black'), guide = 'none')+ 
  facet_nested(yfacetlabel+trainingProportion~TP)+theme_bw()+
  scale_fill_gradient2(limits=c(0,1),low = 'blue',high = 'red', midpoint =.5)+  
  geom_text(aes(label = round(correlation,2)), size = 2.5)+theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))+
  labs(y = 'Family model trained on', x = 'Family Predicted', title = 'Train on one family predict the others') 
dev.off()

# combine by training set 1->3 2->2, 3->1 and examine!
DH20201FamPred3 %>% filter(n_test>4) %>% filter(trainingProportion %in% c(70)) %>% filter(trait =='GI') %>% 
  group_by(trait,trainingProportion, TP, Family, FamilyTrain) %>% summarize(correlation = mean(correlation)) %>%
  mutate(yfacetlabel = 'Training Size', xfacetlabel = 'Correlation within grouping for GI prediction',
         Family = gsub(Family, pattern = '/DH130910', replacement = ''),
         Family = factor(Family, levels = c('Overall',"Flavia", "Scala", "SY_Tepee", "Wintmalt")),
         FamilyTrain = gsub(FamilyTrain, pattern = '/DH130910',replacement = ''),
         Box = ifelse(FamilyTrain == Family, 'Y','N'),
         predictionGroup = "1->3") %>%
  select(Family,FamilyTrain, correlation, predictionGroup, Box,yfacetlabel,trainingProportion,TP) %>%
  rbind(DH20202FamPred2 %>% filter(n_test>4) %>% filter(trainingProportion %in% c(70)) %>% filter(trait == 'GI') %>%
          group_by(trait,trainingProportion, TP, Family, trainingFamilies) %>% summarize(correlation = mean(correlation)) %>%
          mutate(trainingFamilies= gsub(trainingFamilies, pattern = '/DH130910',replacement = '')) %>%
          separate(trainingFamilies, into =c('TF1','TF2'), sep = '\n', remove = F) %>%
          mutate(Family = gsub(Family, pattern = '/DH130910', replacement = ''),
                 Family = factor(Family, levels = c('Overall',"Flavia", "Scala", "SY_Tepee", "Wintmalt")),
                 yfacetlabel = 'Training Size',
                 IntrainingSet = (Family == TF1 | Family==TF2),
                 Box = ifelse(IntrainingSet==TRUE,'Y','N'),
                 FamilyTrain = trainingFamilies,
                 predictionGroup = "2->2") %>%
          select(Family, FamilyTrain,correlation, predictionGroup, Box, yfacetlabel, trainingProportion,TP)) %>%
  rbind(DH20203FamPred1 %>% filter(n_test>7) %>% filter(trainingProportion %in% c(70)) %>% filter(trait =='GI') %>% 
          group_by(trait,trainingProportion, TP, Family, FamilyPredicted) %>% summarize(correlation = mean(correlation)) %>%
          mutate(yfacetlabel = 'Training Size', xfacetlabel = 'Correlation within grouping for GI prediction',
                 Family = gsub(Family, pattern = '/DH130910', replacement = ''),
                 Family = factor(Family, levels = c('Overall',"Flavia", "Scala", "SY_Tepee", "Wintmalt")),
                 FamilyPredicted = gsub(FamilyPredicted, pattern = '/DH130910',replacement = ''),
                 FamilyTrain = mapvalues(FamilyPredicted, from = c('Flavia','Wintmalt','Scala','SY_Tepee'), 
                                         to = c('Wintmalt\nScala\nSY_Tepee', 
                                                'Flavia\nScala\nSY_Tepee',
                                                'Flavia\nWintmalt\nSY_Tepee',
                                                'Flavia\nWintmalt\nScala')),
                 Box = ifelse(FamilyPredicted == Family, 'N','Y'),
                 predictionGroup = '3->1') %>%
          select(Family, FamilyTrain,correlation, predictionGroup, Box, yfacetlabel, trainingProportion,TP)) %>%
  ggplot(aes(x = Family,  y= FamilyTrain, fill = correlation))+
  geom_tile(aes(color = Box), height = .9, width = .9, size = 1)+
  scale_color_manual(values = c('white','black'), guide = 'none')+ 
  facet_nested(yfacetlabel+predictionGroup+trainingProportion~TP, scales = 'free_y', space = 'free')+theme_bw()+
  scale_fill_gradient2(limits=c(0,1),low = 'blue',high = 'red', midpoint =.5)+  
  geom_text(aes(label = round(correlation,2)), size = 2.5)+theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))+
  labs(y = 'Family model trained on', x = 'Family Predicted', title = 'Train on n families (see side facet) predict the others') 


DH20203FamPred1 %>% filter(n_test>7) %>% filter(trainingProportion %in% c(70)) %>% filter(trait =='GI') %>% 
  group_by(trait,trainingProportion, TP, Family, FamilyPredicted) %>% summarize(correlation = mean(correlation)) %>%ungroup() %>%
  mutate(yfacetlabel = 'Training Size', xfacetlabel = 'Correlation within grouping for GI prediction',
         Family = gsub(Family, pattern = '/DH130910', replacement = ''),
         Family = factor(Family, levels = c('Overall',"Flavia", "Scala", "SY_Tepee", "Wintmalt")),
         FamilyPredicted = gsub(FamilyPredicted, pattern = '/DH130910',replacement = ''),
         FamilyTrain = mapvalues(FamilyPredicted, from = c('Flavia','Wintmalt','Scala','SY_Tepee'), 
                                 to = c('Wintmalt\nScala\nSY_Tepee', 
                                        'Flavia\nScala\nSY_Tepee',
                                        'Flavia\nWintmalt\nSY_Tepee',
                                        'Flavia\nWintmalt\nScala')),
         Box = ifelse(FamilyPredicted == Family, 'Y','N'),
         predictionGroup = '3->1') %>% select(FamilyTrain) %>% unique()



# GP: Trait values over time? Need to think about that... ####
# can run similar things to stratified random sampling etc to predict trait values, but not sure here at this point. 
# Also - Should I bring in other traits since this is turning into a genomic prediction paper?
# Plot of Relationship vs the trait value deviance. Also GBLUP vs RRBLUP #####
WinterRelationship = rrBLUP::A.mat(WinterGD[,-1]-1, impute.method = 'EM', return.imputed = F)
Trait_dev_Sample = DH2020Estimates %>% filter(Family != 'Cha') %>% filter(type=='BLUE') %>%
  filter(!(trait =='GE' & TP %in% c('TP4','TP5'))) %>% ungroup() %>%
  filter(TP=='TP2', trait == 'GI')%>% arrange(taxa) %>% ungroup()
# We need one trait to examine trait deviance and the relationship. Use TP2, 2020, GI. 
WinterRelationship[upper.tri(WinterRelationship,diag = TRUE)] <- NA

DHrelation = WinterRelationship %>%  data.frame() %>% rownames_to_column(var = 'taxa') %>%
  pivot_longer(cols = !taxa, names_to = 'taxa2', values_to = 'Kinship') %>% filter(taxa %in% Trait_dev_Sample$taxa) %>%
  filter(taxa2 %in% Trait_dev_Sample$taxa)  %>% 
  filter(!is.na(Kinship)) %>%
  join(Trait_dev_Sample[,c('taxa','value', 'Family')]) %>%
  join(Trait_dev_Sample %>% transmute(taxa2 = taxa, value2 = value, Family2 = Family)) %>% 
  mutate(Family_comparisons = paste0(Family,'\n',Family2),
         trait_dev = value-value2) %>% join(WinterGD[,c('taxa','Qsd1')])%>%
  join(WinterGD[c('taxa','Qsd1')] %>% rename(taxa2 = taxa, Qsd1_2 = Qsd1))
DHrelation %>% select(Family, Family2, Kinship) %>% group_by(Family,Family2) %>%
  summarise(meanK = mean(Kinship))

jpeg(file = 'Presentations/Graphs/TraitdevVkinship.jpg', 1200,900, res = 120)
DHrelation %>% filter(Family %in% c("Flavia/DH130910", "Scala/DH130910", "SY_Tepee/DH130910", "Wintmalt/DH130910"),
                      Family2 %in% c("Flavia/DH130910", "Scala/DH130910", "SY_Tepee/DH130910", "Wintmalt/DH130910"))%>% 
  mutate(Qsd1MonoMorph = Qsd1==Qsd1_2, trait_dev= abs(trait_dev)) %>%
  ggplot(aes(Kinship, y = trait_dev))+geom_point(size = .01)+
  facet_grid(Family~Family2, scales = 'free') +geom_smooth(method = 'lm')
dev.off()

# Lets look at prediction of 20% of the entire Trait_dev_Sample n.
Model.RR = Trait_dev_Sample %>% group_by(year, trait, TP) %>% group_modify(FoldCVGP_tidy_randomTrainS)
Model.RR %>% filter(TrainingPropor==.8) %>% summarise(mean(correlation), sd(correlation))

library(qgg)

K_use = rrBLUP::A.mat((WinterGD %>% filter(taxa%in%Trait_dev_Sample$taxa)%>% select(!taxa))-1, impute.method = 'EM', return.imputed = F)
Trained.Model.GBLUP = mixed.solve(y = (Trait_dev_Sample %>% filter(taxa%in% WinterGD$taxa) %>% select(value) %>% as.matrix())[training_set] ,
                                  K=K_use,
                                  SE=FALSE)

# PHS import and data processing. May need to run with ASREML.R #####
phenotypes
#C:\Users\kars8\git\Cornell-WMB21-selections\data\phenotypes\Germination\2021
getwd()
#33_o

DH_Sprout = read_excel('data/Phenotype_Data/2021Phenotyping/Data/DHs_PHS_2021.xlsx') %>%
  mutate(location = ifelse(PLOT>6999, 'Ketola','McGowan'), year = '2021') %>% 
  # rbind(read_excel('WinterBarley/PhenotypeData/2020Harvest/2020WinterDH_PHS.xlsx')%>%
  #         mutate(location = 'Ketola2020'), year = '2020')  #Not sure if this should be included as these were planed as facultatives. 
  select(!Comment) %>% mutate(taxa = gsub(pattern = '-', replacement = '_',Entry),taxa = gsub(pattern = ' ', replacement = '_',taxa)) %>%
  rename(score = `Sprout Score`) %>% 
  separate(score, into =c('p0','p1','p2','p3','p4','p5'), sep = '') %>% select(!p0) %>% pivot_longer(cols = c(p1,p2,p3,p4,p5)) %>%
  mutate(value = as.numeric(value)) %>%
  group_by(taxa, location, Harv, year) %>% summarise(PHS = mean(value,na.rm = T))
PHS.lmer = lmer(PHS ~ location +(1|location:Harv)+(1|taxa), DH_Sprout) 
VarCorr(PHS.lmer)
BLUPH2(PHS.lmer)
DH.PHS = (ranef(PHS.lmer)$taxa + fixef(PHS.lmer)[1]) %>% as.data.frame() %>% rownames_to_column('taxa') %>% rename(PHS =`(Intercept)`)

# Plot of TimeTo5.0/TimeTo95 vs TP1 GE/GI (with predicted values) and PHS? #####
library(GGally)
GEGIPHSTimeTo2020 = rbind(GE_TP1_2020_predValues,GI_TP1_2020_predValues, 
      DH2020Estimates %>% filter(TP == 'TP1',type == 'BLUE')) %>%
  pivot_wider(values_from = value, names_from = trait) %>%
  join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1 != 1) %>%
  join(DH.PHS) %>%
  join(GEGI2020FPCATimeTo.SplineFit %>% pivot_wider(names_from =c('trait','Param'), values_from = 'value')) %>%
  rename(GE_TP1 = GE, GI_TP1 = GI, fTimeTo5.0 = GI_fTimeTo5.0, fTimeTo95 = GE_fTimeTo95) %>% 
  mutate(Qsd1 = ifelse(Qsd1==2,'Dormant','NonDormant'),
         GE_TP1 = ifelse(GE_TP1<0,0,GE_TP1),
         GI_TP1 = ifelse(GI_TP1<0,0,GI_TP1)) %>% filter(!is.na(fTimeTo5.0), !is.na(fTimeTo95))
  
GEGIPHSTimeTo2020 %>% ggpairs(data = ., 
          columns = c('fTimeTo5.0','fTimeTo95','GE_TP1', 'GI_TP1', 'PHS'),
          ggplot2::aes(color = Qsd1))+theme_bw()+theme(axis.text = element_text(size = 6))

# Genetic Correlations? #####
DH2021Estimates %>% filter(taxa %nin% WinterGD$taxa) %>% select(taxa )%>% unique()
DH2020Estimates%>% filter(taxa %nin% WinterGD$taxa) %>% select(taxa )%>% unique()


















