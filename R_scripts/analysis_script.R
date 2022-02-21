
#Workflow:Winter_DH_data_consolidation.Rmd ->analysis_script.R->WDH_GenotypeData_Processing.Rmd


library(here)
library(lme4)
library(reshape2)
library(sommer);library(arm);library(lme4);library(Hmisc);library(plyr);library(readxl);
library(tibble);library(patchwork);library(ggplot2);library(fda) ; library(magic); 
library(drc);library(rrBLUP);library(tidyr);library(ggh4x);library(dplyr)
load("data/WMB_DH_Germination_data2020_2021.RData")#DH2020_2021
load("data/WMB_DH_Germination_data2020_2021_wide_format.RData")#DH2020_2021_wide

load(here("data/Winter_DH/WMB DH 2020 all data for analysis.RData"))#old file
#linear model analysis for germination traits

#correl
DH2020_2021$Check<-as.factor(DH2020_2021$Check)
table(DH2020_2021$spot_bloch)
#DH 2020 GI
anova(lm(GI~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2020",DH2020_2021$PM_date ==5)))
anova(lm(GI~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2020",DH2020_2021$PM_date ==19)))
#Location significant
anova(lm(GI~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2020",DH2020_2021$PM_date ==47)))
#Location significant
anova(lm(GI~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2020",DH2020_2021$PM_date ==96)))
#replication significant
anova(lm(GI~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2020",DH2020_2021$PM_date ==152)))
#all significant
#DH 2020 GE
anova(lm(GE~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2020",DH2020_2021$PM_date ==5)))
#Location
anova(lm(GE~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2020",DH2020_2021$PM_date ==19)))
#Location significant
anova(lm(GE~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2020",DH2020_2021$PM_date ==47)))
#Location significant and rep
anova(lm(GE~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2020",DH2020_2021$PM_date ==96)))
#location and rep significant
anova(lm(GE~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2020",DH2020_2021$PM_date ==152)))
#field analysis 2020

anova(lm(GE~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2020",DH2020_2021$PM_date ==152)))
#field analysis 2020
DH2020_2021$Maturity_date<-DH2020_2021$Maturity_date+76
colnames(DH2020_2021)
DH2020_2021$fac_field
anova(lm(GI~ taxa +Location+Maturity_date, DH2020_2021%>% filter(DH2020_2021$Year=="2020",DH2020_2021$PM_date ==5)))#Location significant

#Does it make sense that Maturity date is still significant even though we controlled for it?
anova(lm(GI~ taxa +Location+replication+Maturity_date, DH2020_2021%>% filter(DH2020_2021$Year=="2020",DH2020_2021$PM_date ==19)))

#Location significant, winter survival significant
anova(lm(GI~ taxa +Location+replication+Maturity_date, DH2020_2021%>% filter(DH2020_2021$Year=="2020",DH2020_2021$PM_date ==47)))

#Location significant and rep
anova(lm(GI~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2020",DH2020_2021$PM_date ==96)))
#location and rep significant, scald significant for GE, not GI
anova(lm(GI~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2020",DH2020_2021$PM_date ==152)))
anova(lm(Day1Germ~ taxa +Maturity_date:Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2020",DH2020_2021$PM_date ==152)))
#field analysis 2020
#Summation, we need to consider Maturity date in the model, it affects GE at all timepoints and GI for timepoints 1, 2 and 3
#scald affected GE for timepoint 4
##############2021 
#2021 GI
table(DH2020_2021$PM_date,DH2020_2021$Year)#add 12, 33, 68
DH2020_2021$spot_bloch
#scald associated with GI at timepoint 3.5 or so
#scald associated with GE at timepoint 19,33, and 47
anova(lm(GE~ taxa +Maturity_date:Location+scald, DH2020_2021%>% filter(DH2020_2021$Year=="2021",DH2020_2021$PM_date ==19)))
#Location
anova(lm(GI~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2021",DH2020_2021$PM_date ==12)))
#Location
anova(lm(GI~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2021",DH2020_2021$PM_date ==19)))
#Location
anova(lm(GI~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2021",DH2020_2021$PM_date ==33)))
#Location significant
anova(lm(GI~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2021",DH2020_2021$PM_date ==47)))
#neither significant
anova(lm(GI~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2021",DH2020_2021$PM_date ==68)))
#Location significant
anova(lm(GI~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2021",DH2020_2021$PM_date ==96)))
#rep and location significant
anova(lm(GI~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2021",DH2020_2021$PM_date ==152)))
#all significant
#DH 2021 GE

anova(lm(GE~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2021",DH2020_2021$PM_date ==5)))
#Location
anova(lm(GE~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2021",DH2020_2021$PM_date ==12)))
#Location

anova(lm(GE~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2021",DH2020_2021$PM_date ==19)))
#neither
anova(lm(GE~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2021",DH2020_2021$PM_date ==33)))
#rep slightly
anova(lm(GE~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2021",DH2020_2021$PM_date ==47)))
#rep sign
anova(lm(GE~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2021",DH2020_2021$PM_date ==68)))
#Location and rep
anova(lm(GE~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2021",DH2020_2021$PM_date ==96)))
#location
anova(lm(GE~ taxa +Location+replication, DH2020_2021%>% filter(DH2020_2021$Year=="2021",DH2020_2021$PM_date ==152)))
anova(lm(GI~ taxa +Location+replication+spot_bloch,DH2020_2021%>% filter(DH2020_2021$Year=="2021",DH2020_2021$PM_date ==96)))
library(plyr)
DH2020_2021$spot_bloch
lm(GI~taxa + scald,DH2020_2021%>%filter(DH2020_2021$Year=="2021",DH2020_2021$PM_date ==12))

#2020 and 2021 together

table(DH2020_2021$PM_date,DH2020_2021$Year)
anova(lm(GI~ taxa +Year+Year %in% Location+replication+ws_percent+Maturity_date:Year, DH2020_2021%>% filter(DH2020_2021$PM_date ==5)))
#Year significant

#rep not significant
anova(lm(GI~ taxa +replication+Year+Maturity_date%in%Year+Maturity_date, DH2020_2021%>% filter(DH2020_2021$PM_date ==19)))
#Year and year in location significant
anova(lm(GI~ taxa +Year%in%Location+replication+Year, DH2020_2021%>% filter(DH2020_2021$PM_date ==47)))
#Location and year significant
anova(lm(GI~ taxa +Year%in%Location+replication+Year, DH2020_2021%>% filter(DH2020_2021$PM_date ==96)))
#replication, year significant
anova(lm(GI~ taxa +Location%in%Year+replication+Year, DH2020_2021%>% filter(DH2020_2021$PM_date ==152)))
#all significant

#GE for both years
anova(lm(GE~ taxa +Year+Year %in% Location+replication, DH2020_2021%>% filter(DH2020_2021$PM_date ==5)))
#Year significant, Year in locatiobn slightly
anova(lm(GE~ taxa +Year%in%Location+replication+Year, DH2020_2021%>% filter(DH2020_2021$PM_date ==19)))
#Year and year in location significant
anova(lm(GE~ taxa +Year%in%Location+replication+Year, DH2020_2021%>% filter(DH2020_2021$PM_date ==47)))
#Year significant, Year:Location slightly
anova(lm(GI~ taxa +Year%in%Location+replication+Year, DH2020_2021%>% filter(DH2020_2021$PM_date ==96)))
#replication, year significant, Location slightly
anova(lm(GI~ taxa +Location%in%Year+replication+Year, DH2020_2021%>% filter(DH2020_2021$PM_date ==152)))
#all significant
#We should include rep, location within Year, and Year as effects for models


DH2020_2021$TP
library(ggplot2)
library(tidyr)
library(ggh4x)
library(patchwork)
DH2020_2021 %>% pivot_longer(cols = c(GE,GI), names_to = 'trait') %>% filter(trait =='GE') %>%
  ggplot(aes(x = value, fill = TP))+geom_density()+facet_nested(TP ~trait+Year, scales = 'free')+
 theme_bw()+guides(fill = "none")+
  DH2020_2021 %>% pivot_longer(cols = c(GE,GI), names_to = 'trait') %>% filter(trait =='GI') %>%
  ggplot(aes(x = value, fill = TP))+geom_density()+facet_nested(TP ~trait+Year, scales = 'free')+
  theme_bw()+ guides(fill = "none")
  plot_layout(ncol = 2)

#comments about the data table
#each entry is replicated 8 times to cover 4 TPs for each replication
#This means that phenotype data for each entry is replicated that many times
#This shouldn't be a problem if comparing phenotype data within a timepoint,
# but when comparing between ensure that total number of observations are correct

#other notes from the data with field
#scald seems to have a significant effect in the models, this could mean that either scald resistance/susceptibility
  #affects germination OR that the resis/sucespt loci is closely associated with AlaAT
  
#####
  
DH2020_2021=DH2020_2021%>%dplyr::rename(rep=replication,year=Year)
DH2020_2021=DH2020_2021 %>%filter(!Location=="Spring")
BLUPH2 = function(trait.lm) {
  ses<- se.ranef(trait.lm)$'taxa' #where 'm' is your model object from 'lmer' (replace 'genotypes' with whatever you call your individuals in the data)
  v_BLUP<- ses^2
  sigma2_g=VarCorr(trait.lm, comp="Variance")$'taxa'[1]
  Reliability<- 1- v_BLUP/ (2*sigma2_g)  #where sigma2_g is the genetic variance estimated
  H2<- round(mean(Reliability),3)
  return(H2)
}

BlueBlupsH2_Location_rep_taxa <- function(d, groupvars) {

    trait.lmer <- lmer(formula = value ~(1|taxa)+Location+rep, 
                       data = d)
    #in the future we should add maturity date and possibly scald here
      lineEf = (ranef(trait.lmer)$taxa + fixef(trait.lmer)[1])%>% as.data.frame() %>%rownames_to_column('taxa') %>% 
    dplyr::mutate(type = 'BLUP') %>% dplyr::rename(value = '(Intercept)')
      
    trait.lm = broom::tidy(lm(value~ taxa+Location+rep, data=d))
    #in the future we should add maturity date and possibly scald here
    first_taxa = d %>% arrange(taxa) %>% slice_head( n = 1) %>% dplyr::select(taxa) %>% as.character()
    Intercept = trait.lm %>% filter(term == '(Intercept)') %>% dplyr::select(estimate) %>% as.numeric()
    lineBLUE = trait.lm %>% filter(substr(term,1,4)=='taxa') %>% 
      add_row(term = paste0('taxa',first_taxa),
              estimate = 0) %>% mutate(BLUE = estimate + Intercept) %>%
      transmute(taxa = gsub(pattern = 'taxa',replacement = '',x = term),
                value = BLUE, type = 'BLUE')
    H2 = BLUPH2(trait.lmer)
    return(rbind(lineEf, lineBLUE) %>% add_row(value = H2, type = 'H2') %>% arrange(type, taxa) %>%
             join(d %>% dplyr::select(taxa, Family) %>% unique()))
  }
  
d2 = DH2020_2021 %>%dplyr::select(taxa, rep, Location, TP, GE,GI,PM_date,year) %>% 
    mutate(year = factor(year, levels = c('2021','2020'))) %>%
    pivot_longer(cols = c(GE, GI), names_to = 'trait') %>% filter(TP == 'TP2', trait == 'GE')
  
BlueBlupsH2_Year_rep_taxa <- function(d2, groupvars) {
    if (length(unique(d2$year))==2) {
      trait.lmer <- lmer(formula = value ~(1|taxa)+Location + rep, data = d2)
      # fixef(trait.lmer)[2]/need to think about this more. 
      
      Cept = (fixef(trait.lmer)[1]*4+sum(fixef(trait.lmer)[2:4]))/4
      lineEf = (ranef(trait.lmer)$taxa + Cept) %>% as.data.frame() %>% rownames_to_column('taxa') %>% 
        mutate(type = 'BLUP') %>% dplyr::rename(value = '(Intercept)')
      
      trait.lm = broom::tidy(lm(value~ taxa+Location+ rep, data=d2))
      first_taxa = d2 %>% arrange(taxa) %>% slice_head( n = 1) %>% dplyr::select(taxa) %>% as.character()
      Intercept_list = trait.lm %>% filter(term == '(Intercept)'|substr(term,1,8)=='Location')
      Intercept = (Intercept_list$estimate[1]*4+sum(Intercept_list$estimate[2:4]))/4
      
      lineBLUE = trait.lm %>% filter(substr(term,1,4)=='taxa') %>% 
        add_row(term = paste0('taxa',first_taxa),
                estimate = 0) %>% mutate(BLUE = estimate + Intercept) %>%
        transmute(taxa = gsub(pattern = 'taxa',replacement = '',x = term),
                  value = BLUE, type = 'BLUE')
      H2 = BLUPH2(trait.lmer)
      return(rbind(lineEf, lineBLUE) %>% add_row(value = H2, type = 'H2')%>% join(d2 %>% dplyr::select(taxa, Family)%>% unique()))
    }
    else {
      trait.lmer <- lmer(formula = value ~(1|taxa)+Location+rep, 
                         data = d2)
      lineEf = (ranef(trait.lmer)$taxa + fixef(trait.lmer)[1]) %>% as.data.frame() %>% rownames_to_column('taxa') %>% 
        mutate(type = 'BLUP') %>% dplyr::rename(value = '(Intercept)')
      trait.lm = broom::tidy(lm(value~ taxa+Location+rep, data=d2))
      
      first_taxa = d2 %>% arrange(taxa) %>% slice_head( n = 1) %>% dplyr::select(taxa) %>% as.character()
      Intercept = trait.lm %>% filter(term == '(Intercept)') %>% dplyr::select(estimate) %>% as.numeric()
      lineBLUE = trait.lm %>% filter(substr(term,1,4)=='taxa') %>% 
        add_row(term = paste0('taxa',first_taxa),
                estimate = 0) %>% mutate(BLUE = estimate + Intercept) %>%
        transmute(taxa = gsub(pattern = 'taxa',replacement = '',x = term),
                  value = BLUE,type = 'BLUE')
      H2 = BLUPH2(trait.lmer)
      return(rbind(lineEf, lineBLUE) %>% add_row(value = H2, type = 'H2') %>% arrange(type, taxa) %>%
               join(d2 %>% dplyr::select(taxa, Family) %>% unique()))
    }
  }
library(lme4)

DH2020Estimates = DH2020_2021%>%filter(year=="2020") %>% dplyr::select(taxa, rep, Location,TP, GE, GI, PM_date,year,Family) %>% 
    pivot_longer(cols = c(GE, GI), names_to = 'trait') %>%
    group_by(TP, PM_date, trait, year) %>% group_modify(BlueBlupsH2_Location_rep_taxa) %>% ungroup()
DH2021Estimates = DH2020_2021%>%filter(year=="2021") %>% dplyr::select(taxa, rep, Location, TP, GE,GI,PM_date,year,Family) %>% pivot_longer(cols = c(GE, GI), names_to = 'trait') %>%
  group_by(TP,PM_date, trait, year) %>% group_modify(BlueBlupsH2_Location_rep_taxa) %>% ungroup()
#Both
DHCombined =  DH2020_2021 %>% dplyr::select(taxa, rep, Location, TP, GE,GI,PM_date,year,Family) %>% 
  mutate(year = factor(year, levels = c('2021','2020'))) %>%  pivot_longer(cols = c(GE, GI), names_to = 'trait') %>%
  group_by(TP,PM_date, trait) %>% group_modify(BlueBlupsH2_Year_rep_taxa)  %>% mutate(year = '2020/2021')
#correlations
DH2020Estimates %>% join(DH2021Estimates  %>% dplyr::select(!year)%>% dplyr::rename(value2021 = value))  %>%
  filter(!is.na(value2021), type =='BLUE') %>% filter(type !='H2') %>% group_by(type, TP, trait) %>%
  summarise(correlation = cor(value, value2021))
#Heritabilities over both timepoints, both are high
DH2020Estimates %>% rbind(DH2021Estimates, DHCombined) %>%
  filter(type == 'H2') %>% ggplot(aes(x = TP, y = value, fill = trait)) +geom_bar(stat = 'identity', position = 'dodge')+
  facet_wrap(vars(year), ncol = 1)+theme_bw()+labs(title= 'Broard sense heritability\nover time and datasets')

#other ideas to expand here, perhaps calculate narrow sense heritability, since we have the genotype info
#adding genotype data
AllDHBluesPerYear = rbind(DH2020Estimates, DH2021Estimates,DHCombined) %>% filter(type =='BLUE') %>% ungroup()
#winterGD and WinterGM
load('data/Genotype_data/myGM20_LDprune.Rdata')
load('data/Genotype_data/myGD20_LDprune.Rdata')
WinterGM=myGM20_prune;WinterGD=myGD20_prune
rm(myGM20_prune);rm(myGD20_prune)
WinterGM = WinterGM %>% arrange(Chromosome, Position)
WinterGD = WinterGD %>%
  mutate(taxa = gsub(pattern = '-',replacement = '_',taxa),
         taxa = gsub(pattern = ' ', replacement = '_',taxa),
         taxa1 = taxa) %>% remove_rownames()%>% 
  column_to_rownames('taxa1')
sum(WinterGM$SNP == colnames(WinterGD[,-1]))
# Make sure things are in the right order
# Sum should = 8384
dim(WinterGD)
table(substr(WinterGD$taxa,1,5))
WinterRelationship = rrBLUP::A.mat(WinterGD[,-1]-1, impute.method = 'EM', return.imputed = F)
WinterPCA = eigen(WinterRelationship)
WinterPVEPCA = WinterPCA$values/sum(WinterPCA$values)
data.frame(ordinal = 1:10, PVE = WinterPVEPCA[1:10]) %>%plot(., xlab = 'PC', col = 'red') 
winterlinePCAvalues = WinterPCA$vectors %>% data.frame()%>% 
  mutate(family = mapvalues(substr(WinterGD$taxa,1,3), from = c('BS6','BS7','BS8','BS9','DH1','Fla','SY ','Sca','Win'), 
                            to = c('Flavia/DH130910','Scala/DH130910','SY_Tepee/DH130910','Wintmalt/DH130910',
                                   'DH130910','Flavia/DH130910','SY_Tepee/DH130910','Scala/DH130910','Wintmalt/DH130910')),
         taxa = WinterGD$taxa,
         shapes = ifelse(taxa %in% c('DH130910', 'Flavia','SY_Tepee','Wintmalt','Scala'), taxa, 'Lines'),
         size = ifelse(taxa %in% c('DH130910', 'Flavia','SY_Tepee','Wintmalt','Scala'), 3, 2))

winterlinePCAvalues %>% ggplot(aes(x = X1, y = X2, color = family)) + geom_point()+
  winterlinePCAvalues%>% ggplot(aes(x = X1, y = X3,color = family)) + geom_point()+
  winterlinePCAvalues%>% ggplot(aes(x = X2, y = X3,color = family)) + geom_point()

winterlinePCAvalues %>% filter(family != 'Cha') %>%
  ggplot(aes(x = X1, y = X2, color = family, shape = shapes)) + geom_point(aes(size = size))+theme_bw() +guides(size = "none")+
  xlab('PC1')+ylab('PC2')


setwd(rprojroot::find_rstudio_root_file())


#Travis' Time series perspective
DH2020Estimates %>% rbind(DH2021Estimates, DHCombined) %>%  filter(type == 'BLUE') %>% mutate(headerFacet = 'Data source') %>%
  ggplot(aes(x = PM_date, y = value, group = taxa))+
  geom_line()+facet_nested(trait~headerFacet+year, scales = 'free')+
  geom_vline(xintercept = c(12,33,68), color = 'red')+
  labs(title = 'GE and GI BLUEs over time with differeing datasets')+theme_bw()

DHCombined %>% filter(type == 'BLUE' & trait == 'GI') %>% ggplot(aes(x = PM_date, y = value, group = taxa))+geom_line()+
  geom_vline(xintercept = c(12,33,68), color = 'red')

png('plots/Time_series/Sup_BluesByFamilyQsd1AllYears.png', 2400, 1500, res =120)
AllDHBluesPerYear %>%  filter(type == 'BLUE') %>% mutate(year = factor(year, levels = c('2020','2021','2020/2021')))%>%
  join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1!= 1) %>% mutate(Qsd1= as.factor(Qsd1)) %>%
  filter(Family %nin% c('Cha','End','DH130910')) %>% 
  mutate(Qsd1= ifelse(Qsd1==2,'Dormant','Nondormant')) %>%
  filter(!(trait =='GE' &value>1.05)) %>%
  filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  ggplot(aes(x = TP, y = value, fill = Qsd1))+  
  geom_boxplot()+facet_nested(trait~year+Family, scales = 'free', space = 'free_x')
dev.off()

png('WinterBarley/WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/BluesByFamilyQsd12020_2021.png', 1400, 800, res =120)
AllDHBluesPerYear %>%  filter(type == 'BLUE') %>% mutate(year = factor(year, levels = c('2020','2021','2020/2021')))%>%
  join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1!= 1) %>% mutate(Qsd1= ifelse(Qsd1==2,'Dormant','Nondormant')) %>%
  filter(Family %nin% c('Cha','End','DH130910')) %>% filter(year %in% c('2020/2021'))  %>%
  filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  filter(!(trait =='GE' &value>1.05)) %>%
  ggplot(aes(x = TP, y = value, fill = Qsd1))+  
  geom_boxplot()+facet_nested(trait~year+Family, scales = 'free')
dev.off()

AllDHBluesPerYear %>%  filter(type == 'BLUE') %>% mutate(year = factor(year, levels = c('2020','2021','2020/2021')))%>%
  join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1!= 1) %>% mutate(Qsd1= ifelse(Qsd1==2,'Dormant','Nondormant')) %>%
  filter(Family %nin% c('Cha','End','DH130910')) %>% filter(year %in% c('2021'))  %>%
  filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  ggplot(aes(x = TP, y = value, fill = Qsd1))+  
  geom_boxplot()+facet_nested(trait~year+Family, scales = 'free')
# what is the heritability if Qsd1 is accounted for in the model? 

DHs2020 %>% select(taxa, rep, Location,TP, GE, GI,PM_date,year,Family) %>%
  rbind(., DHs2021 %>% select(taxa, rep, Location, TP, GE,GI,PM_date,year,Family)) %>% 
  mutate(year = factor(year, levels = c('2021','2020'))) %>%  pivot_longer(cols = c(GE, GI), names_to = 'trait') %>% 
  filter(Family %nin%  c('Cha','End','DH130910')) %>%
  group_by(TP,PM_date, trait, Family) %>% join(WinterGD[,c('taxa', 'Qsd1')]) %>% filter(Qsd1 != 1) %>%
  group_modify(~{data.frame(H2 = BLUPH2(lmer(value~Location+rep+ Qsd1+(1|taxa),data = .x)))}) %>% ungroup() %>% select(TP,trait,Family,H2) %>%
  pivot_wider(values_from = H2, names_from = c(TP,trait))

######
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

