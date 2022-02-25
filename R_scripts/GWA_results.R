#GWA analysis code

library(sommer);library(arm);library(lme4);library(Hmisc);library(plyr);library(readxl);
library(tibble);library(patchwork);library(ggplot2);library(fda) ; library(magic); 
library(drc);library(rrBLUP);library(tidyr);library(ggh4x);library(dplyr)

load("data/GWA_results/WinterPerTPGWAS.Rdata")
load("data/GWA_results/DH2020.GElogfitGWA_mlmm.Rdata")

load("data/GWA_results/DH20202021_GIlogfitGWA_mlmm.Rdata")



#View(DH2020Estimates %>% rbind(DH2021Estimates, DHCombined) %>%  filter(type == 'BLUE') %>%
   #    ungroup() %>% group_by(year, TP, trait))
table(DHCombined$TP,DHCombined$PM_date)
WinterPerTPGWAS


WinterPerTPGWAS %>% arrange(P.value) %>% ungroup()%>% select(SNP, Chromosome, Position) %>% unique()
WinterPerTPGWAS %>% arrange(P.value) %>% view()
#(P.value<5e-6)
test$PM_date
table(WinterPerTPGWAS$TP)
WinterPerTPGWAS%>%mutate(PM_date=mapvalues(TP,from=c("TP1","TP1.5","TP2","TP2.5","TP3","TP3.5","TP4","TP5"),to=c(5,12,19,33,47,68,96,152)))%>%
  mutate(Time_Point=paste0(TP,"(",PM_date,")"))%>%
  
  ggplot(aes(ordinal, log10PVal, color = Time_Point, shape = year))+geom_point()+
  geom_vline(xintercept = WinterChrLines, color = 'black')+
  geom_vline(xintercept = 4780, color = 'red')+
  annotate(geom= 'text', x = 4780, y = 30, label = 'AlaAT1')+
  geom_vline(xintercept = WinterChrLines)+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = winterOrdinalBreaks)+
  ggtitle("Single Time Point GWA")+
  ylab('-log(p-value)')+xlab('Chromosome')+ geom_hline(yintercept = -log10(5e-5)) +
  facet_grid(rows = vars(trait), scales = 'free_y')+theme_bw()
#Logistics

DH2020.GElogfitGWA.mlmm %>%ggplot(aes(ordinal, log10PVal, color = term))+geom_point()+
  geom_vline(xintercept = WinterChrLines, color = 'black')+
  geom_vline(xintercept = 4780, color = 'red')+
  annotate(geom= 'text', x = 4780, y = 30, label = 'AlaAT1')+
  geom_vline(xintercept = WinterChrLines)+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = winterOrdinalBreaks)+
  ylab('-log(p-value)')+xlab('Chromosome')+ geom_hline(yintercept = -log10(5e-5)) +
  facet_grid(scales = 'free_y')+theme_bw()



