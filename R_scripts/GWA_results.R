#GWA analysis code

library(sommer);library(arm);library(lme4);library(Hmisc);library(plyr);library(readxl);
library(tibble);library(patchwork);library(ggplot2);library(fda) ; library(magic); 
library(drc);library(rrBLUP);library(tidyr);library(ggh4x);library(dplyr)
library(tidyr)
load("data/GWA_results/WinterPerTPGWAS.Rdata")

load("data/GWA_results/DH2020.GElogfitGWA_mlmm.Rdata")
load("data/GWA_results/DH2021.GElogfitGWA_mlmm.Rdata")
load("data/GWA_results/DHCombined.GElogfitGWA_mlmm.Rdata")

load("data/GWA_results/DH2020_GIlogfitGWA_mlmm.Rdata")#
load("data/GWA_results/DH2021_GIlogfitGWA_mlmm.Rdata")
load("data/GWA_results/DH20202021_GIlogfitGWA_mlmm.Rdata")

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



#View(DH2020Estimates %>% rbind(DH2021Estimates, DHCombined) %>%  filter(type == 'BLUE') %>%
   #    ungroup() %>% group_by(year, TP, trait))



library(dplyr)
WinterPerTPGWAS %>% arrange(P.value) %>% ungroup()%>% select(SNP, Chromosome, Position) %>% unique()
WinterPerTPGWAS %>% arrange(P.value) 

SigHitsSingle_TP=WinterPerTPGWAS %>% arrange(P.value) %>% filter(P.value<5e-6) %>% filter(maf>0.05) %>% ungroup()%>% select(SNP, Chromosome,Position,maf,P.value,year,TP,trait)%>%
  mutate(most_significant_model=paste0(year,"_",TP,"_",trait))%>%select(SNP,Chromosome,Position,most_significant_model,year,maf,P.value)%>%group_by(SNP)%>% slice_min(n = 1, P.value)


SigHitsSingle_TP_counts<-WinterPerTPGWAS%>%arrange(P.value)%>%filter(P.value<5e-6)%>%filter(maf>0.05) %>% ungroup()%>% select(SNP, Chromosome,Position,maf,P.value,year,TP,trait)%>%
 mutate(model=paste0(year,"_",TP,"_",trait))%>%select(SNP,Chromosome,Position,model,year,maf,P.value)%>%group_by(SNP,Chromosome,Position)%>%dplyr::summarize(count=n(),.groups="keep")%>%ungroup()%>%select(SNP,count)

SigHitsSingle_TP=SigHitsSingle_TP%>%join(SigHitsSingle_TP_counts,by="SNP")%>%arrange(Chromosome,Position)
SigHitsSingle_TP

SigHitsSingle_TP
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

#Logistic
Log2020_GE=DH2020.GElogfitGWA.mlmm%>%arrange(P.value)%>%filter(P.value<5e-6)%>%filter(maf>0.05) %>% ungroup()%>% select(SNP, Chromosome,Position,maf,P.value,year,term)%>%
mutate(most_significant_model=paste0("Logistic GE ",year,"_",term))
Log2021_GE=DH2021.GElogfitGWA.mlmm%>%arrange(P.value)%>%filter(P.value<5e-6)%>%filter(maf>0.05) %>% ungroup()%>% select(SNP, Chromosome,Position,maf,P.value,year,term)%>%
  mutate(most_significant_model=paste0("Logistic GE ",year,"_",term))
LogCombined_GE=DHCombined.GElogfitGWA.mlmm%>%arrange(P.value)%>%filter(P.value<5e-6)%>%filter(maf>0.05) %>% ungroup()%>% select(SNP, Chromosome,Position,maf,P.value,year,term)%>%
  mutate(most_significant_model=paste0("Logistic GE ",year,"_",term))

Log2020_GI=DH2020_GIlogfitGWA.mlmm%>%arrange(P.value)%>%filter(P.value<5e-6)%>%filter(maf>0.05) %>% ungroup()%>% select(SNP, Chromosome,Position,maf,P.value,year,term)%>%
  mutate(most_significant_model=paste0("Logistic GI ",year,"_",term))
Log2021_GI=DH2021_GIlogfitGWA.mlmm%>%arrange(P.value)%>%filter(P.value<5e-6)%>%filter(maf>0.05) %>% ungroup()%>% select(SNP, Chromosome,Position,maf,P.value,year,term)%>%
  mutate(most_significant_model=paste0("Logistic GI ",year,"_",term))
LogCombined_GI=DH20202021_GIlogfitGWA.mlmm%>%arrange(P.value)%>%filter(P.value<5e-6)%>%filter(maf>0.05) %>% ungroup()%>% select(SNP, Chromosome,Position,maf,P.value,year,term)%>%
  mutate(most_significant_model=paste0("Logistic GI ",year,"_",term))



DH2020.GIlogfitGWA.mlmm%>%arrange(P.value)%>%filter(P.value<5e-6)%>%filter(maf>0.05) %>% ungroup()%>% select(SNP, Chromosome,Position,maf,P.value,year,term)%>%
  mutate(most_significant_model=paste0("Logistic ",year,"_",term))
DH2021.GElogfitGWA.mlmm%>%arrange(P.value)%>%filter(P.value<5e-6)%>%filter(maf>0.05) %>% ungroup()%>% select(SNP, Chromosome,Position,maf,P.value,year,term)%>%
  mutate(most_significant_model=paste0("Logistic ",year,"_",term))

Log_hits_all=rbind(Log2020_GE,Log2021_GE,LogCombined_GE,Log2020_GI,Log2021_GI,LogCombined_GI)%>%select(-term)
Log_hits_filterd=Log_hits_all%>%filter(!most_significant_model=="Upper")%>%group_by(SNP)%>% slice_min(n = 1, P.value)

Log_hits_counts=Log_hits_all%>%select(SNP,Chromosome,Position,most_significant_model,maf,P.value)%>%group_by(SNP,Chromosome,Position)%>%dplyr::summarize(count=n(),.groups="keep")%>%ungroup()%>%select(SNP,count)
Log_hits_filterd=Log_hits_filterd%>%join(Log_hits_counts,by="SNP")%>%arrange(Chromosome,Position)



Total_hits=rbind(SigHitsSingle_TP,Log_hits_filterd)%>%arrange(Chromosome,Position)%>%mutate(maf=round(as.numeric(maf),5))
View(Total_hits)
write.csv(Total_hits,file = "data/GWA_results/Winter_SingleTP_andLogisticHits.csv")
