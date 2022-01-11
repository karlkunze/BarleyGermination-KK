#graphing outputs from WinterTPGWAS
load("data/WinterPerTPGWAS.Rdata")
chrTable = c(785, 1534, 1187,  760, 1348,  823, 1334,16 )
chrLabel = c(1:7, 'UN')
winterOrdinalBreaks = c(chrTable[1]/2)
WinterChrLines = c(785)
for (i in 2:8){
  winterOrdinalBreaks[i] = sum(chrTable[1:i-1])+chrTable[i]/2
  WinterChrLines[i] = sum(chrTable[1:i])
}

detach(asreml)

WinterPerTPGWAS %>% ggplot(aes(ordinal, log10PVal, color = TP, shape = year))+geom_point()+
  geom_vline(xintercept = WinterChrLines, color = 'black')+
  geom_vline(xintercept = 4780, color = 'red',linetype='longdash')+
  annotate(geom= 'text', x = 4780, y = 30, label = 'AlaAT1')+
  geom_vline(xintercept = WinterChrLines)+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = winterOrdinalBreaks)+
  ylab('-log(p-value)')+xlab('Chromosome')+ geom_hline(yintercept = -log10(5e-5)) +
  facet_grid(rows = vars(trait), scales = 'free_y')+theme_bw()
