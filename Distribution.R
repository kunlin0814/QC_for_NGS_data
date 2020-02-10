library(ggplot2)
library(readxl)
library(dplyr)



sample <- read.table("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Randomness/Distribution/Mammary/SRR7780805_DepthofCoverage_Distribution.txt",
                     sep ="",header = F, stringsAsFactors = F)

sample <- sample[-2,]
final <- sample[0:600,]
colnames(sample)=c('Frequency',"Position")
png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Randomness/Distribution/Mammary/Plot/600_cut_distribution.png",width=3000,height=2400,res=300)  
freq = as.numeric(sample$Frequency)
pos = as.numeric(sample$Position)
ggplot()+ 
  geom_line(data = sample,aes(x=pos, y=freq),stat="identity")+
  xlim(0,1000)+
  
dev.off()


final[order(as.numeric(final$Position)), ]

ggplot(sample, aes(V1, colour = V2)) +
  geom_freqpoly()+
  dev.off()  


ggtitle("Mean values of the sequencing randomness")+
  ylab("Mean Value")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(color = "black", size = 20, face = "bold"))