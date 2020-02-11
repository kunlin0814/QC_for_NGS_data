library(ggplot2)
library(readxl)
library(dplyr)



sample <- read.table("C:\\Users\\abc73_000\\Desktop\\Pan_cancer_mapping_result\\Randomness\\Distribution\\Mammary\\SRR7780805_DepthofCoverage_Distribution.txt",
                     sep ="",header = F, stringsAsFactors = F)

sample <- sample[-2,]
final <- sample[0:600,]
colnames(sample)=c('Frequency',"Position")
#png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Randomness/Distribution/Mammary/Plot/600_cut_distribution.png",width=3000,height=2400,res=300)  
freq = as.numeric(sample$Frequency)
pos = as.numeric(sample$Position)
ggplot()+ 
  geom_line(data = sample,aes(x=pos, y=freq),stat="identity")+
  ylim(0,200000)
  
dev.off()

DrawDistribution <- function(file_location, xupper, yupper){
  sample <- read.table(file_location, sep =" ",header =F,
                       stringsAsFactors = F)
  sample <- sample[-2, ]
  final <- sample[0:xupper, ]
  colnames(final)=c('Frequency',"Position")
  freq = as.numeric(final$Frequency)
  pos = as.numeric(final$Position)
  p <- ggplot()+ 
    geom_line(data = final,aes(x=pos, y=freq),stat="identity")+
    ylim(0,yupper)+
    ggtitle("Sequence Randomness")+
    ylab("# of Position")+
    xlab("Cover time")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(color = "black", size = 20, face = "bold"))
  return(p)
  
}

DrawDistribution("C:\\Users\\abc73_000\\Desktop\\Pan_cancer_mapping_result\\Randomness\\Distribution\\Mammary\\SRR7780805_DepthofCoverage_Distribution.txt",
                 600,200000)
