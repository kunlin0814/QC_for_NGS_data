library(ggplot2)
library(readxl)
library(dplyr)

table <- read_excel("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Callable_bases/All_Cancer_callable.xlsx",
                    sheet = 'Sheet1')

png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Callable_Bases.png",width=3000,height=2400,res=300)  
ggplot(table, aes(x= factor(Cancer_type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Bur_Osteo","Lymphoma","Glioma","Unclassified")), 
                               y =  Callable_bases ))+
  geom_boxplot()+
  ylab("Callable bases")+
  #scale_fill_manual(values= c("gray","red"))+
  theme(
    #axis.title=element_text(size=18,face="bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=16),
    axis.text.x = element_text(size=16,hjust = 1, angle = 45),
    axis.text.y = element_text(size=16,hjust = 1))
dev.off()