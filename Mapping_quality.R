library(ggplot2)
library(readxl)
library(dplyr)

total_mapping <- read_excel("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/mapping_quality.xlsx",
                            sheet='Sheet1')

## remove the data frame with na or 0 in each row
total_mapping <- na.omit(total_mapping)
#row_sub = t(apply(total_mapping, 1, function(row) all(row !=0 )))

fina_total_mapping <- total_mapping %>% 
  filter(gt_30_fraction!=0)
  
###
plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/gt30fractions.png",width=3000,height=2400,res=400)  

ggplot(fina_total_mapping, aes(x= factor(Cancer_type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")), 
                          y = gt_30_fraction,fill= Status ))+
  geom_boxplot()+theme_classic()+
  ylab("Fractions of Quality Scores > 30")+
  scale_fill_manual(values= c("gray","red"))+
  theme(
    #axis.title=element_text(size=18,face="bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=16),
    axis.text.x = element_text(size=16,hjust = 1, angle = 45),
    axis.text.y = element_text(size=16,hjust = 1))
dev.off()
plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/gt60fractions.png",width=3000,height=2400,res=400)  

ggplot(fina_total_mapping, aes(x= factor(Cancer_type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")), y = gt_60_fraction,fill= Status ))+
  geom_boxplot()+theme_classic()+
  ylab("Fractions of Quality Scores > 60")+
  scale_fill_manual(values= c("gray","red"))+
  theme(
    #axis.title=element_text(size=18,face="bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=16),
    axis.text.x = element_text(size=16,hjust = 1, angle = 45),
    axis.text.y = element_text(size=16,hjust = 1))
dev.off()
