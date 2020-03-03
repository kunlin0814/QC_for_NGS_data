library(ggplot2)
library(readxl)
library(dplyr)

table <- read_excel("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Callable_bases/All_Cancer_callable.xlsx",
                    sheet = 'Sheet1')
exclude <- table %>% 
  filter(Callable_bases<10000000)
exclude_sample <- read.table("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/total_exclude-samples.txt",
                             header =  T,
                             stringsAsFactors = F,
                             sep = "\t")
finalTable <- table %>% 
  filter(!Sample_Name %in% exclude_sample$Sample_Name) %>% 
  filter(!Cancer_type =="Bur_Osteo")

# write.table(exclude, file= "/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/callable_base_less_than10M.txt", quote = F,
#             sep = "\t", row.names = F)
png("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Exclude-Callable_Bases.png",width=3000,height=2400,res=300)  
ggplot(finalTable, aes(x= factor(Cancer_type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Lymphoma","Glioma","Unclassified")), 
                               y =  Callable_bases ))+
  ylim(0, max(finalTable$Callable_bases))+
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
