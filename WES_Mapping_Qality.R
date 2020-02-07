library(ggplot2)
library(readxl)
library(dplyr)


list <- c('Normal_Mammary_WES_Mapping','Tumor_Mammary_WES_Mapping','Normal_Melanoma_Mapping','Tumor_Melanoma_Mapping',
          'Normal_Osteo_Mapping','Tumor_Osteo_Mapping','Normal_Lym_Mapping','Tumor_Lym_Mapping')
Df  <- "/Users/kun-linho/Desktop/Pan-Cancer-Uniq_exonic_mapped_Distribution.pdf"
pdf(file=Df, w=7, h=5)

for (i in list){
Data <- read_excel("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/V1Pan-Cancer-Supplement_data.xlsx",sheet = i)
Normal <- na.omit(Data)
#Tumor <- read_excel("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/V1Pan-Cancer-Supplement_data.xlsx",sheet = i)


# Meta <- read.table("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Meta/Mammary_Meta.txt",
#                    sep = '\t',
#                    header = T)

#Total <- rbind(Normal, Tumor)

num_sample <-  length(Normal$uniq_CDS_region_paris_rates)
  


plot(ggplot(Normal, aes(x=uniq_CDS_region_paris_rates)) + 
  geom_histogram(bins  = num_sample)+
  xlab("Uniq_exonic_region_mapped_rate")+
  ggtitle(i))
# 
# plot(ggplot(Tumor, aes(x=uniq_CDS_region_paris_rates)) + 
#   geom_histogram(bins  = nrow(Tumor))+
#   xlab("Uniq_exonic_region_mapped_rate")+
#   ggtitle(i))

# plot(ggplot(Total, aes(x = uniq_CDS_region_paris_rates))+
#   geom_histogram(bins  = nrow(Total))+
#   xlab("Uniq_exonic_region_mapped_rate")+
#   ggtitle(i))


}
dev.off()
