library(ggplot2)
library(readxl)
library(dplyr)

Normal <- read_excel("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/V1Pan-Cancer-Supplement_data.xlsx",sheet = 'NormMammary_WES_Mapping_Quality')
  
Tumor <- read_excel("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/V1Pan-Cancer-Supplement_data.xlsx",sheet = 'Tumo_Mammary_WES_Mapping_Qualit')


Meta <- read.table("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Meta/Mammary_Meta.txt",
                   sep = '\t',
                   header = T)

Total <- rbind(Normal, Tumor)

Total %>% 
  filter(uniq_CDS_region_paris_rates >0.65) %>% 
  

Df  <- "/Users/kun-linho/Pan-Cancer-Uniq_exonic_mapped_Distribution.pdf"
pdf(file=Df, w=7, h=5)

plot(ggplot(Normal, aes(x=uniq_CDS_region_paris_rates)) + 
  geom_histogram(bins  = nrow(Normal))+
  xlab("Uniq_exonic_region_mapped_rate")+
  ggtitle("Normal_Samples"))

plot(ggplot(Tumor, aes(x=uniq_CDS_region_paris_rates)) + 
  geom_histogram(bins  = nrow(Tumor))+
  xlab("Uniq_exonic_region_mapped_rate")+
  ggtitle("Tumor_Samples"))

plot(ggplot(Total, aes(x = uniq_CDS_region_paris_rates))+
  geom_histogram(bins  = nrow(Total))+
  xlab("Uniq_exonic_region_mapped_rate")+
  ggtitle("Total_Samples"))

dev.off()
