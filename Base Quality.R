library(tidyverse)
library(readxl)

setwd("C:\\Users\\abc73_000\\Desktop")

total_file <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\V2Supp1_Data.xlsx",
                         sheet ="Total")

PAIR <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\V2Supp1_Data.xlsx",
                   sheet ="PAIRS")  


#Glioma_base <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Base_Quality\\Base_quality.xlsx",
#                          sheet ='After', col_names =F)

Glioma_base <- read.table("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Base_Quality\\Glioma_BaseQuality_total.txt",
                          header = F, stringsAsFactors = F)

Before_Glioma_bases <- read.table("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Base_Quality\\Glioma_data\\Before\\Before_Glioma_total_BasesQuality.txt"
                                  ,header = F, stringsAsFactors =F)
#### extract each id to see if they pass 0.8 #####

gt_0.8 <- total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>%
  filter(Cancer_Type =='Glioma') %>% 
  filter(as.numeric(Uniquely_mapped_rate) >=0.8) %>% 
  select(ID)
#write.table("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Base_Quality\\Glioma_data\\gt0.8_glioma.txt",
#            sep ='\t',row.names = F,quote = F)

lt_0.8 <- total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>%
  filter(Cancer_Type =='Glioma') %>% 
  filter(as.numeric(Uniquely_mapped_rate) <0.8) %>%
  select(ID)
#write.table("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Base_Quality\\Glioma_data\\lt0.8_glioma.txt",
#            sep ='\t',row.names = F,quote = F)


#### before gt lt 0.8 ####

gt0.8_sample <- Before_Glioma_bases %>% 
  filter(V1 %in% gt_0.8$ID) 


lt_0.8_sample <- Before_Glioma_bases %>% 
  filter(V1 %in% lt_0.8$ID) 


gt0.8Normal <- gt0.8_sample %>% 
  filter(V2=='Normal') %>% 
  select(!V1) %>% 
  select(!V2)


gt0.8Tumor <- gt0.8_sample %>% 
  filter(V2=='Tumor') %>% 
  select(!V1) %>% 
  select(!V2)

  
lt0.8Normal <- lt_0.8_sample %>% 
  filter(V2 =='Normal') %>% 
  select(!V1) %>% 
  select(!V2)

lt0.8Tumor <- lt_0.8_sample %>% 
  filter(V2 =='Tumor') %>% 
  select(!V1) %>% 
  select(!V2)


# ### After   ####
# ##### gt 0.8 #####
# 
# gt0.8Normal <- Glioma_base %>% 
#   filter(V1 %in% gt_0.8$ID) %>% 
#   select(!V1) %>%
#   filter(V2=='Normal') %>% 
#   select(!V2) #%>% 
#   #write.table("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Base_Quality\\Glioma_data\\V2gt0.8glioma_normal.txt",
#   #            sep ='\t',row.names = F,col.names = F,quote = F)
# 
# 
# gt0.8Tumor <- Glioma_base %>% 
#   filter(V1 %in% gt_0.8$ID) %>% 
#   select(!V1) %>%
#   filter(V2=='Tumor') %>% 
#   select(!V2) 
#   #%>% 
#   #write.table("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Base_Quality\\Glioma_data\\V2gt0.8glioma_tumor.txt",
#   #            sep ='\t',row.names = F,col.names = F,quote = F)
# 
# ######## lt 0.8 ########
# 
# lt0.8Normal <- Glioma_base %>% 
#   filter(V1 %in% lt_0.8$ID) %>% 
#   select(!V1) %>%
#   filter(V2=='Normal') %>% 
#   select(!V2) 
#   # %>% 
#   # write.table("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Base_Quality\\Glioma_data\\V2lt0.8glioma_normal.txt",
#   #             sep ='\t',row.names = F,col.names = F,quote = F)
# 
# lt0.8Tumor <- Glioma_base %>% 
#   filter(V1 %in% lt_0.8$ID) %>% 
#   select(!V1) %>%
#   filter(V2=='Tumor') %>% 
#   select(!V2) 
# # %>% 
# #   write.table("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Base_Quality\\Glioma_data\\V2lt0.8glioma_tumor.txt",
# #               sep ='\t',row.names = F,col.names = F,quote = F)  

PAIR %>% 
  filter(Cancer_type=='Glioma') %>% 
  filter(Normal %in% lt_0.8$ID | Tumor %in% lt_0.8$ID) %>% 
  select(Cases) %>% 
  write.table("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Base_Quality\\Glioma_data\\lt0.8Cases.txt",
              sep ='\t',row.names = F,col.names = F,quote = F)  


PAIR %>% 
  filter(Cancer_type=='Glioma') %>% 
  filter(! (Normal %in% lt_0.8$ID | Tumor %in% lt_0.8$ID)) %>% 
  select(Cases) %>% 
  write.table("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Base_Quality\\Glioma_data\\gt0.8Cases.txt",
              sep ='\t',row.names = F,col.names = F,quote = F) 


library(ggplot2)


# 
# lym_norm = read.table('C:/Users/abc73_000/Desktop/Base_Quality/Normal_Glioma_BQ15.Baes',header = F)
# lym_tumor = read.table('C:/Users/abc73_000/Desktop/Base_Quality/Tumor_Glioma_BQ15.Baes',header = F)

#gt0.8Normal, gt0.8Tumor, lt0.8Normal, lt0.8Tumor

lym_norm = lt0.8Normal
lym_tumor = lt0.8Tumor


reformat7 = function(lym){
  result = NULL
  for (i in 1:ncol(lym)){
    pos = i
    for (j in 1:nrow(lym)){
      score = as.numeric(as.character(lym[j,i]))
      tmp = c(score,pos)
      result = rbind(result,tmp)
    }
  }
  result = data.frame(result)
  colnames(result)=c('Quality','Position')
  return (result)
}

lym_table = reformat7(lym_norm)

lym_table$Quality = as.numeric(as.character(lym_table$Quality))
lym_table$Position = as.character(lym_table$Position)
lym_table$Position = factor(lym_table$Position,levels = unique(lym_table$Position))
tiff(file = "Before_lt0.8Glioma_WES_baseQuality_norm.tiff", width = 1700, height =1000, units = "px", res = 300)
ggplot(lym_table, aes(x=Position, y=Quality)) +
  ggtitle("Normal")+
  geom_boxplot(outlier.shape = NA) + theme_classic()+ #+ geom_point(position = position_jitterdodge(),size = 1)
  coord_cartesian(ylim=c(22.5,32.5))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=18))
  dev.off()

lym_table = reformat7(lym_tumor)

lym_table$Quality = as.numeric(as.character(lym_table$Quality))
lym_table$Position = as.character(lym_table$Position)
lym_table$Position = factor(lym_table$Position,levels = unique(lym_table$Position))
tiff(file = "Before_lt0.8_Glioma_WES_baseQuality_tumor.tiff", width = 1700, height =1000, units = "px", res = 300)

ggplot(lym_table, aes(x=Position, y=Quality)) +
  ggtitle("Tumor")+
  geom_boxplot(outlier.shape = NA) + theme_classic()+ #+ geom_point(position = position_jitterdodge(),size = 1)
  coord_cartesian(ylim=c(22.5,32.5))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=18))
  dev.off()

  
  
