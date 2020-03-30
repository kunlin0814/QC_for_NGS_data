library(ggplot2)
library(readxl)
library(dplyr)

Mut_table <- read_excel("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Mutation_rate/Mutation_rate.xlsx",sheet ='All_samples_Retro_non_retro')
our_data <- read_excel("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Mutation_rate/Sanger_melanoma.xlsx",sheet ='Our_data_new')
excluded <- read.table("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/excluded_samples.txt",
                       sep='\t',header = T)
sanger_data <- read_excel("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Mutation_rate/Sanger_melanoma.xlsx",
                          sheet ="Sanger_data")

mutation_comp <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Mutation_rate\\Sanger_melanoma.xlsx",
                            sheet ="Double_SNP")
# ele <- Mut_table$Sample
# 
# summary <- list()
# for (i in ele){
#   if (i %in% names(summary)){
#     summary[[i]]=summary[[i]]+1
#   } 
#   else{
#     summary[[i]]=1
#   }
# }
# number <-  c()
# for (i in sort(names(summary))){
#   number <- c(number, summary[[i]])
# }

# total_summary <- data.frame("Sample" = sort(names(summary)), "Mutation_nuber" = number)
# new_total_summary <- total_summary[order(-total_summary$Mutation_rate),]
# write.table(total_summary,file ="C:\\Users\\abc73_000\\Desktop\\Sanger_summary.txt" ,quote = F, sep ="\t", row.names = F)

# png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/sager_mutation.png",width=3000,height=2400,res=300)  
# ggplot(data=new_total_summary, aes(x= factor(new_total_summary$Sample,levels = new_total_summary$Sample), y= Mutation_rate))+
#   geom_bar(stat="identity", fill="steelblue")+
#   xlab("Samples")+
#   ylab("Mutation Number")+
#   theme_minimal()
# 
# dev.off()
# png("/Users/kun-linho/Desktop/Mutation_rate/sanger_mutation_rate.png",width=3000,height=2400,res=300)  
# Mut_table <- Mut_table[order(-Mut_table$Mut_rate),]
# ggplot(data=Mut_table, aes(x= factor(Mut_table$Sample,levels = Mut_table$Sample), y= Mut_rate))+
#   geom_bar(stat="identity", fill="steelblue")+
#   xlab("Samples")+
#   ylab("Mutation Number")+
#   theme_minimal()
# 
# dev.off()
# 
# png("/Users/kun-linho/Desktop/Mutation_rate/sanger_mutation_rate.png",width=3000,height=2400,res=300)  
# sanger_sort_data <- sanger_sort_data[order(-sanger_sort_data),]
# ggplot(data=our_sort_data, aes(x= factor(our_sort_data$file_name,levels = our_sort_data$file_name), y= non_retro_PASS))+
#   geom_bar(stat="identity", fill="steelblue")+
#   xlab("Samples")
#   ylab("Mutation Number")+
#   theme_minimal()+
  

dev.off()

both_table <-  cbind(sanger_data, our_data$Total_Mutation)
both_table$Our_mur_rate= (both_table$`our_data$Total_Mutation`/both_table$`Our call`)*1000000
colnames(both_table) <- c("Samples", "Sanger_mut", "our_callable","Sanger_mut_rate","Our_mut_Total","Our_mut_rate")

clean_table <- Mut_table %>% 
  filter(!file_name %in% excluded$SampleName)
clean_table <- na.omit(clean_table)
png("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Mutation_rate/excluded_rate_include_bur_comparison.png",width=3000,height=2400,res=300)  
ggplot(data = clean_table, aes(x=factor(Cancer_Type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Bur_Osteo_Data","Lymphoma","Glioma","Unclassified")), 
       y=Mutation_Rate,
       color= Status))+
  geom_point(position = position_dodge(width =0.5))+
  theme_classic()+
  theme(

axis.title.x = element_blank(),
axis.title.y = element_text(size=16),
axis.text.x = element_text(size=14,hjust = 1, angle = 45),
axis.text.y = element_text(size=14,hjust = 1))

dev.off()


sig_diff <- both_table$Samples[both_table$Sanger_mut-both_table$Our_mut_Total>20] ## prepare to highlight the difference >20
sig_diff_table <- both_table %>% 
  filter(Samples %in%sig_diff)



#both_table <- both_table[order(both_table$Our_mut_rate),]
png("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Mutation_rate/With_CDS-indel-Mutation_rate_comparison.png",width=3000,height=2400,res=300)  
ggplot(data= both_table, aes(x = Our_mut_Total, y=Sanger_mut))+
  geom_point(shape = 1, size =4)+
  geom_abline(intercept = 0, slope = 1, color="blue", 
              linetype="dashed", size=1.5)+
#xlim(0,)+
#ylim(0,250)+
xlab("Our  Data")+
ylab("Sanger Data")+
  geom_point(data=sig_diff_table, aes(x=Our_mut_Total, y=Sanger_mut),color='red',size=4)+
  theme(
    #axis.title=element_text(size=18,face="bold"),
    axis.title.x = element_text(size=20,vjust = -1),
    axis.title.y = element_text(size=20,vjust = 2),
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16))

  

dev.off()


png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Mutation_comparison.png",width=3000,height=2400,res=300)  
ggplot(data= both_table, aes(x = both_table$Our_mut, y=both_table$Sanger_mut), shape=23)+
  geom_point(shape = 1, size =4)+
  xlab("Our  Data")+
  ylab("Sanger Data")+
  theme(
    #axis.title=element_text(size=18,face="bold"),
    axis.title.x = element_text(size=20,vjust = -1),
    axis.title.y = element_text(size=20,vjust = -1),
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16))

dev.off()
