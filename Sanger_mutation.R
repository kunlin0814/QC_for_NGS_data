library(ggplot2)
library(readxl)
library(dplyr)

Mut_table <- read_excel("G:\\Pan_cancer\\Pan_cancer_mapping_result\\Mutation_rate\\Sanger_melanoma.xlsx",sheet ='Sanger_data')
our_data <- read_excel("G:\\Pan_cancer\\Pan_cancer_mapping_result\\Mutation_rate\\Sanger_melanoma.xlsx",sheet ='Our_data')
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

png("C:\\Users\\abc73_000\\Desktop\\sager_mutation.png",width=3000,height=2400,res=300)  
ggplot(data=new_total_summary, aes(x= factor(new_total_summary$Sample,levels = new_total_summary$Sample), y= Mutation_rate))+
  geom_bar(stat="identity", fill="steelblue")+
  xlab("Samples")+
  ylab("Mutation Number")+
  theme_minimal()

dev.off()
png("/Users/kun-linho/Desktop/Mutation_rate/sanger_mutation_rate.png",width=3000,height=2400,res=300)  
Mut_table <- Mut_table[order(-Mut_table$Mut_rate),]
ggplot(data=Mut_table, aes(x= factor(Mut_table$Sample,levels = Mut_table$Sample), y= Mut_rate))+
  geom_bar(stat="identity", fill="steelblue")+
  xlab("Samples")+
  ylab("Mutation Number")+
  theme_minimal()

dev.off()

png("/Users/kun-linho/Desktop/Mutation_rate/sanger_mutation_rate.png",width=3000,height=2400,res=300)  
sanger_sort_data <- sanger_sort_data[order(-sanger_sort_data),]
ggplot(data=our_sort_data, aes(x= factor(our_sort_data$file_name,levels = our_sort_data$file_name), y= non_retro_PASS))+
  geom_bar(stat="identity", fill="steelblue")+
  xlab("Samples")+
  ylab("Mutation Number")+
  theme_minimal()+
  

dev.off()

both_table <-  cbind(Mut_table, our_data$non_retro_PASS,our_data$non_retro_mutation_rate  )
colnames(both_table) <- c("Samples", "Sanger_mut", "our_callable","Sanger_mut_rate","Our_mut","Our_mut_rate")
both_table <- both_table[order(both_table$Our_mut_rate),]
png("C:\\Users\\abc73_000\\Desktop\\Mutation_rate_comparison.png",width=3000,height=2400,res=300)  
ggplot(data= both_table, aes(x = both_table$Our_mut_rate, y=both_table$Sanger_mut_rate))+
geom_point()+
xlab("Our  Data")+
ylab("Sanger Data")+
  theme(
    #axis.title=element_text(size=18,face="bold"),
    axis.title.x = element_text(size=20,vjust = -1),
    axis.title.y = element_text(size=20,vjust = 2),
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16))

dev.off()


png("C:\\Users\\abc73_000\\Desktop\\Mutation_comparison.png",width=3000,height=2400,res=300)  
ggplot(data= both_table, aes(x = both_table$Our_mut, y=both_table$Sanger_mut))+
  geom_point()+
  xlab("Our  Data")+
  ylab("Sanger Data")+
  theme(
    #axis.title=element_text(size=18,face="bold"),
    axis.title.x = element_text(size=20,vjust = -1),
    axis.title.y = element_text(size=20,vjust = -1),
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16))

dev.off()
