library(ggplot2)
library(readxl)
library(dplyr)


### Without data point ###
Total_table <- read_excel("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Randomness/V2_Randomness.xlsx",sheet ='Total_random')
Total_table <- na.omit(Total_table)  
#Cancer <- factor(Total_table$Cancer_type,levels =Total_table$Cancer_type) 

plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/V2mean.png",width=3000,height=2400,res=400)  
mean <- ggplot(Total_table, aes(x=factor(Cancer_type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")), y=MEAN, fill=Status)) + 
  geom_violin(trim=FALSE)+
  theme_classic()+
  ylab("Coverage Mean")+
  scale_fill_manual(values= c("gray","red"))+
  theme(
    #axis.title=element_text(size=18,face="bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=16,hjust = 1, angle = 45),
    axis.text.y = element_text(size=16,hjust = 1))
#mean
mean
#position_dodge(width = 0.9))
dev.off()

plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/V2STD.png",width=3000,height=2400,res=300)  
std <-  ggplot(Total_table, aes(x=factor(Cancer_type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")), y=STD, fill=Status)) + 
  geom_violin(trim=FALSE)+
  theme_classic()+
  scale_fill_manual(values=c("gray","red"))+
  ylab("Coverage Standard Deviation")+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=16,hjust = 1, angle = 45),
    axis.text.y = element_text(size=16,hjust = 1))
std
dev.off()
#stat_summary(fun.data=mean_sdl,geom="pointrange", color="black",position = position_dodge(width = 0.9))
plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/V2rmse.png",width=3000,height=2400,res=300)  
rmse <-  ggplot(Total_table, aes(x=factor(Cancer_type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")), y=rmse, fill=Status)) + 
  geom_violin(trim=FALSE)+
  theme_classic()+
  scale_fill_manual(values=c("gray","red"))+
  ylab("Coverage RMSE")+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=16,hjust = 1, angle = 45),
    axis.text.y = element_text(size=16,hjust = 1))
rmse
dev.off()

plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/V2sum_of_square_error.png",width=3000,height=2400,res=300)  
sum_of_square_error <-  ggplot(Total_table, aes(x=factor(Cancer_type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")), y=sum_of_square_error, fill=Status)) + 
  geom_violin(trim=FALSE)+
  theme_classic()+
  scale_fill_manual(values=c("gray","red"))+
  ylab("Coverage Sum of Square Error")+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=16,hjust = 1, angle = 45),
    axis.text.y = element_text(size=16,hjust = 1))
sum_of_square_error
dev.off()

plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/V2rmse_count.png",width=3000,height=2400,res=300)  
rmse_count <-  ggplot(Total_table, aes(x=factor(Cancer_type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")), y=rmse_count, fill=Status)) + 
  geom_violin(trim=FALSE)+
  theme_classic()+
  scale_fill_manual(values=c("gray","red"))+
  ylab("Coverage RMSE Counts")+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=16,hjust = 1, angle = 45),
    axis.text.y = element_text(size=16,hjust = 1))
rmse_count
dev.off()

plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/V2sum_of_square_error_count.png",width=3000,height=2400,res=300)  
sum_of_square_error_count <-  ggplot(Total_table, aes(x=factor(Cancer_type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")), y=sum_of_square_error_count, fill=Status)) + 
  geom_violin(trim=FALSE)+
  theme_classic()+
  scale_fill_manual(values=c("gray","red"))+
  ylab("Coverage Sum of Square error Counts")+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=16,hjust = 1, angle = 45),
    axis.text.y = element_text(size=16,hjust = 1))
sum_of_square_error_count
dev.off()




##### Unique CDS and Unique mapping rates #### 
Total_sample <- read_excel("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/V1Pan-Cancer-Supplement_data.xlsx",
                           sheet = "Total_WES_Mapping")
Total_sample <- na.omit(Total_sample)  
plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/v1Uniq_mapping_rate.png",width=3000,height=2400,res=300)
uniq <-  ggplot(Total_sample, aes(x=factor(Cancer_type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma")), y=uniq_mapped_rate*100, fill=Status)) + 
  geom_violin(trim=FALSE)+
  theme_classic()+
  scale_fill_manual(values=c("gray","red"))+
  ylab("Unique mapping Rate %") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=16,hjust = 1, angle = 45),
        axis.text.y = element_text(size=16,hjust = 1))
uniq

dev.off()
plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/v1Uniq_CDS_mapping_rate.png",width=3000,height=2400,res=300)
uniqCDS <-  ggplot(Total_sample, aes(x=factor(Cancer_type,levels =c("Mammary Cancer","Melanoma", "Osteosarcoma",)), y=uniq_CDS_region_paris_rates*100, fill=Status)) + 
  geom_violin(trim=FALSE)+
  theme_classic()+
  scale_fill_manual(values=c("gray","red"))+
  ylab("Unique CDS mapping Value %")+
  ylim(0,100)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=16,hjust = 1, angle = 45),
        axis.text.y = element_text(size=16,hjust = 1))
uniqCDS

dev.off()



########################
## with data point ##
Total_table <- read_excel("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Randomness/Randomness.xlsx",sheet ='Total_random')
Total_table <- na.omit(Total_table)  
Status <- factor(Total_table,levels =as.vector(Total_table$Status)) 
plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/mean.png",width=3000,height=2400,res=300)  
mean <- ggplot(Total_table, aes(x=factor(Cancer_type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassify")), y=MEAN, fill=Status)) + 
    geom_violin(trim=FALSE)+
    theme_classic()+
    ggtitle("Mean values of the sequencing randomness")+
    ylab("Mean Value")+
    xlab("Cancer Type")+
    scale_fill_manual(values= c("gray","red"))+
  theme(
    axis.title=element_text(size=18,face="bold"),
    axis.text.x = element_text(size=16,hjust = 1, angle = 30),
    axis.text.y = element_text(size=16,hjust = 1),
    plot.title = element_text(color = "black", size = 20, face = "bold"))
#mean
mean+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, position = "dodge",colour = "black",
                  binwidth = 3,binpositions="all")
                    #position_dodge(width = 0.9))
dev.off()
plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/STD.png",width=3000,height=2400,res=300)  
std <-  ggplot(Total_table, aes(x=factor(Cancer_type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassify")), y=STD, fill=Status)) + 
  geom_violin(trim=FALSE)+
  theme_classic()+
  scale_fill_manual(values=c("gray","red"))+
  ggtitle("Standard deviation value of the sequencing randomness")+
  ylab("Standard Deviation Value")+
  xlab("Cancer Type")+
  theme(
    axis.title=element_text(size=18,face="bold"),
    axis.text.x = element_text(size=16,hjust = 1, angle = 30),
    axis.text.y = element_text(size=16,hjust = 1),
    plot.title = element_text(color = "black", size = 20, face = "bold"))
std+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, position = "dodge",colour = "black",
                     binwidth = 6,binpositions="all")
dev.off()
#stat_summary(fun.data=mean_sdl,geom="pointrange", color="black",position = position_dodge(width = 0.9))
plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Skewness.png",width=3000,height=2400,res=300)  

skewness <-  ggplot(Total_table, aes(x=factor(Cancer_type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassify")), y=Skewness, fill=Status)) + 
  geom_violin(trim=FALSE)+
  theme_classic()+
  scale_fill_manual(values=c("gray","red"))+
  ggtitle("Skewness value of the sequencing randomness")+
  ylab("Skewness Value")+
  xlab("Cancer Type")+
  theme(
  axis.title=element_text(size=18,face="bold"),
  axis.text.x = element_text(size=16,hjust = 1, angle = 30),
  axis.text.y = element_text(size=16,hjust = 1),
  plot.title = element_text(color = "black", size = 20, face = "bold"))
skewness+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, position = "dodge",colour = "black",
                          binwidth = 0.9,binpositions="all")
dev.off()


Total_sample <- read_excel("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/V1Pan-Cancer-Supplement_data.xlsx",
                           sheet = "Total_WES_Mapping")
Total_sample <- na.omit(Total_sample)  
plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Uniq_mapping_rate.png",width=3000,height=2400,res=300)
uniq <-  ggplot(Total_sample, aes(x=factor(Cancer_type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassify")), y=uniq_mapped_rate*100, fill=Status)) + 
  geom_violin(trim=FALSE)+
  theme_classic()+
  scale_fill_manual(values=c("gray","red"))+
  ggtitle("Total Sample Unique Mapping Rate")+
  ylab("unique mapping Rate %")+
  xlab("Cancer Type")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=16,face="bold"),
        axis.text.x = element_text(size=14,hjust = 1, angle = 30),
        axis.text.y = element_text(size=14,hjust = 1),
        plot.title = element_text(color = "black", size = 20, face = "bold"))
uniq+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3, position = "dodge",colour = "black",
                      binwidth = 0.3,binpositions="all")

dev.off()
plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Uniq_CDS_mapping_rate.png",width=3000,height=2400,res=300)
uniqCDS <-  ggplot(Total_sample, aes(x=factor(Cancer_type,levels =c(Total_table$Cancer_type)), y=uniq_CDS_region_paris_rates*100, fill=Status)) + 
  geom_violin(trim=FALSE)+
  theme_classic()+
  scale_fill_manual(values=c("gray","red"))+
  ggtitle("Total Sample Unique CDS Mapping Rate")+
  ylab("unique CDS mapping Value %")+
  xlab("Cancer Type")+
  ylim(0,100)+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=16,face="bold"),
        axis.text.x = element_text(size=14,hjust = 1, angle = 30),
        axis.text.y = element_text(size=14,hjust = 1),
        plot.title = element_text(color = "black", size = 20, face = "bold"))
uniqCDS+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9, position = "dodge",colour = "black",
                  binwidth = 0.9,binpositions="all")


dev.off()


