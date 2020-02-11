library(ggplot2)
library(readxl)
library(dplyr)


Total_table <- read_excel("C:\\Users\\abc73_000\\Desktop\\Pan_cancer_mapping_result\\Randomness\\Randomness.xlsx",sheet ='Total_random')
Total_table <- na.omit(Total_table)  
 
plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Randomness/mean.png",width=3000,height=2400,res=300)  
mean <- ggplot(Total_table, aes(x=Cancer_type, y=MEAN, fill=Status)) + 
    geom_violin(trim=FALSE)+
    theme_classic()+
    ggtitle("Mean values of the sequencing randomness")+
    ylab("Mean Value")+
    scale_fill_manual(values= c("gray","red"))+
    theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(color = "black", size = 20, face = "bold"))
#mean
mean+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, position = "dodge",colour = "black",
                  binwidth = 2.5,binpositions="all")
                    #position_dodge(width = 0.9))
dev.off()
plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Randomness/STD.png",width=3000,height=2400,res=300)  
std <-  ggplot(Total_table, aes(x=Cancer_type, y=STD, fill=Status)) + 
  geom_violin(trim=FALSE)+
  theme_classic()+
  scale_fill_manual(values=c("gray","red"))+
  ggtitle("Standard deviation value of the sequencing randomness")+
  ylab("Standard Deviation Value")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(color = "black", size = 20, face = "bold"))
std+geom_dotplot(method="histodot",binaxis='y', stackdir='center', dotsize=0.8, position = "dodge",colour = "black",
                     binwidth = 2.5,binpositions="all")
dev.off()
#stat_summary(fun.data=mean_sdl,geom="pointrange", color="black",position = position_dodge(width = 0.9))
plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Randomness/Skewness.png",width=3000,height=2400,res=300)  

skewness <-  ggplot(Total_table, aes(x=Cancer_type, y=Skewness, fill=Status)) + 
  geom_violin(trim=FALSE)+
  theme_classic()+
  scale_fill_manual(values=c("gray","red"))+
  ggtitle("Skewness value of the sequencing_randomness")+
  ylab("Skewness Value")+
  theme(axis.text=element_text(size=10),
  axis.title=element_text(size=14,face="bold"),
  plot.title = element_text(color = "black", size = 20, face = "bold"))
skewness+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, position = "dodge",colour = "black",
                          binwidth = 0.7,binpositions="all")
dev.off()



  # 
  # plot(ggplot(Tumor, aes(x=uniq_CDS_region_paris_rates)) + 
  #   geom_histogram(bins  = nrow(Tumor))+
  #   xlab("Uniq_exonic_region_mapped_rate")+
  #   ggtitle(i))
  
  # plot(ggplot(Total, aes(x = uniq_CDS_region_paris_rates))+
  #   geom_histogram(bins  = nrow(Total))+
  #   xlab("Uniq_exonic_region_mapped_rate")+
  #   ggtitle(i))
  
  






