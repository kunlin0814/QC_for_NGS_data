library(ggplot2)
library(readxl)
library(dplyr)

exclude_sample <- read.table("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/total_exclude-samples.txt",
                             header =  T,
                             stringsAsFactors = F,
                             sep = "\t")
### Without data point ###
table <- read_excel("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Supplement_Figure1/V2Supp1_Data.xlsx",sheet ='%CDS targeting rate')


Total_table  <- table %>% 
  na.omit(table) %>% 
  filter(mean!=0) %>% 
  filter(!Sample_names %in% exclude_sample$Sample_Name)

#### plot the mean value #####
#Cancer <- factor(Total_table$Cancer_type,levels =Total_table$Cancer_type) 

plot_result <- png("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Randomness/CDS-mapping.png",width=3000,height=2400,res=400)  
ggplot(table, aes(x=factor(Cancer_type,levels = c("Mammary_Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")),
                        y=as.numeric(uniq_CDS_region_paris_rates),fill=Status,color=Status)) + 
  geom_point(size=0.01,position = position_jitterdodge(jitter.width = 0.2)) +
  ylab("Mean of the coverage")+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #legend.position = "none",
        legend.title = element_blank(),
        legend.key=element_blank(),
        #legend.text =element_text(color = c("firebrick","black")),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.y = element_text(colour="black",size=18,margin = margin(1,0,0,0)),
        plot.margin = margin(1, 10, 4, 5),
        text = element_text(colour="black",size=14),
        axis.text.x = element_text(colour=c("black"),
                                   size=14,angle=45,vjust=1,hjust = 0.9),
        axis.text.y = element_text(colour="black",size=14),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.1, width = .8,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("firebrick","darkolivegreen"))+
  scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  scale_shape_manual(values = 20)

dev.off()

### plot the mean 1000 ###
plot_result <- png("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Randomness/V7mean.png",width=3000,height=2400,res=400)  
ggplot(Total_table, aes(x=factor(Cancer_type,levels = c("Mammary_Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")),
                        y=Total_table$mean,fill=Status,color=Status)) + 
  geom_point(size=0.01,position = position_jitterdodge(jitter.width = 0.2)) +
  ylab("Mean of the coverage")+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #legend.position = "none",
        legend.title = element_blank(),
        legend.key=element_blank(),
        #legend.text =element_text(color = c("firebrick","black")),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.y = element_text(colour="black",size=18,margin = margin(1,0,0,0)),
        plot.margin = margin(1, 10, 4, 5),
        text = element_text(colour="black",size=14),
        axis.text.x = element_text(colour=c("black"),
                                   size=14,angle=45,vjust=1,hjust = 0.9),
        axis.text.y = element_text(colour="black",size=14),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.1, width = .8,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("firebrick","darkolivegreen"))+
  scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  scale_shape_manual(values = 20)

dev.off()




### plot the rmse ####

plot_result <- png("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Randomness/V7RMSE.png",width=3000,height=2400,res=400)  
ggplot(Total_table, aes(x=factor(Cancer_type,levels = c("Mammary_Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")),
      y=RMSE,fill=Status,color=Status)) + 
  geom_point(size=0.01,position = position_jitterdodge(jitter.width = 0.2)) +
  ylab("RMSE of the coverage")+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #legend.position = "none",
        legend.title = element_blank(),
        legend.key=element_blank(),
        #legend.text =element_text(color = c("firebrick","black")),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.y = element_text(colour="black",size=18,margin = margin(1,0,0,0)),
        plot.margin = margin(1, 10, 4, 5),
        text = element_text(colour="black",size=14),
        axis.text.x = element_text(colour=c("black"),
                                   size=14,angle=45,vjust=1,hjust = 0.9),
        axis.text.y = element_text(colour="black",size=14),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.1, width = .8,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("firebrick","darkolivegreen"))+
  scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  scale_shape_manual(values = 20)

dev.off()

##### plot the rmse 1000 #####
plot_result <- png("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Randomness/V7RMSE.png",width=3000,height=2400,res=400)  
ggplot(Total_table, aes(x=factor(Cancer_type,levels = c("Mammary_Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")),
                        y=RMSE,fill=Status,color=Status)) + 
  geom_point(size=0.01,position = position_jitterdodge(jitter.width = 0.2)) +
  ylab("RMSE of the coverage")+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #legend.position = "none",
        legend.title = element_blank(),
        legend.key=element_blank(),
        #legend.text =element_text(color = c("firebrick","black")),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.y = element_text(colour="black",size=18,margin = margin(1,0,0,0)),
        plot.margin = margin(1, 10, 4, 5),
        text = element_text(colour="black",size=14),
        axis.text.x = element_text(colour=c("black"),
                                   size=14,angle=45,vjust=1,hjust = 0.9),
        axis.text.y = element_text(colour="black",size=14),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.1, width = .8,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("firebrick","darkolivegreen"))+
  scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  scale_shape_manual(values = 20)

dev.off()


###### Violin plot ########


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

plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/V4STD.png",width=3000,height=2400,res=300)  
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
plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/V4rmse.png",width=3000,height=2400,res=300)  
rmse <-  ggplot(Total_table, aes(x=factor(Cancer_type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")), 
                                 y=rmse, fill=Status)) + 
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

plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/V4sum_of_square_error.png",width=3000,height=2400,res=300)  
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

plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/log2-V4rmse_count.png",width=3000,height=2400,res=300)  
rmse_count <-  ggplot(Total_table, aes(x=factor(Cancer_type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")), y=log2(rmse_count), fill=Status)) + 
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

plot_result <- png("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/log2-V4sum_of_square_error_count.png",width=3000,height=2400,res=300)  
sum_of_square_error_count <-  ggplot(Total_table, aes(x=factor(Cancer_type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")), y=log2(sum_of_square_error_count), fill=Status)) + 
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


