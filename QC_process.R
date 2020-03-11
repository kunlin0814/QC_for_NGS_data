library(ggplot2)
library(readxl)
library(dplyr)
# glioma <- read_excel("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Supplement_Figure1/Supp1_Data.xlsx",
#                      sheet= "Sheet1")
# new_glioma <- na.omit(glioma)
# mean(new_glioma$`0...2`[!new_glioma$`0...2`==0])


filter_read <- read_excel("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Supplement_Figure1/V2Supp1_Data.xlsx",
                          sheet = "Total_sequence_pairs")
col_name <- colnames(filter_read)
final_filter_read <- filter_read %>% 
  filter(Total_pairs < 5000000 | Total_pairs==NaN)
rep_time <- nrow(final_filter_read)
reason <- rep("Total_Pair_reads > 5M", rep_time)
final_filter_read <- cbind(final_filter_read,reason)
write.table(final_filter_read,file="/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Exclude_process/seq_read_pairs_lt5M.txt",
            quote = F,sep = "\t",row.names = F)

mapping_quality <- read_excel("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Supplement_Figure1/V2Supp1_Data.xlsx",
                          sheet = "Mapping Quality > 30")
final_filter_mapping <- mapping_quality %>% 
  filter(gt_30_fraction < 0.25 |gt_30_fraction ==NaN)
rep_time <- nrow(final_filter_mapping)
reason <- rep("faction of mapping quality", rep_time)
final_filter_mapping <- cbind(final_filter_mapping,reason)

write.table(final_filter_mapping,file="/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Exclude_process/mapping_qualitylt0.25.txt",
            quote = F,sep = "\t",row.names = F)

CDS_target <-  read_excel("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Supplement_Figure1/V2Supp1_Data.xlsx",
                          sheet = "%CDS targeting rate")
CDS_target %>% 
  filter(Sample_Names=="SAMN03436471")

filter_CDS <- CDS_target %>% 
  filter(as.numeric(uniq_CDS_region_paris_rates) < 0.3|uniq_CDS_region_paris_rates==NaN)
rep_time <- nrow(filter_CDS)
reason <- rep("Unique CDS mapping rate < 0.3", rep_time)
final_filter_CDS <- cbind(filter_CDS,reason)
write.table(final_filter_CDS,file="/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Exclude_process/uniq_CDSlt0.3.txt",
            quote = F,sep = "\t",row.names = F)
  
  # ggplot(aes(x=factor(Cancer_type,levels = c("Mammary_Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")),
  #                         y=as.numeric(uniq_CDS_region_paris_rates),fill=Status,color=Status)) + 
  # geom_point(size=0.01,position = position_jitterdodge(jitter.width = 0.2)) +
  # ylab("Mean of the coverage")+
  # #labs(subtitle = "p<0.01")+
  # #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  # #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  # theme(axis.line = element_line(colour = "black"),
  #       #panel.grid.major = element_blank(),
  #       #panel.grid.minor = element_blank(),
  #       #legend.position = "none",
  #       legend.title = element_blank(),
  #       legend.key=element_blank(),
  #       #legend.text =element_text(color = c("firebrick","black")),
  #       legend.background = element_rect(fill = "transparent"),
  #       panel.border = element_blank(),
  #       axis.title.x = element_blank(),
  #       #axis.ticks.x = element_blank(),
  #       axis.title.y = element_text(colour="black",size=18,margin = margin(1,0,0,0)),
  #       plot.margin = margin(1, 10, 4, 5),
  #       text = element_text(colour="black",size=14),
  #       axis.text.x = element_text(colour=c("black"),
  #                                  size=14,angle=45,vjust=1,hjust = 0.9),
  #       axis.text.y = element_text(colour="black",size=14),
  #       panel.background = element_blank())+
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
  #              geom = "crossbar",size=0.1, width = .8,colour = "black")+
  # #scale_x_discrete(labels=panel2label)+
  # scale_color_manual(values = c("firebrick","darkolivegreen"))+
  # scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  # scale_shape_manual(values = 20)



callable <- read_excel("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Supplement_Figure1/V2Supp1_Data.xlsx",
                       sheet = "Callable_Bases")
callable <- na.omit(callable)

filter_callable <- callable %>%
  filter(Callable_bases<2000000|Callable_bases==NaN) 
rep_time <- nrow(filter_callable)
reason <- rep("callable_bases < 2M", rep_time)  
final_filter_callable <- cbind(filter_callable,reason)
write.table(final_filter_CDS,
            file="/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Exclude_process/callablelt2M.txt",
            quote = F,sep = "\t",row.names = F)

randomness <- read_excel("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Supplement_Figure1/V2Supp1_Data.xlsx",
                       sheet = "Randomness")

filter_randomness <- randomness %>% 
filter(!Mean_of_Coverage<10 |is.nan(Mean_of_Coverage)) 
 
rep_time <- nrow(filter_randomness)
reason <- rep("Mean Coverage < 10", rep_time)  
final_filter_randomness <- cbind(filter_randomness,reason)
write.table(final_filter_randomness,
            file="/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Exclude_process/mean_coveragelt10.txt",
            quote = F,sep = "\t",row.names = F)  

plot_result <- png("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Randomness/V7mean_1000.png",width=3000,height=2400,res=400)  

ggplot(data= filter_randomness,aes(x=factor(Cancer_type,levels = c("Mammary_Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")),
                        y=as.numeric(Mean_of_Coverage),fill=Status,color=Status)) +
geom_point(size=0.01,position = position_jitterdodge(jitter.width = 0.2)) +
ylab("Mean of the coverage to 1000 bases")+
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
plot_result <- png("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Randomness/V7RMSE_1000.png",width=3000,height=2400,res=400)  

ggplot(data= filter_randomness,aes(x=factor(Cancer_type,levels = c("Mammary_Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")),
                                   y=as.numeric(RMSE),fill=Status,color=Status)) +
  geom_point(size=0.01,position = position_jitterdodge(jitter.width = 0.2)) +
  ylab("RMSE of the coverage to 1000 bases")+
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
