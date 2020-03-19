library(readxl)
library(tidyverse)
library(wesanderson)
library(RColorBrewer)
# glioma <- read_excel("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Supplement_Figure1/Supp1_Data.xlsx",
#                      sheet= "Sheet1")
# new_glioma <- na.omit(glioma)
# mean(new_glioma$`0...2`[!new_glioma$`0...2`==0])

total_file <- read_excel("/Volumes/Research_Data/MAC_Research_Data/Pan_cancer/Pan_cancer_mapping_result/Supplement_Figure1/V2Supp1_Data.xlsx",
                         sheet ="Total")
PAIR <- read_excel("/Volumes/Research_Data/MAC_Research_Data/Pan_cancer/Pan_cancer_mapping_result/Supplement_Figure1/V2Supp1_Data.xlsx",
                   sheet ="PAIRS")
callable <- read_excel("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Supplement_Figure1/V2Supp1_Data.xlsx",
                       sheet ="Callable_Bases")
###### Sequence read pairs ######

filtered <- total_file %>% 
  filter(Total_pairs < 5000000 | Total_pairs==NaN) %>% 
  select(ID, Total_pairs, Status,Cancer_Type) 


Cases <- c()

for (i in filtered$ID){
  
  if (i %in% PAIR$Normal){
  Case_name <- PAIR$Cases[PAIR$Normal==i]
  }
  else if ( i %in% PAIR$Tumor){
  Case_name <- PAIR$Cases[PAIR$Tumor==i]  
  }
  
  Cases <- c(Cases,Case_name )
  }
rep_time <- nrow(filtered)
reason <- rep("Total_Pair_reads < 5M", rep_time)
Exclude_Reason <- paste(filtered$Status,filtered$ID,reason,sep = "-")
finalTable <- cbind(Cases,filtered,Exclude_Reason)

write.table(finalTable,file="G:\\Pan_cancer\\Pan_cancer_mapping_result\\Exclude_process\\seq_lt_5M.txt",
            quote = F,sep = "\t",row.names = F)  
###### fraction >30  ######
filtered <- total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN)  %>% 
  filter(!gt_30_fraction < 0.25 |gt_30_fraction ==NaN) %>% 
  #filter(as.numeric(uniq_CDS_region_paris_rates) < 0.3|uniq_CDS_region_paris_rates==NaN) %>% 
  select(ID, gt_30_fraction,Status,Cancer_Type)

############# CDS #############
plot_result <- png("/Volumes/Research_Data/MAC_Research_Data/Pan_cancer/Pan_cancer_mapping_result/Exonic_mapping/Unique_CDS_Mapping_Rate.png",width=3000,height=2400,res=400)  
total_file %>% 
filter(!Total_pairs < 5000000 | Total_pairs==NaN)  %>% 
  filter(!gt_30_fraction < 0.25 |gt_30_fraction ==NaN) %>%
  ggplot(aes(x=factor(Cancer_Type,levels = c("Mammary_Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")),
             y=as.numeric(uniq_CDS_region_paris_rates),fill=Status,color=Status)) +
  geom_point(size=0.1,position = position_jitterdodge(jitter.width = 0.2)) +
  ylab("Unique CDS Mapping Rate")+
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
  scale_shape_manual(values = 20)+
  geom_hline(yintercept=.3, linetype="dashed", color = "red", size =.1)
dev.off()

filter(as.numeric(uniq_CDS_region_paris_rates) < 0.3|uniq_CDS_region_paris_rates==NaN) %>% 
select(ID, uniq_CDS_region_paris_rates,Status,Cancer_Type)

Cases <- c()

for (i in filtered$ID){
  
  if (i %in% PAIR$Normal){
    Case_name <- PAIR$Cases[PAIR$Normal==i]
  }
  else if ( i %in% PAIR$Tumor){
    Case_name <- PAIR$Cases[PAIR$Tumor==i]  
  }
  
  Cases <- c(Cases,Case_name )
}

rep_time <- rep_time <- nrow(filtered)
Exclude_Reason <- rep("Unique CDS mapping rate < 0.3", rep_time)
finalTable <- cbind(Cases,filtered,Exclude_Reason)

write.table(finalTable,file="G:\\Pan_cancer\\Pan_cancer_mapping_result\\Exclude_process\\uniq_CDSlt0.3.txt",
            quote = F,sep = "\t",row.names = F)  

##### Randomness #####

filtered <- total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN)  %>% 
  filter(!gt_30_fraction < 0.25 |gt_30_fraction ==NaN) %>% 
  filter(!as.numeric(uniq_CDS_region_paris_rates) < 0.3|uniq_CDS_region_paris_rates==NaN) %>% 
  filter(mean <10) %>% 
  select(ID,mean,Status,Cancer_Type)

Cases <- c()

for (i in filtered$ID){
  
  if (i %in% PAIR$Normal){
    Case_name <- PAIR$Cases[PAIR$Normal==i]
  }
  else if ( i %in% PAIR$Tumor){
    Case_name <- PAIR$Cases[PAIR$Tumor==i]  
  }
  
  Cases <- c(Cases,Case_name )
}

rep_time <- rep_time <- nrow(filtered)
Exclude_Reason <- rep("Mean Coverage < 10", rep_time)
finalTable <- cbind(Cases,filtered,Exclude_Reason)

write.table(finalTable,
            file="G:\\Pan_cancer\\Pan_cancer_mapping_result\\Exclude_process\\mean_coveragelt10.txt",
            quote = F,sep = "\t",row.names = F)  

### find the run-number need to be excluded ####
exclude <- read_excel("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Supplement_Figure1/V2Supp1_Data.xlsx",
                      sheet = "Total_excluded")

# Total_exclude_run <- PAIR %>% 
#   filter(Cases %in% exclude$Cases)
# 
# write.table(Total_exclude_run,
#             file="/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Exclude_process/Exclude_run_sample.txt",
#             quote = F,sep = "\t",row.names = F)  

plot_result <- png("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/callable_bases.png",width=3000,height=2400,res=400)  

### change the color of the point ###
callable %>% 
  filter(!Sample_Name %in% exclude$Cases) %>%
  ggplot(
  aes(x=factor(Cancer_Type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")),
       y=as.numeric( Callable_bases)))+
  geom_point(size=3,aes(color=Cancer_Type),position = position_jitterdodge(jitter.width = 0.1))+
  ylab("Callable bases")+
  ylim(0, max(callable$Callable_bases))+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        legend.position = "none",
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
  scale_color_brewer(palette="Dark2") ### for scatter plot
  #scale_x_discrete(labels=panel2label)+
  #scale_color_manual(values=wes_palette(n=5, name="Set3")) ## for scatter plot Wes Anderson color palettes
  #scale_fill_manual(values=wes_palette(n=5, name="Royal2")) ## for box plot, bar plot Wes Anderson color palettes
  #scale_shape_manual(values = 20)
dev.off()


plot_result <- png("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/RMSE_coverage.png",width=3000,height=2400,res=400)  

total_file %>% 
  filter(!ID %in% exclude$Normal) %>% 
  filter(!ID %in%exclude$Tumor) %>% 
  ggplot(aes(x=factor(Cancer_Type,levels = c("Mammary_Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")),
  y=as.numeric(RMSE),fill=Status,color=Status)) +
  geom_point(size=0.01,position = position_jitterdodge(jitter.width = 0.2)) +
  ylab("Coverage RMSE")+
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
#### callable bases #####

callable <- read_excel("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Supplement_Figure1/V2Supp1_Data.xlsx",
                       sheet = "Callable_Bases")
callable <- na.omit(callable)

filter_callable <- callable %>%
  filter(!Callable_bases<2000000|Callable_bases==NaN) 
rep_time <- nrow(filter_callable)
reason <- rep("callable_bases < 2M", rep_time)  
final_filter_callable <- cbind(filter_callable,reason)
write.table(final_filter_CDS,
            file="/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Exclude_process/callablelt2M.txt",
            quote = F,sep = "\t",row.names = F)


##### Randomness plot #####

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




