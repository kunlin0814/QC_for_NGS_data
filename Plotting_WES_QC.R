library(ggplot2)
library(readxl)
library(dplyr)
library(wesanderson)
library(RColorBrewer)

total_file <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\V2Supp1_Data.xlsx",
                         sheet ="Total")
PAIR <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\V2Supp1_Data.xlsx",
                   sheet ="PAIRS")
callable <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\V2Supp1_Data.xlsx",
                    sheet ="Callable_Bases")
exclude <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\V2Supp1_Data.xlsx",
                      sheet ="Callable_Bases")


plot_result <- png("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\sequence_pairs.png",
                   width=2800,height=1800,res=450)  
total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>% 
  ggplot(aes(x=factor(Cancer_Type,levels = c("Mammary_Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")),
             y=as.numeric(Total_pairs),fill=Status,color=Status)) +
  geom_point(size=0.01,position = position_jitterdodge(jitter.width = 0.2)) +
  ylab("Sequence Read Pairs")+
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
        axis.title.y = element_text(colour="black",size=26,margin = margin(1,5,.5,0)),
        text = element_text(colour="black",size=14),
        axis.text.x = element_blank(), 
        #axis.text.x = element_text(colour=c("black"),size=12,vjust=1,hjust = 0.9),
        axis.text.y = element_text(colour="black",size=18),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.1, width = .8,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("firebrick","darkolivegreen"))+
  scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  scale_shape_manual(values = 20)+
  coord_cartesian(ylim=c(0,2*10**8))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))
 
dev.off()

### Unique Mapping Rates ###
plot_result <- png("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\unique-mapping_pairs.png",
                   width=2800,height=1800,res=450)
total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>% 
  ggplot(aes(x=factor(Cancer_Type,levels = c("Mammary_Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")),
             y=as.numeric(Uniquely_mapped_rate),fill=Status,color=Status)) +
  geom_point(size=0.01,position = position_jitterdodge(jitter.width = 0.2)) +
  ylab("Unique Mapping Rate")+
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
        axis.title.y = element_text(colour="black",size=26,margin = margin(1,5,.5,0)),
        text = element_text(colour="black",size=14),
        axis.text.x = element_blank(), 
        #axis.text.x = element_text(colour=c("black"),size=12,vjust=1,hjust = 0.9),
        axis.text.y = element_text(colour="black",size=18),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.1, width = .8,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("firebrick","darkolivegreen"))+
  scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  scale_shape_manual(values = 20)+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))+
  coord_cartesian(ylim=c(0.8,1))
  

dev.off()
###### gt 30 #####

plot_result <- png("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\gt-30.png",
                   width=2800,height=1800,res=450)
total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN)  %>% 
  filter(!gt_30_fraction < 0.25 |gt_30_fraction ==NaN) %>% 
  ggplot(aes(x=factor(Cancer_Type,levels = c("Mammary_Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")),
             y=as.numeric(gt_30_fraction),fill=Status,color=Status)) +
  geom_point(size=0.01,position = position_jitterdodge(jitter.width = 0.2)) +
  ylab("Fraction of mapping quality >30")+
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
        axis.title.y = element_text(colour="black",size=18,margin = margin(1,5,.5,0)),
        text = element_text(colour="black",size=14),
        axis.text.x = element_blank(), 
        #axis.text.x = element_text(colour=c("black"),size=12,vjust=1,hjust = 0.9),
        axis.text.y = element_text(colour="black",size=20),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.1, width = .8,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("firebrick","darkolivegreen"))+
  scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  scale_shape_manual(values = 20)+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))+
  coord_cartesian(ylim=c(0.7,1))


dev.off()



#### Unique CDS ####
plot_result <- png("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\unique-CDS-mapping.png",
                   width=2800,height=1800,res=450)
total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>%
  filter(as.numeric(uniq_CDS_region_paris_rates) >=0.3)%>% 
  ggplot(aes(x=factor(Cancer_Type,levels = c("Mammary_Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")),
             y=as.numeric(uniq_CDS_region_paris_rates),fill=Status,color=Status)) +
  geom_point(size=0.01,position = position_jitterdodge(jitter.width = 0.2)) +
  
  ylab("Unique CDS Mapping Rate")+
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
        axis.title.y = element_text(colour="black",size=20,margin = margin(1,5,.5,0)),
        text = element_text(colour="black",size=14),
        axis.text.x = element_blank(), 
        #axis.text.x = element_text(colour=c("black"),size=12,vjust=1,hjust = 0.9),
        axis.text.y = element_text(colour="black",size=20),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.1, width = .8,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("firebrick","darkolivegreen"))+
  scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  scale_shape_manual(values = 20)+
  theme(plot.margin = unit(c(0.5,0.3,1,0.5), "cm"))+
  coord_cartesian(ylim=c(0.2,0.8))


dev.off()

### change the color of the point ###
plot_result <- png("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\callable-bases.png",
                   width=2800,height=1800,res=450)
exclude <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\V2Supp1_Data.xlsx",
                      sheet = "Total_excluded")
callable %>% 
  filter(!Sample_Name %in% exclude$Cases) %>%
  ggplot(
    aes(x=factor(Cancer_Type,levels = c("Mammary Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")),
        y=as.numeric( Callable_bases)))+
  geom_point(size=1,aes(color=Cancer_Type),position = position_jitterdodge(jitter.width = 0.1))+
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
        axis.title.y = element_text(colour="black",size=26,margin = margin(1,5,0.5,0)),
        plot.margin = margin(1, 10, 4, 5),
        text = element_text(colour="black",size=20),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=20),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.1, width = .8,colour = "black")+
  scale_color_brewer(palette="Dark2")+### for scatter plot
  theme(plot.margin = unit(c(0.5,0.3,1,0.5), "cm"))+
  coord_cartesian(ylim=c(0,3*10**7))
#scale_x_discrete(labels=panel2label)+
#scale_color_manual(values=wes_palette(n=5, name="Set3")) ## for scatter plot Wes Anderson color palettes
#scale_fill_manual(values=wes_palette(n=5, name="Royal2")) ## for box plot, bar plot Wes Anderson color palettes
#scale_shape_manual(values = 20)
dev.off()

### Mean Coverage ###
plot_result <- png("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\coverage-Mean.png",
                   width=2800,height=1800,res=450)
total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>%
  filter(as.numeric(uniq_CDS_region_paris_rates) >=0.3)%>%
  filter(mean>10) %>%
  ggplot(aes(x=factor(Cancer_Type,levels = c("Mammary_Cancer","Melanoma", "Osteosarcoma","Lymphoma","Unclassified")),
             y=as.numeric(mean),fill=Status,color=Status)) +
  geom_point(size=0.01,position = position_jitterdodge(jitter.width = 0.2)) +
  
  ylab("Mean Coverage")+
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
        axis.title.y = element_text(colour="black",size=18,margin = margin(1,5,.5,0)),
        text = element_text(colour="black",size=14),
        axis.text.x = element_blank(), 
        #axis.text.x = element_text(colour=c("black"),size=12,vjust=1,hjust = 0.9),
        axis.text.y = element_text(colour="black",size=18),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.1, width = .8,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("firebrick","darkolivegreen"))+
  scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  scale_shape_manual(values = 20)+
  theme(plot.margin = unit(c(0.5,0.3,1,0.5), "cm"))
  #coord_cartesian(ylim=c(0.2,0.8))


dev.off()



#### RMSE ####
plot_result <- png("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\Coverage-RMSE.png",
                   width=2800,height=1800,res=450)
total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>%
  filter(as.numeric(uniq_CDS_region_paris_rates) >=0.3)%>% 
  filter(mean>10) %>% 
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
        legend.position = "none",
        legend.title = element_blank(),
        legend.key=element_blank(),
        #legend.text =element_text(color = c("firebrick","black")),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.y = element_text(colour="black",size=24,margin = margin(1,5,.5,0)),
        text = element_text(colour="black",size=14),
        axis.text.x = element_blank(), 
        #axis.text.x = element_text(colour=c("black"),size=12,vjust=1,hjust = 0.9),
        axis.text.y = element_text(colour="black",size=20),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.1, width = .8,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("firebrick","darkolivegreen"))+
  scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  scale_shape_manual(values = 20)+
  theme(plot.margin = unit(c(0.5,0.3,1,0.5), "cm"))+
  coord_cartesian(ylim=c(0,0.006))


dev.off()
