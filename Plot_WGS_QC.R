library(tidyverse)
library(readxl)
library(wesanderson)
library(RColorBrewer)
library(data.table)


total_file <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\WGS_analysis\\WGS_QC.xlsx",
                         sheet = "Sheet1")

pdf("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\WGS_analysis\\WGS_QC.pdf"
    , height=4.98, width=4.84);

regular.text <- element_text(colour="black",size=20);
file_color <- c("darkblue","red3")


total_file %>% 
  #filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("MT","GLM","OM")),
             y=as.numeric(Total_pairs)/1000000,fill=Status,color=Status)) +
  
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("Sequence read pairs\nin millions")+
  #scale_y_continuous(breaks = c(0,100,200))+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="top", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = file_color)+
  scale_fill_manual(values=file_color)+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  #scale_shape_manual(values = 20)+
  #coord_cartesian(ylim=c(0,200))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))
  #geom_hline(yintercept=5, linetype="longdash", color = "yellow4", size = 0.7)

total_file %>% 
  #filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("MT","GLM","OM")),
             y=as.numeric(Total_pairs)/1000000,fill=Status,color=Status)) +
  
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("Sequence read pairs\nin millions")+
  #scale_y_continuous(breaks = c(0,100,200))+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="top", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = file_color)+
  scale_fill_manual(values=file_color)+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  #scale_shape_manual(values = 20)+
  #coord_cartesian(ylim=c(0,200))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))
#geom_hline(yintercept=5, linetype="longdash", color = "yellow4", size = 0.7)

total_file %>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("MT","GLM","OM")),
             y=as.numeric(total_mean_coverage),fill=Status,color=Status)) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  #scale_y_continuous(breaks = c(0,100,200))+
  ylab("Total Reads Mean coverage")+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="top", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = file_color)+
  scale_fill_manual(values=file_color)+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  #scale_shape_manual(values = 19)+
  #coord_cartesian(ylim=c(0,200))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))
  #geom_hline(yintercept=30, linetype="longdash", color = "yellow4", size = 0.7)

total_file %>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("MT","GLM","OM")),
             y=as.numeric(unique_mapped_mean_coverage),fill=Status,color=Status)) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  #scale_y_continuous(breaks = c(0,100,200))+
  ylab("Total uniquely Reads Mean coverage")+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="top", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = file_color)+
  scale_fill_manual(values=file_color)+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  #scale_shape_manual(values = 19)+
  #coord_cartesian(ylim=c(0,200))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))
  #geom_hline(yintercept=30, linetype="longdash", color = "yellow4", size = 0.7)


total_file %>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("MT","GLM","OM")),
             y=as.numeric(Unique_mapped_rate),fill=Status,color=Status)) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("Uniquely & concordantly\n mapped rate")+
  scale_y_continuous(breaks = c(0,0.4,0.8))+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="top", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = file_color)+
  scale_fill_manual(values=file_color)+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  #scale_shape_manual(values = 20)+
  #coord_cartesian(ylim=c(0,0.8))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))


total_file %>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("MT","GLM","OM")),
             y=as.numeric(RMSE),fill=Status,color=Status)) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  #scale_y_continuous(breaks = c(0,0.005,0.01))+
  ylab("RMSE of \nsequence coverage")+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="top", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.4, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = file_color)+
  scale_fill_manual(values=file_color)+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  #scale_shape_manual(values = 19)+
  #coord_cartesian(ylim=c(0,0.01))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))


dev.off()
