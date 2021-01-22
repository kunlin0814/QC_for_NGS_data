library(tidyverse)
library(readxl)
library(wesanderson)
library(RColorBrewer)
library(data.table)

sep_field = "/"

base_dir <- "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/WGS_analysis"
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/WGS_analysis"
total_file <- fread(paste(base_dir,"WGS_final_table.txt",sep = sep_field))


exclude_sample <- total_file[The_reason_to_exclude!="Pass QC", .(Case_ID)]

summary <- NULL

for (i in 1:nrow(total_file)){
  case_name <- total_file$Case_ID[i]
  status <- total_file$Status[i]
  if (status =="Normal"){
    value = "Normal_sample"
    
  }
  else{
    index = match(case_name,total_callable$V1,nomatch = 0)
    if (index !=0){
      value = as.character(total_callable$V2[index])
      print(value)
    }
  }
  summary <- c(summary,value)
}

total_file$callable_bases <- summary

# fwrite(total_file,file= paste(base_dir,"final_WGS_summary.txt",sep=sep_field),quote = F,
#        col.names = T,sep ="\t")
# 
# # uni <- unique(total_file[Unique_mapped_rate > 0.6 & total_mean_coverage <30,.(Case_ID)])
# 
# a <- total_file[Unique_mapped_rate < 0.6,]

# pdf(paste(base_dir,"WGS_QC.pdf",sep =sep_field),
#      height=5.98, width=6.84);

regular.text <- element_text(colour="black",size=24);
file_color <- c("darkblue","red3")

png(file = paste(base_dir,"Total_Sequence_pair.png",sep =sep_field),
     width = 3500, height =3000, units = "px", res = 500)


total_file %>% 
  #filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("GLM","OM","OSA")),
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
        legend.position="None", 
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

dev.off()
# Gt_30_fraction
png(file = paste(base_dir,"gt30.png",sep =sep_field),
    width = 3500, height =3000, units = "px", res = 500)

total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("GLM","OM","OSA")),
             y=as.numeric(gt30),fill=Status,color=Status)) +
  
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("Fraction of \nmapping quality >30")+
  scale_y_continuous(breaks = c(0,0.5,0.75,1.0))+coord_cartesian(ylim=c(0,1))+
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
        legend.position="None", 
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
dev.off() 

## Unique concordantly map rate
png(file = paste(base_dir,"unique_mapped.png",sep =sep_field),
    width = 3500, height =3000, units = "px", res = 500)

total_file %>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("GLM","OM","OSA")),
             y=as.numeric(total_file$Unique_mapped_rate),fill=Status,color=Status)) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("Uniquely & concordantly\n mapped rate")+
  scale_y_continuous(breaks = c(0,0.5,1.0))+
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
        legend.position="None", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = file_color)+
  scale_fill_manual(values=file_color)+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  #scale_shape_manual(values = 20)+
  coord_cartesian(ylim=c(0,1))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))+
  geom_hline(yintercept=0.6, linetype="longdash", color = "yellow4", size = 0.7)
dev.off()

## mean coverage

png(file = paste(base_dir,"mean_coverage.png",sep =sep_field),
    width = 3500, height =3000, units = "px", res = 500)

total_file %>% 
  filter(Unique_mapped_rate >=0.6) %>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("GLM","OM","OSA")),
             y=as.numeric(total_mean_coverage),fill=Status,color=Status)) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  #scale_y_continuous(breaks = c(0,100,200))+
  ylab("Mean coverage")+
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
        legend.position="None", 
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
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))+
  geom_hline(yintercept=30, linetype="longdash", color = "yellow4", size = 0.7)

dev.off()
# 
# total_file %>% 
#   ggplot(aes(x=factor(Tumor_Type,levels = c("MT","GLM","OM")),
#              y=as.numeric(unique_mapped_mean_coverage),fill=Status,color=Status)) +
#   geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
#   #scale_y_continuous(breaks = c(0,100,200))+
#   ylab("Total uniquely Reads Mean coverage")+
#   #labs(subtitle = "p<0.01")+
#   #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
#   #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
#   theme(axis.text=regular.text, 
#         axis.title.y=regular.text,
#         axis.title.x =element_blank(),
#         axis.text.x = element_text(angle=30, hjust=1), 
#         panel.background=element_blank(), 
#         axis.line=element_line(color="black"),
#         legend.title=regular.text, 
#         legend.position="top", 
#         legend.text=regular.text, 
#         legend.key=element_blank())+
#   stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
#                geom = "crossbar",size=0.5, width = .7,colour = "black")+
#   #scale_x_discrete(labels=panel2label)+
#   scale_color_manual(values = file_color)+
#   scale_fill_manual(values=file_color)+
#   #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
#   #scale_shape_manual(values = 19)+
#   #coord_cartesian(ylim=c(0,200))+
#   theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))
#   #geom_hline(yintercept=30, linetype="longdash", color = "yellow4", size = 0.7)

png(file = paste(base_dir,"RMSE.png",sep =sep_field),
    width = 3500, height =3000, units = "px", res = 500)


total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>% 
  filter(Unique_mapped_rate >=0.6) %>% 
  filter(total_mean_coverage >=30) %>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("GLM","OM","OSA")),
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
        legend.position="None", 
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
## Callable bases

png(file = paste(base_dir,"Callable_bases.png",sep =sep_field),
    width = 3500, height =3000, units = "px", res = 500)

total_file %>% 
  filter(Status=='Tumor') %>% 
  filter(!Case_ID %in% exclude_sample$Case_ID) %>% 
  filter(!callable_bases =='No-Pair')%>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("GLM","OM","OSA")),
             y=as.numeric(callable_bases)/1000000,color='black')) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("Callable bases in millions")+
  scale_y_continuous(breaks = c(0,1000,2000))+
  theme(axis.text=regular.text, 
        axis.title.y=element_text(colour="black",size=20, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="none", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.4, width = .7,colour = "black")+
  #scale_shape_manual(values = 19)+
  scale_color_manual(values = 'black')+
  theme(plot.margin = unit(c(0.5,0.3,1,0.5), "cm"))
  #coord_cartesian(ylim=c(5,30))+
 

dev.off()


