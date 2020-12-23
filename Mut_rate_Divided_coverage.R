library(tidyverse)
library(readxl)
library(wesanderson)
library(RColorBrewer)
library(data.table)
library(ggsignif)

fill_colors <- c("darkblue","red3","deeppink");

category <- c("30<x<50","50<x<100","x>=100")
exclude <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\Original_Data_summary.xlsx",
                      sheet ="Before_Matching_excluded")


total_data <- read_excel("C:\\Users\\abc73_000\\Desktop\\New_WES_QC_dataset.xlsx",
                         sheet = "Mut")
total_data <- setDT(total_data)

total_QC <- read_excel("C:\\Users\\abc73_000\\Desktop\\New_WES_QC_dataset.xlsx",
                       sheet = "whole_new")

total_QC <- setDT(total_QC)

a <- total_QC[,(.N), keyby = .(tumor_coverage,Tumor_Type)]



QC_fill_colors <- c("cyan","darkblue","red3","deeppink");

QC_category <- c("<30","30<x<50","50<x<100","x>=100")

pass <- total_data[!Case_ID %in% exclude$Cases & Status =="Non-Problematic", ]

b <- pass[,(med = median(TMB)),keyby = .(tumor_coverage,Tumor_Type)]

a <- pass[Tumor_Type=="MT"& tumor_coverage =="50<x<100"]
b<- pass[Tumor_Type=="MT"& tumor_coverage =="x>=100"]

wilcox.test(a$TMB,b$TMB)

# creating a random dataset
datasets <- c("MT Korean", "MT SNU","MT UGA","OSA Broad", "OSA Tgen","OSA Sanger",
              "OM Cros.Spcs", "OM Sanger",
              "HSA_47 pairs","HSA Upenn", "GLM Cell","LYM Broad","UCL Broad");

tumor_types <- c("MT","OSA","OM","HSA", "GLM","LYM","UCL")

# plotting
dot_size <- 1.6;
abs_text_size <- 16;
regular.text <- element_text(colour="black",size=abs_text_size);
xangle <- 45;
xaxis_just <- ifelse(xangle > 0, 1, 0.5);
# modify this spacing parameters between groups as desired
group_space <- 0.85;
# x = dataset type

pdf("C:\\Users\\abc73_000\\Desktop\\group_tumor_mut_rates.pdf"
    , height=5.0, width=6.84);

## plot group mut rates
x <- factor(pass$Symbol, datasets);
y <- pass$TMB
fill <- factor(pass$tumor_coverage,category);
group <- factor(pass$Tumor_Type, tumor_types);
data <- data.frame(x=x, y=y, fill=fill, group = group)
p <- ggplot(data, aes(x=x, y=y, fill=fill)) + 
  geom_jitter(aes(color=fill), size=dot_size, shape=20, position=position_jitterdodge())+
  ylab("TMB")
p <- p + facet_grid(. ~ group, scales="free_x", space="free_x");

p <- p + stat_summary(fun=median, fun.min=median, fun.max=median, position="dodge", geom="crossbar", size=0.2, width=0.8, color = "black", show.legend=FALSE);
p <- p + scale_fill_manual(values=fill_colors) 
p <- p + scale_color_manual(values=fill_colors)
p <- p+scale_y_log10(limits=c(0.0,100),breaks=c(0.001,0.01,0.1,1,10,100),labels=c(0,0.01,0.1,1,10,100))
p <- p + theme(
  axis.text=regular.text, axis.title=regular.text, 
  axis.text.x = element_text(angle=xangle, hjust=xaxis_just),
  axis.title.x =element_blank(),
  axis.title.y=regular.text,
  legend.position="top",
  legend.title=regular.text, 
  legend.text=regular.text, 
  legend.key=element_blank(),
  panel.background=element_blank(), 
  axis.line=element_line(color="black"), 
  strip.background=element_rect(color="black", fill="transparent", size=1.5), 
  strip.text = element_text(hjust=0.5, size=12, face="plain",color="black"), 
  panel.spacing=unit(group_space, "lines"))

p <- p + coord_cartesian(ylim=c(0.03,100))
p <- p+labs(title = "Tumor Coverage", fill = " ", color =" ")

p <- p+geom_signif(comparisons = list(c("MT Korean", "MT SUN")), 
                   map_signif_level=TRUE)

print(p)

### group both tumor and normal

x <- factor(pass$Symbol, datasets);
y <- pass$TMB
fill <- factor(pass$both_status,category);
group <- factor(pass$Tumor_Type, tumor_types);
data <- data.frame(x=x, y=y, fill=fill, group = group)
p <- ggplot(data, aes(x=x, y=y, fill=fill)) + 
  geom_jitter(aes(color=fill), size=dot_size, shape=20, position=position_jitterdodge())+
  ylab("TMB")
p <- p + facet_grid(. ~ group, scales="free_x", space="free_x");

p <- p + stat_summary(fun=median, fun.min=median, fun.max=median, position="dodge", geom="crossbar", size=0.2, width=0.8, color = "black", show.legend=FALSE);
p <- p + scale_fill_manual(values=fill_colors) 
p <- p + scale_color_manual(values=fill_colors)
p <- p+scale_y_log10(limits=c(0.0,100),breaks=c(0.001,0.01,0.1,1,10,100),labels=c(0,0.01,0.1,1,10,100))
p <- p + theme(
  axis.text=regular.text, axis.title=regular.text, 
  axis.text.x = element_text(angle=xangle, hjust=xaxis_just),
  axis.title.x =element_blank(),
  axis.title.y=regular.text,
  legend.position="top",
  legend.title=regular.text, 
  legend.text=regular.text, 
  legend.key=element_blank(),
  panel.background=element_blank(), 
  axis.line=element_line(color="black"), 
  strip.background=element_rect(color="black", fill="transparent", size=1.5), 
  strip.text = element_text(hjust=0.5, size=12, face="plain",color="black"), 
  panel.spacing=unit(group_space, "lines"))

p <- p + coord_cartesian(ylim=c(0.03,100))
p <- p+labs(title = "Both Tumor Normal Coverage", fill = " ", color =" ")
print(p)

dev.off()
### Tumor Coverage no group ##

pass <- pass[both_status!="<30",]
x <- factor(pass$Tumor_Type, tumor_types);
y <- pass$TMB
fill <- factor(pass$tumor_coverage,category);
data <- data.frame(x=x, y=y, fill=fill)

p <- ggplot(data, aes(x=x, y=y, fill=fill,color=fill)) + 
  geom_jitter(aes(color=fill), size=dot_size, shape=20, position=position_jitterdodge())+
  ylab("TMB")
#p <- p + facet_grid(. ~ group, scales="free_x", space="free_x");

p <- p + stat_summary(fun=median, fun.min=median, fun.max=median, position="dodge", geom="crossbar", size=0.2, width=0.8, color = "black", show.legend=FALSE);
p <- p + scale_fill_manual(values=fill_colors) 
p <- p + scale_color_manual(values=fill_colors)
p <- p+scale_y_log10(limits=c(0.0,100),breaks=c(0.001,0.01,0.1,1,10,100),labels=c(0,0.01,0.1,1,10,100))
p <- p + theme(axis.text=regular.text, 
               axis.title.y=regular.text,
               axis.title.x =element_blank(),
               axis.text.x = element_text(angle=30, hjust=1), 
               panel.background=element_blank(), 
               axis.line=element_line(color="black"),
               legend.position="top", 
               legend.title=regular.text, 
               legend.text=regular.text, 
               legend.key=element_blank())

p <- p + coord_cartesian(ylim=c(0.03,100))
p <- p+labs(title = "Tumor Coverage", color = " ", fill = " ")
print(p)

##### Both tumor and normal

pass <- pass[both_status!="<30",]
x <- factor(pass$Tumor_Type, tumor_types);
y <- pass$TMB
fill <- factor(pass$both_status,category);
data <- data.frame(x=x, y=y, fill=fill)

p <- ggplot(data, aes(x=x, y=y, fill=fill,color=fill)) + 
  geom_jitter(aes(color=fill), size=dot_size, shape=20, position=position_jitterdodge())+
  ylab("TMB")
#p <- p + facet_grid(. ~ group, scales="free_x", space="free_x");

p <- p + stat_summary(fun=median, fun.min=median, fun.max=median, position="dodge", geom="crossbar", size=0.2, width=0.8, color = "black", show.legend=FALSE);
p <- p + scale_fill_manual(values=fill_colors) 
p <- p + scale_color_manual(values=fill_colors)
p <- p+scale_y_log10(limits=c(0.0,100),breaks=c(0.001,0.01,0.1,1,10,100),labels=c(0,0.01,0.1,1,10,100))
p <- p + theme(axis.text=regular.text, 
               axis.title.y=regular.text,
               axis.title.x =element_blank(),
               axis.text.x = element_text(angle=30, hjust=1), 
               panel.background=element_blank(), 
               axis.line=element_line(color="black"),
               legend.position="top", 
               legend.title=regular.text, 
               legend.text=regular.text, 
               legend.key=element_blank())

p <- p + coord_cartesian(ylim=c(0.03,100))
p <- p+labs(title = "Both Tumor and Normal Coverage", color = " ", fill = " ")
print(p)

dev.off()

#### QC Coverage plot ###

## Tumor Coverage
dot_size <- 1.0;
abs_text_size <- 16;
regular.text <- element_text(colour="black",size=abs_text_size);
xangle <- 45;
xaxis_just <- ifelse(xangle > 0, 1, 0.5);
# modify this spacing parameters between groups as desired
group_space <- 0.85;
# x = dataset type
pdf("C:\\Users\\abc73_000\\Desktop\\tumor_QC.pdf"
    , height=5.0, width=6.84);

total_QC <- total_QC[Total_pairs >= 5000000 & Uniquely_mapped_rate>0.6 & uniq_CDS_region_paris_rates>0.3]
x <- factor(total_QC$Symbol, datasets);
y <- total_QC$mean
fill <- factor(total_QC$Status, c("Normal", "Tumor"));
group <- factor(total_QC$tumor_coverage, QC_category);
data <- data.frame(x=x, y=y, fill=fill, group=group);
p <- ggplot(data, aes(x=x, y=y, fill=fill)) + 
  geom_jitter(aes(color=fill), size=dot_size, shape=20, position=position_jitterdodge())+
  ylab("Mean coverage")
p <- p+ggtitle("Tumor Coverage")
p <- p + facet_grid(. ~ group, scales="free_x", space="free_x");
p <- p + stat_summary(fun=median, fun.min=median, fun.max=median, position="dodge", geom="crossbar", size=0.2, width=0.8, color = "black", show.legend=FALSE);
p <- p + scale_fill_manual(guide=FALSE, values=fill_colors, drop=FALSE) 
p <- p + scale_color_manual(name="", values=fill_colors, drop=FALSE, guide=FALSE)
p <- p+scale_y_continuous(breaks = c(0,100,200))+coord_cartesian(ylim=c(0,200))
p <- p + theme(
  axis.text=regular.text, axis.title=regular.text, 
  axis.text.x = element_text(angle=xangle, hjust=xaxis_just,size = 10),
  axis.title.x =element_blank(),
  axis.title.y=regular.text,
  legend.position="top",
  legend.title=regular.text, 
  legend.text=regular.text, 
  legend.key=element_blank(),
  panel.background=element_blank(), 
  axis.line=element_line(color="black"), 
  strip.background=element_rect(color="black", fill="transparent", size=1.5), 
  strip.text = element_text(hjust=0.5, size=12, face="plain",color="black"), 
  panel.spacing=unit(group_space, "lines"))+
  geom_hline(yintercept=30, linetype="longdash", color = "yellow4", size = 0.7)

p+labs(title = "Tumor Coverage", color = " ", fill = " ")
print(p)

### Both Coverage 

total_QC <- total_QC[Total_pairs >= 5000000 & Uniquely_mapped_rate>0.6 & uniq_CDS_region_paris_rates>0.3]
x <- factor(total_QC$Symbol, datasets);
y <- total_QC$mean
fill <- factor(total_QC$Status, c("Normal", "Tumor"));
group <- factor(total_QC$both_status, QC_category);
data <- data.frame(x=x, y=y, fill=fill, group=group);
p <- ggplot(data, aes(x=x, y=y, fill=fill)) + 
  geom_jitter(aes(color=fill), size=dot_size, shape=20, position=position_jitterdodge())+
  ylab("Mean coverage")
p <- p+ggtitle("Both Tumor and Normal Coverage")
p <- p + facet_grid(. ~ group, scales="free_x", space="free_x");
p <- p + stat_summary(fun=median, fun.min=median, fun.max=median, position="dodge", geom="crossbar", size=0.2, width=0.8, color = "black", show.legend=FALSE);
p <- p + scale_fill_manual(guide=FALSE, values=fill_colors, drop=FALSE) 
p <- p + scale_color_manual(name="", values=fill_colors, drop=FALSE, guide=FALSE)
p <- p+scale_y_continuous(breaks = c(0,100,200))+coord_cartesian(ylim=c(0,200))
p <- p + theme(
  axis.text=regular.text, axis.title=regular.text, 
  axis.text.x = element_text(angle=xangle, hjust=xaxis_just,size = 10),
  axis.title.x =element_blank(),
  axis.title.y=regular.text,
  legend.position="top",
  legend.title=regular.text, 
  legend.text=regular.text, 
  legend.key=element_blank(),
  panel.background=element_blank(), 
  axis.line=element_line(color="black"), 
  strip.background=element_rect(color="black", fill="transparent", size=1.5), 
  strip.text = element_text(hjust=0.5, size=12, face="plain",color="black"), 
  panel.spacing=unit(group_space, "lines"))+
  geom_hline(yintercept=30, linetype="longdash", color = "yellow4", size = 0.7)

print(p)



dev.off()


