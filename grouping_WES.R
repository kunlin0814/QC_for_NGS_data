library(tidyverse)
library(readxl)
library(wesanderson)
library(RColorBrewer)
library(data.table)

base_dir <- "/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Figure1"

total_data <- read_excel(paste(base_dir,"New_WES_QC_dataset.xlsx",sep ="/"),
  #"C:\\Users\\abc73_000\\Desktop\\New_WES_QC_dataset.xlsx",
                         sheet = "Sheet1")
total_data <- setDT(total_data)
total_data <- total_data[Case_ID!="No-Pair"]
# modify this line, put the fill colors for normal and tumor
fill_colors <- c("darkblue","red3");

# creating a random dataset
datasets <- c("MT Korean", "MT SNU","OSA Broad", "OSA Tgen","OSA Sanger",
              "OM Cros.Spcs", "OM Sanger",
              "HSA Broad","HSA UPenn", "GLM Cell","LYM Broad","UCL Broad");

tumor_types <- c("MT", "GLM", "LYM", "OM", "OSA", "HSA" ,"UCL")

# plotting
dot_size <- 1.4;
abs_text_size <- 20;
regular.text <- element_text(colour="black",size=abs_text_size);
xangle <- 45;
xaxis_just <- ifelse(xangle > 0, 1, 0.5);
# modify this spacing parameters between groups as desired
group_space <- 0.85;
# x = dataset type



pdf(paste(base_dir,"V17F1_and_supplementaryF1.pdf",sep ="/")
    , height=5.0, width=6.84);

## plot Total reads
x <- factor(total_data$Symbol, datasets);
y <- total_data$Total_pairs
fill <- factor(total_data$Status, c("Normal", "Tumor"));
# group = cancer type
group <- factor(total_data$Tumor_Type, tumor_types);
data <- data.frame(x=x, y=y, fill=fill, group=group);
p <- ggplot(data, aes(x=x, y=y/1000000, fill=fill)) + 
  geom_jitter(aes(color=fill), size=dot_size, shape=20, position=position_jitterdodge())+
  ylab("Sequence read pairs\nin millions")
p <- p + facet_grid(. ~ group, scales="free_x", space="free_x");
p <- p + stat_summary(fun=median, fun.min=median, fun.max=median, position="dodge", geom="crossbar", size=0.2, width=0.8, color = "black", show.legend=FALSE);
p <- p + scale_fill_manual(guide=FALSE, values=fill_colors, drop=FALSE) 
p <- p + scale_color_manual(name="", values=fill_colors, drop=FALSE, guide=FALSE)
p <- p+scale_y_continuous(breaks = c(0,100,200))+coord_cartesian(ylim=c(0,200))
p <- p + theme(
  axis.text=regular.text, 
  axis.title=regular.text, 
  axis.text.x = element_text(angle=xangle, hjust=xaxis_just),
  axis.title.x =element_blank(),
  axis.title.y= regular.text,
  legend.position="top",
  legend.title=regular.text, 
  legend.text=regular.text, 
  legend.key=element_blank(),
  panel.background=element_blank(), 
  axis.line=element_line(color="black"), 
  strip.background=element_rect(color="black", fill="transparent", size=1.5), 
  strip.text = element_text(hjust=0.5, size=12, face="plain",color="black"), 
  panel.spacing=unit(group_space, "lines"))
p <- p+geom_hline(yintercept=5, linetype="longdash", color = "yellow4", size = 0.7)


print(p)
# Gt_30_fraction
total_data <- total_data[!Total_pairs < 5000000 | Total_pairs==NaN, ]

x <- factor(total_data$Symbol, datasets);
y <- total_data$Gt_30_fraction
fill <- factor(total_data$Status, c("Normal", "Tumor"));
# group = cancer type
group <- factor(total_data$Tumor_Type, tumor_types);
data <- data.frame(x=x, y=y, fill=fill, group=group);
group <- factor(total_data$Tumor_Type, tumor_types);
data <- data.frame(x=x, y=y, fill=fill, group=group);
p <- ggplot(data, aes(x=x, y=y, fill=fill)) + 
  geom_jitter(aes(color=fill), size=dot_size, shape=20, position=position_jitterdodge())+
  ylab("Fraction of \nmapping quality >30")
p <- p + facet_grid(. ~ group, scales="free_x", space="free_x");
p <- p + stat_summary(fun=median, fun.min=median, fun.max=median, position="dodge", geom="crossbar", size=0.2, width=0.8, color = "black", show.legend=FALSE);
p <- p + scale_fill_manual(guide=FALSE, values=fill_colors, drop=FALSE) 
p <- p + scale_color_manual(name="", values=fill_colors, drop=FALSE, guide=FALSE)
p <- p+scale_y_continuous(breaks = c(0,0.5,0.75,1.0))+coord_cartesian(ylim=c(0.3,1))
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


print(p)

# Unique mapped rates

total_data <- total_data[!Total_pairs < 5000000 | Total_pairs==NaN, ]

x <- factor(total_data$Symbol, datasets);
y <- total_data$Uniquely_coordinatly_mapped_rate
fill <- factor(total_data$Status, c("Normal", "Tumor"));
# group = cancer type
group <- factor(total_data$Tumor_Type, tumor_types);
data <- data.frame(x=x, y=y, fill=fill, group=group);
p <- ggplot(data, aes(x=x, y=y, fill=fill)) + 
  geom_jitter(aes(color=fill), size=dot_size, shape=20, position=position_jitterdodge())+
  ylab("Uniquely & concordantly\n mapped rate")
p <- p + facet_grid(. ~ group, scales="free_x", space="free_x");
p <- p + stat_summary(fun=median, fun.min=median, fun.max=median, position="dodge", geom="crossbar", size=0.2, width=0.8, color = "black", show.legend=FALSE);
p <- p + scale_fill_manual(guide=FALSE, values=fill_colors, drop=FALSE) 
p <- p + scale_color_manual(name="", values=fill_colors, drop=FALSE, guide=FALSE)
p <- p+scale_y_continuous(breaks = c(0,0.5,0.75,1.0))+coord_cartesian(ylim=c(0.3,1))
p <- p + theme(
  axis.text=regular.text, 
  axis.title=regular.text, 
  axis.text.x = element_text(angle=xangle, hjust=xaxis_just),
  axis.title.x =element_blank(),
  axis.title.y= regular.text,
  legend.position="top",
  legend.title=regular.text, 
  legend.text=regular.text, 
  legend.key=element_blank(),
  panel.background=element_blank(), 
  axis.line=element_line(color="black"), 
  strip.background=element_rect(color="black", fill="transparent", size=1.5), 
  strip.text = element_text(hjust=0.5, size=12, face="plain",color="black"), 
  panel.spacing=unit(group_space, "lines"))
p <- p+geom_hline(yintercept=0.6, linetype="longdash", color = "yellow4", size = 0.7)


print(p)

## CDS Targeting
total_data <- total_data[Total_pairs >= 5000000 & Uniquely_coordinatly_mapped_rate>0.6, ]

x <- factor(total_data$Symbol, datasets);
y <- total_data$Target_CDS_Mapping_Rates
fill <- factor(total_data$Status, c("Normal", "Tumor"));
group <- factor(total_data$Tumor_Type, tumor_types);
data <- data.frame(x=x, y=y, fill=fill, group=group);
p <- ggplot(data, aes(x=x, y=y, fill=fill)) + 
  geom_jitter(aes(color=fill), size=dot_size, shape=20, position=position_jitterdodge())+
  ylab("CDS targeting rate")
p <- p + facet_grid(. ~ group, scales="free_x", space="free_x");
p <- p + stat_summary(fun=median, fun.min=median, fun.max=median, position="dodge", geom="crossbar", size=0.2, width=0.8, color = "black", show.legend=FALSE);
p <- p + scale_fill_manual(guide=FALSE, values=fill_colors, drop=FALSE) 
p <- p + scale_color_manual(name="", values=fill_colors, drop=FALSE, guide=FALSE)
p <- p+scale_y_continuous(breaks = c(0,0.4,0.8))+coord_cartesian(ylim=c(0,0.8))
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
p <- p+geom_hline(yintercept=0.3, linetype="longdash", color = "yellow4", size = 0.7)


print(p)
## Mean Coverage
total_data <- total_data[Total_pairs >= 5000000 & Uniquely_coordinatly_mapped_rate>0.6 & Target_CDS_Mapping_Rates>0.3]
x <- factor(total_data$Symbol, datasets);
y <- total_data$Mean
fill <- factor(total_data$Status, c("Normal", "Tumor"));
group <- factor(total_data$Tumor_Type, tumor_types);
data <- data.frame(x=x, y=y, fill=fill, group=group);
p <- ggplot(data, aes(x=x, y=y, fill=fill)) + 
  geom_jitter(aes(color=fill), size=dot_size, shape=20, position=position_jitterdodge())+
  ylab("Mean coverage")
p <- p + facet_grid(. ~ group, scales="free_x", space="free_x");
p <- p + stat_summary(fun=median, fun.min=median, fun.max=median, position="dodge", geom="crossbar", size=0.2, width=0.8, color = "black", show.legend=FALSE);
p <- p + scale_fill_manual(guide=FALSE, values=fill_colors, drop=FALSE) 
p <- p + scale_color_manual(name="", values=fill_colors, drop=FALSE, guide=FALSE)
p <- p+scale_y_continuous(breaks = c(0,100,200))+coord_cartesian(ylim=c(0,200))
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
  panel.spacing=unit(group_space, "lines"))+
  geom_hline(yintercept=30, linetype="longdash", color = "yellow4", size = 0.7)

print(p)
### RMSE
total_data <- total_data[Total_pairs >= 5000000 & Uniquely_coordinatly_mapped_rate>=0.6
                         & Target_CDS_Mapping_Rates >=0.3 & Mean >=30] 
                         
x <- factor(total_data$Symbol, datasets);
y <- total_data$RMSE
fill <- factor(total_data$Status, c("Normal", "Tumor"));
group <- factor(total_data$Tumor_Type, tumor_types);
data <- data.frame(x=x, y=y, fill=fill, group=group);

group <- factor(total_data$Tumor_Type, tumor_types);
data <- data.frame(x=x, y=y, fill=fill, group=group);
p <- ggplot(data, aes(x=x, y=y, fill=fill)) + 
  geom_jitter(aes(color=fill), size=dot_size, shape=20, position=position_jitterdodge())+
  ylab("RMSE of \nsequence coverage")
p <- p + facet_grid(. ~ group, scales="free_x", space="free_x");
p <- p + stat_summary(fun=median, fun.min=median, fun.max=median, position="dodge", geom="crossbar", size=0.2, width=0.8, color = "black", show.legend=FALSE);
p <- p + scale_fill_manual(guide=FALSE, values=fill_colors, drop=FALSE) 
p <- p + scale_color_manual(name="", values=fill_colors, drop=FALSE, guide=FALSE)
p <- p+scale_y_continuous(breaks = c(0,0.005,0.01))+coord_cartesian(ylim=c(0,0.01))
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


print(p)

exclude <- read_excel(paste(base_dir,"Original_Data_summary.xlsx",sep="/"),
                      sheet ="Before_Matching_excluded")

total_data <- total_data %>% 
  filter(Case_ID != "No-Pair") %>% 
  filter(Status=='Tumor') %>% 
  filter(!Case_ID %in% exclude$Cases)

## callable
x <- factor(total_data$Symbol, datasets);
y <- total_data$Callable_bases
# group = cancer type
group <- factor(total_data$Tumor_Type, tumor_types);
data <- data.frame(x=x, y=as.numeric(y), group=group);
p <- ggplot(data, aes(x=x, y=y/1000000, color='black')) + 
  geom_jitter(size=1.6, shape=20, position=position_jitterdodge())+
  ylab("Callable bases in millions")
p <- p + facet_grid(. ~ group, scales="free_x", space="free_x");
p <- p + stat_summary(fun=median, fun.min=median, fun.max=median, position="dodge", geom="crossbar", size=0.2, width=0.8, color = "black", show.legend=FALSE);
p <- p + scale_fill_manual(guide=FALSE, values='black', drop=FALSE) 
p <- p + scale_color_manual(name="", values='black', drop=FALSE, guide=FALSE)
p <- p+scale_y_continuous(breaks = c(0,10,20,30))+coord_cartesian(ylim=c(5,30))
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
  panel.spacing=unit(group_space, "lines"))+
  geom_hline(yintercept=10, linetype="longdash", color = "yellow4", size = 0.7)

print(p)
#dev.off()

dev.off()
