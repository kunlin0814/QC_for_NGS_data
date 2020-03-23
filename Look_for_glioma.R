library(tidyverse)
library(readxl)
glioma_meta <- read_excel("/Volumes/Research_Data/MAC_Research_Data/Pan_cancer/Pan_cancer_mapping_result/Mutation_rate/Glioma/ScienceDirect_files_16Mar2020_18-34-11.365/1-s2.0-S153561082030043X-mmc2.xlsx",
                          sheet ="Case_Metadata")
current_cases <- read.table("/Volumes/Research_Data/MAC_Research_Data/Pan_cancer/Pan_cancer_mapping_result/Mutation_rate/Glioma/glioma_mutect_cases.txt",
                            sep = '\n')

glioma_current <- glioma_meta %>% 
  filter(Tumor_Sample_Barcode %in% current_cases$V1) %>% 
  select(Tumor_Sample_Barcode, Dog_Breed) %>% 
  write.table("/Volumes/Research_Data/MAC_Research_Data/Pan_cancer/Pan_cancer_mapping_result/Mutation_rate/Glioma/current_glioma_breed.txt",
              quote = F, sep="\t",
              row.names = F)

glioma_mut_rate <- read_excel("/Volumes/Research_Data/MAC_Research_Data/Pan_cancer/Pan_cancer_mapping_result/Mutation_rate/Mutation_rate.xlsx",
                              sheet = "Glioma")
plot_result <- png("/Volumes/Research_Data/MAC_Research_Data/Pan_cancer/Pan_cancer_mapping_result/Mutation_rate/Glioma/glioma_mut_breeds.png",width=3000,height=2400,res=400)  

glioma_mut_rate %>% 
  na.omit() %>% 
  ggplot(aes(x=Breed,
        y=as.numeric(non_retro_mutation_rate)))+
  geom_point(size=.5,aes(color=Breed),position = position_jitterdodge(jitter.width = 0.1))+
  ylab("Mutation Rates Per-Million Bases")+
  ylim(0, 3)+
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
        text = element_text(colour="black",size=12),
        axis.text.x = element_text(colour=c("black"),
                                   size=12,angle=45,vjust=1,hjust = 0.9),
        axis.text.y = element_text(colour="black",size=14),
        panel.background = element_blank())
+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.1, width = .8,colour = "black")+
  scale_color_brewer(palette="Dark2") ### for scatter plot
#scale_x_discrete(labels=panel2label)+
#scale_color_manual(values=wes_palette(n=5, name="Set3")) ## for scatter plot Wes Anderson color palettes
#scale_fill_manual(values=wes_palette(n=5, name="Royal2")) ## for box plot, bar plot Wes Anderson color palettes
#scale_shape_manual(values = 20)
dev.off()
