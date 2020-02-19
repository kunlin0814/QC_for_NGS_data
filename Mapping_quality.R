library(ggplot2)
mapping <- read.table("/Users/kun-linho/Desktop/SRR7781089_mapping_quality", header = F)
new_mapping <- mapping$V1[1:10000]
new_mapping <- as.data.frame(new_mapping)
png("SRR7781089_mapping_quality.png",width=3000,height=2400,res=400)
ggplot(data=mapping, aes(x= mapping$V1))+
  geom_histogram(bins= 100, na.rm = T)+
  xlab("Mapping Quality Scores")+
  ylab("Frequency")+
  xlim(0,70)+
  theme(
    #axis.title=element_text(size=18,face="bold"),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=16,hjust = 1),
    axis.text.y = element_text(size=16,hjust = 1))
dev.off()
#mapping[1:10,]
