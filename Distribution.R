DrawDistribution <- function(file_location, withlimit, xupper,file_name){
  library(ggplot2)
  library(readxl)
  library(dplyr)
  title <- strsplit(file_name,split ="_")[[1]][1]
  if (withlimit){
  sample <- read.table(file_location, sep =" ",header =F,
                       stringsAsFactors = F)
  sample <- sample[-2, ]
  final <- sample[0:xupper, ]
  yupper = as.numeric(max(final[1]))
  colnames(final)=c('Frequency',"Position")
  freq = as.numeric(final$Frequency)
  pos = as.numeric(final$Position)
  p <- ggplot()+ 
    geom_line(data = final,aes(x=pos, y=freq),stat="identity")+
    ylim(0,yupper)+
    ggtitle(paste(title,"Sequence Randomness",sep="_"))+
    ylab("# of Position")+
    xlab("The frequency of being cover")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(color = "black", size = 20, face = "bold"))
  return(plot(p))
  }
  else{
    sample <- read.table(file_location, sep =" ",header =F,
                         stringsAsFactors = F)
    sample <- sample[-2, ]
    final <- sample
    colnames(final)=c('Frequency',"Position")
    freq = as.numeric(final$Frequency)
    pos = as.numeric(final$Position)
    p <- ggplot()+ 
      geom_line(data = final,aes(x=pos, y=freq),stat="identity")+
      ggtitle(paste(title,"Sequence Randomness",sep="_"))+
      ylab("# of Position")+
      xlab("Cover time")+
      theme(axis.text=element_text(size=10),
            axis.title=element_text(size=14,face="bold"),
            plot.title = element_text(color = "black", size = 16, face = "bold"))
    return(plot(p))
    }
}



# DrawDistribution("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Distribution/Mammary/Normal/SRR7780741_DepthofCoverage_Distribution.txt",
#                  T,600)
# 
# dev.off()
Df  <- "/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Distribution/Lym/Lymphoma_Tumor_Sequence_distribution.pdf"
pdf(file=Df, w=7, h=5)
file_list <- read.table("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Distribution/Lym/Tumor/Tumor_list",
                        sep = "", header = F, stringsAsFactors = F)
for (i in as.vector(file_list$V1)){
  
DrawDistribution(paste("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Distribution/Lym/Tumor/",i,sep = ""),
                   T, 600,i)
  #print(paste("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Distribution/Mammary/Normal/",i,sep = ""))
 
}
dev.off()
