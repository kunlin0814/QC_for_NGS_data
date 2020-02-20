data <- read.table("G:\\Pan_cancer\\Pan_cancer_mapping_result\\Distribution\\Mammary\\Normal\\SRR7780741_DepthofCoverage_Distribution.txt",
                   sep = " ",
                   header = F)
data <- data[-2,]


last_value <- as.numeric(as.vector(tail(data$V2, n= 1)))

for (i in 0: last_value){
  if(as.numeric(as.vector(data$V2[i+1]))== i){
    print("same")
  }
  else{
    print(i)
    data <- rbind(data, c(0, i))
  }

}
data <- cbind(data, as.data.frame(c(1,2)))


