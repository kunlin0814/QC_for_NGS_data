library(ggplot2)
library(readxl)
library(dplyr)
library(wesanderson)
library(RColorBrewer)
# glioma <- read_excel("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Supplement_Figure1/Supp1_Data.xlsx",
#                      sheet= "Sheet1")
# new_glioma <- na.omit(glioma)
# mean(new_glioma$`0...2`[!new_glioma$`0...2`==0])

total_file <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\V2Supp1_Data.xlsx",
                         sheet ="Total")
PAIR <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\V2Supp1_Data.xlsx",
                   sheet ="PAIRS")
#callable <- read_excel("/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Supplement_Figure1/V2Supp1_Data.xlsx",
#                       sheet ="Callable_Bases")
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
  filter(gt_30_fraction < 0.25) %>% 
  #filter(as.numeric(uniq_CDS_region_paris_rates) < 0.3|uniq_CDS_region_paris_rates==NaN) %>% 
  select(ID, gt_30_fraction,Status,Cancer_Type)

############# CDS #############

filtered <- total_file %>% 
filter(!Total_pairs < 5000000 | Total_pairs==NaN)  %>% 
filter(!gt_30_fraction < 0.25 |gt_30_fraction ==NaN) %>% 
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




