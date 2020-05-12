library(tidyverse)
library(readxl)
library(wesanderson)
library(RColorBrewer)

which(Final_table$Case_ID=='BluKim')


Final_table <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\Summary_of_public_data.xlsx",
                          sheet ='NGS_Data_summary')

Mut_data <- read_excel("G:/MAC_Research_Data/Pan_cancer/Pan_cancer_mapping_result/Mutation_rate/Mutation_rate.xlsx",
                       sheet ='All_samples_Retro+indel' )

total_file <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\Data_summary.xlsx",
                         sheet ='Total')

PAIR <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\Data_summary.xlsx",
                   sheet ='PAIRS')

callable <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\Data_summary.xlsx",
                       sheet ='non-retro-mut')

Lymphoma_subtype <- read.table("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Meta\\Pancancer_metadata_05_11_2020.txt",
                               sep ='\t',header =T, stringsAsFactors =F)

match_data <- read.table("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Meta\\Tumor_normal_pair_metadata.txt",
                         sep ='\t',header =T, stringsAsFactors =F)

exclude <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\V2_Data_summary.xlsx",
                      sheet ="Total_excluded")

Lym <- Lymphoma_subtype %>%
  filter(CancerType=='Lymphoma') %>% 
  select(SampleName,Sample_id,DiseaseAcr)



Lym_final <- Final_table %>% 
  filter(Cancer_Type=='Lymphoma') %>% 
  select(Sample_ID)

Lym_final <- Lym_final$Sample_ID

# we can use which function to look at the position and we can use match
# Match will be faster match('CMT-100',non_retro$file_name)

check_nonretro_mut <- function(x){
  non_retro <- Mut_data %>% 
    filter(status=="Non-Retro")
    value <- which(non_retro$file_name==x)
  if (length(value)>0){
  return (non_retro$Mutation_rate[value])
  }
    else{
      return ("Non-Pair")
    }
}

Non_retro <- sapply(Final_table$Case_ID,check_nonretro_mut)



non_retro <- Mut_data %>% 
  filter(status=="Non-Retro")


a <- which(non_retro$file_name=='CMT-100')

match('DD0001',non_retro$file_name)

length(a)


check_subtype <- function(x){
  total_row <- nrow(Lym)
  for (i in 1:total_row){
    if (x %in% Lym[i,]){
      return (Lym$DiseaseAcr[i])
    }
  }
  return ("No-Pair")
}

check_new_subtype <- function(x){
  total_row <- nrow(Lym)
  for (i in 1:total_row){
    if (x %in% Lym[i,]){
      return (Lym$DiseaseAcr[i])
    }
  }
  return ("No-Pair")
}



Lym_subtype <- sapply(Lym_final,check_subtype )

Lym_final <- Final_table %>%
  filter(Cancer_Type=="Lymphoma") %>% 
  mutate(Cancer_Type=Lym_subtype)
other_final <- Final_table %>% 
  filter(!Cancer_Type=='Lymphoma')

Final_table <- rbind(other_final,Lym_final)

check_sef_match <- function(x){
  total_row <- nrow(match_data)
  for (i in 1:total_row){
    if (x %in% match_data[i,]){
      return (match_data$SelfMatch[i])
    }
  }
  return ("NaN")
}

check_diff_best <- function(x){
  total_row <- nrow(match_data)
  for (i in 1:total_row){
    if (x %in% match_data[i,]){
      return (match_data$Diff.from.best[i])
    }
  }
  return ("NaN")
}


exclude_check <- function(x){
  total_row <- nrow(exclude)
  for (i in 1:total_row){
    if (x %in% exclude[i,]){
      return (exclude$Exclude_Reason[i])
    }
  }
  return ("Pass QC")
}


exclude_reason <- sapply(as.vector(Final_table$Case_ID),exclude_check)
Self_Match <- sapply(as.vector(Final_table$Case_ID), check_sef_match)
diff_best <- sapply(as.vector(Final_table$Case_ID), check_diff_best)

Final_Table <- Final_table %>% 
  mutate(Self_Match=Self_Match) %>% 
  mutate(Diff.from.best=diff_best) %>% 
  mutate(Exclude_reason = exclude_reason)


Final_Table %>% 
  filter(!Case_ID %in% exclude$Cases) %>% 
  write.table("C:\\Users\\abc73_000\\Desktop\\PASS_Total_table.txt",
              sep ='\t',row.names = F,quote = F)

Final_Table %>% 
  filter(Case_ID %in% exclude$Cases) %>% 
  write.table("C:\\Users\\abc73_000\\Desktop\\Fail_Total_table.txt",
              sep ='\t',row.names = F,quote = F)
  
  


###### arrange table ######

Total_sample <-  total_file$Sample_ID

check_inside <- function(x){
  #PAIR <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\V2Supp1_Data.xlsx",
  #                     sheet ='PAIRS')
  
  total_row <- nrow(PAIR)
  for (i in 1:total_row){
    if (x %in% PAIR[i, ]){
      return (PAIR$Cases_ID[i])
    }
  }
  return ("No-Pair")
}

Total_case <- sapply(Total_sample,check_inside) 

total_file <- total_file %>% 
  mutate(Case_ID = Total_case)


check_Non_Retro_mut_inside <- function(x){
  #Mut_data <- read_excel("G:/MAC_Research_Data/Pan_cancer/Pan_cancer_mapping_result/Mutation_rate/Mutation_rate.xlsx",
  #                       sheet ='All_samples_Retro+indel' )
  
  Non_Retro_Mut_data <- Mut_data %>% 
    filter(Status == 'Non-Retro')
  
  total_Non_row <- nrow(Non_Retro_Mut_data)
  
  
  for (i in 1:total_Non_row){
    if (x %in% Non_Retro_Mut_data[i, ]){
      return (Non_Retro_Mut_data$Mutation_rate[i])
    }
  }
  return ("No-Pair")
}

check_Retro_mut_inside <- function(x){
  #Mut_data <-read_excel("G:/MAC_Research_Data/Pan_cancer/Pan_cancer_mapping_result/Mutation_rate/Mutation_rate.xlsx",
  #                      sheet ='All_samples_Retro+indel' )
  
  
  Retro_Mut_data <- Mut_data %>% 
    filter(Status == 'Retro')
  
  total_Retro_row <- nrow(Retro_Mut_data)
  
  for (i in 1:total_Retro_row){
    if (x %in% Retro_Mut_data[i, ]){
      return (Retro_Mut_data$Mutation_rate[i])
    }
  }
  return ("No-Pair")
}


Non_mutation <- sapply(as.vector(Total_case),check_Non_Retro_mut_inside)
Retro_mutation <- sapply(as.vector(Total_case),check_Retro_mut_inside)


Total_Table<- total_file %>% 
  mutate(Non_Retro_Mutation_Rate = Non_mutation) %>% 
  mutate(Retro_Mutation_Rate = Retro_mutation)


a <- Total_Table %>% 
  filter(Status=='Normal') %>% 
  mutate(Non_Retro_Mutation_Rate = 'Normal_sample') %>%
  mutate(Retro_Mutation_Rate = 'Normal_sample')

b <- Total_Table %>% 
  filter(Status=='Tumor')

Final_Total <- rbind(a,b)
Final_Total <- Final_Total[order(Final_Total$Cancer_Type), ] 

Final_Total[,c(10,1,2,3,4,5,6,7,8,9,11,12)] %>% 
  write.table("C:\\Users\\abc73_000\\Desktop\\Final_Total_table.txt",
              sep ='\t',row.names = F,quote = F)




