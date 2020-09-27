library(tidyverse)
library(readxl)
#library(wesanderson)
library(RColorBrewer)
library(data.table)


total_file <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan-Cancer-Manuscript\\Methods_legends_tables\\TableS1_8-28-20.xlsx",
                         sheet ='DiscoveryWESQCdata', skip =1)
location <- "G:\\MAC_Research_Data\\MHC_ANN\\Somatic_Germline"
#"C:\\Users\\abc73_000\\Desktop\\Bioproject_check"
#"G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Grant_table\\Cancer_Sample\\MC\\WXS\\"
file <- "CMT-100_sINDEL_VEP.vcf"
#"PRJNA489159.txt"

#data <- fread(paste(location,file,sep = "\\"))

data <- fread("C:/Users/abc73_000/Desktop/Meta/PRJNA247493_LYM_Unclasified.txt")

data[,c(Age)]

strsplit(data$`Library Name`,"_")
#1.
glm <- glm_Meta_data[,c("Sample Name","Age")]
mc <- data[,c("Sample_ID","Age")]
hsa <- data[,c("Run","Age")]

rbind(glm,mc, fill = T)



Mut_data <- read_excel("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate/Mutation_rate.xlsx",
                       sheet ='All_samples_Retro+indel' )

callable <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\Original_Data_summary.xlsx",
                       sheet ='non-retro-mut')

exclude <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\Original_Data_summary.xlsx",
                      sheet ="Total_excluded")


library(readxl)
total_breeds <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan-Cancer-Manuscript\\Methods_legends_tables\\Combined_Summary_Public_Data.xlsx",
                    sheet = "Sheet1")  

sort(unique(total_breeds$Breed)) %>% 
  write.table("C:\\Users\\abc73_000\\Desktop\\Pan_Cancer_Breed.txt",
              sep ='\t',row.names = F,quote = F)



breeds <- fread("C:\\Users\\abc73_000\\Desktop\\Breeds_info.txt")
breeds_info <- breeds[, c("Sample_id","SampleName","Breed","QC_result","BreedCluster")]


## This function is used to identify the case name
## the input file needs case normal-tsumor file
check_sample_name <- function(x){
  total_row <- nrow(sample_pair)
  for (i in 1:total_row){
    if (x %in% sample_pair[i,]){
      return (sample_pair$Cases[i])
    }
  }
  return ("NaN")
}

# we can use which function to look at the position and we can use match
# Match will be faster , but
# eg. match('CMT-100',non_retro$file_name)

check_origin_breeds <- function(x){
    value <- match(x, breeds$Sample_id,nomatch=0)
  if (value!=0){
  return (breeds[value, Breed])
  }
    else{
      return ("NaN")
    }
}

check_breeds_QC <- function(x){
  value <- match(x, breeds$Sample_id,nomatch=0)
  if (value!=0){
    return (breeds[value, BreedQC])
  }
  else{
    return ("NaN")
  }
}


check_breed_cluster <- function(x){
  value <- match(x, breeds$Sample_id,nomatch=0)
  if (value!=0){
    return (breeds[value, BreedCluster])
  }
  else{
    return ("NaN")
  }
}

check_breed_table <- function(x){
  value <- match(x, breeds_info$Sample_id,nomatch=0)
  if (value!=0){
    return (breeds_info[value, ])
  }
  else{
    return ("NaN")
  }
}


origin_breeds <- sapply(total_file$Sample_ID,check_CaseID)
breeds_cluster <- sapply(total_file$Sample_ID,check_breed_cluster)
breeds_QC <- sapply(total_file$Sample_ID,check_breeds_QC)
total_file$Breeds= origin_breeds
total_file$Breeds_cluster <- breeds_cluster
total_file$Breeds_QC <- breeds_QC

test <- lapply(total_file$Sample_ID,check_breed_table)

name = names(test[[1000]])

data.frame(matrix(unlist(test),nrow=length(test),byrow=TRUE))

a = data.frame(matrix(unlist(test),nrow=length(test),byrow=TRUE))

b <- matrix(unlist(test), nrow = length(test), byrow = T)


length(unlist(test))

a = as.data.frame(test, col.names = name)

mylist = do.call(rbind,test)

write.table(total_file,"C:\\Users\\abc73_000\\Desktop\\final_table.txt",
            quote = F,sep = '\t',col.names = T, row.names = F)




a <- which(non_retro$file_name=='CMT-100')

#match('DD0001',total_file$Case_ID)

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




