######## FA length #################################################

# Data:2_Lipids_Imputed_CVfiltered (Output_Data)
#      normalized, imputed, CV_filtered

# Approach: 1) Fatty acid abundance in all lipids (individual fatty acid concentration matches the corresponding lipid). Odd chain fatty acids and fatty acids with d, O and t modifications are not analyzed.
#           2) Categorize in length
#           3) Fold change and stats 
#           4) Plot length


# Goal: Plot length of fatty acid sidechains in all lipids


# Output: Csv: 10_FA_length_all lipids.csv   



#Packages
library(tidyverse)
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))# Set wd to the current file location
rm(list=ls())

load(file = "../Output_Data/2_Lipids_Imputed_CVfiltered.Rdata")
final_table<-df.impute_CVfiltered


########Fatty acid abundance in all lipids #################################################

########FA count. List with individual fatty acids listed with the corresponding lipid concentration in one table

FA_table<-final_table[,c(2,6:ncol(final_table))] #make a table with all FA1

FA_temp<-final_table[,c(3,6:ncol(final_table))] ##make a table with all FA2

colnames(FA_temp)[1]<-c("FA1") #adjust colnames so that they match for bind_rows

FA_table<-bind_rows(FA_table,FA_temp) #merge

FA_temp<-final_table[,c(4,6:ncol(final_table))] #repeat step above with FA3 

colnames(FA_temp)[1]<-c("FA1")

FA_table<-bind_rows(FA_table,FA_temp)

FA_temp<-final_table[,c(5,6:ncol(final_table))] #repeat step above with FA4

colnames(FA_temp)[1]<-c("FA1")

FA_table<-bind_rows(FA_table,FA_temp)

FA_table_1<-bind_rows(FA_table,FA_temp)

colnames(FA_table_1)[1]<-"FA"


#Filter out odd chain fatty acids, O,t,d modifications 

FA_table_1<- filter(FA_table_1, !grepl("t",FA))

FA_table_1<- filter(FA_table_1, !grepl("d",FA))

FA_table_1<- filter(FA_table_1, !grepl("O",FA))

Length_FA_table<- FA_table_1 %>% 
  select(1:25) %>%
  filter(!grepl("6:0",FA)) %>%
  filter(!grepl("1:",FA)) %>%
  filter(!grepl("3:",FA)) %>%
  filter(!grepl("5:",FA)) %>%
  filter(!grepl("7:",FA)) %>%
  filter(!grepl("9:",FA))


#Table with the number of unsaturations (ether fatty acids are kept)

Length_FA_table<-Length_FA_table %>%
  mutate(FA = str_replace_all(FA, "e", "")) %>% #replace all the "e"s with spaces
  mutate(FA = sub(":.*", "", FA)) %>% #delete the info after the ":"
  select(FA,2:25)


#Aggregate

Length_FA_table<-aggregate(Length_FA_table[,-1], by = list(Length_FA_table$FA), FUN = sum) #summarize intensity of the same fatty acids
Length_FA_table<-Length_FA_table[-1,]#delete the first row with the empty FA loadings 
Length_FA_table <-as.data.frame(lapply(Length_FA_table,as.numeric)) #make numeric

rownames(Length_FA_table)<-Length_FA_table$Group.1
Length_FA_table<-Length_FA_table[,-1]

########4) Fold change and stats #################################################
#Mean, Coefficient of Variantion and Fold Change

Length_FA_table<-Length_FA_table %>% 
  mutate(ASH = apply(Length_FA_table %>% select(1:6), 1, mean ))%>% 
  mutate(EV = apply(Length_FA_table %>% select(7:12), 1, mean ))%>%
  mutate(PA = apply(Length_FA_table %>% select(13:18), 1, mean ))%>% 
  mutate(PRX = apply(Length_FA_table %>% select(19:24), 1, mean ))%>%
  mutate(cvASH = apply(Length_FA_table %>% select(1:6), 1, sd )/ASH)%>% 
  mutate(cvEV = apply(Length_FA_table %>% select(7:12), 1, sd )/EV)%>%
  mutate(cvPA = apply(Length_FA_table %>% select(13:18), 1, sd )/PA)%>% 
  mutate(cvPRX = apply(Length_FA_table %>% select(19:24), 1, sd )/PRX)%>%
  mutate(FC_EVvsASH = (ASH/EV))%>%
  mutate(FC_PRXvsPA = (PRX/PA))%>%
  mutate(FC_EVvsPRX = (PRX/EV))%>%
  mutate(FC_EVvsPA = (PA/EV))

#Wilcoxon

for (i in 1:nrow(Length_FA_table)) {
  
  Length_FA_table$pval_EVvsASH[i]<-wilcox.test(as.numeric(Length_FA_table[i,1:6]),as.numeric(Length_FA_table[i,7:12]))$p.value #to compare ASH vs EV
  Length_FA_table$pval_PRXvsPA[i]<-wilcox.test(as.numeric(Length_FA_table[i,19:24]),as.numeric(Length_FA_table[i,13:18]))$p.value #to compare PRX-ASH with PRX
  Length_FA_table$pval_EVvsPA[i]<-wilcox.test(as.numeric(Length_FA_table[i,7:12]),as.numeric(Length_FA_table[i,13:18]))$p.value #to compare EV with PRX-ASH
  Length_FA_table$pval_EVvsPRX[i]<-wilcox.test(as.numeric(Length_FA_table[i,7:12]),as.numeric(Length_FA_table[i,19:24]))$p.value #to compare EV with PRX
}

#Adjust p-val

Length_FA_table<- Length_FA_table %>%
  mutate(FDR_EVvsASH = p.adjust(pval_EVvsASH,method = "BH"))%>%
  mutate(FDR_PRXvsPA = p.adjust(pval_PRXvsPA,method = "BH")) %>%
  mutate(FDR_EVvsPA = p.adjust(pval_EVvsPA,method = "BH"))%>%
  mutate(FDR_EVvsPRX = p.adjust(pval_EVvsPRX,method = "BH"))


#Add names

Length_FA_table$Length = rownames(Length_FA_table)

Length_FA_table <- Length_FA_table %>%
  select(Length, everything()) #move length to the front



#Save
write.csv(Length_FA_table, '../Output_Data/10_FA_Length_All_Lipids.csv')

# R version 3.6.3 (2020-02-29)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Catalina 10.15.7
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggpubr_0.4.0.999  data.table_1.13.6 forcats_0.5.1     stringr_1.4.0     dplyr_1.0.5       purrr_0.3.4      
# [7] readr_1.4.0       tidyr_1.1.3       tibble_3.1.1      ggplot2_3.3.3     tidyverse_1.3.1  
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.1.0  haven_2.3.1       carData_3.0-4     colorspace_1.4-1  vctrs_0.3.7       generics_0.1.0   
# [7] utf8_1.1.4        rlang_0.4.10      pillar_1.6.0      foreign_0.8-76    glue_1.4.2        withr_2.4.2      
# [13] DBI_1.1.0         dbplyr_2.1.1      sessioninfo_1.1.1 modelr_0.1.8      readxl_1.3.1      lifecycle_1.0.0  
# [19] munsell_0.5.0     ggsignif_0.6.1    gtable_0.3.0      cellranger_1.1.0  zip_2.1.1         rvest_1.0.0      
# [25] rio_0.5.26        labeling_0.4.2    curl_4.3          fansi_0.4.1       broom_0.7.6       Rcpp_1.0.5       
# [31] backports_1.2.0   scales_1.1.1      jsonlite_1.7.2    abind_1.4-5       farver_2.0.3      fs_1.5.0         
# [37] hms_1.0.0         digest_0.6.27     openxlsx_4.2.3    stringi_1.5.3     rstatix_0.7.0     cowplot_1.1.1    
# [43] grid_3.6.3        cli_2.5.0         tools_3.6.3       magrittr_2.0.1    car_3.0-10        crayon_1.4.1     
# [49] pkgconfig_2.0.3   ellipsis_0.3.2    MASS_7.3-53       xml2_1.3.2        reprex_2.0.0      lubridate_1.7.10 
# [55] assertthat_0.2.1  httr_1.4.2        rstudioapi_0.13   R6_2.5.0          compiler_3.6.3
