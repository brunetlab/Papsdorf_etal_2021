######## Lipid droplet number and peroxisome number screen statistics #################################################

# Data: (1) Normalized lipid droplet and peroxisome numbers that are used to generate the Manhattan plot Fig 3i-l. 
#           R is used to calculate exact p-values and p-adjust 

#Approach: 1) Perform lipid droplet number statistics
#          2) Perform peroxisome number statistics
#          3) Merge and save

# Goal: Apply Wilcoxon rank test and Benjamini Hochberg to correct for multiple hypotheses

# Output: Csv: 13_Lipid_droplet_Peroxisome_Screen_Stats.csv 



rm(list=ls())

options(stringsAsFactors = FALSE)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))# Set wd to the current file locati



########1) Lipid droplet number statistics #################################################

dfLD<- read.csv ('../Raw_data/Normalized LD numbers.csv', stringsAsFactors = FALSE)
dfLD<- as.data.frame(dfLD)


######## Calculate exact p-values for lipid droplet number of RNAi conditions compared to control 

x=as.numeric (dfLD[,"Control"]) #to compare each condition to the control 

p_val_LD<- NULL

for (i in 1:ncol(dfLD))
{
  p_value<-wilcox.test(x,dfLD[,i])[["p.value"]] 
  {
    p_val_LD[[length(p_val_LD) + 1]] <- p_value
  }
}

p_val_LD<-as.data.frame(p_val_LD)

row.names(p_val_LD)<-colnames(dfLD)

rm(i)


########## Calculate p_adjust for lipid droplet number

p_adj_LD<- NULL

for (i in 1:ncol(p_val_LD))
{
  p_adj<-p.adjust(p_val_LD$p_val_LD,method = "BH")
  {
    p_adj_LD[[length(p_adj_LD) + 1]] <- p_adj
  }
}

rm(i)

p_adj_LD<-as.data.frame(p_adj_LD)

row.names(p_adj_LD)<-colnames(dfLD)

colnames(p_adj_LD)<-"p_adj_LD"


######### Calculate fold change for lipid droplet number of RNAi conditions compared to control 

FC_LD = colMeans(dfLD, na.rm = T)/ mean(dfLD[,1], na.rm = T)

FC_LD<-as.data.frame(FC_LD)


######### Bind tables together

LD_stats <-cbind(p_val_LD, p_adj_LD, FC_LD)




######## 2) Peroxisome number statistics #################################################

dfP<- read.csv ('../Raw_data/Normalized P numbers.csv', stringsAsFactors = FALSE)

dfP<- as.data.frame(dfP)


########## Calculate exact p-values for peroxisome number of RNAi conditions compared to control 

x=as.numeric (dfP[,"Control"]) #to compare each condition to the control 

p_val_P<- NULL

for (i in 1:ncol(dfP))
{
  p_value<-wilcox.test(x,dfP[,i])[["p.value"]] 
  {
    p_val_P[[length(p_val_P) + 1]] <- p_value
  }
}

rm(i)

p_val_P<-as.data.frame(p_val_P)

row.names(p_val_P)<-colnames(dfP)


######### Calculate p_adjust for peroxisome number

p_adj_P<- NULL

for (i in 1:ncol(p_val_P))
{
  p_adj<-p.adjust(p_val_P$p_val_P,method = "BH")
  {
    p_adj_P[[length(p_adj_P) + 1]] <- p_adj
  }
}

rm(i)

p_adj_P<-as.data.frame(p_adj_P)

row.names(p_adj_P)<-colnames(dfP)

colnames(p_adj_P)<-"p_adj_P"


######### Calculate fold change for peroxisome number of RNAi conditions compared to control 

FC_P = colMeans(dfP, na.rm = T)/ mean(dfP[,1], na.rm = T)
          
FC_P<-as.data.frame(FC_P)


####### Bind tables tables together

P_stats <-cbind(p_val_P, p_adj_P, FC_P)



######## 3) Sort and merge data #################################################
LD_stats [nrow(LD_stats ) + 1,] = c("1","1", "NA") #add one row with values for dhs-3
row.names(LD_stats)[53] <- 'dhs.3'

#P_stats<-P_stats[row.names(P_stats) != "dhs.3", , drop = FALSE] #delete dhs-3 because it does not occur in both datasets

LD_stats <- LD_stats[ order(row.names(LD_stats)), ]

P_stats <- P_stats[ order(row.names(P_stats)), ]


#double check if rownames are the same

identical(rownames(LD_stats), rownames(P_stats)) #true


######## 3) Merge and save
Combined_Screen <-as.data.frame(cbind(LD_stats, P_stats)) #make one table

write.csv (Combined_Screen, "../Output_Data/13_Lipid_droplet_Peroxisome_Screen_Stats.csv")



# sessionInfo(package=NULL)
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
#   [1] forcats_0.5.1   stringr_1.4.0   dplyr_1.0.5     purrr_0.3.4     readr_1.4.0     tidyr_1.1.3     tibble_3.1.1    ggplot2_3.3.3  
# [9] tidyverse_1.3.1 pheatmap_1.0.12
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.5         cellranger_1.1.0   pillar_1.6.0       compiler_3.6.3     RColorBrewer_1.1-2 dbplyr_2.1.1       tools_3.6.3       
# [8] lubridate_1.7.10   jsonlite_1.7.2     lifecycle_1.0.0    gtable_0.3.0       pkgconfig_2.0.3    rlang_0.4.10       reprex_2.0.0      
# [15] cli_2.5.0          DBI_1.1.0          rstudioapi_0.13    haven_2.3.1        xml2_1.3.2         withr_2.4.2        httr_1.4.2        
# [22] fs_1.5.0           generics_0.1.0     vctrs_0.3.7        hms_1.0.0          grid_3.6.3         tidyselect_1.1.0   glue_1.4.2        
# [29] R6_2.5.0           fansi_0.4.1        readxl_1.3.1       modelr_0.1.8       magrittr_2.0.1     scales_1.1.1       backports_1.2.0   
# [36] ellipsis_0.3.2     rvest_1.0.0        assertthat_0.2.1   colorspace_1.4-1   utf8_1.1.4         stringi_1.5.3      munsell_0.5.0     
# [43] broom_0.7.6        crayon_1.4.1

