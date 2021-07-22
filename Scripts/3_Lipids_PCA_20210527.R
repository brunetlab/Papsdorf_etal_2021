######## PCA #################################################

# Data:2_Lipids_Imputed_CVfiltered (Output_Data)

# Approach: 1) PCA with log scaled lipids
#           2) Identify top 50 lipids that drive PCA in each direction. Plots are generated to better visualize. 
#              The highest lipids are considered the drivers and listed in the PCA Fig4c. 
#              Fatty acids are analyed to classify ether lipids. Otherwise fatty acid composition is not further evaluated.

# Output:     Pdf: 3_PCA.pdf
#             Csv: 3_PCA.csv
#                  3_Lipids_PC_loading.csv


#Packages
library(tidyverse)
library(data.table)
library(dplyr)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))# Set wd to the current file location
rm(list=ls())

load(file = "../Output_Data/2_Lipids_Imputed_CVfiltered.Rdata")
final_table<-df.impute_CVfiltered

################### 1)PCA         ############################################

#Logscale data and merge name
pca_table<-final_table%>%select(6:ncol(final_table)) #select rows that contain lipid values
pca_table <-  as.data.frame(log2(pca_table)) #logscale 
row.names(pca_table)<-paste0(final_table[,1],"-",final_table[,2],"-",final_table[,3],"-",final_table[,4],"-",final_table[,5]) #merges lipid features into one column and set as rownames

# Calculate PCs
data.pca <- prcomp(t(pca_table),center=T, scale=T)
#summary(data.pca)
sum.pca<-summary(data.pca)
imp.pca<-sum.pca$importance;
std.pca<-imp.pca[1,] # standard deviation
var.pca<-imp.pca[2,] # variance explained by each PC
cum.pca<-imp.pca[3,] # cummulated variance explained


# PCA plot PC1 and PC2

Comp <- c(1,2) # Choose PC to plot
ggdata <- data.frame(data.pca$x[,Comp[1]],data.pca$x[,Comp[2]])

# Column category to apply a color scheme. 
ggdata$Category[c(1:6)] <- "ASH"
ggdata$Category[c(7:12)] <- "EV"
ggdata$Category[c(13:18)] <- "PRX_ASH"
ggdata$Category[c(19:24)] <- "PRX"

#Order Category
ggdata$Category <- factor(ggdata$Category , levels=c("EV", "ASH", "PRX", "PRX_ASH")) 

colnames(ggdata) <- c(paste("PC", Comp, sep=""),"Category")

################### Plot

plot<-ggplot(ggdata, aes(x=PC1,y=PC2,fill=factor(Category), label=rownames(ggdata))) +
  stat_ellipse(geom="polygon", level=0.95, alpha=0.25) +
  geom_point(size=7, shape=21, color="black", alpha = 0.85) +
  scale_fill_manual(values = c("EV" = "grey55",
                               "ASH" = "red1",
                               "PRX" = "blue",
                               "PRX_ASH" = "cyan2")) +
  #geom_text(aes(label = rownames(ggdata)) ,hjust=0, vjust=0) +
  guides(color=guide_legend("Category"),fill=guide_legend("Category")) +
  scale_x_continuous(paste("PC1 (", round(100*var.pca[Comp[1]],1),"%)",sep="")) +
  scale_y_continuous(paste("PC2 (", round(100*var.pca[Comp[2]],1),"%)",sep="")) +
  theme_bw() + theme(legend.justification=c(1,1), legend.position=c(1,1)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.title = element_text(face="bold", size=16)) + theme(legend.title = element_text(colour="white"),legend.position= c(0.5,0.95), plot.title = element_text(hjust = 0.5, face="bold", size=16))+
  ggtitle("PCA_Median quantitation")


plot(plot)

pdf("../Output_Figures/3_PCA.pdf")
print(plot)
dev.off() 


################### 2) Identify lipids that drive PCA#################################

Lipids_PCA<-data.pca$rotation #Loadings for lipids
Lipids_PCA<-as.data.frame(Lipids_PCA)

write.csv(Lipids_PCA, '../Output_Data/3_Lipids_PC_loading.csv')

#              Plots are generated to better visualize. 
#              The highest lipids are considered the drivers and listed in the PCA Fig4c. 
#              Fatty acids are analyzed to classify ether lipids. Otherwise fatty acid composition is not further evaluated.

###############PC1

Lipids_PC1<-Lipids_PCA[order(Lipids_PCA$PC1),]        #this orders them by PC1 value
Lipids_PC1<-setDT(Lipids_PC1, keep.rownames = TRUE)   #convert to rownames to separate column

Lipids_PC1<-Lipids_PC1%>%
          separate(rn, c("Class", "FA1", "FA2", "FA3", "FA4"), sep = "([-])") #separate names

#Bottom PC1 Class = TG, PE
Bottom_PC1 <- as.data.frame(Lipids_PC1[1:50,]) 
plot(as.factor(Bottom_PC1[,1])) 

#Bottom PC1 FA = ether-fatty acids and 18:1
Bottom_PC1_FA <- Bottom_PC1[,2:5] #select only FA columns
Bottom_PC1_FA <-Bottom_PC1_FA  %>% gather(key = FA, value = Value, FA1:FA4) #gather in one column
plot(as.factor(Bottom_PC1_FA[,2]), las=2) 

#Ether lipids have a ether bound fatty acid only in the sn1(FA1) position. 38% of all lipids
plot(as.factor(Bottom_PC1[,2]), las=2) # an ether-lipid has an ether-bound fatty acid in the sn1 position (FA1). Only FA1!
x=100/nrow(Bottom_PC1)*nrow(Bottom_PC1%>%filter(grepl("e", FA1))) 

#Top PC1 Class = PE
Top_PC1 <- as.data.frame(Lipids_PC1[448:498,]) 
plot(as.factor(Top_PC1[,1]))

#Top PC1 FA = diverse. Not specified in Figure.
Top_PC1_FA <- Top_PC1[,2:5]
Top_PC1_FA <-Top_PC1_FA  %>% gather(key = FA, value = Value, FA1:FA4)
plot(as.factor(Top_PC1_FA[,2]), las=2)  



###############PC2

Lipids_PC2<-Lipids_PCA[order(Lipids_PCA$PC2),]        #this orders them by PC2 value
Lipids_PC2<-setDT(Lipids_PC2, keep.rownames = TRUE)   #convert to rownames to separate column
Lipids_PC2<-Lipids_PC2%>%
  separate(rn, c("Class", "FA1", "FA2", "FA3", "FA4"), sep = "([-])")

#Bottom PC2 Class = TG
Bottom_PC2 <- as.data.frame(Lipids_PC2[1:50,]) 
plot(as.factor(Bottom_PC2[,1])) 

#Bottom PCA2 FA =  20:5, 17:1. 
Bottom_PC2_FA <- Bottom_PC2[,2:5] #select only FA columns
Bottom_PC2_FA <-Bottom_PC2_FA  %>% gather(key = FA, value = Value, FA1:FA4) #gather them all in one column
plot(as.factor(Bottom_PC2_FA[,2]), las=2)  

#Top PC2 Class = CL, PS
Top_PC2 <- as.data.frame(Lipids_PC2[448:498,]) 
plot(as.factor(Top_PC2[,1]))

#Top PCA2 FA = diverse but mostly highly unsaturated. Not specified in Figure.
Top_PC2_FA <- Top_PC2[,2:5] #select only FA columns
Top_PC2_FA <-Top_PC2_FA  %>% gather(key = FA, value = Value, FA1:FA4) #gather them all in one column
plot(as.factor(Top_PC2_FA[,2]), las=2)  



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
#   [1] data.table_1.13.6 forcats_0.5.1     stringr_1.4.0     dplyr_1.0.5       purrr_0.3.4       readr_1.4.0      
# [7] tidyr_1.1.3       tibble_3.1.1      ggplot2_3.3.3     tidyverse_1.3.1  
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.5       cellranger_1.1.0 pillar_1.6.0     compiler_3.6.3   dbplyr_2.1.1     tools_3.6.3      digest_0.6.27   
# [8] lubridate_1.7.10 jsonlite_1.7.2   lifecycle_1.0.0  gtable_0.3.0     pkgconfig_2.0.3  rlang_0.4.10     reprex_2.0.0    
# [15] DBI_1.1.0        cli_2.5.0        rstudioapi_0.13  haven_2.3.1      xml2_1.3.2       withr_2.4.2      httr_1.4.2      
# [22] fs_1.5.0         generics_0.1.0   vctrs_0.3.7      hms_1.0.0        grid_3.6.3       tidyselect_1.1.0 glue_1.4.2      
# [29] R6_2.5.0         fansi_0.4.1      readxl_1.3.1     farver_2.0.3     modelr_0.1.8     magrittr_2.0.1   MASS_7.3-53     
# [36] scales_1.1.1     backports_1.2.0  ellipsis_0.3.2   rvest_1.0.0      assertthat_0.2.1 colorspace_1.4-1 labeling_0.4.2  
# [43] utf8_1.1.4       stringi_1.5.3    munsell_0.5.0    broom_0.7.6      crayon_1.4.1  



