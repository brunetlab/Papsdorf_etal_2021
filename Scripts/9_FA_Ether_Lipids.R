######## Individual fatty acids in Ether lipids #################################################

# Data:2_Lipids_Imputed_CVfiltered (Output_Data)
#      normalized, imputed, CV_filtered

# Approach: 1) Fatty acid abundance in ether lipids (individual fatty acid concentration matches the corresponding lipid). Odd chain fatty acids and fatty acids with d, O and t modifications are not analyzed.
#           2) Normalize by abundance to get relative fatty acid composition
#           3) Fold change and stats 
#           4) Plot individual fatty acids with an intensity > 0.0045

# Goal: Plot individual fatty acids

# Output:  Pdf: 9_FA_Ether_Lipids.pdf
#          Csv: 9_FA_Ether_Lipids.csv 

#Packages
library(tidyverse)
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))# Set wd to the current file location
rm(list=ls())

load(file = "../Output_Data/2_Lipids_Imputed_CVfiltered.Rdata")
final_table<-df.impute_CVfiltered


######## Filter ether-lipids

final_table_e<-final_table%>%filter(grepl("e",FA1)) #44 lipids with ether fatty acids identified


######### Normalize for composition. Silence the next two lines if absolute concentration should be analyzed

total<-as.list(apply(final_table_e[,6:ncol(final_table_e)],2,sum))

final_table_e[,6:ncol(final_table_e)]<-final_table_e[,6:ncol(final_table_e)]/total



########Fatty acid abundance in all lipids #################################################

########FA count. List with individual fatty acids listed with the corresponding lipid concentration in one table

FA_table<-final_table_e[,c(2,6:ncol(final_table_e))] #make a table with all FA1

FA_temp<-final_table_e[,c(3,6:ncol(final_table_e))] ##make a table with all FA2

colnames(FA_temp)[1]<-c("FA1") #adjust colnames so that they match for bind_rows

FA_table<-bind_rows(FA_table,FA_temp) #merge

FA_temp<-final_table_e[,c(4,6:ncol(final_table_e))] #repeat step above with FA3 

colnames(FA_temp)[1]<-c("FA1")

FA_table<-bind_rows(FA_table,FA_temp)

FA_temp<-final_table_e[,c(5,6:ncol(final_table_e))] #repeat step above with FA4

colnames(FA_temp)[1]<-c("FA1")

FA_table<-bind_rows(FA_table,FA_temp)

FA_table_1<-bind_rows(FA_table,FA_temp)

colnames(FA_table_1)[1]<-"FA"


#Filter out odd chain fatty acids, O,t,d modifications 

FA_table_1<- filter(FA_table_1, !grepl("t",FA))

FA_table_1<- filter(FA_table_1, !grepl("d",FA))

FA_table_1<- filter(FA_table_1, !grepl("O",FA))


FA_table_1<- FA_table_1 %>% 
  select(1:25) %>%
  filter(!grepl("6:0",FA)) %>%
  filter(!grepl("1:",FA)) %>%
  filter(!grepl("3:",FA)) %>%
  filter(!grepl("5:",FA)) %>%
  filter(!grepl("7:",FA)) %>%
  filter(!grepl("9:",FA)) 


#Aggregate for fatty acid plot

FA_table_2<-aggregate(FA_table_1[,-1], by = list(FA_table_1$FA), FUN = sum) #summarize intensity of the same fatty acids

FA_table_2<-FA_table_2[-1,]#delete the first row with the empty FA loadings 

FA_table_2 [,2:25] <-as.data.frame(lapply(FA_table_2[,2:25],as.numeric)) #make numeric




########4) Fold change and stats #################################################
#Mean, Coefficient of Variantion and Fold Change

FA_table_2<-FA_table_2 %>% 
  mutate(ASH = apply(FA_table_2 %>% select(2:7), 1, mean ))%>% 
  mutate(EV = apply(FA_table_2 %>% select(8:13), 1, mean ))%>%
  mutate(PA = apply(FA_table_2 %>% select(14:19), 1, mean ))%>% 
  mutate(PRX = apply(FA_table_2 %>% select(20:25), 1, mean ))%>%
  mutate(cvASH = apply(FA_table_2 %>% select(2:7), 1, sd )/ASH)%>% 
  mutate(cvEV = apply(FA_table_2 %>% select(8:13), 1, sd )/EV)%>%
  mutate(cvPA = apply(FA_table_2 %>% select(14:19), 1, sd )/PA)%>% 
  mutate(cvPRX = apply(FA_table_2 %>% select(20:25), 1, sd )/PRX)%>%
  mutate(FC_EVvsASH = (ASH/EV))%>%
  mutate(FC_PRXvsPA = (PRX/PA))

#Wilcoxon

for (i in 1:nrow(FA_table_2)) {
  
  FA_table_2$pval_EVvsASH[i]<-wilcox.test(as.numeric(FA_table_2[i,2:7]),as.numeric(FA_table_2[i,8:13]))$p.value #to compare ASH vs EV
  FA_table_2$pval_PRXvsPA[i]<-wilcox.test(as.numeric(FA_table_2[i,20:25]),as.numeric(FA_table_2[i,14:19]))$p.value #to compare PRX-ASH with PRX
  FA_table_2$pval_EVvsPA[i]<-wilcox.test(as.numeric(FA_table_2[i,8:13]),as.numeric(FA_table_2[i,14:19]))$p.value #to compare EV with PRX-ASH
  FA_table_2$pval_EVvsPRX[i]<-wilcox.test(as.numeric(FA_table_2[i,8:13]),as.numeric(FA_table_2[i,20:25]))$p.value #to compare EV with PRX
}

rm(i)

#Adjust p-val

FA_table_2<- FA_table_2 %>%
  mutate(FDR_EVvsASH = p.adjust(pval_EVvsASH,method = "BH"))%>%
  mutate(FDR_PRXvsPA = p.adjust(pval_PRXvsPA,method = "BH")) %>%
  mutate(FDR_EVvsPA = p.adjust(pval_EVvsPA,method = "BH"))%>%
  mutate(FDR_EVvsPRX = p.adjust(pval_EVvsPRX,method = "BH"))


#Add names
rownames(FA_table_2)<-FA_table_2$Group.1

#Sort
uns_FA_table<-FA_table_2
uns_FA_table$Uns<-row.names(FA_table_2) #make row with unsaturations

uns_FA_table_2<-uns_FA_table %>%
  mutate(Uns = sub(".*:", "", Uns))%>% #make unsaturations one number
  mutate(Uns = sub("e", "", Uns))  #delete "e"

uns_FA_table_2$Uns<-as.numeric(uns_FA_table_2$Uns) #numeric values in uns column
uns_FA_table_2<-uns_FA_table_2[order(uns_FA_table_2$Uns),] #order

FA_table_2<-uns_FA_table_2[,-44]#delete uns column



#Save
write.csv(FA_table_2, '../Output_Data/9_FA_Ether_Lipids.csv')


########Organize data for plot

tFA_table_2<-as.data.frame(t(FA_table_2[,2:25])) %>%
  mutate(Group = NA)

tFA_table_2[c(1:6),ncol(tFA_table_2)]<- "ASH"

tFA_table_2[c(7:12),ncol(tFA_table_2)]<- "EV"

tFA_table_2[c(13:18),ncol(tFA_table_2)]<- "PRX_ASH"

tFA_table_2[c(19:24),ncol(tFA_table_2)]<- "PRX"

tFA_table_2$Group <- factor(tFA_table_2$Group , levels=c("EV", "ASH", "PRX", "PRX_ASH")) #order the groups

ptSMP_table<-pivot_longer(tFA_table_2, cols = c(1:(ncol(tFA_table_2)-1)), names_to = "Var", values_to = "Val")

ptSMP_table$Var<- factor(ptSMP_table$Var , levels=c("12:0e", "18:0", "18:0e", "12:1e", "14:1", "16:1", "18:1", "18:1e", "18:2","18:3e", "20:2", "20:2e", "20:3", "20:3e", "20:4", "20:5")) #order by unsaturation

ptSMP_table<-ptSMP_table%>%filter(Val>0.0045) #only keep values above 0.0045 to have show higher abundant FAs

######## Plot

plot<-ggplot(ptSMP_table,aes(x = Var, y = Val, fill = Group)) +
  geom_boxplot(outlier.alpha = 0.99, outlier.size=1.5)+
  xlab("")+
  ylab("Relative fatty acid composition in Ether lipids")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 50, vjust = 0.5), plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size=15))+
  theme(legend.position = "none")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_fill_manual(values = c("EV" = "grey55",
                               "ASH" = "red1",
                               "PRX" = "blue",
                               "PRX_ASH" = "cyan2"))

plot(plot)

ggsave("../Output_Figures/9_FA_Ether_Lipids.pdf", plot, width=4.5, height=2, units="in", scale=3) #change filename 

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
