######## Ether lipid abundance #################################################

# Data:2_Lipids_Imputed_CVfiltered (Output_Data)
#      Lipids are normalized, imputed, CV_filtered

# Approach: 1) Filter ether lipids and analyze class abundance for ether lipids
#           2) Fold change and stats 
#           3) Plot all ether lipids together
#           4) Plot PE, TG and PS together. Cut off value is class abundance >250


# Goal: Plot ether lipid abundance changes across all lipids and in TG, PE


# Output:  Pdf: 7_Ether_Lipids.pdf
#               7_Ether_Lipids_PE_TG_PS
#          Csv: 7_Ether_Lipids.csv   

#Packages
library(tidyverse)
library(ggpubr)

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))# Set wd to the current file location

load(file = "../Output_Data/2_Lipids_Imputed_CVfiltered.Rdata")
final_table<-df.impute_CVfiltered


######## 1) Class abundance Ether Lipids  #################################################

########Filter ether-lipids

final_table_e<-final_table%>%filter(grepl("e",FA1)) #44 lipids with ether fatty acids identified

class_table<-final_table_e%>%
  select(c(1,6:ncol(final_table_e))) #ignore fatty acid composition in columns

class_table<-aggregate(class_table[,-1], by = list(class_table$Class), FUN = sum) #aggregate all the classes
sum<-colSums(class_table [,2:25], na.rm= F) #summarize all ether-lipid species
class_table[nrow(class_table) + 1,2:25] = c(sum) #add to class table
class_table[6,1]<-"All" #rename



######## 2) Fold change and stats #################################################

#Mean, Coefficient of Variantion and Fold Change

class_table<-class_table %>% 
  mutate(ASH = apply(class_table %>% select(2:7), 1, mean ))%>% 
  mutate(EV = apply(class_table %>% select(8:13), 1, mean ))%>%
  mutate(PA = apply(class_table %>% select(14:19), 1, mean ))%>% 
  mutate(PRX = apply(class_table %>% select(20:25), 1, mean ))%>%
  mutate(cvASH = apply(class_table %>% select(2:7), 1, sd )/ASH)%>% 
  mutate(cvEV = apply(class_table %>% select(8:13), 1, sd )/EV)%>%
  mutate(cvPA = apply(class_table %>% select(14:19), 1, sd )/PA)%>% 
  mutate(cvPRX = apply(class_table %>% select(20:25), 1, sd )/PRX)%>%
  mutate(FC_ASHvsEV = (ASH/EV))%>%
  mutate(FC_PRXvsEV = (PRX/EV))%>%
  mutate(FC_PAvsEV = (PA/EV))%>%
  mutate(FC_PRXvsPA = (PRX/PA))


#Wilcoxon

for (i in 1:nrow(class_table)) {
  
  class_table$pval_ASHvsEV[i]<-wilcox.test(as.numeric(class_table[i,2:7]),as.numeric(class_table[i,8:13]))$p.value #comparing ASH vs EV
  class_table$pval_PAvsEV[i]<-wilcox.test(as.numeric(class_table[i,14:19]),as.numeric(class_table[i,8:13]))$p.value #comparing PRX-5_ASH vs EV
  class_table$pval_PRXvsEV[i]<-wilcox.test(as.numeric(class_table[i,20:25]),as.numeric(class_table[i,8:13]))$p.value #comparing PRX-5 vs EV
  class_table$pval_PRXvsPA[i]<-wilcox.test(as.numeric(class_table[i,20:25]),as.numeric(class_table[i,14:19]))$p.value #comparing PRX-5 vs PRX_ASH
}

rm(i)

#Adjust p-val

class_table<- class_table %>%
  mutate(FDR_ASHvsEV = p.adjust(pval_ASHvsEV,method = "BH"))%>%
  mutate(FDR_PAvsEV = p.adjust(pval_PAvsEV,method = "BH"))%>%
  mutate(FDR_PRXvsEV = p.adjust(pval_PRXvsEV,method = "BH"))%>%
  mutate(FDR_PRXvsPA = p.adjust(pval_PRXvsPA,method = "BH"))


#Add class names

colnames(class_table)[1]<-"Class"
rownames(class_table)<-class_table$Class

#Save
write.csv(class_table, '../Output_Data/7_Ether_Lipids.csv')


######## Organize data for plot

tclass_table<-as.data.frame(t(class_table[,2:25])) %>%
  mutate(Group = NA)

tclass_table[c(1:6),ncol(tclass_table)]<- "ASH"

tclass_table[c(7:12),ncol(tclass_table)]<- "EV"

tclass_table[c(13:18),ncol(tclass_table)]<- "PRX_ASH"

tclass_table[c(19:24),ncol(tclass_table)]<- "PRX"

tclass_table$Group <- factor(tclass_table$Group , levels=c("EV", "ASH", "PRX", "PRX_ASH")) #order the groups

ptclass_table<-pivot_longer(tclass_table, cols = c(1:(ncol(tclass_table)-1)), names_to = "Var", values_to = "Val")



######## 3) Plot all ether lipids together     #################################################

########Plot all 

plot<- filter(ptclass_table, grepl("All",Var)) 

Classes<-ggplot(plot,aes(x = Group, y = Val, fill = Group)) +
  geom_boxplot(outlier.alpha = 0.99, outlier.size=1.5)+
  xlab("")+
  ylab("Ether lipid abundance (ng/mg_prot)")+  
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

plot(Classes)

########Save Documents

ggsave("../Output_Figures/7_Ether_Lipids.pdf", Classes, width=0.8, height=2, units="in", scale=3) 



######## 4) Plot PE, TG and PS  together     #################################################


########Plot PE

PE<- filter(ptclass_table, grepl("PE",Var))
PE<- filter(PE, !grepl("LPE",Var)) #unselect lyso species

PE<-ggplot(PE,aes(x = Group, y = Val, fill = Group)) +
  geom_boxplot(outlier.alpha = 0.99, outlier.size=1.5)+
  xlab("")+
  ylab("Ether lipid PE abundance (ng/mg_prot)")+  
  ylim(0,23000)+
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


########Plot TG

TG<- filter(ptclass_table, grepl("TG",Var)) 

TG<-ggplot(TG,aes(x = Group, y = Val, fill = Group)) +
  geom_boxplot(outlier.alpha = 0.99, outlier.size=1.5)+
  xlab("")+
  ylab("Ether lipid TG abundance  (ng/mg_prot)")+ 
  ylim(0,23000)+
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


########Plot PS

PS<- filter(ptclass_table, grepl("PS",Var)) 

PS<-ggplot(PS,aes(x = Group, y = Val, fill = Group)) +
  geom_boxplot(outlier.alpha = 0.99, outlier.size=1.5)+
  xlab("")+
  ylab("Ether lipid PS abundance  (ng/mg_prot)")+ 
  ylim(0,23000)+
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

######## Plot all together  

PE_TG_PS<-ggarrange(PE, PS, TG, 
               ncol = 3, nrow = 1)

plot(PE_TG_PS)


########Save Documents

ggsave("../Output_Figures/7_Ether_Lipids_PE_PS_TG.pdf", PE_TG_PS, width=2.4, height=2, units="in", scale=3) 

# sessionInfo(package=NULL)
# R version 3.6.1 (2019-07-05)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14.6
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
#   [1] ggpubr_0.4.0    forcats_0.5.1   stringr_1.4.0   dplyr_1.0.6     purrr_0.3.4     readr_1.4.0     tidyr_1.1.3     tibble_3.1.2    ggplot2_3.3.3   tidyverse_1.3.1
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.1.1  haven_2.4.1       carData_3.0-4     colorspace_2.0-1  vctrs_0.3.8       generics_0.1.0    utf8_1.2.1        rlang_0.4.11      pillar_1.6.1     
# [10] foreign_0.8-71    glue_1.4.2        withr_2.4.2       DBI_1.1.1         dbplyr_2.1.1      modelr_0.1.8      readxl_1.3.1      lifecycle_1.0.0   munsell_0.5.0    
# [19] ggsignif_0.6.1    gtable_0.3.0      cellranger_1.1.0  zip_2.1.1         rvest_1.0.0       labeling_0.4.2    rio_0.5.26        curl_4.3.1        fansi_0.4.2      
# [28] broom_0.7.6       Rcpp_1.0.6        backports_1.2.1   scales_1.1.1      jsonlite_1.7.2    abind_1.4-5       farver_2.1.0      fs_1.5.0          digest_0.6.27    
# [37] hms_1.1.0         stringi_1.6.2     openxlsx_4.2.3    rstatix_0.7.0     cowplot_1.1.1     grid_3.6.1        cli_2.5.0         tools_3.6.1       magrittr_2.0.1   
# [46] crayon_1.4.1      car_3.0-10        pkgconfig_2.0.3   ellipsis_0.3.2    data.table_1.13.6 xml2_1.3.2        reprex_2.0.0      lubridate_1.7.10  assertthat_0.2.1 
# [55] httr_1.4.2        rstudioapi_0.13   R6_2.5.0          compiler_3.6.1 



