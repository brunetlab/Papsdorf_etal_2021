######## Shared GO terms between Greer et al. Nature 2010 (EG) and Han et al. Nature 2017 (SH)  #################################################

# Data:   (1) Greer et al. Nature 2010 whole worm microarray of adult day 6 and (2) Han et al. Nature 2017 worm intestine of adult day 1 upon ash-2 RNAi.
#         Genes with a log2FC bigger than 1 and a p-adjust < 0.05 were used to input in EnrichR
#         All GO terms from EnrichR were used as an input with a combined score > 4
 
#Approach: 1) Dataframe with shared GO terms. Save shared GO-terms
#          2) Manually annotate shared GO-terms are loaded (this is to break up the GO-terms into larger categories)
#          3) Make heatmap

# Goal: Plot shared GO-terms in a heatmap

# Output:  Pdf: 12_Shared_Go_term_heatmap.pdf
#          Csv: Shared_GO_terms.csv.csv 

#Of note: In the final paper figure the individual GO-terms were not plotted

rm(list=ls())
library(dplyr)
library(pheatmap)
library(RColorBrewer)


options(stringsAsFactors = FALSE)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))# Set wd to the current file location


########Organizing and Cleaning Data  #################################################

# Loading data
GOterms_SH <- read.delim("../Raw_data/GO_terms_SH.txt", header= T)

# Making sure these are factors
GOterms_SH$Term_SH = as.factor(GOterms_SH$Term_SH)

# Loading data
GOterms_EG <- read.delim('../Raw_data/GO_terms_EG.txt', header= T)

# Making sure these are factors
GOterms_EG$Term_EG = as.factor(GOterms_EG$Term_EG)

########Merging two tables  #################################################

x<- GOterms_EG
y<- GOterms_SH

########1) Dataframe with shared GO terms  #################################################

#make a list with the shared GO terms and the corresponding combined scores for SH and EG
df1 <- merge(x, y, by = intersect(names(x), names(y)),
                          by.x = "Term_EG", by.y = "Term_SH", all = FALSE,
                          sort = TRUE, suffixes = c(".x",".y"), no.dups = TRUE,
                          incomparables = NULL)

#rename first column
colnames (df1) <- c("Terms", "Combined_score_EG","Combined_score_SH")

#delete duplicates. Duplicates are deleted based on a higher combined score in the Han et al. Nature 2017 (SH) dataset.
df1 = df1[order(df1[,'Terms'],-df1[,'Combined_score_SH']),]
df2 = df1[!duplicated(df1$Term),]

#make "terms" to rownames
df3 <- df2[,-1]
rownames(df3) <- df2[,1]

 
write.csv(df3, file='../Output_Data/12_Shared_GO_terms.csv', row.names = T)

########2) Heatmap of shared GO terms  #################################################

#manual annotation of the GO-terms in excel. This is to summarize the GO-terms into larger categories for visibility.

########3) Heatmap of shared GO terms  #################################################

df4<- read.csv('../Raw_data/Shared_GO_terms_manual_annotation.csv', header= T) #manually load sorted GO terms

#clean data
rownames(df4)<-df4[,1]
df5<-df4[,-c(1,4:6)]

#double check that no mistake happened during manual annotation
sortdf3<-sort(df3[,1]) #shared GO-terms sorted by column 1
sortdf5<-sort(df5[,1]) #manual annotated shared GO-terms sorted by column 1

identical(sortdf3,sortdf5) #true

#plot
heatmap<-pheatmap(df5, 
         cluster_rows = F,
         fontsize_row=4,
         show_rownames = T,
         color = colorRampPalette(c("seashell", "chocolate1", "red3"))(100)) 

#save
pdf("../Output_Figures/12_Shared_GO_term_heatmap.pdf")
print(heatmap)
dev.off() 


# sessionInfo(package = NULL)
# 
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
#   [1] RColorBrewer_1.1-2 pheatmap_1.0.12    dplyr_1.0.6       
# 
# loaded via a namespace (and not attached):
#   [1] rstudioapi_0.13  magrittr_2.0.1   tidyselect_1.1.1 munsell_0.5.0    colorspace_2.0-1 R6_2.5.0         rlang_0.4.11     fansi_0.4.2      tools_3.6.1     
# [10] grid_3.6.1       gtable_0.3.0     utf8_1.2.1       cli_2.5.0        DBI_1.1.1        ellipsis_0.3.2   assertthat_0.2.1 tibble_3.1.2     lifecycle_1.0.0 
# [19] crayon_1.4.1     purrr_0.3.4      ggplot2_3.3.3    vctrs_0.3.8      glue_1.4.2       compiler_3.6.1   pillar_1.6.1     generics_0.1.0   scales_1.1.1    
# [28] pkgconfig_2.0.3
