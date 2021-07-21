# All code associated with Papsdorf et. al. 2021
## Scripts:
Scripts 1-10 are code analyzing the lipidomics data and generate the corresponding figures and tables.

Script 12 generates the GO-term heatmap

Script 13 calculates the lipid droplet peroxisome screen statistics

## Raw data:
Contains all the raw data that is necessary to run the scripts

## Table of scripts and corresponding figures
If a script relies on output from previous scripts, this is indicated in the table below. The corresponding figure panels are listed below.

|SCRIPT NAME|REQUIRED OUTPUTS FROM OTHER SCRIPTS TO RUN| OUTPUT|FIGURE PANELS/ TABLES|
---|---|---|---
|1_LIPIDS_PREPROCESSING_IS_NORMALIZATION_20210527|no|1_Lipids_Preprocessed_Normalized.csv<br>1_Lipids_Preprocessed_Normalized.Rdata|none|
|2_LIPIDS_IMPUTATION_CV_FILTERING|1_Lipids_Preprocessed_Normalized.Rdata|2_Lipids_Imputed_CVfiltered.csv <br>2_Lipids_Imputed_CVfiltered.Rdata|Extended Data Table 6 <br> Lipidomic Data.csv|
|3_LIPIDS_PCA_20210527|2_Lipids_Imputed_CVfiltered.Rdata|3_PCA.pdf<br> 3_Lipids_PC_loading.csv|Figure 4c<br>Extended Data Table 3 PCA loadings|
|4_CLASSES_20210527|	2_Lipids_Imputed_CVfiltered.Rdata	|4_Classes.csv<br>	4_Class_TG.pdf|Source Data Extended Data Figure 4a <br>4_Class_TG.pdf<br>	Figure S4a|
|5_SFAS_MUFAS_PUFAS_ALL_LIPIDS	|2_Lipids_Imputed_CVfiltered.Rdata	|5_SFAs_MUFAs_PUFAs_All_Lipids.csv<br>5_SFAs_MUFAs_PUFAs_All_Lipids.pdf|	Source Data Figure 4b<br>	Figure 4b|		
|6_SFAS_MUFAS_PUFAS_ETHER_LIPIDS	|2_Lipids_Imputed_CVfiltered.Rdata	|6_SFAs_MUFAs_PUFAs_Ether_Lipids.csv<br>6_SFAs_MUFAs_PUFAs_Ether_Lipids.pdf<br>6_MUFAtoPUFA_Ether_Lipids.pdf|	Source Data Figure 4e,f <br>Figure 4e<br>Figure 4f left panel|
|7_ETHER_LIPID_ABUNDANCE_20210528|	2_Lipids_Imputed_CVfiltered.Rdata|	7_Ether_Lipids.csv	<br>7_Ether_Lipids.pdf<br> 7_Ether_Lipid_TG_PE_PS.pdf| Source Data Extended Data Figure 4b, d <br>Figure S4b<br>Figure S4d|
|8_SFA_MUFA_PUFA_PE_ETHER_LIPIDS|	2_Lipids_Imputed_CVfiltered.Rdata	|8_MUFAtoPUFA_Ether_Lipids_PE_PS_TG.csv<br>8_MUFAtoPUFA_Ether_Lipids_PE_PS_TG.pdf|	Source Data Extended Data Figure 4e,f right panel<br> Figure S4e, Figure 4f right panel
|9_FA_ETHER_LIPIDS|	2_Lipids_Imputed_CVfiltered.Rdata|	9_FA_Ether_Lipids.csv	<br>9_FA_Ether_Lipids.pdf| Source Data Extended Data Figure 4c <br>Figure S4c|
|10_FA_LENGTH_ALL_LIPIDS	|2_Lipids_Imputed_CVfiltered.Rdata	|10_FA_Lenght_All_Lipids.csv|	Extended Data Table 7 Fatty acid length|
|12 SHARED_GO_TERMS_UPREGULATED_GREER_HAN|	no	|12_Shared_GO_terms.csv	<br> 12_Shared_Go_term_heatmap|Extended Data Table 2 GO Terms<br> Figure 3b|
|13_LIPID_DROPLET_PEROXISOME_SCREEN_STATS|	no	|13_Lipid_droplet_Peroxisome_Screen_Stats.csv	|Source Data Figure 3i, k|


