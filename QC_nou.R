

###################################################
###################################################
################## MAIN      ######################
###################################################
###################################################


setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/")


# VERBAL

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/pdp_verm_m100sd15_MAIN/INMA_202200601_pdp_verm_m100sd15_MAIN_allanalyses.RData")
ls()
summary(alldataout)
verbal<-as.data.frame(alldataout$AdjustedwithCellType)
colnames(verbal)<-c("probeID","N", "BETA", "SE", "zval", "P_VAL", "warnings")
head(verbal)
#probeID   N          BETA           SE       zval       P_VAL warnings
#1 cg14817997 255  5.523218e-04 1.824671e-04  3.0269668 0.002470211     none
#2 cg26928153 255  7.465912e-05 1.499919e-04  0.4977545 0.618657084     none
#3 cg16269199 255  5.007095e-05 1.786920e-04  0.2802081 0.779317881     none
#4 cg13869341 255 -5.517950e-05 1.335543e-04 -0.4131616 0.679488209     none
#5 cg14008030 255 -2.558632e-05 1.844984e-04 -0.1386804 0.889702673     none
#6 cg12045430 255  3.558604e-05 9.504152e-05  0.3744262 0.708087240     none


write.table(verbal, "verbalMAIN.txt", col.names=TRUE)

# to get lambdas

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/pdp_verm_m100sd15_MAIN/INMA_202200601_pdp_verm_m100sd15_MAIN_lambdas.RData")
ls()
#[1] "alldataout" "alllambda"  "verbal"     "verbalsig"

alllambda$Adjusted
alllambda$AdjustedwithCellType #0.9711832

# save or directly write the numbers in an excel file

# NONVERBAL

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/pdp_perm_m100sd15_MAIN/INMA_202200601_pdp_perm_m100sd15_MAIN_allanalyses.RData")
ls()
summary(alldataout)
nonverbal<-as.data.frame(alldataout$AdjustedwithCellType)
head(nonverbal)
#CpG   n          coef           se       zval       pval warnings
#1 cg14817997 255  3.409748e-04 1.840935e-04  1.8521831 0.06399954     none
#2 cg26928153 255 -4.277895e-05 1.447365e-04 -0.2955644 0.76756276     none
#3 cg16269199 255  5.600827e-05 1.883860e-04  0.2973058 0.76623304     none
#4 cg13869341 255  1.986705e-05 1.496629e-04  0.1327453 0.89439481     none
#5 cg14008030 255 -4.762588e-05 1.947259e-04 -0.2445790 0.80678242     none
#6 cg12045430 255 -9.340348e-05 8.504743e-05 -1.0982516 0.27209463     none

colnames(nonverbal)<-c("probeID","N", "BETA", "SE", "zval", "P_VAL", "warnings")
write.table(nonverbal, "nonverbalMAIN.txt", col.names=TRUE)

# to get lambdas

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/pdp_perm_m100sd15_MAIN/INMA_202200601_pdp_perm_m100sd15_MAIN_lambdas.RData")
ls()
#[1] "alldataout" "alllambda"  "verbal"     "verbalsig"

alllambda$Adjusted

alllambda$AdjustedwithCellType 
#0.9155405


# save or directly write the numbers in an excel file

# GENERAL 

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/pdp_genm_m100sd15_MAIN/INMA_202200601_pdp_genm_m100sd15_MAIN_allanalyses.RData")
ls()
summary(alldataout)
general<-as.data.frame(alldataout$AdjustedwithCellType)
head(general)
#CpG   n          coef           se         zval        pval warnings
#1 cg14817997 255  5.696398e-04 1.886366e-04  3.019773004 0.002529642     none
#2 cg26928153 255 -2.633021e-05 1.514363e-04 -0.173869846 0.861967751     none
#3 cg16269199 255  4.636216e-06 1.808294e-04  0.025638625 0.979545578     none
#4 cg13869341 255  1.081859e-06 1.451411e-04  0.007453844 0.994052748     none
#5 cg14008030 255 -5.119515e-05 1.865306e-04 -0.274459809 0.783731286     none
#6 cg12045430 255  7.991962e-06 9.072282e-05  0.088092086 0.929803487     none


colnames(general)<-c("probeID","N", "BETA", "SE", "zval", "P_VAL", "warnings")
write.table(general, "generalMAIN.txt", col.names=TRUE)

# to get lambdas

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/pdp_genm_m100sd15_MAIN/INMA_202200601_pdp_genm_m100sd15_MAIN_lambdas.RData")
ls()
#[1] "alldataout" "alllambda"  "verbal"     "verbalsig"

alllambda$Adjusted

alllambda$AdjustedwithCellType ## 0.8915743


# save or directly write the numbers in an excel file



##### QC extra amb paquet EASIER (to exclude CpGs, and to get FDR and Bonferroni significance as extra columns in the dataframe)
## -------------------------------------
##  Install EASIER Package Code
## -------------------------------------
##
##  Uncomment this code to install EASIER package
#
# # Install devtools
#install.packages("devtools")
#
# # Install required packages
#devtools::source_url("https://raw.githubusercontent.com/isglobal-brge/EASIER/HEAD/installer.R")

# # Install EASIER package
#devtools::install_github("isglobal-brge/EASIER@HEAD")

##  END -  Install EASIER Package Code
## -------------------------------------


library(EASIER)



########## ----------  VARIABLES DEFINED BY USER  ----------  ##########

# Set working directory to new folder
setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/")


# Files used in QC, needed in meta-analysis to plot ForestPlot (these files are the ones we have just saved)
files <- c('verbalMAIN.txt',
           'nonverbalMAIN.txt',
           'generalMAIN.txt')

# Result folder
results_folder <- 'QC_MAIN_NOU'

# Prefixes for each file
prefixes <- c('verbal', 'nonverbal', 'general')


# Exclude - MASK snp5
ethnic <- c('EUR', 'EUR', 'EUR')

# Array type, used : EPIC or 450K
artype <- c('EPIC', 'EPIC', 'EPIC')

# Parameters to exclude CpGs
exclude <- c( 'control_probes',
              'noncpg_probes',
              'Sex',
              'MASK_mapping',
              'MASK_sub30_copy',
              'MASK_extBase',
              'MASK_typeINextBaseSwitch',
              'MASK_snp5_ethnic',
              'Unrel_450_EPIC_blood')



N <- c(255, 255, 255)
n <- c(NA, NA, NA)

# Minimum sample representation percentage required for CpGs
# Filter minimum percentage of missings for each CpG in cohort
# We need to define two parameters,
#  - colname_NforProbe:
#        Column name with Number of individuals per probe, this variable only needs
#           to be defined if you want to filter CpGs with low representation.
#         If defined value in colname_NforProbe not exists, no filter will be applied
#  - pcMissingSamples :
#        Máximum percent of missing samples allowed,

colname_NforProbe <- 'N_for_probe'
pcMissingSamples <- 0.9


########## ----------  END VARIABLES DEFINED BY USER  ----------  ########## a partir d'aquí ja no s'ha de canviar res



## ###################### ##
##  QC - Quality Control  ##
## ###################### ##

# Variable declaration to perform precision plot
medianSE <- numeric(length(files))
value_N <- numeric(length(files))

if(length(n) == length(N))
  value_n <- numeric(length(files))

cohort_label <- character(length(files))

# Prepare output folder for results (create if not exists)
if(!dir.exists(file.path(getwd(), results_folder )))
  suppressWarnings(dir.create(file.path(getwd(), results_folder)))

## Remove duplicates, Exclude CpGs and adjust data (BN and FDR)
for ( i in 1:length(files) )
{
  
  # Prepare output subfolder for cohort-model results (create if not exists)
  if(!dir.exists(file.path(getwd(), results_folder, prefixes[i] )))
    suppressWarnings(dir.create(file.path(getwd(), results_folder, prefixes[i])))
  
  # Creates an empty file to resume all data if an old file exist  removes
  # the file and creates a new one
  fResumeName <- paste0( file.path(getwd(), results_folder, prefixes[i]),"/",prefixes[i], "_descriptives.txt")
  if ( file.exists(fResumeName) ) {
    file.remove(fResumeName)
  }
  file.create(fResumeName)
  
  # Read data.
  cohort <- read.table(files[i], header = TRUE, as.is = TRUE)
  print(paste0("Cohort file : ",files[i]," - readed OK", sep = " "))
  
  # Remove rows with NA from data
  cohort <- clean_NA_from_data(cohort)
  
  # Descriptives - Before CpGs deletion #
  descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], before = TRUE)
  
  # Remove duplicates
  # cohort <- remove_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i], '/',prefixes[i],'_descriptives_duplic.txt'), paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
  
  test_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
  
  # Remove cpGs with low representation
  # first, we test if colname_NforProbe and pcMissingSampes are defined
  if( !exists("colname_NforProbe") ) { colname_NforProbe <- NULL }
  if( !exists("pcMissingSamples") ) { pcMissingSamples <- NULL }
  
  cohort <- filterLowRepresentedCpGsinCohort(cohort, colname_NforProbe, pcMissingSamples, N[i], fileresume = fResumeName )
  
  # Exclude CpGs not meet conditions
  if("MASK_snp5_ethnic" %in% exclude ){
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = ethnic[i], filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
  } else {
    #..# if( !is.null(exclude) && exclude!='') {
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = "", filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
    #..# }
  }
  
  # Descriptives - After CpGs deletion #
  descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], before = FALSE )
  
  # Adjust data by Bonferroni and FDR
  cohort <- adjust_data(cohort, "P_VAL", bn=TRUE, fdr=TRUE, fResumeName, N[i]  )
  
  # Write QC complete data to external file
  write_QCData(cohort, paste0(results_folder, '/',prefixes[i], '/',prefixes[i]))
  
  ## Visualization - Plots
  #. Problems in some workstations and servers.# rasterpdf::raster_pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'), res = 300)
  #..# Problems in some cases --> Get png plots : pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'))
  
  # Distribution plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_SE_plot.png'))
  plot_distribution(cohort$SE, main = paste('Standard Errors of', prefixes[i]), xlab = 'SE')
  dev.off()
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_pvals_plot.png'))
  plot_distribution(cohort$P_VAL, main = paste('p-values of', prefixes[i]), xlab = 'p-value')
  dev.off()
  
  # QQ plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_QQ_plot.png'))
  qqman::qq(cohort$P_VAL, main = sprintf('QQ plot of %s (lambda = %f)', prefixes[i], lambda=get_lambda(cohort,"P_VAL")))
  dev.off()
  
  # Volcano plot.
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_Volcano_plot.png'))
  plot_volcano(cohort, "BETA", "P_VAL", main =paste('Volcano plot of', prefixes[i]) )
  dev.off()
  
  # Add mandatory data for precisionplot
  medianSE[i] <-  median(cohort$SE)
  value_N[i] <- N[i]
  cohort_label[i] <-  prefixes[i]
  
  # if n is defined for dichotomic condition :
  if(length(n) == length(N))  value_n[i] <- n[i]
  
  # Store data for Beta Box-Plot
  if( i == 1)
    betas.data <- list()
  betas.data[[prefixes[i]]] <- cohort[,"BETA"]
  
}




#########################################################
#########################################################
######################MATERNAL IQ #######################
#########################################################
#########################################################


setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220601NOUIQ_Output/Resultats_matIQ/")


# VERBAL

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220601NOUIQ_Output/pdp_verm_m100sd15_adj_maternal_IQNOU/INMA_20220601NOUIQ_pdp_verm_m100sd15_adj_maternal_IQNOU_allanalyses.RData")
ls()
summary(alldataout)
verbal<-as.data.frame(alldataout$AdjustedwithCellType)
colnames(verbal)<-c("probeID","N", "BETA", "SE", "zval", "P_VAL", "warnings")
write.table(verbal, "verbalmaternalIQ.txt", col.names=TRUE)

# to get lambdas
load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220601NOUIQ_Output/pdp_verm_m100sd15_adj_maternal_IQNOU/INMA_20220601NOUIQ_pdp_verm_m100sd15_adj_maternal_IQNOU_lambdas.RData")
ls()
#[1] "alldataout" "alllambda"  "verbal"     "verbalsig"

alllambda$Adjusted
alllambda$AdjustedwithCellType ## la buena! anar apuntant!

#0.9657854
# save or directly write the numbers in an excel file

# NONVERBAL

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220601NOUIQ_Output/pdp_perm_m100sd15_adj_maternal_IQNOU/INMA_20220601NOUIQ_pdp_perm_m100sd15_adj_maternal_IQNOU_allanalyses.RData")
ls()
summary(alldataout)
nonverbal<-as.data.frame(alldataout$AdjustedwithCellType)
colnames(nonverbal)<-c("probeID","N", "BETA", "SE", "zval", "P_VAL", "warnings")
write.table(nonverbal, "nonverbalIQ.txt", col.names=TRUE)

# to get lambdas

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220601NOUIQ_Output/pdp_perm_m100sd15_adj_maternal_IQNOU/INMA_20220601NOUIQ_pdp_perm_m100sd15_adj_maternal_IQNOU_lambdas.RData")
ls()
#[1] "alldataout" "alllambda"  "verbal"     "verbalsig"

alllambda$Adjusted

alllambda$AdjustedwithCellType ## la buena! anar apuntant!
#0.9110402

# save or directly write the numbers in an excel file

# GENERAL 

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220601NOUIQ_Output/pdp_genm_m100sd15_adj_maternal_IQNOU/INMA_20220601NOUIQ_pdp_genm_m100sd15_adj_maternal_IQNOU_allanalyses.RData")
ls()
summary(alldataout)
general<-as.data.frame(alldataout$AdjustedwithCellType)
colnames(general)<-c("probeID","N", "BETA", "SE", "zval", "P_VAL", "warnings")
write.table(general, "generalmaternalIQ.txt", col.names=TRUE)

# to get lambdas

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220601NOUIQ_Output/pdp_genm_m100sd15_adj_maternal_IQNOU/INMA_20220601NOUIQ_pdp_genm_m100sd15_adj_maternal_IQNOU_lambdas.RData")
ls()
#[1] "alldataout" "alllambda"  "verbal"     "verbalsig"

alllambda$Adjusted

alllambda$AdjustedwithCellType ## la buena! anar apuntant!
#0.8761689

# save or directly write the numbers in an excel file



##### QC extra amb paquet EASIER (to exclude CpGs, and to get FDR and Bonferroni significance as extra columns in the dataframe)


library(EASIER)



########## ----------  VARIABLES DEFINED BY USER  ----------  ##########

# Set working directory to new folder

setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220601NOUIQ_Output/Resultats_matIQ/")


# Files used in QC, needed in meta-analysis to plot ForestPlot (these files are the ones we have just saved)
files <- c('verbalmaternalIQ.txt',
           'nonverbalIQ.txt',
           'generalmaternalIQ.txt')

# Result folder
results_folder <- 'QC_ResultsmaternalIQ'

# Prefixes for each file
prefixes <- c('verbal', 'nonverbal', 'general')


# Exclude - MASK snp5
ethnic <- c('EUR', 'EUR', 'EUR')

# Array type, used : EPIC or 450K
artype <- c('EPIC', 'EPIC', 'EPIC')

# Parameters to exclude CpGs
exclude <- c( 'control_probes',
              'noncpg_probes',
              'Sex',
              'MASK_mapping',
              'MASK_sub30_copy',
              'MASK_extBase',
              'MASK_typeINextBaseSwitch',
              'MASK_snp5_ethnic',
              'Unrel_450_EPIC_blood')



N <- c(255, 255, 255)
n <- c(NA, NA, NA)

# Minimum sample representation percentage required for CpGs
# Filter minimum percentage of missings for each CpG in cohort
# We need to define two parameters,
#  - colname_NforProbe:
#        Column name with Number of individuals per probe, this variable only needs
#           to be defined if you want to filter CpGs with low representation.
#         If defined value in colname_NforProbe not exists, no filter will be applied
#  - pcMissingSamples :
#        Máximum percent of missing samples allowed,

colname_NforProbe <- 'N_for_probe'
pcMissingSamples <- 0.9


########## ----------  END VARIABLES DEFINED BY USER  ----------  ########## a partir d'aquí ja no s'ha de canviar res



## ###################### ##
##  QC - Quality Control  ##
## ###################### ##

# Variable declaration to perform precision plot
medianSE <- numeric(length(files))
value_N <- numeric(length(files))

if(length(n) == length(N))
  value_n <- numeric(length(files))

cohort_label <- character(length(files))

# Prepare output folder for results (create if not exists)
if(!dir.exists(file.path(getwd(), results_folder )))
  suppressWarnings(dir.create(file.path(getwd(), results_folder)))

## Remove duplicates, Exclude CpGs and adjust data (BN and FDR)
for ( i in 1:length(files) )
{
  
  # Prepare output subfolder for cohort-model results (create if not exists)
  if(!dir.exists(file.path(getwd(), results_folder, prefixes[i] )))
    suppressWarnings(dir.create(file.path(getwd(), results_folder, prefixes[i])))
  
  # Creates an empty file to resume all data if an old file exist  removes
  # the file and creates a new one
  fResumeName <- paste0( file.path(getwd(), results_folder, prefixes[i]),"/",prefixes[i], "_descriptives.txt")
  if ( file.exists(fResumeName) ) {
    file.remove(fResumeName)
  }
  file.create(fResumeName)
  
  # Read data.
  cohort <- read.table(files[i], header = TRUE, as.is = TRUE)
  print(paste0("Cohort file : ",files[i]," - readed OK", sep = " "))
  
  # Remove rows with NA from data
  cohort <- clean_NA_from_data(cohort)
  
  # Descriptives - Before CpGs deletion #
  descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], before = TRUE)
  
  # Remove duplicates
  # cohort <- remove_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i], '/',prefixes[i],'_descriptives_duplic.txt'), paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
  
  test_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
  
  # Remove cpGs with low representation
  # first, we test if colname_NforProbe and pcMissingSampes are defined
  if( !exists("colname_NforProbe") ) { colname_NforProbe <- NULL }
  if( !exists("pcMissingSamples") ) { pcMissingSamples <- NULL }
  
  cohort <- filterLowRepresentedCpGsinCohort(cohort, colname_NforProbe, pcMissingSamples, N[i], fileresume = fResumeName )
  
  # Exclude CpGs not meet conditions
  if("MASK_snp5_ethnic" %in% exclude ){
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = ethnic[i], filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
  } else {
    #..# if( !is.null(exclude) && exclude!='') {
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = "", filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
    #..# }
  }
  
  # Descriptives - After CpGs deletion #
  descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], before = FALSE )
  
  # Adjust data by Bonferroni and FDR
  cohort <- adjust_data(cohort, "P_VAL", bn=TRUE, fdr=TRUE, fResumeName, N[i]  )
  
  # Write QC complete data to external file
  write_QCData(cohort, paste0(results_folder, '/',prefixes[i], '/',prefixes[i]))
  
  ## Visualization - Plots
  #. Problems in some workstations and servers.# rasterpdf::raster_pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'), res = 300)
  #..# Problems in some cases --> Get png plots : pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'))
  
  # Distribution plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_SE_plot.png'))
  plot_distribution(cohort$SE, main = paste('Standard Errors of', prefixes[i]), xlab = 'SE')
  dev.off()
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_pvals_plot.png'))
  plot_distribution(cohort$P_VAL, main = paste('p-values of', prefixes[i]), xlab = 'p-value')
  dev.off()
  
  # QQ plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_QQ_plot.png'))
  qqman::qq(cohort$P_VAL, main = sprintf('QQ plot of %s (lambda = %f)', prefixes[i], lambda=get_lambda(cohort,"P_VAL")))
  dev.off()
  
  # Volcano plot.
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_Volcano_plot.png'))
  plot_volcano(cohort, "BETA", "P_VAL", main =paste('Volcano plot of', prefixes[i]) )
  dev.off()
  
  # Add mandatory data for precisionplot
  medianSE[i] <-  median(cohort$SE)
  value_N[i] <- N[i]
  cohort_label[i] <-  prefixes[i]
  
  # if n is defined for dichotomic condition :
  if(length(n) == length(N))  value_n[i] <- n[i]
  
  # Store data for Beta Box-Plot
  if( i == 1)
    betas.data <- list()
  betas.data[[prefixes[i]]] <- cohort[,"BETA"]
  
}




#########################################################
#########################################################
######################SENSITIVITY #######################
#########################################################
#########################################################


setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220707_Output/SENSITIVITY_NOU/")


# VERBAL

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220707_Output/pdp_verm_m100sd15_sensitivity/INMA_20220707_pdp_verm_m100sd15_sensitivity_allanalyses.RData")
ls()
summary(alldataout)
verbal<-as.data.frame(alldataout$AdjustedwithCellType)
head(verbal)
#         CpG   n          coef           se         zval        pval warnings
#1 cg14817997 210  5.752185e-04 1.930221e-04  2.980065877 0.002881864     none
#2 cg26928153 210 -4.631046e-07 1.595916e-04 -0.002901811 0.997684693     none
#3 cg16269199 210  1.022585e-05 1.781021e-04  0.057415697 0.954214059     none
#4 cg13869341 210  3.704182e-06 1.346273e-04  0.027514354 0.978049491     none
#5 cg14008030 210 -6.921928e-05 1.940098e-04 -0.356782495 0.721254642     none
#6 cg12045430 210 -9.362871e-06 9.950132e-05 -0.094097961 0.925031340     none
colnames(verbal)<-c("probeID","N", "BETA", "SE", "zval", "P_VAL", "warnings")
write.table(verbal, "verbalsens.txt", col.names=TRUE)

# to get lambdas
load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220707_Output/pdp_verm_m100sd15_sensitivity/INMA_20220707_pdp_verm_m100sd15_sensitivity_lambdas.RData")
ls()
#[1] "alldataout" "alllambda"  "verbal"     "verbalsig"

alllambda$Adjusted
alllambda$AdjustedwithCellType ## la buena! anar apuntant!

#1.10016
# save or directly write the numbers in an excel file

# NONVERBAL

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220707_Output/pdp_perm_m100sd15_sensitivity/INMA_20220707_pdp_perm_m100sd15_sensitivity_allanalyses.RData")
ls()
summary(alldataout)
nonverbal<-as.data.frame(alldataout$AdjustedwithCellType)
head(nonverbal)
#   CpG   n          coef           se       zval      pval warnings
#1 cg14817997 210  2.744763e-04 1.924782e-04  1.4260128 0.1538646     none
#2 cg26928153 210 -7.309167e-05 1.526811e-04 -0.4787213 0.6321369     none
#3 cg16269199 210 -2.279057e-05 1.967121e-04 -0.1158575 0.9077655     none
#4 cg13869341 210 -1.819256e-05 1.667732e-04 -0.1090856 0.9131346     none
#5 cg14008030 210 -2.721984e-05 1.800202e-04 -0.1512044 0.8798145     none
#6 cg12045430 210 -1.133677e-04 8.874105e-05 -1.2775117 0.2014217     none

colnames(nonverbal)<-c("probeID","N", "BETA", "SE", "zval", "P_VAL", "warnings")
write.table(nonverbal, "nonverbalsens.txt", col.names=TRUE)

# to get lambdas

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220707_Output/pdp_perm_m100sd15_sensitivity/INMA_20220707_pdp_perm_m100sd15_sensitivity_lambdas.RData")
ls()
#[1] "alldataout" "alllambda"  "verbal"     "verbalsig"

alllambda$Adjusted

alllambda$AdjustedwithCellType ## la buena! anar apuntant!
#0.9961188


# save or directly write the numbers in an excel file

# GENERAL 

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220707_Output/pdp_genm_m100sd15_sensitivity/INMA_20220707_pdp_genm_m100sd15_sensitivity_allanalyses.RData")
ls()
summary(alldataout)
general<-as.data.frame(alldataout$AdjustedwithCellType)
head(general)
#     probeID   N          BETA           SE       zval       P_VAL warnings
#1 cg14817997 210  5.430941e-04 2.026557e-04  2.6798854 0.007364738     none
#2 cg26928153 210 -1.003440e-04 1.562202e-04 -0.6423239 0.520662872     none
#3 cg16269199 210 -6.429859e-05 1.926202e-04 -0.3338103 0.738522741     none
#4 cg13869341 210  2.086720e-05 1.509645e-04  0.1382258 0.890061935     none
#5 cg14008030 210 -5.488648e-05 1.872499e-04 -0.2931189 0.769431253     none
#6 cg12045430 210 -2.566694e-05 9.152707e-05 -0.2804300 0.779147587     none
colnames(general)<-c("probeID","N", "BETA", "SE", "zval", "P_VAL", "warnings")
write.table(general, "generalsens.txt", col.names=TRUE)

# to get lambdas

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220707_Output/pdp_genm_m100sd15_sensitivity/INMA_20220707_pdp_genm_m100sd15_sensitivity_lambdas.RData")
ls()
#[1] "alldataout" "alllambda"  "verbal"     "verbalsig"

alllambda$Adjusted

alllambda$AdjustedwithCellType ## la buena! anar apuntant!
#1.017348


# save or directly write the numbers in an excel file



##### QC extra amb paquet EASIER (to exclude CpGs, and to get FDR and Bonferroni significance as extra columns in the dataframe)


library(EASIER)



########## ----------  VARIABLES DEFINED BY USER  ----------  ##########

# Set working directory to new folder

setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220707_Output/SENSITIVITY_NOU/")


# Files used in QC, needed in meta-analysis to plot ForestPlot (these files are the ones we have just saved)
files <- c('verbalsens.txt',
           'nonverbalsens.txt',
           'generalsens.txt')

# Result folder
results_folder <- 'QC_ResultsSENSITIVITY'

# Prefixes for each file
prefixes <- c('verbal', 'nonverbal', 'general')


# Exclude - MASK snp5
ethnic <- c('EUR', 'EUR', 'EUR')

# Array type, used : EPIC or 450K
artype <- c('EPIC', 'EPIC', 'EPIC')

# Parameters to exclude CpGs
exclude <- c( 'control_probes',
              'noncpg_probes',
              'Sex',
              'MASK_mapping',
              'MASK_sub30_copy',
              'MASK_extBase',
              'MASK_typeINextBaseSwitch',
              'MASK_snp5_ethnic',
              'Unrel_450_EPIC_blood')



N <- c(255, 255, 255)
n <- c(NA, NA, NA)

# Minimum sample representation percentage required for CpGs
# Filter minimum percentage of missings for each CpG in cohort
# We need to define two parameters,
#  - colname_NforProbe:
#        Column name with Number of individuals per probe, this variable only needs
#           to be defined if you want to filter CpGs with low representation.
#         If defined value in colname_NforProbe not exists, no filter will be applied
#  - pcMissingSamples :
#        Máximum percent of missing samples allowed,

colname_NforProbe <- 'N_for_probe'
pcMissingSamples <- 0.9


########## ----------  END VARIABLES DEFINED BY USER  ----------  ########## a partir d'aquí ja no s'ha de canviar res



## ###################### ##
##  QC - Quality Control  ##
## ###################### ##

# Variable declaration to perform precision plot
medianSE <- numeric(length(files))
value_N <- numeric(length(files))

if(length(n) == length(N))
  value_n <- numeric(length(files))

cohort_label <- character(length(files))

# Prepare output folder for results (create if not exists)
if(!dir.exists(file.path(getwd(), results_folder )))
  suppressWarnings(dir.create(file.path(getwd(), results_folder)))

## Remove duplicates, Exclude CpGs and adjust data (BN and FDR)
for ( i in 1:length(files) )
{
  
  # Prepare output subfolder for cohort-model results (create if not exists)
  if(!dir.exists(file.path(getwd(), results_folder, prefixes[i] )))
    suppressWarnings(dir.create(file.path(getwd(), results_folder, prefixes[i])))
  
  # Creates an empty file to resume all data if an old file exist  removes
  # the file and creates a new one
  fResumeName <- paste0( file.path(getwd(), results_folder, prefixes[i]),"/",prefixes[i], "_descriptives.txt")
  if ( file.exists(fResumeName) ) {
    file.remove(fResumeName)
  }
  file.create(fResumeName)
  
  # Read data.
  cohort <- read.table(files[i], header = TRUE, as.is = TRUE)
  print(paste0("Cohort file : ",files[i]," - readed OK", sep = " "))
  
  # Remove rows with NA from data
  cohort <- clean_NA_from_data(cohort)
  
  # Descriptives - Before CpGs deletion #
  descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], before = TRUE)
  
  # Remove duplicates
  # cohort <- remove_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i], '/',prefixes[i],'_descriptives_duplic.txt'), paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
  
  test_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
  
  # Remove cpGs with low representation
  # first, we test if colname_NforProbe and pcMissingSampes are defined
  if( !exists("colname_NforProbe") ) { colname_NforProbe <- NULL }
  if( !exists("pcMissingSamples") ) { pcMissingSamples <- NULL }
  
  cohort <- filterLowRepresentedCpGsinCohort(cohort, colname_NforProbe, pcMissingSamples, N[i], fileresume = fResumeName )
  
  # Exclude CpGs not meet conditions
  if("MASK_snp5_ethnic" %in% exclude ){
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = ethnic[i], filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
  } else {
    #..# if( !is.null(exclude) && exclude!='') {
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = "", filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
    #..# }
  }
  
  # Descriptives - After CpGs deletion #
  descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], before = FALSE )
  
  # Adjust data by Bonferroni and FDR
  cohort <- adjust_data(cohort, "P_VAL", bn=TRUE, fdr=TRUE, fResumeName, N[i]  )
  
  # Write QC complete data to external file
  write_QCData(cohort, paste0(results_folder, '/',prefixes[i], '/',prefixes[i]))
  
  ## Visualization - Plots
  #. Problems in some workstations and servers.# rasterpdf::raster_pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'), res = 300)
  #..# Problems in some cases --> Get png plots : pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'))
  
  # Distribution plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_SE_plot.png'))
  plot_distribution(cohort$SE, main = paste('Standard Errors of', prefixes[i]), xlab = 'SE')
  dev.off()
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_pvals_plot.png'))
  plot_distribution(cohort$P_VAL, main = paste('p-values of', prefixes[i]), xlab = 'p-value')
  dev.off()
  
  # QQ plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_QQ_plot.png'))
  qqman::qq(cohort$P_VAL, main = sprintf('QQ plot of %s (lambda = %f)', prefixes[i], lambda=get_lambda(cohort,"P_VAL")))
  dev.off()
  
  # Volcano plot.
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_Volcano_plot.png'))
  plot_volcano(cohort, "BETA", "P_VAL", main =paste('Volcano plot of', prefixes[i]) )
  dev.off()
  
  # Add mandatory data for precisionplot
  medianSE[i] <-  median(cohort$SE)
  value_N[i] <- N[i]
  cohort_label[i] <-  prefixes[i]
  
  # if n is defined for dichotomic condition :
  if(length(n) == length(N))  value_n[i] <- n[i]
  
  # Store data for Beta Box-Plot
  if( i == 1)
    betas.data <- list()
  betas.data[[prefixes[i]]] <- cohort[,"BETA"]
  
}


#########################################################
#########################################################
####################### SGES ############################
#########################################################
#########################################################


setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220722NOU_Output/resultats_nous_sges/")


# VERBAL

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220722NOU_Output/pdp_verm_m100sd15_not_adj_sges/INMA_20220722NOU_pdp_verm_m100sd15_not_adj_sges_allanalyses.RData")
ls()
summary(alldataout)
verbal<-as.data.frame(alldataout$AdjustedwithCellType)
head(verbal)
#       CpG   n          coef           se       zval       pval warnings
#1 cg14817997 255  5.487403e-04 1.821659e-04  3.0123108 0.00259267     none
#2 cg26928153 255  7.713974e-05 1.461374e-04  0.5278577 0.59759813     none
#3 cg16269199 255  4.991648e-05 1.764986e-04  0.2828151 0.77731857     none
#4 cg13869341 255 -5.243179e-05 1.345949e-04 -0.3895526 0.69686743     none
#5 cg14008030 255 -2.628334e-05 1.840860e-04 -0.1427775 0.88646590     none
#6 cg12045430 255  3.713033e-05 9.738041e-05  0.3812916 0.70298689     none


colnames(verbal)<-c("probeID","N", "BETA", "SE", "zval", "P_VAL", "warnings")
write.table(verbal, "verbalsges.txt", col.names=TRUE)

# to get lambdas
load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220722NOU_Output/pdp_verm_m100sd15_not_adj_sges/INMA_20220722NOU_pdp_verm_m100sd15_not_adj_sges_lambdas.RData")
ls()
#[1] "alldataout" "alllambda"  "verbal"     "verbalsig"

alllambda$Adjusted
alllambda$AdjustedwithCellType ## la buena! anar apuntant!

#0.9687708


# save or directly write the numbers in an excel file

# NONVERBAL

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220722NOU_Output/pdp_perm_m100sd15_not_adj_sges/INMA_20220722NOU_pdp_perm_m100sd15_not_adj_sges_allanalyses.RData")
ls()
summary(alldataout)
nonverbal<-as.data.frame(alldataout$AdjustedwithCellType)
head(nonverbal)
#            CpG   n          coef           se       zval       pval warnings
#1 cg14817997 255  3.133304e-04 1.874411e-04  1.6716211 0.09459906     none
#2 cg26928153 255 -5.870685e-05 1.418764e-04 -0.4137888 0.67902876     none
#3 cg16269199 255  5.395176e-05 1.887039e-04  0.2859069 0.77494942     none
#4 cg13869341 255  2.720427e-05 1.531782e-04  0.1775988 0.85903803     none
#5 cg14008030 255 -4.458685e-05 1.940430e-04 -0.2297782 0.81826414     none
#6 cg12045430 255 -9.699235e-05 8.255347e-05 -1.1749033 0.24003339     none

colnames(nonverbal)<-c("probeID","N", "BETA", "SE", "zval", "P_VAL", "warnings")
write.table(nonverbal, "nonverbalsges.txt", col.names=TRUE)

# to get lambdas

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220722NOU_Output/pdp_perm_m100sd15_not_adj_sges/INMA_20220722NOU_pdp_perm_m100sd15_not_adj_sges_lambdas.RData")
ls()
#[1] "alldataout" "alllambda"  "verbal"     "verbalsig"

alllambda$Adjusted

alllambda$AdjustedwithCellType ## la buena! anar apuntant!
# 0.9302214



# save or directly write the numbers in an excel file

# GENERAL 

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220722NOU_Output/pdp_genm_m100sd15_not_adj_sges/INMA_20220722NOU_pdp_genm_m100sd15_not_adj_sges_allanalyses.RData")
ls()
summary(alldataout)
general<-as.data.frame(alldataout$AdjustedwithCellType)
head(general)
# CpG   n          coef           se        zval        pval warnings
#1 cg14817997 255  5.599544e-04 1.900398e-04  2.94651121 0.003213808     none
#2 cg26928153 255 -3.360315e-05 1.479860e-04 -0.22706976 0.820369497     none
#3 cg16269199 255  3.593525e-06 1.827743e-04  0.01966100 0.984313803     none
#4 cg13869341 255  6.539855e-06 1.436726e-04  0.04551915 0.963693512     none
#5 cg14008030 255 -5.032207e-05 1.868704e-04 -0.26928860 0.787707603     none
#6 cg12045430 255  6.182691e-06 9.120609e-05  0.06778814 0.945954282     none


colnames(general)<-c("probeID","N", "BETA", "SE", "zval", "P_VAL", "warnings")
write.table(general, "generalsges.txt", col.names=TRUE)

# to get lambdas

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220722NOU_Output/pdp_genm_m100sd15_not_adj_sges/INMA_20220722NOU_pdp_genm_m100sd15_not_adj_sges_lambdas.RData")
ls()
#[1] "alldataout" "alllambda"  "verbal"     "verbalsig"

alllambda$Adjusted

alllambda$AdjustedwithCellType ## la buena! anar apuntant!
#0.8933315



# save or directly write the numbers in an excel file



##### QC extra amb paquet EASIER (to exclude CpGs, and to get FDR and Bonferroni significance as extra columns in the dataframe)


library(EASIER)



########## ----------  VARIABLES DEFINED BY USER  ----------  ##########

# Set working directory to new folder

setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220722NOU_Output/resultats_nous_sges/")


# Files used in QC, needed in meta-analysis to plot ForestPlot (these files are the ones we have just saved)
files <- c('verbalsges.txt',
           'nonverbalsges.txt',
           'generalsges.txt')

# Result folder
results_folder <- 'QC_ResultsSGES'

# Prefixes for each file
prefixes <- c('verbal', 'nonverbal', 'general')


# Exclude - MASK snp5
ethnic <- c('EUR', 'EUR', 'EUR')

# Array type, used : EPIC or 450K
artype <- c('EPIC', 'EPIC', 'EPIC')

# Parameters to exclude CpGs
exclude <- c( 'control_probes',
              'noncpg_probes',
              'Sex',
              'MASK_mapping',
              'MASK_sub30_copy',
              'MASK_extBase',
              'MASK_typeINextBaseSwitch',
              'MASK_snp5_ethnic',
              'Unrel_450_EPIC_blood')



N <- c(255, 255, 255)
n <- c(NA, NA, NA)

# Minimum sample representation percentage required for CpGs
# Filter minimum percentage of missings for each CpG in cohort
# We need to define two parameters,
#  - colname_NforProbe:
#        Column name with Number of individuals per probe, this variable only needs
#           to be defined if you want to filter CpGs with low representation.
#         If defined value in colname_NforProbe not exists, no filter will be applied
#  - pcMissingSamples :
#        Máximum percent of missing samples allowed,

colname_NforProbe <- 'N_for_probe'
pcMissingSamples <- 0.9


########## ----------  END VARIABLES DEFINED BY USER  ----------  ########## a partir d'aquí ja no s'ha de canviar res



## ###################### ##
##  QC - Quality Control  ##
## ###################### ##

# Variable declaration to perform precision plot
medianSE <- numeric(length(files))
value_N <- numeric(length(files))

if(length(n) == length(N))
  value_n <- numeric(length(files))

cohort_label <- character(length(files))

# Prepare output folder for results (create if not exists)
if(!dir.exists(file.path(getwd(), results_folder )))
  suppressWarnings(dir.create(file.path(getwd(), results_folder)))

## Remove duplicates, Exclude CpGs and adjust data (BN and FDR)
for ( i in 1:length(files) )
{
  
  # Prepare output subfolder for cohort-model results (create if not exists)
  if(!dir.exists(file.path(getwd(), results_folder, prefixes[i] )))
    suppressWarnings(dir.create(file.path(getwd(), results_folder, prefixes[i])))
  
  # Creates an empty file to resume all data if an old file exist  removes
  # the file and creates a new one
  fResumeName <- paste0( file.path(getwd(), results_folder, prefixes[i]),"/",prefixes[i], "_descriptives.txt")
  if ( file.exists(fResumeName) ) {
    file.remove(fResumeName)
  }
  file.create(fResumeName)
  
  # Read data.
  cohort <- read.table(files[i], header = TRUE, as.is = TRUE)
  print(paste0("Cohort file : ",files[i]," - readed OK", sep = " "))
  
  # Remove rows with NA from data
  cohort <- clean_NA_from_data(cohort)
  
  # Descriptives - Before CpGs deletion #
  descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], before = TRUE)
  
  # Remove duplicates
  # cohort <- remove_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i], '/',prefixes[i],'_descriptives_duplic.txt'), paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
  
  test_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
  
  # Remove cpGs with low representation
  # first, we test if colname_NforProbe and pcMissingSampes are defined
  if( !exists("colname_NforProbe") ) { colname_NforProbe <- NULL }
  if( !exists("pcMissingSamples") ) { pcMissingSamples <- NULL }
  
  cohort <- filterLowRepresentedCpGsinCohort(cohort, colname_NforProbe, pcMissingSamples, N[i], fileresume = fResumeName )
  
  # Exclude CpGs not meet conditions
  if("MASK_snp5_ethnic" %in% exclude ){
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = ethnic[i], filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
  } else {
    #..# if( !is.null(exclude) && exclude!='') {
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = "", filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
    #..# }
  }
  
  # Descriptives - After CpGs deletion #
  descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], before = FALSE )
  
  # Adjust data by Bonferroni and FDR
  cohort <- adjust_data(cohort, "P_VAL", bn=TRUE, fdr=TRUE, fResumeName, N[i]  )
  
  # Write QC complete data to external file
  write_QCData(cohort, paste0(results_folder, '/',prefixes[i], '/',prefixes[i]))
  
  ## Visualization - Plots
  #. Problems in some workstations and servers.# rasterpdf::raster_pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'), res = 300)
  #..# Problems in some cases --> Get png plots : pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'))
  
  # Distribution plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_SE_plot.png'))
  plot_distribution(cohort$SE, main = paste('Standard Errors of', prefixes[i]), xlab = 'SE')
  dev.off()
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_pvals_plot.png'))
  plot_distribution(cohort$P_VAL, main = paste('p-values of', prefixes[i]), xlab = 'p-value')
  dev.off()
  
  # QQ plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_QQ_plot.png'))
  qqman::qq(cohort$P_VAL, main = sprintf('QQ plot of %s (lambda = %f)', prefixes[i], lambda=get_lambda(cohort,"P_VAL")))
  dev.off()
  
  # Volcano plot.
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_Volcano_plot.png'))
  plot_volcano(cohort, "BETA", "P_VAL", main =paste('Volcano plot of', prefixes[i]) )
  dev.off()
  
  # Add mandatory data for precisionplot
  medianSE[i] <-  median(cohort$SE)
  value_N[i] <- N[i]
  cohort_label[i] <-  prefixes[i]
  
  # if n is defined for dichotomic condition :
  if(length(n) == length(N))  value_n[i] <- n[i]
  
  # Store data for Beta Box-Plot
  if( i == 1)
    betas.data <- list()
  betas.data[[prefixes[i]]] <- cohort[,"BETA"]
  
}






#########################################################
#########################################################
####################### CELLTYPES########################
#########################################################
#########################################################

#######GENERAL

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/pdp_genm_m100sd15_MAINcelltypes/INMA_202200715_pdp_genm_m100sd15_MAINcelltypes_allanalyses.RData")
summary(alldataout$CellInteraction$coe)
#Length Class      Mode
#Trophoblasts        5      data.frame list
#Hofbauer            5      data.frame list
#Endothelial         5      data.frame list
#nRBC                5      data.frame list
#Syncytiotrophoblast 5      data.frame list
trophoblast<- alldataout$CellInteraction$coe$Trophoblasts
hofbauer<- alldataout$CellInteraction$coe$Hofbauer
endothelial<-alldataout$CellInteraction$coe$Endothelial
nRBC<- alldataout$CellInteraction$coe$nRBC
syncytiotrophoblast<-alldataout$CellInteraction$coe$Syncytiotrophoblast

#Make sure all of them are data.frames. 
class(trophoblast)
#[1] "data.frame"
class(hofbauer)
class(endothelial)
class(nRBC)
class(syncytiotrophoblast)

write.csv(trophoblast,"generaltrophoblast.csv")
write.csv(hofbauer,"generalhofbauer.csv")
write.csv(endothelial,"generalendothelial.csv")
write.csv(nRBC,"generalnRBC.csv")
write.csv(syncytiotrophoblast,"generalsyncytiotrophoblast.csv")


# to get lambdas

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/pdp_genm_m100sd15_MAINcelltypes/INMA_202200715_pdp_genm_m100sd15_MAINcelltypes_lambdas.RData")
ls()
#[1] "alldataout" "alllambda"  "verbal"     "verbalsig"

alllambda$Adjusted
#0.7913024
alllambda$AdjustedwithCellType 
#0.8915743

#####################################################
#
# NON-VERBAL
load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/pdp_perm_m100sd15_MAINcelltypes/INMA_202200715_pdp_perm_m100sd15_MAINcelltypes_allanalyses.RData")
summary(alldataout$CellInteraction$coe)
trophoblast<- alldataout$CellInteraction$coe$Trophoblasts
hofbauer<- alldataout$CellInteraction$coe$Hofbauer
endothelial<-alldataout$CellInteraction$coe$Endothelial
nRBC<- alldataout$CellInteraction$coe$nRBC
syncytiotrophoblast<-alldataout$CellInteraction$coe$Syncytiotrophoblast

#Make sure all of them are data.frames. 
class(trophoblast)
#[1] "data.frame"
class(hofbauer)
class(endothelial)
class(nRBC)
class(syncytiotrophoblast)

write.csv(trophoblast,"nonverbaltrophoblast.csv")
write.csv(hofbauer,"nonverbalhofbauer.csv")
write.csv(endothelial,"nonverbalendothelial.csv")
write.csv(nRBC,"nonverbalnRBC.csv")
write.csv(syncytiotrophoblast,"nonverbalsyncytiotrophoblast.csv")

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/pdp_perm_m100sd15_MAINcelltypes/INMA_202200715_pdp_perm_m100sd15_MAINcelltypes_lambdas.RData")
ls()
#[1] "alldataout" "alllambda"  "verbal"     "verbalsig"

alllambda$Adjusted
#0.9232329
alllambda$AdjustedwithCellType 
#0.9155405

#####################################################
#
#VERBAL
load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/pdp_verm_m100sd15_MAINcelltypes/INMA_202200715_pdp_verm_m100sd15_MAINcelltypes_allanalyses.RData")

summary(alldataout$CellInteraction$coe)
trophoblast<- alldataout$CellInteraction$coe$Trophoblasts
hofbauer<- alldataout$CellInteraction$coe$Hofbauer
endothelial<-alldataout$CellInteraction$coe$Endothelial
nRBC<- alldataout$CellInteraction$coe$nRBC
syncytiotrophoblast<-alldataout$CellInteraction$coe$Syncytiotrophoblast

#Make sure all of them are data.frames. 
class(trophoblast)
#[1] "data.frame"
class(hofbauer)
class(endothelial)
class(nRBC)
class(syncytiotrophoblast)

write.csv(trophoblast,"verbaltrophoblast.csv")
write.csv(hofbauer,"verbalhofbauer.csv")
write.csv(endothelial,"verbalendothelial.csv")
write.csv(nRBC,"verbalnRBC.csv")
write.csv(syncytiotrophoblast,"verbalsyncytiotrophoblast.csv")

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/pdp_verm_m100sd15_MAINcelltypes/INMA_202200715_pdp_verm_m100sd15_MAINcelltypes_lambdas.RData")
ls()
#[1] "alldataout" "alllambda"  "verbal"     "verbalsig"

alllambda$Adjusted
#0.8632186
alllambda$AdjustedwithCellType 
#0.9711832


########################LET's DO THE QC for each one.



#######################
#ENDOTHELIAL
#######################

generalendothelial <- read.csv("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/generalendothelial.csv")
head(generalendothelial)
#   X     Estimate          SE          t         p      adjP
#1 cg14817997  0.004795756 0.008860875  0.5412282 0.5888696 0.9999763
#2 cg26928153  0.013292849 0.006333929  2.0986736 0.0369272 0.9705185
#3 cg16269199  0.002894123 0.008685489  0.3332136 0.7392740 0.9999763
#4 cg13869341 -0.005114711 0.005495237 -0.9307536 0.3529485 0.9999763
#5 cg14008030  0.002348124 0.007284765  0.3223335 0.7474903 0.9999763
#6 cg12045430  0.002466542 0.003492715  0.7061961 0.4807748 0.9999763



N<- rep(255, times= 811990)
N<- as.vector(N)
generalendothelial1<- data.frame(probeID=generalendothelial$X, N= N, BETA= generalendothelial$Estimate, SE= generalendothelial$SE, zval= generalendothelial$t, P_VAL=generalendothelial$p )
head(generalendothelial1)
#   probeID   N         BETA          SE       zval     P_VAL
#1 cg14817997 255  0.004795756 0.008860875  0.5412282 0.5888696
#2 cg26928153 255  0.013292849 0.006333929  2.0986736 0.0369272
#3 cg16269199 255  0.002894123 0.008685489  0.3332136 0.7392740
#4 cg13869341 255 -0.005114711 0.005495237 -0.9307536 0.3529485
#5 cg14008030 255  0.002348124 0.007284765  0.3223335 0.7474903
#6 cg12045430 255  0.002466542 0.003492715  0.7061961 0.4807748

write.table(generalendothelial1, "generalendothelial1.txt", col.names=TRUE)

###################non verbal

nonverbalendothelial<- read.csv("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/nonverbalendothelial.csv")
head(nonverbalendothelial)
#       X      Estimate          SE          t         p      adjP
#1 cg14817997  0.0078774578 0.008724324  0.9029304 0.3674994 0.8900535
#2 cg26928153  0.0024323496 0.006250781  0.3891273 0.6975389 0.9603024
#3 cg16269199  0.0004183192 0.008556843  0.0488871 0.9610513 0.9958441
#4 cg13869341 -0.0019151362 0.005379655 -0.3559961 0.7221670 0.9643202
#5 cg14008030  0.0033923849 0.007101628  0.4776912 0.6333198 0.9493441
#6 cg12045430 -0.0017124591 0.003403142 -0.5031994 0.6153015 0.9461732

N<- rep(255, times= 811990)
N<- as.vector(N)
nonverbalendothelial1<- data.frame(probeID=nonverbalendothelial$X, N= N, BETA= nonverbalendothelial$Estimate, SE= nonverbalendothelial$SE, zval= nonverbalendothelial$t, P_VAL=nonverbalendothelial$p )
head(nonverbalendothelial1)
#   probeID   N          BETA          SE       zval     P_VAL
#1 cg14817997 255  0.0078774578 0.008724324  0.9029304 0.3674994
#2 cg26928153 255  0.0024323496 0.006250781  0.3891273 0.6975389
#3 cg16269199 255  0.0004183192 0.008556843  0.0488871 0.9610513
#4 cg13869341 255 -0.0019151362 0.005379655 -0.3559961 0.7221670
#5 cg14008030 255  0.0033923849 0.007101628  0.4776912 0.6333198
#6 cg12045430 255 -0.0017124591 0.003403142 -0.5031994 0.6153015

write.table(nonverbalendothelial1, "nonverbalendothelial1.txt", col.names=TRUE)


######################verbal

verbalendothelial<- read.csv("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/verbalsyncytiotrophoblast.csv")
head(verbalendothelial)
#  X     Estimate          SE          t           p      adjP
#1 cg14817997  0.002257937 0.009347970  0.2415431 0.809347603 0.9642086
#2 cg26928153  0.021419618 0.006626903  3.2322212 0.001406739 0.4601419
#3 cg16269199  0.004381953 0.009117651  0.4806011 0.631252902 0.9198484
#4 cg13869341 -0.004715146 0.005822187 -0.8098583 0.418851501 0.8460703
#5 cg14008030  0.001427534 0.007713306  0.1850742 0.853332561 0.9735388
#6 cg12045430  0.003809615 0.003693448  1.0314521 0.303402924 0.7932151

N<- rep(255, times= 811990)
N<- as.vector(N)
verbalendothelial1<- data.frame(probeID=verbalendothelial$X, N= N, BETA= verbalendothelial$Estimate, SE= verbalendothelial$SE, zval= verbalendothelial$t, P_VAL=verbalendothelial$p )
head(verbalendothelial1)
#probeID   N         BETA          SE       zval       P_VAL
#1 cg14817997 255  0.002257937 0.009347970  0.2415431 0.809347603
#2 cg26928153 255  0.021419618 0.006626903  3.2322212 0.001406739
#3 cg16269199 255  0.004381953 0.009117651  0.4806011 0.631252902
#4 cg13869341 255 -0.004715146 0.005822187 -0.8098583 0.418851501
#5 cg14008030 255  0.001427534 0.007713306  0.1850742 0.853332561
#6 cg12045430 255  0.003809615 0.003693448  1.0314521 0.303402924

write.table(verbalendothelial1, "verbalendothelial1.txt", col.names=TRUE)






library(EASIER)



########## ----------  VARIABLES DEFINED BY USER  ----------  ##########

# Set working directory to new folder
setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/")


# Files used in QC, needed in meta-analysis to plot ForestPlot (these files are the ones we have just saved)
files <- c('generalendothelial1.txt',
           'nonverbalendothelial1.txt',
           'verbalendothelial1.txt')

# Result folder
results_folder <- 'QC_ResultsENDOTHELIAL'

# Prefixes for each file
prefixes <- c('verbal', 'nonverbal', 'general')


# Exclude - MASK snp5
ethnic <- c('EUR', 'EUR', 'EUR')

# Array type, used : EPIC or 450K
artype <- c('EPIC', 'EPIC', 'EPIC')

# Parameters to exclude CpGs
exclude <- c( 'control_probes',
              'noncpg_probes',
              'Sex',
              'MASK_mapping',
              'MASK_sub30_copy',
              'MASK_extBase',
              'MASK_typeINextBaseSwitch',
              'MASK_snp5_ethnic',
              'Unrel_450_EPIC_blood')



N <- c(255, 255, 255)
n <- c(NA, NA, NA)

# Minimum sample representation percentage required for CpGs
# Filter minimum percentage of missings for each CpG in cohort
# We need to define two parameters,
#  - colname_NforProbe:
#        Column name with Number of individuals per probe, this variable only needs
#           to be defined if you want to filter CpGs with low representation.
#         If defined value in colname_NforProbe not exists, no filter will be applied
#  - pcMissingSamples :
#        Máximum percent of missing samples allowed,

colname_NforProbe <- 'N_for_probe'
pcMissingSamples <- 0.9


########## ----------  END VARIABLES DEFINED BY USER  ----------  ##########



## ###################### ##
##  QC - Quality Control  ##
## ###################### ##

# Variable declaration to perform precision plot
medianSE <- numeric(length(files))
value_N <- numeric(length(files))

if(length(n) == length(N))
  value_n <- numeric(length(files))

cohort_label <- character(length(files))

# Prepare output folder for results (create if not exists)
if(!dir.exists(file.path(getwd(), results_folder )))
  suppressWarnings(dir.create(file.path(getwd(), results_folder)))

## Remove duplicates, Exclude CpGs and adjust data (BN and FDR)
for ( i in 1:length(files) )
{
  
  # Prepare output subfolder for cohort-model results (create if not exists)
  if(!dir.exists(file.path(getwd(), results_folder, prefixes[i] )))
    suppressWarnings(dir.create(file.path(getwd(), results_folder, prefixes[i])))
  
  # Creates an empty file to resume all data if an old file exist  removes
  # the file and creates a new one
  fResumeName <- paste0( file.path(getwd(), results_folder, prefixes[i]),"/",prefixes[i], "_descriptives.txt")
  if ( file.exists(fResumeName) ) {
    file.remove(fResumeName)
  }
  file.create(fResumeName)
  
  # Read data.
  cohort <- read.table(files[i], header = TRUE, as.is = TRUE)
  print(paste0("Cohort file : ",files[i]," - readed OK", sep = " "))
  
  # Remove rows with NA from data
  cohort <- clean_NA_from_data(cohort)
  
  # Descriptives - Before CpGs deletion #
  descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], before = TRUE)
  
  # Remove duplicates
  # cohort <- remove_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i], '/',prefixes[i],'_descriptives_duplic.txt'), paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
  
  test_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
  
  # Remove cpGs with low representation
  # first, we test if colname_NforProbe and pcMissingSampes are defined
  if( !exists("colname_NforProbe") ) { colname_NforProbe <- NULL }
  if( !exists("pcMissingSamples") ) { pcMissingSamples <- NULL }
  
  cohort <- filterLowRepresentedCpGsinCohort(cohort, colname_NforProbe, pcMissingSamples, N[i], fileresume = fResumeName )
  
  # Exclude CpGs not meet conditions
  if("MASK_snp5_ethnic" %in% exclude ){
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = ethnic[i], filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
  } else {
    #..# if( !is.null(exclude) && exclude!='') {
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = "", filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
    #..# }
  }
  
  # Descriptives - After CpGs deletion #
  descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], before = FALSE )
  
  # Adjust data by Bonferroni and FDR
  cohort <- adjust_data(cohort, "P_VAL", bn=TRUE, fdr=TRUE, fResumeName, N[i]  )
  
  # Write QC complete data to external file
  write_QCData(cohort, paste0(results_folder, '/',prefixes[i], '/',prefixes[i]))
  
  ## Visualization - Plots
  #. Problems in some workstations and servers.# rasterpdf::raster_pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'), res = 300)
  #..# Problems in some cases --> Get png plots : pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'))
  
  # Distribution plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_SE_plot.png'))
  plot_distribution(cohort$SE, main = paste('Standard Errors of', prefixes[i]), xlab = 'SE')
  dev.off()
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_pvals_plot.png'))
  plot_distribution(cohort$P_VAL, main = paste('p-values of', prefixes[i]), xlab = 'p-value')
  dev.off()
  
  # QQ plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_QQ_plot.png'))
  qqman::qq(cohort$P_VAL, main = sprintf('QQ plot of %s (lambda = %f)', prefixes[i], lambda=get_lambda(cohort,"P_VAL")))
  dev.off()
  
  # Volcano plot.
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_Volcano_plot.png'))
  plot_volcano(cohort, "BETA", "P_VAL", main =paste('Volcano plot of', prefixes[i]) )
  dev.off()
  
  # Add mandatory data for precisionplot
  medianSE[i] <-  median(cohort$SE)
  value_N[i] <- N[i]
  cohort_label[i] <-  prefixes[i]
  
  # if n is defined for dichotomic condition :
  if(length(n) == length(N))  value_n[i] <- n[i]
  
  # Store data for Beta Box-Plot
  if( i == 1)
    betas.data <- list()
  betas.data[[prefixes[i]]] <- cohort[,"BETA"]
  
}



###################################
##### syncytiotrophoblast
###################################


generalsyncytiotrophoblast <- read.csv("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/generalsyncytiotrophoblast.csv")
head(generalsyncytiotrophoblast)
#X      Estimate           SE          t         p      adjP
#1 cg14817997 -0.0001375936 0.0009683181 -0.1420954 0.8871279 0.9872763
#2 cg26928153 -0.0008820622 0.0006921729 -1.2743380 0.2038186 0.8429583
#3 cg16269199 -0.0009567483 0.0009491518 -1.0080034 0.3145030 0.8757481
#4 cg13869341  0.0000672576 0.0006005205  0.1119988 0.9109212 0.9899296
#5 cg14008030 -0.0003942486 0.0007960805 -0.4952371 0.6209015 0.9465795
#6 cg12045430  0.0001438201 0.0003816845  0.3768035 0.7066641 0.9614575



N<- rep(255, times= 811990)
N<- as.vector(N)
generalsyncytiotrophoblast1<- data.frame(probeID=generalsyncytiotrophoblast$X, N= N, BETA= generalsyncytiotrophoblast$Estimate, SE= generalsyncytiotrophoblast$SE, zval= generalsyncytiotrophoblast$t, P_VAL=generalsyncytiotrophoblast$p )
head(generalsyncytiotrophoblast1)
#probeID   N          BETA           SE       zval     P_VAL
#1 cg14817997 255 -0.0001375936 0.0009683181 -0.1420954 0.8871279
#2 cg26928153 255 -0.0008820622 0.0006921729 -1.2743380 0.2038186
#3 cg16269199 255 -0.0009567483 0.0009491518 -1.0080034 0.3145030
#4 cg13869341 255  0.0000672576 0.0006005205  0.1119988 0.9109212
#5 cg14008030 255 -0.0003942486 0.0007960805 -0.4952371 0.6209015
#6 cg12045430 255  0.0001438201 0.0003816845  0.3768035 0.7066641


write.table(generalsyncytiotrophoblast1, "generalsyncytiotrophoblast1.txt", col.names=TRUE)

###################non verbal

nonverbalsyncytiotrophoblast<- read.csv("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/nonverbalsyncytiotrophoblast.csv")
head(nonverbalsyncytiotrophoblast)
#       X      Estimate          SE          t         p      adjP
#1 cg14817997  0.0078774578 0.008724324  0.9029304 0.3674994 0.8900535
#2 cg26928153  0.0024323496 0.006250781  0.3891273 0.6975389 0.9603024
#3 cg16269199  0.0004183192 0.008556843  0.0488871 0.9610513 0.9958441
#4 cg13869341 -0.0019151362 0.005379655 -0.3559961 0.7221670 0.9643202
#5 cg14008030  0.0033923849 0.007101628  0.4776912 0.6333198 0.9493441
#6 cg12045430 -0.0017124591 0.003403142 -0.5031994 0.6153015 0.9461732

N<- rep(255, times= 811990)
N<- as.vector(N)
nonverbalsyncytiotrophoblast1<- data.frame(probeID=nonverbalsyncytiotrophoblast$X, N= N, BETA= nonverbalsyncytiotrophoblast$Estimate, SE= nonverbalsyncytiotrophoblast$SE, zval= nonverbalsyncytiotrophoblast$t, P_VAL=nonverbalsyncytiotrophoblast$p )
head(nonverbalsyncytiotrophoblast1)
#   probeID   N          BETA           SE        zval     P_VAL
#1 cg14817997 255 -7.058139e-04 0.0010338343 -0.68271469 0.4954683
#2 cg26928153 255 -3.481090e-05 0.0007407189 -0.04699610 0.9625568
#3 cg16269199 255  1.994502e-04 0.0010139877  0.19669881 0.8442354
#4 cg13869341 255 -2.385189e-05 0.0006374903 -0.03741529 0.9701861
#5 cg14008030 255 -5.990947e-04 0.0008415444 -0.71189908 0.4772425
#6 cg12045430 255  4.169252e-04 0.0004032731  1.03385324 0.3022813



write.table(nonverbalsyncytiotrophoblast1, "nonverbalsyncytiotrophoblast1.txt", col.names=TRUE)


######################verbal

verbalsyncytiotrophoblast<- read.csv("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/verbalsyncytiotrophoblast.csv")
head(verbalsyncytiotrophoblast)
#            X      Estimate           SE           t          p      adjP
#1 cg14817997  1.524792e-04 0.0009993216  0.15258274 0.87885996 0.9735446
#2 cg26928153 -1.452687e-03 0.0007084327 -2.05056487 0.04143321 0.4432414
#3 cg16269199 -1.496303e-03 0.0009746999 -1.53514196 0.12611156 0.5853344
#4 cg13869341  3.658387e-06 0.0006224065  0.00587781 0.99531526 0.9989955
#5 cg14008030  1.523682e-05 0.0008245719  0.01847846 0.98527304 0.9970341
#6 cg12045430 -3.378513e-05 0.0003948389 -0.08556686 0.93188454 0.9855186


N<- rep(255, times= 811990)
N<- as.vector(N)
verbalsyncytiotrophoblast1<- data.frame(probeID=verbalsyncytiotrophoblast$X, N= N, BETA= verbalsyncytiotrophoblast$Estimate, SE= verbalsyncytiotrophoblast$SE, zval= verbalsyncytiotrophoblast$t, P_VAL=verbalsyncytiotrophoblast$p )
head(verbalsyncytiotrophoblast1)
#probeID   N          BETA           SE        zval      P_VAL
#1 cg14817997 255  1.524792e-04 0.0009993216  0.15258274 0.87885996
#2 cg26928153 255 -1.452687e-03 0.0007084327 -2.05056487 0.04143321
#3 cg16269199 255 -1.496303e-03 0.0009746999 -1.53514196 0.12611156
#4 cg13869341 255  3.658387e-06 0.0006224065  0.00587781 0.99531526
#5 cg14008030 255  1.523682e-05 0.0008245719  0.01847846 0.98527304
#6 cg12045430 255 -3.378513e-05 0.0003948389 -0.08556686 0.93188454

write.table(verbalsyncytiotrophoblast1, "verbalsyncytiotrophoblast1.txt", col.names=TRUE)






library(EASIER)



########## ----------  VARIABLES DEFINED BY USER  ----------  ##########

# Set working directory to new folder
setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/")


# Files used in QC, needed in meta-analysis to plot ForestPlot (these files are the ones we have just saved)
files <- c('generalsyncytiotrophoblast1.txt',
           'nonverbalsyncytiotrophoblast1.txt',
           'verbalsyncytiotrophoblast1.txt')

# Result folder
results_folder <- 'QC_Resultssyncytiotrophoblast'

# Prefixes for each file
prefixes <- c('verbal', 'nonverbal', 'general')


# Exclude - MASK snp5
ethnic <- c('EUR', 'EUR', 'EUR')

# Array type, used : EPIC or 450K
artype <- c('EPIC', 'EPIC', 'EPIC')

# Parameters to exclude CpGs
exclude <- c( 'control_probes',
              'noncpg_probes',
              'Sex',
              'MASK_mapping',
              'MASK_sub30_copy',
              'MASK_extBase',
              'MASK_typeINextBaseSwitch',
              'MASK_snp5_ethnic',
              'Unrel_450_EPIC_blood')



N <- c(255, 255, 255)
n <- c(NA, NA, NA)

# Minimum sample representation percentage required for CpGs
# Filter minimum percentage of missings for each CpG in cohort
# We need to define two parameters,
#  - colname_NforProbe:
#        Column name with Number of individuals per probe, this variable only needs
#           to be defined if you want to filter CpGs with low representation.
#         If defined value in colname_NforProbe not exists, no filter will be applied
#  - pcMissingSamples :
#        Máximum percent of missing samples allowed,

colname_NforProbe <- 'N_for_probe'
pcMissingSamples <- 0.9


########## ----------  END VARIABLES DEFINED BY USER  ----------  ##########



## ###################### ##
##  QC - Quality Control  ##
## ###################### ##

# Variable declaration to perform precision plot
medianSE <- numeric(length(files))
value_N <- numeric(length(files))

if(length(n) == length(N))
  value_n <- numeric(length(files))

cohort_label <- character(length(files))

# Prepare output folder for results (create if not exists)
if(!dir.exists(file.path(getwd(), results_folder )))
  suppressWarnings(dir.create(file.path(getwd(), results_folder)))

## Remove duplicates, Exclude CpGs and adjust data (BN and FDR)
for ( i in 1:length(files) )
{
  
  # Prepare output subfolder for cohort-model results (create if not exists)
  if(!dir.exists(file.path(getwd(), results_folder, prefixes[i] )))
    suppressWarnings(dir.create(file.path(getwd(), results_folder, prefixes[i])))
  
  # Creates an empty file to resume all data if an old file exist  removes
  # the file and creates a new one
  fResumeName <- paste0( file.path(getwd(), results_folder, prefixes[i]),"/",prefixes[i], "_descriptives.txt")
  if ( file.exists(fResumeName) ) {
    file.remove(fResumeName)
  }
  file.create(fResumeName)
  
  # Read data.
  cohort <- read.table(files[i], header = TRUE, as.is = TRUE)
  print(paste0("Cohort file : ",files[i]," - readed OK", sep = " "))
  
  # Remove rows with NA from data
  cohort <- clean_NA_from_data(cohort)
  
  # Descriptives - Before CpGs deletion #
  descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], before = TRUE)
  
  # Remove duplicates
  # cohort <- remove_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i], '/',prefixes[i],'_descriptives_duplic.txt'), paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
  
  test_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
  
  # Remove cpGs with low representation
  # first, we test if colname_NforProbe and pcMissingSampes are defined
  if( !exists("colname_NforProbe") ) { colname_NforProbe <- NULL }
  if( !exists("pcMissingSamples") ) { pcMissingSamples <- NULL }
  
  cohort <- filterLowRepresentedCpGsinCohort(cohort, colname_NforProbe, pcMissingSamples, N[i], fileresume = fResumeName )
  
  # Exclude CpGs not meet conditions
  if("MASK_snp5_ethnic" %in% exclude ){
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = ethnic[i], filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
  } else {
    #..# if( !is.null(exclude) && exclude!='') {
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = "", filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
    #..# }
  }
  
  # Descriptives - After CpGs deletion #
  descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], before = FALSE )
  
  # Adjust data by Bonferroni and FDR
  cohort <- adjust_data(cohort, "P_VAL", bn=TRUE, fdr=TRUE, fResumeName, N[i]  )
  
  # Write QC complete data to external file
  write_QCData(cohort, paste0(results_folder, '/',prefixes[i], '/',prefixes[i]))
  
  ## Visualization - Plots
  #. Problems in some workstations and servers.# rasterpdf::raster_pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'), res = 300)
  #..# Problems in some cases --> Get png plots : pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'))
  
  # Distribution plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_SE_plot.png'))
  plot_distribution(cohort$SE, main = paste('Standard Errors of', prefixes[i]), xlab = 'SE')
  dev.off()
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_pvals_plot.png'))
  plot_distribution(cohort$P_VAL, main = paste('p-values of', prefixes[i]), xlab = 'p-value')
  dev.off()
  
  # QQ plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_QQ_plot.png'))
  qqman::qq(cohort$P_VAL, main = sprintf('QQ plot of %s (lambda = %f)', prefixes[i], lambda=get_lambda(cohort,"P_VAL")))
  dev.off()
  
  # Volcano plot.
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_Volcano_plot.png'))
  plot_volcano(cohort, "BETA", "P_VAL", main =paste('Volcano plot of', prefixes[i]) )
  dev.off()
  
  # Add mandatory data for precisionplot
  medianSE[i] <-  median(cohort$SE)
  value_N[i] <- N[i]
  cohort_label[i] <-  prefixes[i]
  
  # if n is defined for dichotomic condition :
  if(length(n) == length(N))  value_n[i] <- n[i]
  
  # Store data for Beta Box-Plot
  if( i == 1)
    betas.data <- list()
  betas.data[[prefixes[i]]] <- cohort[,"BETA"]
  
}


###################################
##### nRBC
###################################


generalnRBC <- read.csv("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/generalnRBC.csv")
head(generalnRBC)
#          X      Estimate         SE           t         p      adjP
#1 cg14817997  0.0066374864 0.03253899  0.20398565 0.8385435 0.9999979
#2 cg26928153  0.0166323213 0.02325951  0.71507612 0.4752810 0.9999979
#3 cg16269199 -0.0068996410 0.03189493 -0.21632406 0.8289252 0.9999979
#4 cg13869341  0.0005720549 0.02017966  0.02834809 0.9774089 0.9999979
#5 cg14008030  0.0307055863 0.02675118  1.14782173 0.2522244 0.9999979
#6 cg12045430 -0.0034973575 0.01282598 -0.27267761 0.7853436 0.9999979




N<- rep(255, times= 811990)
N<- as.vector(N)
generalnRBC1<- data.frame(probeID=generalnRBC$X, N= N, BETA= generalnRBC$Estimate, SE= generalnRBC$SE, zval= generalnRBC$t, P_VAL=generalnRBC$p )
head(generalnRBC1)
# probeID   N          BETA         SE        zval     P_VAL
#1 cg14817997 255  0.0066374864 0.03253899  0.20398565 0.8385435
#2 cg26928153 255  0.0166323213 0.02325951  0.71507612 0.4752810
#3 cg16269199 255 -0.0068996410 0.03189493 -0.21632406 0.8289252
#4 cg13869341 255  0.0005720549 0.02017966  0.02834809 0.9774089
#5 cg14008030 255  0.0307055863 0.02675118  1.14782173 0.2522244
#6 cg12045430 255 -0.0034973575 0.01282598 -0.27267761 0.7853436



write.table(generalnRBC1, "generalnRBC1.txt", col.names=TRUE)

###################non verbal

nonverbalnRBC<- read.csv("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/nonverbalnRBC.csv")
head(nonverbalnRBC)
#            X     Estimate         SE          t          p      adjP
#1 cg14817997 -0.031837624 0.03115262 -1.0219885 0.30785083 0.9999756
#2 cg26928153 -0.009741159 0.02232015 -0.4364289 0.66293143 0.9999756
#3 cg16269199 -0.020287224 0.03055458 -0.6639666 0.50737097 0.9999756
#4 cg13869341  0.005823732 0.01920955  0.3031685 0.76203311 0.9999756
#5 cg14008030  0.048648733 0.02535834  1.9184513 0.05628118 0.9999756
#6 cg12045430 -0.009169007 0.01215187 -0.7545350 0.45129327 0.9999756


N<- rep(255, times= 811990)
N<- as.vector(N)
nonverbalnRBC1<- data.frame(probeID=nonverbalnRBC$X, N= N, BETA= nonverbalnRBC$Estimate, SE= nonverbalnRBC$SE, zval= nonverbalnRBC$t, P_VAL=nonverbalnRBC$p )
head(nonverbalnRBC1)
#  probeID   N         BETA         SE       zval      P_VAL
#1 cg14817997 255 -0.031837624 0.03115262 -1.0219885 0.30785083
#2 cg26928153 255 -0.009741159 0.02232015 -0.4364289 0.66293143
#3 cg16269199 255 -0.020287224 0.03055458 -0.6639666 0.50737097
#4 cg13869341 255  0.005823732 0.01920955  0.3031685 0.76203311
#5 cg14008030 255  0.048648733 0.02535834  1.9184513 0.05628118
#6 cg12045430 255 -0.009169007 0.01215187 -0.7545350 0.45129327




write.table(nonverbalnRBC1, "nonverbalnRBC1.txt", col.names=TRUE)


######################verbal

verbalnRBC<- read.csv("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/verbalnRBC.csv")
head(verbalnRBC)
#            X      Estimate         SE           t         p     adjP
#1 cg14817997  0.0124720207 0.03013649  0.41385117 0.6793654 0.999988
#2 cg26928153  0.0304366222 0.02136417  1.42465761 0.1556000 0.999988
#3 cg16269199  0.0116020793 0.02939397  0.39470947 0.6934199 0.999988
#4 cg13869341 -0.0009087383 0.01876988 -0.04841472 0.9614274 0.999988
#5 cg14008030  0.0025787439 0.02486657  0.10370325 0.9174944 0.999988
#6 cg12045430  0.0008499437 0.01190714  0.07138104 0.9431560 0.999988



N<- rep(255, times= 811990)
N<- as.vector(N)
verbalnRBC1<- data.frame(probeID=verbalnRBC$X, N= N, BETA= verbalnRBC$Estimate, SE= verbalnRBC$SE, zval= verbalnRBC$t, P_VAL=verbalnRBC$p )
head(verbalnRBC1)
#probeID   N          BETA         SE        zval     P_VAL
#1 cg14817997 255  0.0124720207 0.03013649  0.41385117 0.6793654
#2 cg26928153 255  0.0304366222 0.02136417  1.42465761 0.1556000
#3 cg16269199 255  0.0116020793 0.02939397  0.39470947 0.6934199
#4 cg13869341 255 -0.0009087383 0.01876988 -0.04841472 0.9614274
#5 cg14008030 255  0.0025787439 0.02486657  0.10370325 0.9174944
#6 cg12045430 255  0.0008499437 0.01190714  0.07138104 0.9431560


write.table(verbalnRBC1, "verbalnRBC1.txt", col.names=TRUE)






library(EASIER)



########## ----------  VARIABLES DEFINED BY USER  ----------  ##########

# Set working directory to new folder
setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/")


# Files used in QC, needed in meta-analysis to plot ForestPlot (these files are the ones we have just saved)
files <- c('generalnRBC1.txt',
           'nonverbalnRBC1.txt',
           'verbalnRBC1.txt')

# Result folder
results_folder <- 'QC_ResultsnRBC'

# Prefixes for each file
prefixes <- c('verbal', 'nonverbal', 'general')


# Exclude - MASK snp5
ethnic <- c('EUR', 'EUR', 'EUR')

# Array type, used : EPIC or 450K
artype <- c('EPIC', 'EPIC', 'EPIC')

# Parameters to exclude CpGs
exclude <- c( 'control_probes',
              'noncpg_probes',
              'Sex',
              'MASK_mapping',
              'MASK_sub30_copy',
              'MASK_extBase',
              'MASK_typeINextBaseSwitch',
              'MASK_snp5_ethnic',
              'Unrel_450_EPIC_blood')



N <- c(255, 255, 255)
n <- c(NA, NA, NA)

# Minimum sample representation percentage required for CpGs
# Filter minimum percentage of missings for each CpG in cohort
# We need to define two parameters,
#  - colname_NforProbe:
#        Column name with Number of individuals per probe, this variable only needs
#           to be defined if you want to filter CpGs with low representation.
#         If defined value in colname_NforProbe not exists, no filter will be applied
#  - pcMissingSamples :
#        Máximum percent of missing samples allowed,

colname_NforProbe <- 'N_for_probe'
pcMissingSamples <- 0.9


########## ----------  END VARIABLES DEFINED BY USER  ----------  ##########



## ###################### ##
##  QC - Quality Control  ##
## ###################### ##

# Variable declaration to perform precision plot
medianSE <- numeric(length(files))
value_N <- numeric(length(files))

if(length(n) == length(N))
  value_n <- numeric(length(files))

cohort_label <- character(length(files))

# Prepare output folder for results (create if not exists)
if(!dir.exists(file.path(getwd(), results_folder )))
  suppressWarnings(dir.create(file.path(getwd(), results_folder)))

## Remove duplicates, Exclude CpGs and adjust data (BN and FDR)
for ( i in 1:length(files) )
{
  
  # Prepare output subfolder for cohort-model results (create if not exists)
  if(!dir.exists(file.path(getwd(), results_folder, prefixes[i] )))
    suppressWarnings(dir.create(file.path(getwd(), results_folder, prefixes[i])))
  
  # Creates an empty file to resume all data if an old file exist  removes
  # the file and creates a new one
  fResumeName <- paste0( file.path(getwd(), results_folder, prefixes[i]),"/",prefixes[i], "_descriptives.txt")
  if ( file.exists(fResumeName) ) {
    file.remove(fResumeName)
  }
  file.create(fResumeName)
  
  # Read data.
  cohort <- read.table(files[i], header = TRUE, as.is = TRUE)
  print(paste0("Cohort file : ",files[i]," - readed OK", sep = " "))
  
  # Remove rows with NA from data
  cohort <- clean_NA_from_data(cohort)
  
  # Descriptives - Before CpGs deletion #
  descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], before = TRUE)
  
  # Remove duplicates
  # cohort <- remove_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i], '/',prefixes[i],'_descriptives_duplic.txt'), paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
  
  test_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
  
  # Remove cpGs with low representation
  # first, we test if colname_NforProbe and pcMissingSampes are defined
  if( !exists("colname_NforProbe") ) { colname_NforProbe <- NULL }
  if( !exists("pcMissingSamples") ) { pcMissingSamples <- NULL }
  
  cohort <- filterLowRepresentedCpGsinCohort(cohort, colname_NforProbe, pcMissingSamples, N[i], fileresume = fResumeName )
  
  # Exclude CpGs not meet conditions
  if("MASK_snp5_ethnic" %in% exclude ){
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = ethnic[i], filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
  } else {
    #..# if( !is.null(exclude) && exclude!='') {
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = "", filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
    #..# }
  }
  
  # Descriptives - After CpGs deletion #
  descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], before = FALSE )
  
  # Adjust data by Bonferroni and FDR
  cohort <- adjust_data(cohort, "P_VAL", bn=TRUE, fdr=TRUE, fResumeName, N[i]  )
  
  # Write QC complete data to external file
  write_QCData(cohort, paste0(results_folder, '/',prefixes[i], '/',prefixes[i]))
  
  ## Visualization - Plots
  #. Problems in some workstations and servers.# rasterpdf::raster_pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'), res = 300)
  #..# Problems in some cases --> Get png plots : pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'))
  
  # Distribution plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_SE_plot.png'))
  plot_distribution(cohort$SE, main = paste('Standard Errors of', prefixes[i]), xlab = 'SE')
  dev.off()
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_pvals_plot.png'))
  plot_distribution(cohort$P_VAL, main = paste('p-values of', prefixes[i]), xlab = 'p-value')
  dev.off()
  
  # QQ plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_QQ_plot.png'))
  qqman::qq(cohort$P_VAL, main = sprintf('QQ plot of %s (lambda = %f)', prefixes[i], lambda=get_lambda(cohort,"P_VAL")))
  dev.off()
  
  # Volcano plot.
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_Volcano_plot.png'))
  plot_volcano(cohort, "BETA", "P_VAL", main =paste('Volcano plot of', prefixes[i]) )
  dev.off()
  
  # Add mandatory data for precisionplot
  medianSE[i] <-  median(cohort$SE)
  value_N[i] <- N[i]
  cohort_label[i] <-  prefixes[i]
  
  # if n is defined for dichotomic condition :
  if(length(n) == length(N))  value_n[i] <- n[i]
  
  # Store data for Beta Box-Plot
  if( i == 1)
    betas.data <- list()
  betas.data[[prefixes[i]]] <- cohort[,"BETA"]
  
}


###################################
##### hofbauer
###################################


generalhofbauer <- read.csv("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/generalhofbauer.csv")
head(generalhofbauer)
#         X      Estimate         SE           t          p      adjP
#1 cg14817997  0.0007295108 0.01801145  0.04050262 0.96772725 0.9999886
#2 cg26928153 -0.0271397280 0.01287494 -2.10795010 0.03610843 0.9999886
#3 cg16269199 -0.0136184303 0.01765494 -0.77136655 0.44127462 0.9999886
#4 cg13869341 -0.0089165242 0.01117013 -0.79824681 0.42554362 0.9999886
#5 cg14008030  0.0064704613 0.01480770  0.43696607 0.66254239 0.9999886
#6 cg12045430 -0.0039976786 0.00709962 -0.56308344 0.57392184 0.9999886



N<- rep(255, times= 811990)
N<- as.vector(N)
generalhofbauer1<- data.frame(probeID=generalhofbauer$X, N= N, BETA= generalhofbauer$Estimate, SE= generalhofbauer$SE, zval= generalhofbauer$t, P_VAL=generalhofbauer$p )
head(generalhofbauer1)
# probeID   N          BETA         SE        zval      P_VAL
#1 cg14817997 255  0.0007295108 0.01801145  0.04050262 0.96772725
#2 cg26928153 255 -0.0271397280 0.01287494 -2.10795010 0.03610843
#3 cg16269199 255 -0.0136184303 0.01765494 -0.77136655 0.44127462
#4 cg13869341 255 -0.0089165242 0.01117013 -0.79824681 0.42554362
#5 cg14008030 255  0.0064704613 0.01480770  0.43696607 0.66254239
#6 cg12045430 255 -0.0039976786 0.00709962 -0.56308344 0.57392184




write.table(generalhofbauer1, "generalhofbauer1.txt", col.names=TRUE)

###################non verbal

nonverbalhofbauer<- read.csv("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/nonverbalhofbauer.csv")
head(nonverbalhofbauer)
#           X      Estimate          SE           t          p      adjP
#1 cg14817997  0.0165291592 0.017924738  0.92214229 0.35741233 0.9637142
#2 cg26928153 -0.0238431618 0.012842670 -1.85655807 0.06464176 0.9567057
#3 cg16269199 -0.0160761000 0.017580636 -0.91442087 0.36144512 0.9639961
#4 cg13869341 -0.0132479101 0.011052880 -1.19859348 0.23190895 0.9572796
#5 cg14008030 -0.0010540894 0.014590794 -0.07224346 0.94247040 0.9979123
#6 cg12045430 -0.0006804047 0.006991995 -0.09731195 0.92256268 0.9972445


N<- rep(255, times= 811990)
N<- as.vector(N)
nonverbalhofbauer1<- data.frame(probeID=nonverbalhofbauer$X, N= N, BETA= nonverbalhofbauer$Estimate, SE= nonverbalhofbauer$SE, zval= nonverbalhofbauer$t, P_VAL=nonverbalhofbauer$p )
head(nonverbalhofbauer1)
#  probeID   N          BETA          SE        zval      P_VAL
#1 cg14817997 255  0.0165291592 0.017924738  0.92214229 0.35741233
#2 cg26928153 255 -0.0238431618 0.012842670 -1.85655807 0.06464176
#3 cg16269199 255 -0.0160761000 0.017580636 -0.91442087 0.36144512
#4 cg13869341 255 -0.0132479101 0.011052880 -1.19859348 0.23190895
#5 cg14008030 255 -0.0010540894 0.014590794 -0.07224346 0.94247040
#6 cg12045430 255 -0.0006804047 0.006991995 -0.09731195 0.92256268



write.table(nonverbalhofbauer1, "nonverbalhofbauer1.txt", col.names=TRUE)


######################verbal

verbalhofbauer<- read.csv("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/verbalhofbauer.csv")
head(verbalhofbauer)
#      X     Estimate          SE          t          p      adjP
# 1 cg14817997 -0.010298442 0.017943848 -0.5739261 0.56657381 0.9999953
#2 cg26928153 -0.026216082 0.012720638 -2.0609093 0.04042663 0.9999953
#3 cg16269199 -0.014242219 0.017501740 -0.8137601 0.41661678 0.9999953
#4 cg13869341 -0.003273080 0.011175949 -0.2928682 0.76988467 0.9999953
#5 cg14008030  0.011653346 0.014806036  0.7870673 0.43204578 0.9999953
#6 cg12045430 -0.006001018 0.007089739 -0.8464370 0.39818125 0.9999953




N<- rep(255, times= 811990)
N<- as.vector(N)
verbalhofbauer1<- data.frame(probeID=verbalhofbauer$X, N= N, BETA= verbalhofbauer$Estimate, SE= verbalhofbauer$SE, zval= verbalhofbauer$t, P_VAL=verbalhofbauer$p )
head(verbalhofbauer1)
#probeID   N         BETA          SE       zval      P_VAL
#1 cg14817997 255 -0.010298442 0.017943848 -0.5739261 0.56657381
#2 cg26928153 255 -0.026216082 0.012720638 -2.0609093 0.04042663
#3 cg16269199 255 -0.014242219 0.017501740 -0.8137601 0.41661678
#4 cg13869341 255 -0.003273080 0.011175949 -0.2928682 0.76988467
#5 cg14008030 255  0.011653346 0.014806036  0.7870673 0.43204578
#6 cg12045430 255 -0.006001018 0.007089739 -0.8464370 0.39818125



write.table(verbalhofbauer1, "verbalhofbauer1.txt", col.names=TRUE)






library(EASIER)



########## ----------  VARIABLES DEFINED BY USER  ----------  ##########

# Set working directory to new folder
setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/")


# Files used in QC, needed in meta-analysis to plot ForestPlot (these files are the ones we have just saved)
files <- c('generalhofbauer1.txt',
           'nonverbalhofbauer1.txt',
           'verbalhofbauer1.txt')

# Result folder
results_folder <- 'QC_Resultshofbauer'

# Prefixes for each file
prefixes <- c('verbal', 'nonverbal', 'general')


# Exclude - MASK snp5
ethnic <- c('EUR', 'EUR', 'EUR')

# Array type, used : EPIC or 450K
artype <- c('EPIC', 'EPIC', 'EPIC')

# Parameters to exclude CpGs
exclude <- c( 'control_probes',
              'noncpg_probes',
              'Sex',
              'MASK_mapping',
              'MASK_sub30_copy',
              'MASK_extBase',
              'MASK_typeINextBaseSwitch',
              'MASK_snp5_ethnic',
              'Unrel_450_EPIC_blood')



N <- c(255, 255, 255)
n <- c(NA, NA, NA)

# Minimum sample representation percentage required for CpGs
# Filter minimum percentage of missings for each CpG in cohort
# We need to define two parameters,
#  - colname_NforProbe:
#        Column name with Number of individuals per probe, this variable only needs
#           to be defined if you want to filter CpGs with low representation.
#         If defined value in colname_NforProbe not exists, no filter will be applied
#  - pcMissingSamples :
#        Máximum percent of missing samples allowed,

colname_NforProbe <- 'N_for_probe'
pcMissingSamples <- 0.9


########## ----------  END VARIABLES DEFINED BY USER  ----------  ##########



## ###################### ##
##  QC - Quality Control  ##
## ###################### ##

# Variable declaration to perform precision plot
medianSE <- numeric(length(files))
value_N <- numeric(length(files))

if(length(n) == length(N))
  value_n <- numeric(length(files))

cohort_label <- character(length(files))

# Prepare output folder for results (create if not exists)
if(!dir.exists(file.path(getwd(), results_folder )))
  suppressWarnings(dir.create(file.path(getwd(), results_folder)))

## Remove duplicates, Exclude CpGs and adjust data (BN and FDR)
for ( i in 1:length(files) )
{
  
  # Prepare output subfolder for cohort-model results (create if not exists)
  if(!dir.exists(file.path(getwd(), results_folder, prefixes[i] )))
    suppressWarnings(dir.create(file.path(getwd(), results_folder, prefixes[i])))
  
  # Creates an empty file to resume all data if an old file exist  removes
  # the file and creates a new one
  fResumeName <- paste0( file.path(getwd(), results_folder, prefixes[i]),"/",prefixes[i], "_descriptives.txt")
  if ( file.exists(fResumeName) ) {
    file.remove(fResumeName)
  }
  file.create(fResumeName)
  
  # Read data.
  cohort <- read.table(files[i], header = TRUE, as.is = TRUE)
  print(paste0("Cohort file : ",files[i]," - readed OK", sep = " "))
  
  # Remove rows with NA from data
  cohort <- clean_NA_from_data(cohort)
  
  # Descriptives - Before CpGs deletion #
  descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], before = TRUE)
  
  # Remove duplicates
  # cohort <- remove_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i], '/',prefixes[i],'_descriptives_duplic.txt'), paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
  
  test_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
  
  # Remove cpGs with low representation
  # first, we test if colname_NforProbe and pcMissingSampes are defined
  if( !exists("colname_NforProbe") ) { colname_NforProbe <- NULL }
  if( !exists("pcMissingSamples") ) { pcMissingSamples <- NULL }
  
  cohort <- filterLowRepresentedCpGsinCohort(cohort, colname_NforProbe, pcMissingSamples, N[i], fileresume = fResumeName )
  
  # Exclude CpGs not meet conditions
  if("MASK_snp5_ethnic" %in% exclude ){
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = ethnic[i], filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
  } else {
    #..# if( !is.null(exclude) && exclude!='') {
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = "", filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
    #..# }
  }
  
  # Descriptives - After CpGs deletion #
  descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], before = FALSE )
  
  # Adjust data by Bonferroni and FDR
  cohort <- adjust_data(cohort, "P_VAL", bn=TRUE, fdr=TRUE, fResumeName, N[i]  )
  
  # Write QC complete data to external file
  write_QCData(cohort, paste0(results_folder, '/',prefixes[i], '/',prefixes[i]))
  
  ## Visualization - Plots
  #. Problems in some workstations and servers.# rasterpdf::raster_pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'), res = 300)
  #..# Problems in some cases --> Get png plots : pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'))
  
  # Distribution plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_SE_plot.png'))
  plot_distribution(cohort$SE, main = paste('Standard Errors of', prefixes[i]), xlab = 'SE')
  dev.off()
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_pvals_plot.png'))
  plot_distribution(cohort$P_VAL, main = paste('p-values of', prefixes[i]), xlab = 'p-value')
  dev.off()
  
  # QQ plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_QQ_plot.png'))
  qqman::qq(cohort$P_VAL, main = sprintf('QQ plot of %s (lambda = %f)', prefixes[i], lambda=get_lambda(cohort,"P_VAL")))
  dev.off()
  
  # Volcano plot.
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_Volcano_plot.png'))
  plot_volcano(cohort, "BETA", "P_VAL", main =paste('Volcano plot of', prefixes[i]) )
  dev.off()
  
  # Add mandatory data for precisionplot
  medianSE[i] <-  median(cohort$SE)
  value_N[i] <- N[i]
  cohort_label[i] <-  prefixes[i]
  
  # if n is defined for dichotomic condition :
  if(length(n) == length(N))  value_n[i] <- n[i]
  
  # Store data for Beta Box-Plot
  if( i == 1)
    betas.data <- list()
  betas.data[[prefixes[i]]] <- cohort[,"BETA"]
  
}




###################################
##### trophoblast
###################################


generaltrophoblast <- read.csv("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/generaltrophoblast.csv")
head(generaltrophoblast)
#         X     Estimate          SE          t          p      adjP
#1 cg14817997  0.001256973 0.005415477  0.2321076 0.81665905 0.9900776
#2 cg26928153 -0.004801480 0.003871090 -1.2403432 0.21610142 0.9590119
#3 cg16269199  0.005893166 0.005308286  1.1101825 0.26806975 0.9594701
#4 cg13869341  0.005733363 0.003358509  1.7071155 0.08913848 0.9590119
#5 cg14008030 -0.002708689 0.004452210 -0.6083921 0.54352261 0.9723143
#6 cg12045430 -0.002893755 0.002134633 -1.3556216 0.17653777 0.9590119


N<- rep(255, times= 811990)
N<- as.vector(N)
generaltrophoblast1<- data.frame(probeID=generaltrophoblast$X, N= N, BETA= generaltrophoblast$Estimate, SE= generaltrophoblast$SE, zval= generaltrophoblast$t, P_VAL=generaltrophoblast$p )
head(generaltrophoblast1)
# probeID   N         BETA          SE       zval      P_VAL
#1 cg14817997 255  0.001256973 0.005415477  0.2321076 0.81665905
#2 cg26928153 255 -0.004801480 0.003871090 -1.2403432 0.21610142
#3 cg16269199 255  0.005893166 0.005308286  1.1101825 0.26806975
#4 cg13869341 255  0.005733363 0.003358509  1.7071155 0.08913848
#5 cg14008030 255 -0.002708689 0.004452210 -0.6083921 0.54352261
#6 cg12045430 255 -0.002893755 0.002134633 -1.3556216 0.17653777



write.table(generaltrophoblast1, "generaltrophoblast1.txt", col.names=TRUE)

###################non verbal

nonverbaltrophoblast<- read.csv("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/nonverbaltrophoblast.csv")
head(nonverbaltrophoblast)
#          X     Estimate          SE          t         p      adjP
#1 cg14817997 -0.000631565 0.005452473 -0.1158309 0.9078867 0.9851867
#2 cg26928153  0.002005661 0.003906573  0.5134067 0.6081554 0.9179930
#3 cg16269199  0.002668448 0.005347801  0.4989804 0.6182660 0.9210808
#4 cg13869341  0.004021281 0.003362143  1.1960473 0.2328991 0.7598056
#5 cg14008030 -0.002574219 0.004438330 -0.5799972 0.5624794 0.9048072
#6 cg12045430 -0.001558430 0.002126874 -0.7327324 0.4644613 0.8720458



N<- rep(255, times= 811990)
N<- as.vector(N)
nonverbaltrophoblast1<- data.frame(probeID=nonverbaltrophoblast$X, N= N, BETA= nonverbaltrophoblast$Estimate, SE= nonverbaltrophoblast$SE, zval= nonverbaltrophoblast$t, P_VAL=nonverbaltrophoblast$p )
head(nonverbaltrophoblast1)
#  probeID   N         BETA          SE       zval     P_VAL
#1 cg14817997 255 -0.000631565 0.005452473 -0.1158309 0.9078867
#2 cg26928153 255  0.002005661 0.003906573  0.5134067 0.6081554
#3 cg16269199 255  0.002668448 0.005347801  0.4989804 0.6182660
#4 cg13869341 255  0.004021281 0.003362143  1.1960473 0.2328991
#5 cg14008030 255 -0.002574219 0.004438330 -0.5799972 0.5624794
#6 cg12045430 255 -0.001558430 0.002126874 -0.7327324 0.4644613




write.table(nonverbaltrophoblast1, "nonverbaltrophoblast1.txt", col.names=TRUE)


######################verbal

verbaltrophoblast<- read.csv("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/verbaltrophoblast.csv")
head(verbaltrophoblast)
#     X     Estimate          SE          t          p      adjP
#1 cg14817997  0.003110506 0.005619938  0.5534770 0.58046980 0.9999909
#2 cg26928153 -0.009377692 0.003984050 -2.3538088 0.01941738 0.9999909
#3 cg16269199  0.007419973 0.005481471  1.3536463 0.17716653 0.9999909
#4 cg13869341  0.004687447 0.003500260  1.3391712 0.18182549 0.9999909
#5 cg14008030 -0.003457358 0.004637188 -0.7455721 0.45668066 0.9999909
#6 cg12045430 -0.002866823 0.002220476 -1.2910845 0.19795934 0.9999909



N<- rep(255, times= 811990)
N<- as.vector(N)
verbaltrophoblast1<- data.frame(probeID=verbaltrophoblast$X, N= N, BETA= verbaltrophoblast$Estimate, SE= verbaltrophoblast$SE, zval= verbaltrophoblast$t, P_VAL=verbaltrophoblast$p )
head(verbaltrophoblast1)
# probeID   N         BETA          SE       zval      P_VAL
# 1 cg14817997 255  0.003110506 0.005619938  0.5534770 0.58046980
#2 cg26928153 255 -0.009377692 0.003984050 -2.3538088 0.01941738
#3 cg16269199 255  0.007419973 0.005481471  1.3536463 0.17716653
#4 cg13869341 255  0.004687447 0.003500260  1.3391712 0.18182549
#5 cg14008030 255 -0.003457358 0.004637188 -0.7455721 0.45668066
#6 cg12045430 255 -0.002866823 0.002220476 -1.2910845 0.19795934




write.table(verbaltrophoblast1, "verbaltrophoblast1.txt", col.names=TRUE)






library(EASIER)



########## ----------  VARIABLES DEFINED BY USER  ----------  ##########

# Set working directory to new folder
setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/")


# Files used in QC, needed in meta-analysis to plot ForestPlot (these files are the ones we have just saved)
files <- c('generaltrophoblast1.txt',
           'nonverbaltrophoblast1.txt',
           'verbaltrophoblast1.txt')

# Result folder
results_folder <- 'QC_Resultstrophoblast'

# Prefixes for each file
prefixes <- c('verbal', 'nonverbal', 'general')


# Exclude - MASK snp5
ethnic <- c('EUR', 'EUR', 'EUR')

# Array type, used : EPIC or 450K
artype <- c('EPIC', 'EPIC', 'EPIC')

# Parameters to exclude CpGs
exclude <- c( 'control_probes',
              'noncpg_probes',
              'Sex',
              'MASK_mapping',
              'MASK_sub30_copy',
              'MASK_extBase',
              'MASK_typeINextBaseSwitch',
              'MASK_snp5_ethnic',
              'Unrel_450_EPIC_blood')



N <- c(255, 255, 255)
n <- c(NA, NA, NA)

# Minimum sample representation percentage required for CpGs
# Filter minimum percentage of missings for each CpG in cohort
# We need to define two parameters,
#  - colname_NforProbe:
#        Column name with Number of individuals per probe, this variable only needs
#           to be defined if you want to filter CpGs with low representation.
#         If defined value in colname_NforProbe not exists, no filter will be applied
#  - pcMissingSamples :
#        Máximum percent of missing samples allowed,

colname_NforProbe <- 'N_for_probe'
pcMissingSamples <- 0.9


########## ----------  END VARIABLES DEFINED BY USER  ----------  ##########



## ###################### ##
##  QC - Quality Control  ##
## ###################### ##

# Variable declaration to perform precision plot
medianSE <- numeric(length(files))
value_N <- numeric(length(files))

if(length(n) == length(N))
  value_n <- numeric(length(files))

cohort_label <- character(length(files))

# Prepare output folder for results (create if not exists)
if(!dir.exists(file.path(getwd(), results_folder )))
  suppressWarnings(dir.create(file.path(getwd(), results_folder)))

## Remove duplicates, Exclude CpGs and adjust data (BN and FDR)
for ( i in 1:length(files) )
{
  
  # Prepare output subfolder for cohort-model results (create if not exists)
  if(!dir.exists(file.path(getwd(), results_folder, prefixes[i] )))
    suppressWarnings(dir.create(file.path(getwd(), results_folder, prefixes[i])))
  
  # Creates an empty file to resume all data if an old file exist  removes
  # the file and creates a new one
  fResumeName <- paste0( file.path(getwd(), results_folder, prefixes[i]),"/",prefixes[i], "_descriptives.txt")
  if ( file.exists(fResumeName) ) {
    file.remove(fResumeName)
  }
  file.create(fResumeName)
  
  # Read data.
  cohort <- read.table(files[i], header = TRUE, as.is = TRUE)
  print(paste0("Cohort file : ",files[i]," - readed OK", sep = " "))
  
  # Remove rows with NA from data
  cohort <- clean_NA_from_data(cohort)
  
  # Descriptives - Before CpGs deletion #
  descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], before = TRUE)
  
  # Remove duplicates
  # cohort <- remove_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i], '/',prefixes[i],'_descriptives_duplic.txt'), paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
  
  test_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
  
  # Remove cpGs with low representation
  # first, we test if colname_NforProbe and pcMissingSampes are defined
  if( !exists("colname_NforProbe") ) { colname_NforProbe <- NULL }
  if( !exists("pcMissingSamples") ) { pcMissingSamples <- NULL }
  
  cohort <- filterLowRepresentedCpGsinCohort(cohort, colname_NforProbe, pcMissingSamples, N[i], fileresume = fResumeName )
  
  # Exclude CpGs not meet conditions
  if("MASK_snp5_ethnic" %in% exclude ){
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = ethnic[i], filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
  } else {
    #..# if( !is.null(exclude) && exclude!='') {
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = "", filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
    #..# }
  }
  
  # Descriptives - After CpGs deletion #
  descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], before = FALSE )
  
  # Adjust data by Bonferroni and FDR
  cohort <- adjust_data(cohort, "P_VAL", bn=TRUE, fdr=TRUE, fResumeName, N[i]  )
  
  # Write QC complete data to external file
  write_QCData(cohort, paste0(results_folder, '/',prefixes[i], '/',prefixes[i]))
  
  ## Visualization - Plots
  #. Problems in some workstations and servers.# rasterpdf::raster_pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'), res = 300)
  #..# Problems in some cases --> Get png plots : pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'))
  
  # Distribution plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_SE_plot.png'))
  plot_distribution(cohort$SE, main = paste('Standard Errors of', prefixes[i]), xlab = 'SE')
  dev.off()
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_pvals_plot.png'))
  plot_distribution(cohort$P_VAL, main = paste('p-values of', prefixes[i]), xlab = 'p-value')
  dev.off()
  
  # QQ plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_QQ_plot.png'))
  qqman::qq(cohort$P_VAL, main = sprintf('QQ plot of %s (lambda = %f)', prefixes[i], lambda=get_lambda(cohort,"P_VAL")))
  dev.off()
  
  # Volcano plot.
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_Volcano_plot.png'))
  plot_volcano(cohort, "BETA", "P_VAL", main =paste('Volcano plot of', prefixes[i]) )
  dev.off()
  
  # Add mandatory data for precisionplot
  medianSE[i] <-  median(cohort$SE)
  value_N[i] <- N[i]
  cohort_label[i] <-  prefixes[i]
  
  # if n is defined for dichotomic condition :
  if(length(n) == length(N))  value_n[i] <- n[i]
  
  # Store data for Beta Box-Plot
  if( i == 1)
    betas.data <- list()
  betas.data[[prefixes[i]]] <- cohort[,"BETA"]
  
}


