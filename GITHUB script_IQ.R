##########################################################################################
#### PACE placenta IQ CODE                                                            ####
#### Marta Cosin (marta.cosin@isglobal.org)                                           ####
#### 20220322                                                                         ####
#### Laia Diez Ahijado                                                                ####
##########################################################################################

#########################################################################################
# The following R code will allow you to complete all the EWAS requested in the PACE IQ and placenta DNA methylation analysis plan.
# The code also produces files summarising the variables included in the EWAS.
# You shouldn't have to rewrite or add to the following code, unless otherwise stated.
# There are just two inputs required for this analysis:
# 1) phenodataframe: a dataframe containing all the "phenotype" data needed for this project. 
#    Each row is a sample (individual) and each column is a different variable: 
#    - Outcomes: pdp_perm_m100sd15" (Perceptive performance IQ estimated using the MSCA, continuous), "pdp_verm_m100sd15" #(Verbal IQ estimated using the MSCA, continuous),
#                "pdp_genm_m100sd15" (general cognitive score estimated using the MSCA, continuous) 
#    - Main covariates: "Basename", "Sex","edMcCarthy","estudios3c", "edadm", "parity2c", "Smoke", 
#                       "sges", "Meanlog2oddsContamination"(contamination score), "cohort" (mother-child from Sabadell, Guipuzcoa or Valencia)
#                       
# 2) Betasnooutliers: a matrix of methylation illumina beta values. Each column is a sample and each
#row is a probe on the array (450k or EPIC). 
#####################################################################################################################################################################

##########################################################################################
### Go to the directory where to save results (change to your own working directory)
setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/") 

##########################################################################################
### Load required package (PACEanalysis package)
#if this is not already installed, you will have to install it as the first step --> follow the instructions from here: https://www.epicenteredresearch.com/pace/birthsize/

### Load required package (for correlations)
library(DescTools)
#install.packages("psych")   #crec que es: install.packages("psych")
library(psych)# for correlations
library(magrittr) #for piping
#install.packages("tidyverse")
library(tidyverse)
library(PACEanalysis)

### Data loading and preprocessing. 
#If your data have been preprocessed according to PACEanalysis package instructions, then proceed to the next step (Data preparation). 
#Otherwise, please follow the preprocessing steps found in https://www.epicenteredresearch.com/pace/birthsize/

##########################################################################################

### Data Preparation

#GETTING STARTED WITH YOUR PHENOTYPE DATA
#Check your data

#If you encounter any issues, please check out our troubleshooting guide to see if there is any guidance that may help: https://www.epicenteredresearch.com/pace/troubleshooting/
#If you closed prior R session (step 1), you can load list of pre-processed objects that is automatically saved by the preprocessingofData function, e.g.

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/EPIC_placenta/QC_PACEanalysis_20220508_slide/paceanalysis_017_3/INMA_20220426_Output/INMA_20220426_PreprocessedBetas_nooutliers.RData")

#First goal in this step is to create a phenodataframe that contains all the necessary variables and only complete cases. You can get the phenodata from the preprocessed data object or load your phenodataframe previously prepared:
load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/db/phenodataframe07062022.Rdata")

#per comprobar que tot estigui be.mirar que el sex sigui character.
str(phenodataframe)
  
#rownames(phenodataframe)<-as.character(phenodataframe$Basename)

# load celltypes
Celltypes<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/db/celltypes.txt", header=TRUE)
Celltypes<-as.matrix(Celltypes)
# Match celltypes with phenodataframe: 

phenodataframe<-phenodataframe[match(rownames(Celltypes),phenodataframe$Basename),]
all(phenodataframe$Basename == rownames(Celltypes))
#TRUE


# Calculate descriptives for your outcome variables, and then transform and standardize them as follows (if you haven't done it yet): 

mean_verbal <- mean(phenodataframe$pdp_verm_m100sd15) 
#100.4699
sd_verbal <-  sd(phenodataframe$pdp_verm_m100sd15)
#15.10878
min_verbal <- min(phenodataframe$pdp_verm_m100sd15)
#53.92244
max_verbal <- max(phenodataframe$pdp_verm_m100sd15)
#140.8315
mean_ppIQ <- mean(phenodataframe$pdp_perm_m100sd15) 
#100.6825
sd_ppIQ<-  sd(phenodataframe$pdp_perm_m100sd15)
#14.99538
min_ppIQ <- min(phenodataframe$pdp_perm_m100sd15)
#52.82082
max_ppIQ <- max(phenodataframe$pdp_perm_m100sd15)
#129.7736
mean_general <- mean(phenodataframe$pdp_genm_m100sd15) 
#100.8223
sd_general <-  sd(phenodataframe$pdp_genm_m100sd15)
#15.12766
min_general <- min(phenodataframe$pdp_genm_m100sd15)
#51.9369
max_general <- max(phenodataframe$pdp_genm_m100sd15)
#136.1481

Tabledescriptives <- rbind(mean_verbal, sd_verbal, min_verbal, max_verbal, mean_ppIQ,
                           sd_ppIQ, min_ppIQ, max_ppIQ, mean_general, sd_general, 
                           min_general, max_general)
Tabledescriptives <- round(Tabledescriptives,2)
Tabledescriptives 
#mean_verbal  100.47
#sd_verbal     15.11
#min_verbal    53.92
#max_verbal   140.83
#mean_ppIQ    100.68
#sd_ppIQ       15.00
#min_ppIQ      52.82
#max_ppIQ     129.77
#mean_general 100.82
#sd_general    15.13
#min_general   51.94
#max_general  136.15



#BiocManager::install("FDb.InfiniumMethylation.hg19")
library(FDb.InfiniumMethylation.hg19)


# Save results
setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/")
cohort<- "INMA" #CHANGE IT TO YOUR COHORT NAME
write.csv(Tabledescriptives, file=paste0("PACE_Pla_IQ", 
                                         cohort,"_Descriptives.Outcomes_transformed_",
                                         format(Sys.Date(), "%d%m%Y"), ".csv"))
# optional histograms before transforming outcomes

pdf(file=paste0(cohort,"_VERBAL_Histogram_", format(Sys.Date(), "%d%m%Y"),".pdf"))
hist(phenodataframe$pdp_verm_m100sd15)
dev.off()

pdf(file=paste0(cohort,"_ppIQ_Histogram_", format(Sys.Date(), "%d%m%Y"),".pdf"))
hist(phenodataframe$pdp_perm_m100sd15)
dev.off()

pdf(file=paste0(cohort,"_general_Histogram_", format(Sys.Date(), "%d%m%Y"),".pdf"))
hist(phenodataframe$pdp_genm_m100sd15)
dev.off()


# You can see the amount of complete cases
#please make sure your phenodataframe only includes complete cases

sum(complete.cases(phenodataframe)==TRUE)
# 228
dim(phenodataframe)
# 255 22



# does the amount look right? Remember that some samples dropped in the methylation QC.

# you can check that the amount does not change in subset of only necessary variables 
#(no need if you know you don't have unnecessary or optional variables such as maternal_IQ with NAs)

pheno_subset <- subset(phenodataframe, select=c(pdp_perm_m100sd15, pdp_verm_m100sd15, pdp_genm_m100sd15, Basename, Sex, sges, cohort, edadm, Smoke, parity2c,
                                                estudios3c, edMcCarthy, Meanlog2oddsContamination)) 


sum(complete.cases(pheno_subset)==TRUE)		
#255

#############################

# check that the phenodataframe has Basename as rownames: 
rownames(phenodataframe)<-as.character(phenodataframe$Basename)

#Make sure all categorical adjustment variables are coded as factors or characters

rownames(phenodataframe)<-as.character(phenodataframe$Basename)
phenodataframe$Sex<-as.character(phenodataframe$Sex)
phenodataframe$estudios3c<-as.factor(phenodataframe$estudios3c)
phenodataframe$Smoke<-as.factor(phenodataframe$Smoke)
phenodataframe$maternal_IQ<-as.numeric(as.character(phenodataframe$maternal_IQ))
phenodataframe$parity2c<-as.factor(phenodataframe$parity2c)
phenodataframe$cesarean<-as.factor(phenodataframe$cesarean)
phenodataframe$preterm<-as.factor(phenodataframe$preterm)
phenodataframe$labor_iniciation<-as.factor(phenodataframe$labor_iniciation)
phenodataframe$PregnancyComplications<-as.factor(phenodataframe$PregnancyComplications)
phenodataframe$cohort<-as.factor(phenodataframe$cohort)
phenodataframe$sga<-as.factor(phenodataframe$sga)
phenodataframe$Sex<-as.character(phenodataframe$Sex)

save(phenodataframe,file="phenofinal_IQ.Rdata")

######################################################
######################################################


### Calculate CORRELATIONS between main variables and save them: 

#For data sets with continuous, polytomous and dichotomous variables.

#Select relevant variables *NoTE: verbal, ppIQ and general scores must be the last columns

relevant_variables <- c("edadm","estudios3c","Smoke",
                        "edMcCarthy", "parity2c", "sges", "cohort", "pdp_perm_m100sd15", "pdp_verm_m100sd15", "pdp_genm_m100sd15", "maternal_IQ", "PregnancyComplications", "preterm", "cesarean", "labor_iniciation", "sga") 


column_positions <- match(relevant_variables, colnames(phenodataframe))
corrdata <- phenodataframe[column_positions]
colnames(corrdata)

#Check data type and coerce to numeric (only for variables that not numeric and continuous as indicated)
#NOTE: do not call as.numeric() directly since as.numeric() gives the internal codes 
str(corrdata)
#######################

corrdata$estudios3c <- as.numeric(as.character(corrdata$estudios3c))
corrdata$Smoke <- as.numeric(as.character(corrdata$Smoke))
corrdata$parity2c <- as.numeric(as.character(corrdata$parity2c))
corrdata$cohort <- as.numeric(as.character(corrdata$cohort))
corrdata$PregnancyComplications <- as.numeric(as.character(corrdata$PregnancyComplications))
corrdata$cesarean <- as.numeric(as.character(corrdata$cesarean))
corrdata$preterm <- as.numeric(as.character(corrdata$preterm))
corrdata$labor_iniciation <- as.numeric(as.character(corrdata$labor_iniciation))
corrdata$sga <- as.numeric(as.character(corrdata$labor_iniciation))

#Calculate correlations  

source("mycor.ci.r") #modified version of cor.ci from psych package, supplied with this example code.

continuous_variables <- c("edadm","sges","edMcCarthy","pdp_verm_m100sd15","pdp_perm_m100sd15",
                          "pdp_genm_m100sd15", "maternal_IQ") 


polychoric_variables <- c("estudios3c","cohort","labor_iniciation", "Smoke") # edit according to your binary or polychoric covariates (for example some cohort might have Smoke as binary)
dichotomous_variables <- c("cesarean", "parity2c","PregnancyComplications", "preterm", "sga"
) # as before, edit accordingly

c_vars_positions <- match(continuous_variables, colnames(corrdata))
p_vars_positions <- match(polychoric_variables, colnames(corrdata))
d_vars_positions <- match(dichotomous_variables, colnames(corrdata))

corrdata<-as.matrix(corrdata)

correlations <- mycor.ci(x=corrdata, keys = NULL, n.iter = 1000, p = 0.05, overlap = FALSE, 
                         poly = FALSE, method = "pearson", plot = FALSE, minlength = 5, 
                         cvars= c_vars_positions, 
                         pvars= p_vars_positions, 
                         dvars= d_vars_positions)
correlations

#Get table of correlation coefficients
corr_table <- correlations["rho"]
corr_table

#Get table of confidence intervals and only keep those with p<.0.05
sig_table <- data.frame(correlations["ci"]) %>% mutate(Correlation=row.names(.)) %>% filter(ci.p<.05)
sig_table

#Write results files and save as .csv file

write.csv(corr_table, file=paste0(cohort,"_IQ_correlations_r_", 
                                  format(Sys.Date(), "%d%m%Y"),".csv"), col.names=FALSE, row.names=TRUE, quote=FALSE)

write.csv(sig_table, file=paste0(cohort,"_IQ_correlations_ci_", 
                                 format(Sys.Date(), "%d%m%Y"),".csv"), col.names=FALSE, row.names=TRUE, quote=FALSE)


#Now you can move on to the next step (Data Analysis)



####################################################
### DATA ANALISYS
####################################################

##Site-Specific Data Analysis

#This stage of the analysis is specific to the chosen exposure/outcome and the specified adjustment variables. Below is the code for all of the analyses to run for the mental health project. Please be sure to update the cohort and date information in the below code for your analysis, as well as the destination path.
#Finally, be sure to update the column names of the exposure/outcome(s) of interest, the adjustment variables, and the table 1 variables. These should correspond to column names in the dataframe specified in the phenofinal argument of the dataAnalysis function.
#Quick check to make sure the function runs in your cohort
#Given the modeling approaches used, the dataAnalysis function requires a good deal of time to run. We recommend first checking whether the function runs on a relatively small subset of sites (i.e. 100 CpG loci). If you encounter any issues, please let us know. If not, proceed to the next step.

#if needed reload Preprocessed data of step 1 and the final phenodataframe from step 2
library(PACEanalysis)

load("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/EPIC_placenta/QC_PACEanalysis_20220508_slide/paceanalysis_017_3/INMA_20220426_Output/INMA_20220426_PreprocessedBetas_nooutliers.RData")



## Running the final models

#Now running the models for all CpG loci. Please note that the R package automatically saves the results in the end.


## Models M1-M3 (Main models)

phenodataframe<-phenodataframe2 # to recover original dimension of the phenodataframe (all cases present even if they have missings in optional variables) 

modelstorun<-data.frame(varofinterest=c("pdp_verm_m100sd15","pdp_perm_m100sd15",
                                        "pdp_genm_m100sd15"))
modelstorun$varofinterest<-as.character(modelstorun$varofinterest)
modelstorun$vartype<-"ExposureCont"
Betasnooutliers<-betafinal.nooutlier ## Change betafinal.nooutlier to the name that appears when you load the methylation data object
Betas<-Betasnooutliers[,colnames(Betasnooutliers)%in%phenodataframe$Basename]
dim(Betas)
Betasnooutliers<-Betas  ## Change betafinal.nooutlier to the name that appears when you load the methylation data object

## check that betas and pheno data are matched: 
table(ifelse(phenodataframe$Basename) ==colnames(Betasnooutliers),"Matched","--NOT MATCHED--")

# if not matched, then match them as follows: 

phenodataframe<-phenodataframe[match(colnames(Betasnooutliers),phenodataframe$Basename),]
all(phenodataframe$Basename == colnames(Betasnooutliers))
# TRUE


for (i in 1:nrow(modelstorun)){
  
  cat("Exposure:",modelstorun$varofinterest[i],"\n")
  
  tempresults<-dataAnalysis(phenofinal=phenodataframe,
                            betafinal=Betasnooutliers, 
                            array="EPIC", ## edit if you have 450K data
                            maxit=100,
                            robust=TRUE,
                            Omega=Celltypes,
                            vartype=modelstorun$vartype[i],
                            varofinterest=modelstorun$varofinterest[i],
                            Table1vars=c("pdp_verm_m100sd15","pdp_perm_m100sd15",
                                         "pdp_genm_m100sd15", "Sex", "edMcCarthy", "estudios3c", "edadm", 
                                         "parity2c", "Smoke", "sges", "Meanlog2oddsContamination", "cohort"),
                            StratifyTable1=FALSE,
                            StratifyTable1var=NULL,             
                            adjustmentvariables=c("Sex", "edMcCarthy", "estudios3c", "edadm", 
                                                  "parity2c", "Smoke", "sges", "Meanlog2oddsContamination", "cohort"),
                            RunUnadjusted=FALSE,
                            RunAdjusted=TRUE,
                            RunCellTypeAdjusted=TRUE,
                            RunSexSpecific=FALSE,
                            RestrictToSubset = FALSE,
                            RestrictionVar = NULL,
                            RestrictToIndicator = NULL,
                            RunCellTypeInteract = FALSE,
                            destinationfolder="/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/",
                            savelog=TRUE,
                            cohort="INMA",
                            analysisdate = "202200601",
                            analysisname = "MAIN")
}

#RUN MAIN WITH THE CELLTYPE INTERACTION:


modelstorun<-data.frame(varofinterest=c("pdp_verm_m100sd15","pdp_perm_m100sd15",
                                        "pdp_genm_m100sd15"))
modelstorun$varofinterest<-as.character(modelstorun$varofinterest)
modelstorun$vartype<-"ExposureCont"
Betasnooutliers<-betafinal.nooutlier ## Change betafinal.nooutlier to the name that appears when you load the methylation data object
Betas<-Betasnooutliers[,colnames(Betasnooutliers)%in%phenodataframe$Basename]
dim(Betas)
Betasnooutliers<-Betas
#811990    255



# if not matched, then match them as follows: 

phenodataframe<-phenodataframe[match(colnames(Betasnooutliers),phenodataframe$Basename),]
all(phenodataframe$Basename == colnames(Betasnooutliers))
# TRUE

dim(phenodataframe)
#[1] 255  22

prova<-celltypes_new2[rownames(celltypes_new2)%in%phenodataframe$Basename,]

phenodataframe<-phenodataframe[match(rownames(prova),phenodataframe$Basename),]
all(phenodataframe$Basename == rownames(prova))
#TRUE


celltypes_newprova<-as.matrix(prova)
dim(celltypes_newprova)
#[1] 255   5


for (i in 1:nrow(modelstorun)){
  
  cat("Exposure:",modelstorun$varofinterest[i],"\n")
  
  tempresults<-dataAnalysis(phenofinal=phenodataframe,
                            betafinal=Betasnooutliers, 
                            array="EPIC", ## edit if you have 450K data
                            maxit=100,
                            robust=TRUE,
                            Omega=celltypes_newprova,
                            vartype=modelstorun$vartype[i],
                            varofinterest=modelstorun$varofinterest[i],
                            Table1vars=c("pdp_verm_m100sd15","pdp_perm_m100sd15",
                                         "pdp_genm_m100sd15", "Sex", "edMcCarthy", "estudios3c", "edadm", 
                                         "parity2c", "Smoke", "sges", "Meanlog2oddsContamination", "cohort"),
                            StratifyTable1=FALSE,
                            StratifyTable1var=NULL,             
                            adjustmentvariables=c("Sex", "edMcCarthy", "estudios3c", "edadm", 
                                                  "parity2c", "Smoke", "sges", "Meanlog2oddsContamination", "cohort"),
                            RunUnadjusted=FALSE,
                            RunAdjusted=TRUE,
                            RunCellTypeAdjusted=TRUE,
                            RunSexSpecific=FALSE,
                            RestrictToSubset = FALSE,
                            RestrictionVar = NULL,
                            RestrictToIndicator = NULL,
                            RunCellTypeInteract = TRUE,
                            destinationfolder="/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results",
                            savelog=TRUE,
                            cohort="INMA",
                            analysisdate = "202200715",
                            analysisname = "MAINcelltypes")
}


## Models M4-M9 (Main models stratified by sex)



modelstorun<-data.frame(varofinterest=c("pdp_verm_m100sd15","pdp_perm_m100sd15",
                                        "pdp_genm_m100sd15"))
modelstorun$varofinterest<-as.character(modelstorun$varofinterest)
modelstorun$vartype<-"ExposureCont"
Betasnooutliers<-Betas ## Change betafinal.nooutlier to the name that appears when you load the methylation data object

## check that betas and pheno data are matched: 
table(ifelse(phenodataframe$Basename) ==colnames(Betasnooutliers),"Matched","--NOT MATCHED--")

# if not matched, then match them as follows: 

phenodataframe<-phenodataframe[match(colnames(Betasnooutliers),phenodataframe$Basename),]
all(phenodataframe$Basename == colnames(Betasnooutliers))
# TRUE

for (i in 1:nrow(modelstorun)){
  
  cat("Exposure:",modelstorun$varofinterest[i],"\n")
  
  tempresults<-dataAnalysis(phenofinal=phenodataframe,
                            betafinal=Betasnooutliers, 
                            array="EPIC", ## edit if you have 450K data
                            maxit=100,
                            robust=TRUE,
                            Omega=Celltypes,
                            vartype=modelstorun$vartype[i],
                            varofinterest=modelstorun$varofinterest[i],
                            Table1vars=c("pdp_verm_m100sd15","pdp_perm_m100sd15",
                                         "pdp_genm_m100sd15", "edMcCarthy", "estudios3c", "edadm", 
                                         "parity2c", "Smoke", "sges", "Meanlog2oddsContamination", "cohort"),
                            StratifyTable1=TRUE,
                            StratifyTable1var=FALSE,             
                            adjustmentvariables=c("edMcCarthy", "estudios3c", "edadm", 
                                                  "parity2c", "Smoke", "sges", "Meanlog2oddsContamination", "cohort"),
                            RunUnadjusted=FALSE,
                            RunAdjusted=TRUE,
                            RunCellTypeAdjusted=TRUE,
                            RunSexSpecific=TRUE,
                            RestrictToSubset = FALSE,
                            RestrictionVar = NULL,
                            RestrictToIndicator = NULL,
                            RunCellTypeInteract = FALSE,
                            destinationfolder="/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/",
                            savelog=TRUE,
                            cohort="INMA",
                            analysisdate = "20220601",
                            analysisname = "MAIN_sex")
}

## Models M10-M12 (models adjusted by maternal IQ)

#please make sure your phenodataframe only includes complete cases in the necessary variables (for these models maternal_IQ becomes necessary)

# to not erase the phenodataframe used in the models with only the main covariates, save it with another name to recover it later in the script
phenodataframe2<- phenodataframe

sum(complete.cases(phenodataframe)==TRUE)   # does the amount look right? Remember that some samples dropped in the methylation QC.

# you can check that the amount does not change in subset of only necessary variables (no need if you know you don't have unnecessary or optional variables such as maternal_IQ with NAs)

phenodataframe<-phenodataframe[!is.na(phenodataframe$maternal_IQ),]
dim(phenodataframe)
#243  22


pheno_subset <- 
  subset(phenodataframe, select=c(pdp_perm_m100sd15, pdp_verm_m100sd15, pdp_genm_m100sd15, Basename, Sex, sges, cohort, edadm, Smoke, parity2c,
                                  estudios3c, edMcCarthy, maternal_IQ, Meanlog2oddsContamination))   # EDIT if you don't e.g have any optional variable or have different amount of PCAs 

sum(complete.cases(pheno_subset)==TRUE)
#243
sum(complete.cases(pheno_subset)==FALSE)  # e.g. this should be 0.
#0

# finally, check that you have done the possible exclusions - check Analysis plan or above

#If you don't have all necessary variables in phenodataframe, you can just load another phenodata (pheno2) and merge them by ID-variable as shown before

modelstorun<-data.frame(varofinterest=c("pdp_verm_m100sd15","pdp_perm_m100sd15",
                                        "pdp_genm_m100sd15"))
modelstorun$varofinterest<-as.character(modelstorun$varofinterest)
modelstorun$vartype<-"ExposureCont"
Betasnooutliers<-Betas

## check that betas and pheno data are matched: 
table(ifelse(phenodataframe$Basename) ==colnames(Betasnooutliers),"Matched","--NOT MATCHED--")

# if not matched, then match them as follows: 

phenodataframe<-phenodataframe[match(colnames(Betasnooutliers),phenodataframe$Basename),]
all(phenodataframe$Basename == colnames(Betasnooutliers))
# TRUE

for (i in 1:nrow(modelstorun)){
  
  cat("Exposure:",modelstorun$varofinterest[i],"\n")
  
  tempresults<-dataAnalysis(phenofinal=phenodataframe,
                            betafinal=Betasnooutliers, 
                            array="EPIC", ## edit if you have 450K data
                            maxit=100,
                            robust=TRUE,
                            Omega=Celltypes,
                            vartype=modelstorun$vartype[i],
                            varofinterest=modelstorun$varofinterest[i],
                            Table1vars=c("pdp_verm_m100sd15","pdp_perm_m100sd15",
                                         "pdp_genm_m100sd15", "Sex", "edMcCarthy", "estudios3c", "edadm", 
                                         "parity2c", "Smoke", "sges", "Meanlog2oddsContamination", "maternal_IQ", "cohort"),
                            StratifyTable1=FALSE,
                            StratifyTable1var=NULL,             
                            adjustmentvariables=c("Sex", "edMcCarthy", "estudios3c", "edadm", 
                                                  "parity2c", "Smoke", "sges", "Meanlog2oddsContamination", "maternal_IQ", "cohort"),
                            RunUnadjusted=FALSE,
                            RunAdjusted=TRUE,
                            RunCellTypeAdjusted=TRUE,
                            RunSexSpecific=FALSE,
                            RestrictToSubset = FALSE,
                            RestrictionVar = NULL,
                            RestrictToIndicator = NULL,
                            RunCellTypeInteract = FALSE,
                            destinationfolder="/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/",
                            savelog=TRUE,
                            cohort="INMA",
                            analysisdate = "20220601NOUIQ",
                            analysisname = "adj_maternal_IQNOU")
  
}


  
  
## Models M13-M15 (Models not adjusted by gestational age (because it could be a mediator and it is interesting to compare adjusted and not adjusted))
  
  # If you have reduced your phenodataframe general n in the previous model adjusting by maternal mental health, then recover the original phenodataframe that was saved: 
modelstorun<-data.frame(varofinterest=c("pdp_verm_m100sd15","pdp_perm_m100sd15",
                                        "pdp_genm_m100sd15"))
modelstorun$varofinterest<-as.character(modelstorun$varofinterest)
modelstorun$vartype<-"ExposureCont"
Betasnooutliers<-Betas  


  for (i in 1:nrow(modelstorun)){
    
    cat("Exposure:",modelstorun$varofinterest[i],"\n")
    
    
    tempresults<-dataAnalysis(phenofinal=phenodataframe,
                              betafinal=Betasnooutliers,
                              array="EPIC", ## edit if you have 450K data
                              maxit=100,
                              robust=TRUE,
                              Omega=Celltypes,
                              vartype=modelstorun$vartype[i],
                              varofinterest=modelstorun$varofinterest[i],
                              Table1vars=c("pdp_verm_m100sd15","pdp_perm_m100sd15",
                                           "pdp_genm_m100sd15", "Sex", "edMcCarthy", "estudios3c", "edadm", 
                                           "parity2c", "Smoke", "Meanlog2oddsContamination", "cohort"),
                              StratifyTable1=FALSE,
                              StratifyTable1var=NULL,             
                              adjustmentvariables=c("Sex", "edMcCarthy", "estudios3c", "edadm", 
                                                    "parity2c", "Smoke", "Meanlog2oddsContamination", "cohort"),
                              RunUnadjusted=FALSE,
                              RunAdjusted=TRUE,
                              RunCellTypeAdjusted=TRUE,
                              RunSexSpecific=FALSE,
                              RestrictToSubset = FALSE,
                              RestrictionVar = NULL,
                              RestrictToIndicator = NULL,
                              RunCellTypeInteract = FALSE,
                              destinationfolder="/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/",
                              savelog=TRUE,
                              cohort="INMA",
                              analysisdate = "20220722NOU",
                              analysisname = "not_adj_sges")
  }
  
  
## Sensitivity analysis ((M16, 17, 18) excluding preterm neonates and participant mothers who experienced pregnancy complications (preeclampsia, diabetes 2, others considered by cohorts that can imply pregnancy complications and an impact on neuro and molecular landscape of placenta))
  
# If you have reduced your phenodataframe general n in the previous model adjusting by maternal mental health, then recover the original phenodataframe that was saved: 
  

phenodataframe3 <- phenodataframe[phenodataframe$PregnancyComplications=="0",]

  dim(phenodataframe3) ## check that the number makes sense
  phenodataframe<-phenodataframe3
  phenodataframe4 <- phenodataframe[phenodataframe$preterm=="0",]    
  dim(phenodataframe4) ## check that the number makes sense
  phenodataframe<-phenodataframe4
  phenodataframe5 <- phenodataframe[phenodataframe$sga=="0",]    
  dim(phenodataframe5) ## check that the number makes sense
  
  phenodataframe<-phenodataframe5
  #211  22
  
  

  modelstorun<-data.frame(varofinterest=c("pdp_verm_m100sd15","pdp_perm_m100sd15",
                                          "pdp_genm_m100sd15"))
  modelstorun$varofinterest<-as.character(modelstorun$varofinterest)
  modelstorun$vartype<-"ExposureCont"
  Betasnooutliers<-Betas  
  
  phenodataframe<-phenodataframe[match(colnames(Betasnooutliers),phenodataframe$Basename),]
  all(phenodataframe$Basename == colnames(Betasnooutliers))
  #TRUE
  
  
  
  for (i in 1:nrow(modelstorun)){
    
    cat("Exposure:",modelstorun$varofinterest[i],"\n")
    
    
    tempresults<-dataAnalysis(phenofinal=phenodataframe,
                              betafinal=Betasnooutliers,
                              array="EPIC", ## edit if you have 450K data
                              maxit=100,
                              robust=TRUE,
                              Omega=Celltypes,
                              vartype=modelstorun$vartype[i],
                              varofinterest=modelstorun$varofinterest[i],
                              Table1vars=c("Sex", "edMcCarthy", "estudios3c", "edadm", "parity2c", 
                                           "sges", "Smoke",  "cohort","Meanlog2oddsContamination"),
                              StratifyTable1=FALSE,
                              StratifyTable1var=NULL,             
                              adjustmentvariables=c("Sex", "edMcCarthy", "estudios3c", "edadm", "parity2c", 
                                                    "sges", "Smoke", "cohort","Meanlog2oddsContamination"),
                              RunUnadjusted=FALSE,
                              RunAdjusted=TRUE,
                              RunCellTypeAdjusted=TRUE,
                              RunSexSpecific=FALSE,
                              RestrictToSubset = FALSE,
                              RestrictionVar = NULL,
                              RestrictToIndicator = NULL,
                              RunCellTypeInteract = FALSE,
                              destinationfolder="/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/",
                              savelog=TRUE,
                              cohort="INMA",
                              analysisdate = "20220707",
                              analysisname = "sensitivity")
  }
  