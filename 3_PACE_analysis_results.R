###################################################
################## MAIN ###########################
###################################################

 
setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/")


# VERBAL

verbalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/verbal/verbal_QCData.txt", header=TRUE)
dim(verbalqc)
#[1] 708105      9     

verbalqcBonf<-verbalqc[verbalqc$padj.bonf=="yes",]
dim(verbalqcBonf)
# 0

verbalqcfdr<-verbalqc[verbalqc$padj.fdr<0.05,]
dim(verbalqcfdr)
# 0

verbalqcsig5<-verbalqc[verbalqc$P_VAL<0.00001,] #  (p_val 10^-5)
dim(verbalqcsig5)
# 29
verbalqcsig4<-verbalqc[verbalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(verbalqcsig4)
# 172

write.table(verbalqcsig5, "verbal_sig_5.txt", col.names=TRUE)
write.table(verbalqcsig4, "verbal_sig_4.txt", col.names=TRUE)


# NON-VERBAL

nonverbalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/nonverbal/nonverbal_QCData.txt", header=TRUE)
nonverbalqcfdr<-nonverbalqc[nonverbalqc$padj.fdr<0.05,]
dim(nonverbalqcfdr)
# 0 
nonverbalqcBonf<-nonverbalqc[nonverbalqc$padj.bonf=="yes",]
dim(nonverbalqcBonf)
# 0
nonverbalqcsig5<-nonverbalqc[nonverbalqc$P_VAL<0.00001,]
dim(nonverbalqcsig5)
# 17
nonverbalqcsig4<-nonverbalqc[nonverbalqc$P_VAL<0.0001,]
dim(nonverbalqcsig4)
# 112

write.table(nonverbalqcsig5, "nonverbal_sig_5.txt", col.names=TRUE)
write.table(nonverbalqcsig4, "nonverbal_sig_4.txt", col.names=TRUE)

 # GENERAL


generalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/general/general_QCData.txt", header=TRUE)
generalqcfdr<-generalqc[generalqc$padj.fdr<0.05,]
dim(generalqcfdr)
# 4 
generalqcBonf<-generalqc[generalqc$padj.bonf=="yes",]
dim(generalqcBonf)
# 1 
generalqcsig5<-generalqc[generalqc$P_VAL<0.00001,]
dim(generalqcsig5)
# 16
generalqcsig4<-generalqc[generalqc$P_VAL<0.0001,]
dim(generalqcsig4)
#103

write.table(generalqcsig5, "general_sig_5.txt", col.names=TRUE)
write.table(generalqcsig4, "general_sig_4.txt", col.names=TRUE)


###############################
### Add ANNOTATION     MAIN ###
###############################


# Install package for annotation

#if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

library(minfi)

data(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

annotation.table = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

annotation.table2<-annotation.table[,c(1, 2, 3, 12, 13, 14, 15, 18, 19, 22, 23, 24)]

# VERBAL
verbalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/verbal/verbal_QCData.txt", header=TRUE)
verbalqcsig4<-verbalqc[verbalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(verbalqcsig4) #172   9
verbalqcsig4ann<-merge(verbalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
dim(verbalqcsig4ann) #[1] 172  21
write.table(verbalqcsig4ann, "verbal_sig_4ann.txt", col.names=TRUE) # save annotated dataframe


# NONVERBAL
nonverbalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/nonverbal/nonverbal_QCData.txt", header=TRUE)
nonverbalqcsig4<-nonverbalqc[nonverbalqc$P_VAL<0.0001,]
dim(nonverbalqcsig4) #112   9
nonverbalqcsig4ann<-merge(nonverbalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
dim(nonverbalqcsig4ann)#[1] 112  21
write.table(nonverbalqcsig4ann, "nonverbal_sig_4ann.txt", col.names=TRUE)


# GENERAL

generalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/general/general_QCData.txt", header=TRUE)
generalqcsig4<-generalqc[generalqc$P_VAL<0.0001,]
dim(generalqcsig4) #103   9
generalqcsig4ann<-merge(generalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
dim(generalqcsig4ann)#103  21
write.table(generalqcsig4ann, "general_sig_4ann.txt", col.names=TRUE)

##########################
####MATERNAL IQ###########
##########################

setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220601NOUIQ_Output/Resultats_matIQ/")

 

# VERBAL

verbalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220601NOUIQ_Output/Resultats_matIQ/QC_ResultsmaternalIQ/verbal/verbal_QCData.txt", header=TRUE)
dim(verbalqc)
#[1] 708105      9     

verbalqcBonf<-verbalqc[verbalqc$padj.bonf=="yes",]
dim(verbalqcBonf)
# 0

verbalqcfdr<-verbalqc[verbalqc$padj.fdr<0.05,]
dim(verbalqcfdr)
# 0
 

verbalqcsig5<-verbalqc[verbalqc$P_VAL<0.00001,] #  (p_val 10^-5)
dim(verbalqcsig5)
#  24  9
verbalqcsig4<-verbalqc[verbalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(verbalqcsig4)
# [1] 164   9


write.table(verbalqcsig5, "verbal_sig_5.txt", col.names=TRUE)
write.table(verbalqcsig4, "verbal_sig_4.txt", col.names=TRUE)

# NON-VERBAL

nonverbalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220601NOUIQ_Output/Resultats_matIQ/QC_ResultsmaternalIQ/nonverbal/nonverbal_QCData.txt", header=TRUE)
nonverbalqcfdr<-nonverbalqc[nonverbalqc$padj.fdr<0.05,]
dim(nonverbalqcfdr)
#2 9
nonverbalqcBonf<-nonverbalqc[nonverbalqc$padj.bonf=="yes",]
dim(nonverbalqcBonf)
# 1 9
 
nonverbalqcsig5<-nonverbalqc[nonverbalqc$P_VAL<0.00001,]
dim(nonverbalqcsig5)
# 17  9
nonverbalqcsig4<-nonverbalqc[nonverbalqc$P_VAL<0.0001,]
dim(nonverbalqcsig4)
# 105   9

write.table(nonverbalqcsig5, "nonverbal_sig_5.txt", col.names=TRUE)
write.table(nonverbalqcsig4, "nonverbal_sig_4.txt", col.names=TRUE)

# GENERAL


generalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220601NOUIQ_Output/Resultats_matIQ/QC_ResultsmaternalIQ/general/general_QCData.txt", header=TRUE)
generalqcfdr<-generalqc[generalqc$padj.fdr<0.05,]
dim(generalqcfdr)
# 0 9
generalqcBonf<-generalqc[generalqc$padj.bonf=="yes",]
dim(generalqcBonf)
# 0 9
# tenim significatives per bonferroni i fdr però per tenir més cpg's per "treballar" agafem pval menys restrictiu i guardem taula. 
generalqcsig5<-generalqc[generalqc$P_VAL<0.00001,]
dim(generalqcsig5)
# 15  9
generalqcsig4<-generalqc[generalqc$P_VAL<0.0001,]
dim(generalqcsig4)
#91  9

write.table(generalqcsig5, "general_sig_5.txt", col.names=TRUE)
write.table(generalqcsig4, "general_sig_4.txt", col.names=TRUE)


###############################
### Add ANNOTATION          ###
###############################

# Install package for annotation

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

library(minfi)

data(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

annotation.table = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

annotation.table2<-annotation.table[,c(1, 2, 3, 12, 13, 14, 15, 18, 19, 22, 23, 24)]

# VERBAL

verbalqcsig4<-verbalqc[verbalqc$P_VAL<0.0001,]
dim(verbalqcsig4) #164   9 if you don't have this file loaded from before, then reload it with read.table, etc...)
verbalqcsig4ann<-merge(verbalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
dim(verbalqcsig4ann) #164  21
write.table(verbalqcsig4ann, "verbal_sig_4ann.txt", col.names=TRUE) # save annotated dataframe

# NONVERBAL
nonverbalqcsig4<-nonverbalqc[nonverbalqc$P_VAL<0.0001,]
dim(nonverbalqcsig4) #105   9
nonverbalqcsig4ann<-merge(nonverbalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
dim(nonverbalqcsig4ann)#105  21
write.table(nonverbalqcsig4ann, "nonverbal_sig_4ann.txt", col.names=TRUE)

# GENERAL

generalqcsig4<-generalqc[generalqc$P_VAL<0.0001,]
dim(generalqcsig4) #91  9
generalqcsig4ann<-merge(generalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
dim(generalqcsig4ann)#91 21
write.table(generalqcsig4ann, "general_sig_4ann.txt", col.names=TRUE)


###############################
####### SENSITIVITY ###########
###############################

setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/")

# VERBAL


verbalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/verbal/verbal_QCData.txt", header=TRUE)
dim(verbalqc)
#[1] 708105      9  

verbalqcBonf<-verbalqc[verbalqc$padj.bonf=="yes",]
dim(verbalqcBonf)
# 0 9
verbalqcfdr<-verbalqc[verbalqc$padj.fdr<0.05,]
dim(verbalqcfdr)
# 0 9

verbalqcsig5<-verbalqc[verbalqc$P_VAL<0.00001,] #  (p_val 10^-5)
dim(verbalqcsig5)
# 29  9
verbalqcsig4<-verbalqc[verbalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(verbalqcsig4)
#  172   9

write.table(verbalqcsig5, "verbal_sig_5.txt", col.names=TRUE)
write.table(verbalqcsig4, "verbal_sig_4.txt", col.names=TRUE)

# NON-VERBAL


nonverbalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/nonverbal/nonverbal_QCData.txt", header=TRUE)
nonverbalqcfdr<-nonverbalqc[nonverbalqc$padj.fdr<0.05,]
dim(nonverbalqcfdr)
#0 9
nonverbalqcBonf<-nonverbalqc[nonverbalqc$padj.bonf=="yes",]
dim(nonverbalqcBonf)
# 0 9

nonverbalqcsig5<-nonverbalqc[nonverbalqc$P_VAL<0.00001,]
dim(nonverbalqcsig5)
# 17  9
nonverbalqcsig4<-nonverbalqc[nonverbalqc$P_VAL<0.0001,]
dim(nonverbalqcsig4)
# 112   9
write.table(nonverbalqcsig5, "nonverbal_sig_5.txt", col.names=TRUE)
write.table(nonverbalqcsig4, "nonverbal_sig_4.txt", col.names=TRUE)

# GENERAL


generalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/general/general_QCData.txt", header=TRUE)
generalqcfdr<-generalqc[generalqc$padj.fdr<0.05,]
dim(generalqcfdr)
# 4 9


generalqcBonf<-generalqc[generalqc$padj.bonf=="yes",]
dim(generalqcBonf)
# 1] 1 9

generalqcsig5<-generalqc[generalqc$P_VAL<0.00001,]
dim(generalqcsig5)
# [1] 16  9
generalqcsig4<-generalqc[generalqc$P_VAL<0.0001,]
dim(generalqcsig4)
#[1] 103   9

write.table(generalqcsig5, "general_sig_5.txt", col.names=TRUE)
write.table(generalqcsig4, "general_sig_4.txt", col.names=TRUE)



###############################
### Add ANNOTATION          ###
###############################


# Install package for annotation

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

library(minfi)

data(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

annotation.table = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

annotation.table2<-annotation.table[,c(1, 2, 3, 12, 13, 14, 15, 18, 19, 22, 23, 24)]

# VERBAL

verbalqcsig4<-verbalqc[verbalqc$P_VAL<0.0001,]
dim(verbalqcsig4) # 172   9 
verbalqcsig4ann<-merge(verbalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
dim(verbalqcsig4ann) #172  21

write.table(verbalqcsig4ann, "verbal_sig_4ann.txt", col.names=TRUE) # save annotated dataframe


# NONVERBAL
nonverbalqcsig4<-nonverbalqc[nonverbalqc$P_VAL<0.0001,]
dim(nonverbalqcsig4) #112   9

nonverbalqcsig4ann<-merge(nonverbalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
dim(nonverbalqcsig4ann)#112  21

write.table(nonverbalqcsig4ann, "nonverbal_sig_4ann.txt", col.names=TRUE)



# GENERAL

generalqcsig4<-generalqc[generalqc$P_VAL<0.0001,]
dim(generalqcsig4) #[1] 103   9

generalqcsig4ann<-merge(generalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
dim(generalqcsig4ann)#[1] 103  21

write.table(generalqcsig4ann, "general_sig_4ann.txt", col.names=TRUE)


###############################
########### SGES ##############
###############################

setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220722NOU_Output/resultats_nous_sges/")

# VERBAL


verbalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220722NOU_Output/resultats_nous_sges/QC_ResultsSGES/verbal/verbal_QCData.txt", header=TRUE)
dim(verbalqc)
#[1] 708105      9

verbalqcBonf<-verbalqc[verbalqc$padj.bonf=="yes",]
dim(verbalqcBonf)
#  0 9


verbalqcfdr<-verbalqc[verbalqc$padj.fdr<0.05,]
dim(verbalqcfdr)
# 0 9

verbalqcsig5<-verbalqc[verbalqc$P_VAL<0.00001,] #  (p_val 10^-5)
dim(verbalqcsig5)
# 27 9

verbalqcsig4<-verbalqc[verbalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(verbalqcsig4)
#189


write.table(verbalqcsig5, "verbal_sig_5.txt", col.names=TRUE)
write.table(verbalqcsig4, "verbal_sig_4.txt", col.names=TRUE)

# NON-VERBAL


nonverbalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220722NOU_Output/resultats_nous_sges/QC_ResultsSGES/nonverbal/nonverbal_QCData.txt", header=TRUE)
nonverbalqcfdr<-nonverbalqc[nonverbalqc$padj.fdr<0.05,]
dim(nonverbalqcfdr)
#0 9
nonverbalqcBonf<-nonverbalqc[nonverbalqc$padj.bonf=="yes",]
dim(nonverbalqcBonf)
# 0 9
nonverbalqcsig5<-nonverbalqc[nonverbalqc$P_VAL<0.00001,]
dim(nonverbalqcsig5)
# 19  9
nonverbalqcsig4<-nonverbalqc[nonverbalqc$P_VAL<0.0001,]
dim(nonverbalqcsig4)
# 115   9


write.table(nonverbalqcsig5, "nonverbal_sig_5.txt", col.names=TRUE)
write.table(nonverbalqcsig4, "nonverbal_sig_4.txt", col.names=TRUE)

# GENERAL


generalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220722NOU_Output/resultats_nous_sges/QC_ResultsSGES/general/general_QCData.txt", header=TRUE)
generalqcfdr<-generalqc[generalqc$padj.fdr<0.05,]
dim(generalqcfdr)
# 1 9
generalqcBonf<-generalqc[generalqc$padj.bonf=="yes",]
dim(generalqcBonf)
# 1 9
generalqcsig5<-generalqc[generalqc$P_VAL<0.00001,]
dim(generalqcsig5)
# 17  9
generalqcsig4<-generalqc[generalqc$P_VAL<0.0001,]
dim(generalqcsig4)
#102 9

write.table(generalqcsig5, "general_sig_5.txt", col.names=TRUE)
write.table(generalqcsig4, "general_sig_4.txt", col.names=TRUE)


###############################
### Add ANNOTATION     sges ###
###############################


# Install package for annotation

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

library(minfi)

data(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

annotation.table = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

annotation.table2<-annotation.table[,c(1, 2, 3, 12, 13, 14, 15, 18, 19, 22, 23, 24)]

# VERBAL

verbalqcsig4<-verbalqc[verbalqc$P_VAL<0.0001,]
dim(verbalqcsig4) # 172   9
verbalqcsig4ann<-merge(verbalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
dim(verbalqcsig4ann) #172  21
write.table(verbalqcsig4ann, "verbal_sig_4ann.txt", col.names=TRUE) # save annotated dataframe


# NONVERBAL
nonverbalqcsig4<-nonverbalqc[nonverbalqc$P_VAL<0.0001,]
dim(nonverbalqcsig4) #112   9
nonverbalqcsig4ann<-merge(nonverbalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
dim(nonverbalqcsig4ann)#112  21
write.table(nonverbalqcsig4ann, "nonverbal_sig_4ann.txt", col.names=TRUE)



# GENERAL

generalqcsig4<-generalqc[generalqc$P_VAL<0.0001,]
dim(generalqcsig4) #103   9
generalqcsig4ann<-merge(generalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
dim(generalqcsig4ann)#103  21
write.table(generalqcsig4ann, "general_sig_4ann.txt", col.names=TRUE)


###########################################
####CELLTYPES results PACE. and anotation
##############################################

# Install package for annotation

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

library(minfi)

data(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

annotation.table = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

annotation.table2<-annotation.table[,c(1, 2, 3, 12, 13, 14, 15, 18, 19, 22, 23, 24)]


#######################################ENDOTHELIAL 

setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/")


# VERBAL

verbalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/QC_ResultsENDOTHELIAL/verbal/verbal_QCData.txt", header=TRUE)
dim(verbalqc)
#708105      9    

verbalqcBonf<-verbalqc[verbalqc$padj.bonf=="yes",]
dim(verbalqcBonf)
#  0 9

verbalqcfdr<-verbalqc[verbalqc$padj.fdr<0.05,]
dim(verbalqcfdr)
#  0 9

verbalqcsig6<-verbalqc[verbalqc$P_VAL<0.000001,] #  (p_val 10^-6)
dim(verbalqcsig6) # 0 9

verbalqcsig5<-verbalqc[verbalqc$P_VAL<0.00001,] #  (p_val 10^-5)
dim(verbalqcsig5) #13  9

verbalqcsig4<-verbalqc[verbalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(verbalqcsig4) #117   9

write.table(verbalqcsig5, "verbal_sig_5ENDOTHELIAL.txt", col.names=TRUE)
write.table(verbalqcsig4, "verbal_sig_4ENDOTHELIAL.txt", col.names=TRUE)


verbalqcsig4ann<-merge(verbalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
write.table(verbalqcsig4ann, "verbal_sig_4ENDOTHELIALann.txt", col.names=TRUE) # save annotated dataframe

verbalqcsig5ann<-merge(verbalqcsig5, annotation.table2, by.x="probeID", by.y="row.names")
write.table(verbalqcsig5ann, "verbal_sig_5ENDOTHELIALann.txt", col.names=TRUE) # save annotated dataframe


# NON-VERBAL

nonverbalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/QC_ResultsENDOTHELIAL/nonverbal/nonverbal_QCData.txt", header=TRUE)
dim(nonverbalqc)
#[1] 708105      9


nonverbalqcfdr<-nonverbalqc[nonverbalqc$padj.fdr<0.05,]
dim(nonverbalqcfdr)
# 0 
nonverbalqcBonf<-nonverbalqc[nonverbalqc$padj.bonf=="yes",]
dim(nonverbalqcBonf)
# 0
nonverbalqcsig6<-nonverbalqc[nonverbalqc$P_VAL<0.000001,] #  (p_val 10^-6)
dim(nonverbalqcsig6) #1

nonverbalqcsig5<-nonverbalqc[nonverbalqc$P_VAL<0.00001,] #  (p_val 10^-5)
dim(nonverbalqcsig5) # 8

nonverbalqcsig4<-nonverbalqc[nonverbalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(nonverbalqcsig4) #107

write.table(nonverbalqcsig5, "nonverbal_sig_5ENDOTHELIAL.txt", col.names=TRUE)
write.table(nonverbalqcsig4, "nonverbal_sig_4ENDOTHELIAL.txt", col.names=TRUE)

nonverbalqcsig4ann<-merge(nonverbalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
write.table(nonverbalqcsig4ann, "nonverbal_sig_4ENDOTHELIALann.txt", col.names=TRUE) # save annotated dataframe

nonverbalqcsig5ann<-merge(nonverbalqcsig5, annotation.table2, by.x="probeID", by.y="row.names")
write.table(nonverbalqcsig5ann, "nonverbal_sig_5ENDOTHELIALann.txt", col.names=TRUE) # save annotated dataframe


# GENERAL


generalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/QC_ResultsENDOTHELIAL/general/general_QCData.txt", header=TRUE)
dim(nonverbalqc)
#[1] 708105      9
generalqcfdr<-generalqc[generalqc$padj.fdr<0.05,]
dim(generalqcfdr)
# 0 
generalqcBonf<-generalqc[generalqc$padj.bonf=="yes",]
dim(generalqcBonf)
# 0 

generalqcsig6<-generalqc[generalqc$P_VAL<0.000001,] #  (p_val 10^-5)
dim(generalqcsig6) #1

generalqcsig5<-generalqc[generalqc$P_VAL<0.00001,] #  (p_val 10^-5)
dim(generalqcsig5) #16

generalqcsig4<-generalqc[generalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(generalqcsig4) #172

write.table(generalqcsig5, "generalqc_sig_5ENDOTHELIAL.txt", col.names=TRUE)
write.table(generalqcsig4, "generalqc_sig_4ENDOTHELIAL.txt", col.names=TRUE)


generalqcsig4ann<-merge(generalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
write.table(generalqcsig4ann, "generalqc_sig_4ENDOTHELIALann.txt", col.names=TRUE) # save annotated dataframe

generalqcsig5ann<-merge(generalqcsig5, annotation.table2, by.x="probeID", by.y="row.names")
write.table(generalqcsig5ann, "generalqc_sig_5ENDOTHELIALann.txt", col.names=TRUE) # save annotated dataframe



#######################################syncytiotrophoblast

 
setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/")


# VERBAL

verbalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/QC_Resultssyncytiotrophoblast/verbal/verbal_QCData.txt", header=TRUE)
dim(verbalqc)
#708105      9    

verbalqcBonf<-verbalqc[verbalqc$padj.bonf=="yes",]
dim(verbalqcBonf)
#  0 9

verbalqcfdr<-verbalqc[verbalqc$padj.fdr<0.05,]
dim(verbalqcfdr)
#  0 9
 

verbalqcsig6<-verbalqc[verbalqc$P_VAL<0.000001,] #  (p_val 10^-6)
dim(verbalqcsig6) # 3 9

verbalqcsig5<-verbalqc[verbalqc$P_VAL<0.00001,] #  (p_val 10^-5)
dim(verbalqcsig5) #12

verbalqcsig4<-verbalqc[verbalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(verbalqcsig4) #79

write.table(verbalqcsig5, "verbal_sig_5txtsyncytiotrophoblast.txt", col.names=TRUE)
write.table(verbalqcsig4, "verbal_sig_4txtsyncytiotrophoblast.txt", col.names=TRUE)


verbalqcsig4ann<-merge(verbalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
write.table(verbalqcsig4ann, "verbal_sig_4syncytiotrophoblastann.txt", col.names=TRUE) # save annotated dataframe

verbalqcsig5ann<-merge(verbalqcsig5, annotation.table2, by.x="probeID", by.y="row.names")
write.table(verbalqcsig5ann, "verbal_sig_5syncytiotrophoblastann.txt", col.names=TRUE) # save annotated dataframe


# NON-VERBAL

nonverbalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/QC_Resultssyncytiotrophoblast/nonverbal/nonverbal_QCData.txt", header=TRUE)
dim(nonverbalqc)
#[1] 708105      9


nonverbalqcfdr<-nonverbalqc[nonverbalqc$padj.fdr<0.05,]
dim(nonverbalqcfdr)
# 0 
nonverbalqcBonf<-nonverbalqc[nonverbalqc$padj.bonf=="yes",]
dim(nonverbalqcBonf)
# 0
nonverbalqcsig6<-nonverbalqc[nonverbalqc$P_VAL<0.000001,] #  (p_val 10^-6)
dim(nonverbalqcsig6) #0

nonverbalqcsig5<-nonverbalqc[nonverbalqc$P_VAL<0.00001,] #  (p_val 10^-5)
dim(nonverbalqcsig5) # 6

nonverbalqcsig4<-nonverbalqc[nonverbalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(nonverbalqcsig4) #66

write.table(nonverbalqcsig5, "nonverbal_sig_5syncytiotrophoblast.txt", col.names=TRUE)
write.table(nonverbalqcsig4, "nonverbal_sig_4syncytiotrophoblast.txt", col.names=TRUE)

nonverbalqcsig4ann<-merge(nonverbalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
write.table(nonverbalqcsig4ann, "nonverbal_sig_4syncytiotrophoblastann.txt", col.names=TRUE) # save annotated dataframe

nonverbalqcsig5ann<-merge(nonverbalqcsig5, annotation.table2, by.x="probeID", by.y="row.names")
write.table(nonverbalqcsig5ann, "nonverbal_sig_5syncytiotrophoblastann.txt", col.names=TRUE) # save annotated dataframe

# GENERAL


generalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/QC_Resultssyncytiotrophoblast/general/general_QCData.txt", header=TRUE)
dim(nonverbalqc)
#[1] 708105      9
generalqcfdr<-generalqc[generalqc$padj.fdr<0.05,]
dim(generalqcfdr)
# 0 
generalqcBonf<-generalqc[generalqc$padj.bonf=="yes",]
dim(generalqcBonf)
# 0 

generalqcsig6<-generalqc[generalqc$P_VAL<0.000001,] #  (p_val 10^-5)
dim(generalqcsig6) #2

generalqcsig5<-generalqc[generalqc$P_VAL<0.00001,] #  (p_val 10^-5)
dim(generalqcsig5) #17

generalqcsig4<-generalqc[generalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(generalqcsig4) #214

write.table(generalqcsig5, "generalqc_sig_5syncytiotrophoblast.txt", col.names=TRUE)
write.table(generalqcsig4, "generalqc_sig_4syncytiotrophoblast.txt", col.names=TRUE)

generalqcsig4ann<-merge(generalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
write.table(generalqcsig4ann, "generalqc_sig_4syncytiotrophoblastann.txt", col.names=TRUE) # save annotated dataframe

generalqcsig5ann<-merge(generalqcsig5, annotation.table2, by.x="probeID", by.y="row.names")
write.table(generalqcsig5ann, "generalqc_sig_5syncytiotrophoblastann.txt", col.names=TRUE) # save annotated dataframe


#######################################nRBC

 
setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/")


# VERBAL

verbalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/QC_ResultsnRBC/verbal/verbal_QCData.txt", header=TRUE)
dim(verbalqc)
#708105      9    

verbalqcBonf<-verbalqc[verbalqc$padj.bonf=="yes",]
dim(verbalqcBonf)
#  0 9

verbalqcfdr<-verbalqc[verbalqc$padj.fdr<0.05,]
dim(verbalqcfdr)
#  0 9
 

verbalqcsig6<-verbalqc[verbalqc$P_VAL<0.000001,] #  (p_val 10^-6)
dim(verbalqcsig6) # 2

verbalqcsig5<-verbalqc[verbalqc$P_VAL<0.00001,] #  (p_val 10^-5)
dim(verbalqcsig5) #8

verbalqcsig4<-verbalqc[verbalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(verbalqcsig4) #35

write.table(verbalqcsig5, "verbal_sig_5nRBC.txt", col.names=TRUE)
write.table(verbalqcsig4, "verbal_sig_4nRBC.txt", col.names=TRUE)


verbalqcsig4ann<-merge(verbalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
write.table(verbalqcsig4ann, "verbal_sig_4nRBCann.txt", col.names=TRUE) # save annotated dataframe

verbalqcsig5ann<-merge(verbalqcsig5, annotation.table2, by.x="probeID", by.y="row.names")
write.table(verbalqcsig5ann, "verbal_sig_5nRBCann.txt", col.names=TRUE) # save annotated dataframe


# NON-VERBAL

nonverbalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/QC_ResultsnRBC/nonverbal/nonverbal_QCData.txt", header=TRUE)
dim(nonverbalqc)
#[1] 708105      9


nonverbalqcfdr<-nonverbalqc[nonverbalqc$padj.fdr<0.05,]
dim(nonverbalqcfdr)
# 0 
nonverbalqcBonf<-nonverbalqc[nonverbalqc$padj.bonf=="yes",]
dim(nonverbalqcBonf)
# 0
nonverbalqcsig6<-nonverbalqc[nonverbalqc$P_VAL<0.000001,] #  (p_val 10^-6)
dim(nonverbalqcsig6) #0

nonverbalqcsig5<-nonverbalqc[nonverbalqc$P_VAL<0.00001,] #  (p_val 10^-5)
dim(nonverbalqcsig5) # 2

nonverbalqcsig4<-nonverbalqc[nonverbalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(nonverbalqcsig4) #25

write.table(nonverbalqcsig5, "nonverbal_sig_5nRBC.txt", col.names=TRUE)
write.table(nonverbalqcsig4, "nonverbal_sig_4nRBC.txt", col.names=TRUE)

nonverbalqcsig4ann<-merge(nonverbalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
write.table(nonverbalqcsig4ann, "nonverbal_sig_4nRBCann.txt", col.names=TRUE) # save annotated dataframe

nonverbalqcsig5ann<-merge(nonverbalqcsig5, annotation.table2, by.x="probeID", by.y="row.names")
write.table(nonverbalqcsig5ann, "nonverbal_sig_5nRBCann.txt", col.names=TRUE) # save annotated dataframe

# GENERAL


generalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/QC_ResultsnRBC/general/general_QCData.txt", header=TRUE)
dim(generalqc)
#[1] 708105      9
generalqcfdr<-generalqc[generalqc$padj.fdr<0.05,]
dim(generalqcfdr)
# 0 
generalqcBonf<-generalqc[generalqc$padj.bonf=="yes",]
dim(generalqcBonf)
# 0 

generalqcsig6<-generalqc[generalqc$P_VAL<0.000001,] #  (p_val 10^-5)
dim(generalqcsig6) #1

generalqcsig5<-generalqc[generalqc$P_VAL<0.00001,] #  (p_val 10^-5)
dim(generalqcsig5) #2

generalqcsig4<-generalqc[generalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(generalqcsig4) #23

write.table(generalqcsig5, "generalqc_sig_5nRBC.txt", col.names=TRUE)
write.table(generalqcsig4, "generalqc_sig_4nRBC.txt", col.names=TRUE)

generalqcsig4ann<-merge(generalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
write.table(generalqcsig4ann, "generalqc_sig_4nRBCann.txt", col.names=TRUE) # save annotated dataframe

generalqcsig5ann<-merge(generalqcsig5, annotation.table2, by.x="probeID", by.y="row.names")
write.table(generalqcsig5ann, "generalqc_sig_5nRBCann.txt", col.names=TRUE) # save annotated dataframe


#######################################Hofbauer

 
setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/")


# VERBAL

verbalqc<- read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/QC_Resultshofbauer/verbal/verbal_QCData.txt", header=TRUE)
dim(verbalqc)
#708105      9    

verbalqcBonf<-verbalqc[verbalqc$padj.bonf=="yes",]
dim(verbalqcBonf)
#  0 9

verbalqcfdr<-verbalqc[verbalqc$padj.fdr<0.05,]
dim(verbalqcfdr)
#  0 9
 

verbalqcsig6<-verbalqc[verbalqc$P_VAL<0.000001,] #  (p_val 10^-6)
dim(verbalqcsig6) # 1 9

verbalqcsig5<-verbalqc[verbalqc$P_VAL<0.00001,] #  (p_val 10^-5)
dim(verbalqcsig5) #3

verbalqcsig4<-verbalqc[verbalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(verbalqcsig4) #42

write.table(verbalqcsig5, "verbal_sig_5Hofbauer.txt", col.names=TRUE)
write.table(verbalqcsig4, "verbal_sig_4Hofbauer.txt", col.names=TRUE)


verbalqcsig4ann<-merge(verbalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
write.table(verbalqcsig4ann, "verbal_sig_4Hofbauerann.txt", col.names=TRUE) # save annotated dataframe

verbalqcsig5ann<-merge(verbalqcsig5, annotation.table2, by.x="probeID", by.y="row.names")
write.table(verbalqcsig5ann, "verbal_sig_5Hofbauerann.txt", col.names=TRUE) # save annotated dataframe


# NON-VERBAL

nonverbalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/QC_Resultshofbauer/nonverbal/nonverbal_QCData.txt", header=TRUE)
dim(nonverbalqc)
#[1] 708105      9


nonverbalqcfdr<-nonverbalqc[nonverbalqc$padj.fdr<0.05,]
dim(nonverbalqcfdr)
# 0 
nonverbalqcBonf<-nonverbalqc[nonverbalqc$padj.bonf=="yes",]
dim(nonverbalqcBonf)
# 0
nonverbalqcsig6<-nonverbalqc[nonverbalqc$P_VAL<0.000001,] #  (p_val 10^-6)
dim(nonverbalqcsig6) #0

nonverbalqcsig5<-nonverbalqc[nonverbalqc$P_VAL<0.00001,] #  (p_val 10^-5)
dim(nonverbalqcsig5) # 5

nonverbalqcsig4<-nonverbalqc[nonverbalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(nonverbalqcsig4) #41

write.table(nonverbalqcsig5, "nonverbal_sig_5Hofbauer.txt", col.names=TRUE)
write.table(nonverbalqcsig4, "nonverbal_sig_4Hofbauer.txt", col.names=TRUE)

nonverbalqcsig4ann<-merge(nonverbalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
write.table(nonverbalqcsig4ann, "nonverbal_sig_4Hofbauerann.txt", col.names=TRUE) # save annotated dataframe

nonverbalqcsig5ann<-merge(nonverbalqcsig5, annotation.table2, by.x="probeID", by.y="row.names")
write.table(nonverbalqcsig5ann, "nonverbal_sig_5Hofbauerann.txt", col.names=TRUE) # save annotated dataframe

# GENERAL


generalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/QC_Resultshofbauer/general/general_QCData.txt", header=TRUE)
dim(generalqc)
#[1] 708105      9
generalqcfdr<-generalqc[generalqc$padj.fdr<0.05,]
dim(generalqcfdr)
# 0 
generalqcBonf<-generalqc[generalqc$padj.bonf=="yes",]
dim(generalqcBonf)
# 0 

generalqcsig6<-generalqc[generalqc$P_VAL<0.000001,] #  (p_val 10^-5)
dim(generalqcsig6) #1

generalqcsig5<-generalqc[generalqc$P_VAL<0.00001,] #  (p_val 10^-5)
dim(generalqcsig5) #9

generalqcsig4<-generalqc[generalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(generalqcsig4) #57

write.table(generalqcsig5, "generalqc_sig_5Hofbauer.txt", col.names=TRUE)
write.table(generalqcsig4, "generalqc_sig_4Hofbauer.txt", col.names=TRUE)

generalqcsig4ann<-merge(generalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
write.table(generalqcsig4ann, "generalqc_sig_4Hofbauerann.txt", col.names=TRUE) # save annotated dataframe

generalqcsig5ann<-merge(generalqcsig5, annotation.table2, by.x="probeID", by.y="row.names")
write.table(generalqcsig5ann, "generalqc_sig_5Hofbauerann.txt", col.names=TRUE) # save annotated dataframe



#######################################Trophoblast.

 
setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/")


# VERBAL

verbalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/QC_Resultstrophoblast/verbal/verbal_QCData.txt", header=TRUE)
dim(verbalqc)
#708105      9    

verbalqcBonf<-verbalqc[verbalqc$padj.bonf=="yes",]
dim(verbalqcBonf)
#  0 9

verbalqcfdr<-verbalqc[verbalqc$padj.fdr<0.05,]
dim(verbalqcfdr)
#  0 9
 

verbalqcsig6<-verbalqc[verbalqc$P_VAL<0.000001,] #  (p_val 10^-6)
dim(verbalqcsig6) # 0 9

verbalqcsig5<-verbalqc[verbalqc$P_VAL<0.00001,] #  (p_val 10^-5)
dim(verbalqcsig5) #2

verbalqcsig4<-verbalqc[verbalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(verbalqcsig4) #62

write.table(verbalqcsig5, "verbal_sig_5Trophoblast.txt", col.names=TRUE)
write.table(verbalqcsig4, "verbal_sig_4Trophoblast.txt", col.names=TRUE)


verbalqcsig4ann<-merge(verbalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
write.table(verbalqcsig4ann, "verbal_sig_4Trophoblastann.txt", col.names=TRUE) # save annotated dataframe

verbalqcsig5ann<-merge(verbalqcsig5, annotation.table2, by.x="probeID", by.y="row.names")
write.table(verbalqcsig5ann, "verbal_sig_5Trophoblastann.txt", col.names=TRUE) # save annotated dataframe

# NON-VERBAL

nonverbalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/QC_Resultstrophoblast/nonverbal/nonverbal_QCData.txt", header=TRUE)
dim(nonverbalqc)
#[1] 708105      9


nonverbalqcfdr<-nonverbalqc[nonverbalqc$padj.fdr<0.05,]
dim(nonverbalqcfdr)
# 0 
nonverbalqcBonf<-nonverbalqc[nonverbalqc$padj.bonf=="yes",]
dim(nonverbalqcBonf)
# 0
nonverbalqcsig6<-nonverbalqc[nonverbalqc$P_VAL<0.000001,] #  (p_val 10^-6)
dim(nonverbalqcsig6) #1

nonverbalqcsig5<-nonverbalqc[nonverbalqc$P_VAL<0.00001,] #  (p_val 10^-5)
dim(nonverbalqcsig5) # 14

nonverbalqcsig4<-nonverbalqc[nonverbalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(nonverbalqcsig4) #144

write.table(nonverbalqcsig5, "nonverbal_sig_5Trophoblast.txt", col.names=TRUE)
write.table(nonverbalqcsig4, "nonverbal_sig_4Trophoblast.txt", col.names=TRUE)

nonverbalqcsig4ann<-merge(nonverbalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
write.table(nonverbalqcsig4ann, "nonverbal_sig_4Trophoblastann.txt", col.names=TRUE) # save annotated dataframe

nonverbalqcsig5ann<-merge(nonverbalqcsig5, annotation.table2, by.x="probeID", by.y="row.names")
write.table(nonverbalqcsig5ann, "nonverbal_sig_5Trophoblastann.txt", col.names=TRUE) # save annotated dataframe



# GENERAL


generalqc<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200715_Output/RESULTATScelltypes/QC_Resultstrophoblast/general/general_QCData.txt", header=TRUE)
dim(generalqc)
#[1] 708105      9
generalqcfdr<-generalqc[generalqc$padj.fdr<0.05,]
dim(generalqcfdr)
# 0 
generalqcBonf<-generalqc[generalqc$padj.bonf=="yes",]
dim(generalqcBonf)
# 0 

generalqcsig6<-generalqc[generalqc$P_VAL<0.000001,] #  (p_val 10^-5)
dim(generalqcsig6) #1

generalqcsig5<-generalqc[generalqc$P_VAL<0.00001,] #  (p_val 10^-5)
dim(generalqcsig5) #7

generalqcsig4<-generalqc[generalqc$P_VAL<0.0001,] #  (p_val 10^-5)
dim(generalqcsig4) #63

write.table(generalqcsig5, "generalqc_sig_5Trophoblast.txt", col.names=TRUE)
write.table(generalqcsig4, "generalqc_sig_4Trophoblast.txt", col.names=TRUE)

generalqcsig4ann<-merge(generalqcsig4, annotation.table2, by.x="probeID", by.y="row.names")
write.table(generalqcsig4ann, "generalqc_sig_4Trophoblastann.txt", col.names=TRUE) # save annotated dataframe

generalqcsig5ann<-merge(generalqcsig5, annotation.table2, by.x="probeID", by.y="row.names")
write.table(generalqcsig5ann, "generalqc_sig_5Trophoblastann.txt", col.names=TRUE) # save annotated dataframe
