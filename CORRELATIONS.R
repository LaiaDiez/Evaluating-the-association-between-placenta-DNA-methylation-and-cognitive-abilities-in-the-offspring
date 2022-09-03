
############
###GENERAL MAIN vs SENSITIVITY###
### ALL

setwd("/home/isglobal.lan/mcosin/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/COMPARE_MODELS/marta_28072022/")

gen_sens<-read.table("/home/isglobal.lan/mcosin/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220707_Output/SENSITIVITY_NOU/QC_ResultsSENSITIVITY/general/general_QCData.txt", header=TRUE)
head(gen_sens)
#probeID          BETA           SE        P_VAL    padj.fdr padj.bonf
#1 cg13073105  4.908501e-04 8.195346e-05 2.106474e-09 0.001491605       yes
#2 cg05409601 -6.419241e-05 1.236778e-05 2.099631e-07 0.074337972        no
#3 cg25185985  2.386892e-04 4.775170e-05 5.776363e-07 0.096862474        no
#4 cg16470259 -1.495487e-03 3.003753e-04 6.400278e-07 0.096862474        no
#5 cg00866476  5.486221e-04 1.104787e-04 6.839556e-07 0.096862474        no
#6 cg01282822 -7.876320e-04 1.618632e-04 1.138598e-06 0.110607358        no

#CpG_chrm   CpG_beg   CpG_end
#1     chr1 201355349 201355351
#2     chr4  77156758  77156760
#3    chr14  50911590  50911592
#4    chr16  22824946  22824948
#5     chr2  95074936  95074938
#6    chr18  25164493  25164495

gen_all<-read.table("/home/isglobal.lan/mcosin/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/general/general_QCData.txt", header=TRUE)

head(gen_all)
#probeID          BETA           SE        P_VAL   padj.fdr padj.bonf
#1 cg15480200 -1.792379e-03 3.322577e-04 6.869765e-08 0.03927355       yes
#2 cg02986379 -6.112358e-05 1.151569e-05 1.109258e-07 0.03927355        no
#3 cg00866476  5.480824e-04 1.058830e-04 2.263267e-07 0.04876029        no
#4 cg14113931 -4.040245e-04 7.861133e-05 2.754410e-07 0.04876029        no
#5 cg14613533 -9.486538e-04 1.886799e-04 4.960130e-07 0.07024586        no
#6 cg22965956 -1.107075e-03 2.220194e-04 6.151800e-07 0.07260200        no
#CpG_chrm   CpG_beg   CpG_end
#1    chr10 132817110 132817112
#2     chr5  39425327  39425329
#3     chr2  95074936  95074938
#4    chr18  12701670  12701672
#5    chr17   5192674   5192676
#6    chr18  74212614  74212616

probes<-gen_all$probeID #probes from the general model that are sig by pval*10^-4

sensord<-gen_sens[order(gen_sens$probeID),]
genallord<-gen_all[order(gen_all$probeID),]

sensordb<- sensord[,c("probeID","BETA")]
genallordb<- genallord[,c("probeID","BETA")]

colnames(sensordb)<-c("probeID","BETA_gen_sensitivity")
colnames(genallordb)<-c("probeID","BETA_main_gen_all")

corr_matrix<-merge(sensordb, genallordb,by="probeID") 

cor.test(corr_matrix$BETA_gen_sensitivity, corr_matrix$BETA_main_gen_all)

#Pearson's product-moment correlation

#data:  corr_matrix$BETA_gen_sensitivity and corr_matrix$BETA_main_gen_all
#t = 1949.3, df = 708103, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9177361 0.9184679
#sample estimates:
#      cor
#0.9181028

library(GGally)

corr_matrix2<-corr_matrix[,2:3]

plot.cor.BETA<-ggpairs(corr_matrix2, title="Correlation between main general all vs general sensitivity all ")
file.create(paste0("QCCorrelation_general_sens_all.pdf"))
path1<-file.path(paste0("QCCorrelation_general_sens_all.pdf"))
pdf(path1)
print(plot.cor.BETA)
dev.off()

#  SAME CORRELATION given by the THE COR.TEST FUNCTION THAN IN THE PDF RESULTING FROM GGALLY package

############
###GENERAL MAIN vs SENSITIVITY###
### CpGs that are SIG with PVAL< 10-4 (IN THE MAIN MODEL)

setwd("/home/isglobal.lan/mcosin/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/COMPARE_MODELS/marta_28072022/")

gen_sens<-read.table("/home/isglobal.lan/mcosin/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220707_Output/SENSITIVITY_NOU/QC_ResultsSENSITIVITY/general/general_QCData.txt", header=TRUE)
head(gen_sens)

gen_all<-read.table("/home/isglobal.lan/mcosin/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/general/general_QCData.txt", header=TRUE)

head(gen_all)

gen_all4<-gen_all[gen_all$P_VAL<0.0001,]

probes<-gen_all4$probeID #probes from the general model that are sig by pval*10^-4

gen_sens<-gen_sens[gen_sens$probeID %in% probes,] #(has seleccionat les cpg's dins el gen_sens que donen sig a la menys 4 al MAIN)

dim(gen_sens)
#103

sensord<-gen_sens[order(gen_sens$probeID),]
genallord<-gen_all4[order(gen_all4$probeID),]

sensordb<- sensord[,c("probeID","BETA")]
genallordb<- genallord[,c("probeID","BETA")]

colnames(sensordb)<-c("probeID","BETA_gen_sensitivity")
colnames(genallordb)<-c("probeID","BETA_main_gen_all")

corr_matrix<-merge(sensordb, genallordb,by="probeID") 

cor.test(corr_matrix$BETA_gen_sensitivity, corr_matrix$BETA_main_gen_all)

# 0.9915


library(GGally)

corr_matrix2<-corr_matrix[,2:3]

plot.cor.BETA<-ggpairs(corr_matrix2, title="Correlation between main general sig-4 vs general sensitivity sig-4 ")
file.create(paste0("QCCorrelation_general_sens_sig4.pdf"))
path1<-file.path(paste0("QCCorrelation_general_sens_sig4.pdf"))
pdf(path1)
print(plot.cor.BETA)
dev.off()

#  SAME CORRELATION given by the THE COR.TEST FUNCTION THAN IN THE PDF RESULTING FROM GGALLY package


############
###GENERAL MAIN vs Mat_iq and model without sgest###
### ALL

setwd("/home/isglobal.lan/mcosin/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/COMPARE_MODELS/marta_28072022/")

gen_mat<-read.table("/home/isglobal.lan/mcosin/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/INMA_20220702_Output/RESULTATS_iQ/QC_ResultsmaternalIQ/general/general_QCData.txt", header=TRUE)
head(gen_mat)
gen_sges<- read.table("/home/isglobal.lan/mcosin/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220627_Output/RESULTATS_sges_NOU/QC_ResultsSGES/general/general_QCData.txt", header=TRUE)
head(gen_sges)

gen_all<-read.table("/home/isglobal.lan/mcosin/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/general/general_QCData.txt", header=TRUE)

head(gen_all)

probes<-gen_all$probeID #probes from the general model that are sig by pval*10^-4

genmatord<-gen_mat[order(gen_mat$probeID),]
gensgesord<-gen_sges[order(gen_sges$probeID),]
genallord<-gen_all[order(gen_all$probeID),]

genmatordb<- genmatord[,c("probeID","BETA")]
gensgesordb<- gensgesord[,c("probeID","BETA")]
genallordb<- genallord[,c("probeID","BETA")]

colnames(genmatordb)<-c("probeID","BETA_gen_maternal_iq")
colnames(gensgesordb)<-c("probeID","BETA_gen_wo_sges")
colnames(genallordb)<-c("probeID","BETA_main_gen_all")

corr_matrix<-merge(genmatordb, genallordb,by="probeID") 
corr_matrix2<-merge(gensgesordb, corr_matrix,by="probeID") 
head(corr_matrix2)

cor.test(corr_matrix2$BETA_gen_maternal_iq, corr_matrix2$BETA_main_gen_all)

#Pearson's product-moment correlation

#0.965

cor.test(corr_matrix2$BETA_gen_wo_sges, corr_matrix2$BETA_main_gen_all)

#Pearson's product-moment correlation

#1 ## la comparació MAIN amb model wo sges DÓNA 1 PERQUÈ S'HA DE REFER EL QC DEL MODEL SENSE GESTATIONAL AGE I CARREGAR DE NOU EN AQUEST CODI L'ARXIU ACTUALITZAT

library(GGally)

corr_matrix2<-corr_matrix2[,2:4]

plot.cor.BETA<-ggpairs(corr_matrix2, title="Correlation between main general all vs general maternal iq all and wo gestational age model ")
file.create(paste0("QCCorrelation_general_matiq_sges_all.pdf"))
path1<-file.path(paste0("QCCorrelation_general_matiq_sges_all.pdf"))
pdf(path1)
print(plot.cor.BETA)
dev.off()

#  SAME CORRELATION given by the THE COR.TEST FUNCTION THAN IN THE PDF RESULTING FROM GGALLY package

############
###GENERAL MAIN vs Mat_iq and model without sgest###
### CpGs that are SIG with PVAL< 10-4 (IN THE MAIN MODEL)

setwd("/home/isglobal.lan/mcosin/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/COMPARE_MODELS/marta_28072022/")

gen_all<-read.table("/home/isglobal.lan/mcosin/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/general/general_QCData.txt", header=TRUE)

head(gen_all)

gen_all4<-gen_all[gen_all$P_VAL<0.0001,]

probes<-gen_all4$probeID #probes from the general model that are sig by pval*10^-4

gen_mat4<-gen_mat[gen_mat$probeID %in% probes,] #(has seleccionat les cpg's dins el gen_sens que donen sig a la menys 4 al MAIN)
gen_sges4<-gen_sges[gen_sges$probeID %in% probes,]

dim(gen_sges4)
#103
dim(gen_mat4)
#103
dim(gen_all4)
#103

genmatord<-gen_mat4[order(gen_mat4$probeID),]
gensgesord<-gen_sges4[order(gen_sges4$probeID),]
genallord<-gen_all4[order(gen_all4$probeID),]

genmatordb<- genmatord[,c("probeID","BETA")]
gensgesordb<- gensgesord[,c("probeID","BETA")]
genallordb<- genallord[,c("probeID","BETA")]

colnames(genmatordb)<-c("probeID","BETA_gen_maternal_iq_4")
colnames(gensgesordb)<-c("probeID","BETA_gen_wo_sges_4")
colnames(genallordb)<-c("probeID","BETA_main_gen_4")

corr_matrix<-merge(genmatordb, genallordb,by="probeID") 
corr_matrix2<-merge(gensgesordb, corr_matrix,by="probeID") 
head(corr_matrix2)

cor.test(corr_matrix2$BETA_gen_maternal_iq_4, corr_matrix2$BETA_main_gen_4)

#Pearson's product-moment correlation

#0.995

cor.test(corr_matrix2$BETA_gen_wo_sges_4, corr_matrix2$BETA_main_gen_4)

#Pearson's product-moment correlation

#1 ## la comparació MAIN amb model wo sges DÓNA 1 PERQUÈ S'HA DE REFER EL QC DEL MODEL SENSE GESTATIONAL AGE I CARREGAR DE NOU EN AQUEST CODI L'ARXIU ACTUALITZAT

library(GGally)

corr_matrix2<-corr_matrix2[,2:4]

plot.cor.BETA<-ggpairs(corr_matrix2, title="Correlation between main general sig -4 vs general maternal iq and wo gestational age model ")
file.create(paste0("QCCorrelation_general_matiq_sges_sig4.pdf"))
path1<-file.path(paste0("QCCorrelation_general_matiq_sges_sig4.pdf"))
pdf(path1)
print(plot.cor.BETA)
dev.off()

#  SAME CORRELATION given by the THE COR.TEST FUNCTION THAN IN THE PDF RESULTING FROM GGALLY package




