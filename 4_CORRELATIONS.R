
#################################
###GENERAL MAIN vs SENSITIVITY###
#################################


### ALL

setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/COMPARE_MODELS/net_NEW/")

gen_sens<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220707_Output/SENSITIVITY_NOU/QC_ResultsSENSITIVITY/general/general_QCData.txt", header=TRUE)
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

gen_all<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/general/general_QCData.txt", header=TRUE)

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

probes<-gen_all$probeID 

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

#################################
###GENERAL MAIN vs SENSITIVITY###
#################################

### CpGs that are SIG with PVAL< 10-4 (IN THE MAIN MODEL)

setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/COMPARE_MODELS/net_NEW/")

gen_sens<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220707_Output/SENSITIVITY_NOU/QC_ResultsSENSITIVITY/general/general_QCData.txt", header=TRUE)
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


gen_all<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/general/general_QCData.txt", header=TRUE)

gen_all4<-gen_all[gen_all$P_VAL<0.0001,]

probes<-gen_all4$probeID #probes from the general model that are sig by pval*10^-4

gen_sens<-gen_sens[gen_sens$probeID %in% probes,]

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

# Pearson's product-moment correlation

#data:  corr_matrix$BETA_gen_sensitivity and corr_matrix$BETA_main_gen_all
#t = 76.95, df = 101, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.9875627 0.9943021
#sample estimates:
#  cor
# 0.991579



library(GGally)

corr_matrix2<-corr_matrix[,2:3]

plot.cor.BETA<-ggpairs(corr_matrix2, title="Correlation between main general sig-4 vs general sensitivity sig-4 ")
file.create(paste0("QCCorrelation_general_sens_sig4.pdf"))
path1<-file.path(paste0("QCCorrelation_general_sens_sig4.pdf"))
pdf(path1)
print(plot.cor.BETA)
dev.off()

#  SAME CORRELATION given by the THE COR.TEST FUNCTION THAN IN THE PDF RESULTING FROM GGALLY package


####################################################
###GENERAL MAIN vs Mat_iq and model without sgest###
####################################################

### ALL

setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/COMPARE_MODELS/net_NEW/")

gen_mat<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220601NOUIQ_Output/Resultats_matIQ/QC_ResultsmaternalIQ/general/general_QCData.txt", header=TRUE)
head(gen_mat)

#probeID          BETA           SE        P_VAL  padj.fdr padj.bonf
#1 cg02986379 -6.149673e-05 1.203959e-05 3.258022e-07 0.2307022        no
#2 cg15480200 -1.777487e-03 3.598521e-04 7.832551e-07 0.2773134        no
#3 cg13550802 -1.496209e-03 3.139744e-04 1.884928e-06 0.3342624        no
#4 cg22965956 -1.092937e-03 2.295088e-04 1.916178e-06 0.3342624        no
#5 cg06057918 -1.047108e-04 2.218522e-05 2.360260e-06 0.3342624        no
#6 cg21294960 -1.044911e-03 2.244991e-04 3.249050e-06 0.3834448        no
#CpG_chrm   CpG_beg   CpG_end
#1     chr5  39425327  39425329
#2    chr10 132817110 132817112
#3    chr10 133122349 133122351
#4    chr18  74212614  74212616
#5     chr6  22069864  22069866
#6     chr2 207124227 207124229



gen_sges<- read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220722NOU_Output/resultats_nous_sges/QC_ResultsSGES/general/general_QCData.txt", header=TRUE)
head(gen_sges)
#probeID          BETA           SE        P_VAL   padj.fdr padj.bonf
#1 cg02986379 -6.106364e-05 0.0000115379 1.206965e-07 0.08546577        no
#2 cg00866476  5.444129e-04 0.0001072825 3.883641e-07 0.09982773        no
#3 cg15480200 -1.744898e-03 0.0003479216 5.297693e-07 0.09982773        no
#4 cg14613533 -9.384119e-04 0.0001878776 5.889548e-07 0.09982773        no
#5 cg22965956 -1.103393e-03 0.0002224577 7.048936e-07 0.09982773        no
#6 cg13550802 -1.442138e-03 0.0002933113 8.799283e-07 0.10384694        no
#CpG_chrm   CpG_beg   CpG_end
#1     chr5  39425327  39425329
#2     chr2  95074936  95074938
#3    chr10 132817110 132817112
#4    chr17   5192674   5192676
#5    chr18  74212614  74212616
#6    chr10 133122349 133122351

gen_all<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/general/general_QCData.txt", header=TRUE)

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
#probeID BETA_gen_wo_sges BETA_gen_maternal_iq BETA_main_gen_all
#1 cg00000029    -3.809294e-05        -1.027170e-05     -3.609285e-05
#2 cg00000103    -5.930111e-05        -7.388988e-05     -1.236557e-04
#3 cg00000109     1.074353e-06        -1.276734e-05      1.837904e-06
#4 cg00000155    -1.015392e-05        -6.105371e-06     -9.723200e-06
#5 cg00000158    -1.760491e-05        -1.860305e-05     -1.773702e-05
#6 cg00000165    -3.934910e-04        -2.028160e-04     -3.959174e-04



cor.test(corr_matrix2$BETA_gen_maternal_iq, corr_matrix2$BETA_main_gen_all)

#Pearson's product-moment correlation

#data:  corr_matrix2$BETA_gen_maternal_iq and corr_matrix2$BETA_main_gen_all
#t = 3079.4, df = 708103, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9644711 0.9647948
#sample estimates:
#      cor
#0.9646333


cor.test(corr_matrix2$BETA_gen_wo_sges, corr_matrix2$BETA_main_gen_all)

#Pearson's product-moment correlation

#data:  corr_matrix2$BETA_gen_wo_sges and corr_matrix2$BETA_main_gen_all
#t = 9761, df = 708103, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9962873 0.9963217
#sample estimates:
#      cor
#0.9963045

library(GGally)

corr_matrix2<-corr_matrix2[,2:4]

plot.cor.BETA<-ggpairs(corr_matrix2, title="Correlation between main general all vs general maternal iq all and wo gestational age model ")
file.create(paste0("QCCorrelation_general_matiq_sges_all.pdf"))
path1<-file.path(paste0("QCCorrelation_general_matiq_sges_all.pdf"))
pdf(path1)
print(plot.cor.BETA)
dev.off()

#  SAME CORRELATION given by the THE COR.TEST FUNCTION THAN IN THE PDF RESULTING FROM GGALLY package

####################################################
###GENERAL MAIN vs Mat_iq and model without sgest###
####################################################

### CpGs that are SIG with PVAL< 10-4 (IN THE MAIN MODEL)

setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/COMPARE_MODELS/net_NEW/")

gen_all<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/general/general_QCData.txt", header=TRUE)

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


gen_all4<-gen_all[gen_all$P_VAL<0.0001,]

probes<-gen_all4$probeID #probes from the general model that are sig by pval*10^-4

gen_mat4<-gen_mat[gen_mat$probeID %in% probes,]
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
#probeID BETA_gen_wo_sges_4 BETA_gen_maternal_iq_4 BETA_main_gen_4
#1 cg00169184      -0.0001024811          -1.104908e-04   -1.034933e-04
#2 cg00647233      -0.0010772402          -9.888902e-04   -1.087627e-03
#3 cg00734567      -0.0004655130          -4.058831e-04   -4.997743e-04
#4 cg00866476       0.0005444129           4.929689e-04    5.480824e-04
#5 cg01153132       0.0000925057           9.788733e-05    9.243933e-05
#6 cg01157470      -0.0008234726          -7.892696e-04   -8.345778e-04



cor.test(corr_matrix2$BETA_gen_maternal_iq_4, corr_matrix2$BETA_main_gen_4)

#Pearson's product-moment correlation

#data:  corr_matrix2$BETA_gen_maternal_iq_4 and corr_matrix2$BETA_main_gen_4
#t = 106.76, df = 101, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9934934 0.9970239
#sample estimates:
#      cor
#0.9955988



cor.test(corr_matrix2$BETA_gen_wo_sges_4, corr_matrix2$BETA_main_gen_4)
#Pearson's product-moment correlation

#data:  corr_matrix2$BETA_gen_wo_sges_4 and corr_matrix2$BETA_main_gen_4
#t = 427.59, df = 101, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9995914 0.9998134
#sample estimates:
#      cor
#0.9997239




library(GGally)

corr_matrix2<-corr_matrix2[,2:4]

plot.cor.BETA<-ggpairs(corr_matrix2, title="Correlation between main general sig -4 vs general maternal iq and wo gestational age model ")
file.create(paste0("QCCorrelation_general_matiq_sges_sig4.pdf"))
path1<-file.path(paste0("QCCorrelation_general_matiq_sges_sig4.pdf"))
pdf(path1)
print(plot.cor.BETA)
dev.off()

#  SAME CORRELATION given by the THE COR.TEST FUNCTION THAN IN THE PDF RESULTING FROM GGALLY package




################################
###verbal MAIN vs SENSITIVITY###
################################

### ALL

setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/COMPARE_MODELS/net_NEW/")

verb_sens<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220707_Output/SENSITIVITY_NOU/QC_ResultsSENSITIVITY/verbal/verbal_QCData.txt", header=TRUE)
head(verb_sens)
#probeID          BETA           SE        P_VAL     padj.fdr padj.bonf
#1 cg01566028 -3.669924e-04 5.884537e-05 4.473082e-10 0.0003167412       yes
#2 cg25818763 -1.180464e-03 2.124903e-04 2.770110e-08 0.0098076450       yes
#3 cg00856839  3.281992e-05 6.113878e-06 7.956938e-08 0.0143656525        no
#4 cg13073105  4.494284e-04 8.392733e-05 8.557004e-08 0.0143656525        no
#5 cg26106417 -2.861690e-04 5.374943e-05 1.014373e-07 0.0143656525        no
#6 cg15480200 -1.785643e-03 3.395937e-04 1.454925e-07 0.0171706589        no
#CpG_chrm   CpG_beg   CpG_end
#1     chr9 120578988 120578990
#2    chr18  52345643  52345645
#3     chr9 109780162 109780164
#4     chr1 201355349 201355351
#5     chr4   3423653   3423655
#6    chr10 132817110 132817112


verb_all<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/verbal/verbal_QCData.txt", header=TRUE)

head(verb_all)
#probeID          BETA           SE        P_VAL  padj.fdr padj.bonf
#1 cg15480200 -0.0017058791 3.282710e-04 2.030122e-07 0.1025309        no
#2 cg01566028 -0.0002960664 5.772491e-05 2.914105e-07 0.1025309        no
#3 cg16897193 -0.0015097311 2.987653e-04 4.343885e-07 0.1025309        no
#4 cg01572157 -0.0022572901 4.576186e-04 8.110547e-07 0.1241466        no
#5 cg22965956 -0.0010666854 2.180143e-04 9.944656e-07 0.1241466        no
#6 cg01725682 -0.0002951650 6.077573e-05 1.194028e-06 0.1241466        no
#CpG_chrm   CpG_beg   CpG_end
#1    chr10 132817110 132817112
#2     chr9 120578988 120578990
#3    chr19  45940542  45940544
#4     chr7 125438150 125438152
#5    chr18  74212614  74212616
#6    chr17  50359485  50359487



probes<-verb_all$probeID #probes from the verbal model that are sig by pval*10^-4

sensord<-verb_sens[order(verb_sens$probeID),]
verballord<-verb_all[order(verb_all$probeID),]

sensordb<- sensord[,c("probeID","BETA")]
verballordb<- verballord[,c("probeID","BETA")]

colnames(sensordb)<-c("probeID","BETA_verb_sensitivity")
colnames(verballordb)<-c("probeID","BETA_main_verb_all")

corr_matrix<-merge(sensordb, verballordb,by="probeID")
head(corr_matrix)
#probeID BETA_verb_sensitivity BETA_main_verb_all
#1 cg00000029         -4.963895e-05      -6.043817e-05
#2 cg00000103         -3.788326e-04      -1.913580e-04
#3 cg00000109         -8.885526e-06       1.477226e-05
#4 cg00000155         -1.363407e-05      -7.754261e-06
#5 cg00000158         -1.641456e-05      -2.657303e-06
#6 cg00000165         -3.740744e-04      -1.172157e-04



cor.test(corr_matrix$BETA_verb_sensitivity, corr_matrix$BETA_main_verb_all)
#Pearson's product-moment correlation

#data:  corr_matrix$BETA_verb_sensitivity and corr_matrix$BETA_main_verb_all
#t = 1963.2, df = 708103, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9187646 0.9194876
#sample estimates:
#      cor
#0.9191269


library(GGally)

corr_matrix2<-corr_matrix[,2:3]

plot.cor.BETA<-ggpairs(corr_matrix2, title="Correlation between main verbal all vs verbal sensitivity all ")
file.create(paste0("QCCorrelation_verbal_sens_all.pdf"))
path1<-file.path(paste0("QCCorrelation_verbal_sens_all.pdf"))
pdf(path1)
print(plot.cor.BETA)
dev.off()

#  SAME CORRELATION given by the THE COR.TEST FUNCTION THAN IN THE PDF RESULTING FROM GGALLY package

############
###verbal MAIN vs SENSITIVITY###
### CpGs that are SIG with PVAL< 10-4 (IN THE MAIN MODEL)

setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/COMPARE_MODELS/net_NEW/")

verb_sens<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220707_Output/SENSITIVITY_NOU/QC_ResultsSENSITIVITY/verbal/verbal_QCData.txt", header=TRUE)
head(verb_sens)
#probeID          BETA           SE        P_VAL     padj.fdr padj.bonf
#1 cg01566028 -3.669924e-04 5.884537e-05 4.473082e-10 0.0003167412       yes
#2 cg25818763 -1.180464e-03 2.124903e-04 2.770110e-08 0.0098076450       yes
#3 cg00856839  3.281992e-05 6.113878e-06 7.956938e-08 0.0143656525        no
#4 cg13073105  4.494284e-04 8.392733e-05 8.557004e-08 0.0143656525        no
#5 cg26106417 -2.861690e-04 5.374943e-05 1.014373e-07 0.0143656525        no
#6 cg15480200 -1.785643e-03 3.395937e-04 1.454925e-07 0.0171706589        no
#CpG_chrm   CpG_beg   CpG_end
#1     chr9 120578988 120578990
#2    chr18  52345643  52345645
#3     chr9 109780162 109780164
#4     chr1 201355349 201355351
#5     chr4   3423653   3423655
#6    chr10 132817110 132817112



verb_all<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/verbal/verbal_QCData.txt", header=TRUE)

head(verb_all)
#probeID          BETA           SE        P_VAL  padj.fdr padj.bonf
#1 cg15480200 -0.0017058791 3.282710e-04 2.030122e-07 0.1025309        no
#2 cg01566028 -0.0002960664 5.772491e-05 2.914105e-07 0.1025309        no
#3 cg16897193 -0.0015097311 2.987653e-04 4.343885e-07 0.1025309        no
#4 cg01572157 -0.0022572901 4.576186e-04 8.110547e-07 0.1241466        no
#5 cg22965956 -0.0010666854 2.180143e-04 9.944656e-07 0.1241466        no
#6 cg01725682 -0.0002951650 6.077573e-05 1.194028e-06 0.1241466        no
#CpG_chrm   CpG_beg   CpG_end
#1    chr10 132817110 132817112
#2     chr9 120578988 120578990
#3    chr19  45940542  45940544
#4     chr7 125438150 125438152
#5    chr18  74212614  74212616
#6    chr17  50359485  50359487


verb_all4<-verb_all[verb_all$P_VAL<0.0001,]

probes<-verb_all4$probeID #probes from the verbal model that are sig by pval*10^-4

verb_sens<-verb_sens[verb_sens$probeID %in% probes,]

dim(verb_sens)
#103

sensord<-verb_sens[order(verb_sens$probeID),]
verballord<-verb_all4[order(verb_all4$probeID),]

sensordb<- sensord[,c("probeID","BETA")]
verballordb<- verballord[,c("probeID","BETA")]

colnames(sensordb)<-c("probeID","BETA_verb_sensitivity")
colnames(verballordb)<-c("probeID","BETA_main_verb_all")

corr_matrix<-merge(sensordb, verballordb,by="probeID") 
head(corr_matrix)
#probeID BETA_verb_sensitivity BETA_main_verb_all
#1 cg00295161         -1.412009e-04      -1.584223e-04
#2 cg00431874          3.247535e-05       3.250708e-05
#3 cg00470352         -5.087538e-05      -5.977362e-05
#4 cg00590134         -5.497319e-04      -4.556521e-04
#5 cg00618746         -4.196034e-04      -4.144913e-04
#6 cg00856839          3.281992e-05       2.917391e-05


cor.test(corr_matrix$BETA_verb_sensitivity, corr_matrix$BETA_main_verb_all)

# Pearson's product-moment correlation

#data:  corr_matrix$BETA_verb_sensitivity and corr_matrix$BETA_main_verb_all
#t = 84.693, df = 170, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.9842908 0.9913743
#sample estimates:
#  cor
#0.9883563



library(GGally)

corr_matrix2<-corr_matrix[,2:3]

plot.cor.BETA<-ggpairs(corr_matrix2, title="Correlation between main verbal sig-4 vs verbal sensitivity sig-4 ")
file.create(paste0("QCCorrelation_verbal_sens_sig4.pdf"))
path1<-file.path(paste0("QCCorrelation_verbal_sens_sig4.pdf"))
pdf(path1)
print(plot.cor.BETA)
dev.off()

#  SAME CORRELATION given by the THE COR.TEST FUNCTION THAN IN THE PDF RESULTING FROM GGALLY package


##################################################
###verbal MAIN vs Mat_iq and model without sges###
##################################################

### ALL

setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/COMPARE_MODELS/net_NEW/")

verb_mat<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220601NOUIQ_Output/Resultats_matIQ/QC_ResultsmaternalIQ/verbal/verbal_QCData.txt", header=TRUE)
head(verb_mat)
#1 cg00856839  3.247666e-05 6.118705e-06 1.109762e-07 0.07258772        no
#2 cg16897193 -1.443719e-03 2.779201e-04 2.050196e-07 0.07258772        no
#3 cg10152475 -1.628104e-03 3.329455e-04 1.008352e-06 0.14335210        no
#4 cg15480200 -1.658546e-03 3.392604e-04 1.014989e-06 0.14335210        no
#5 cg25677261 -3.450604e-04 7.127096e-05 1.288440e-06 0.14335210        no
#6 cg03226208 -8.289710e-04 1.712437e-04 1.292616e-06 0.14335210        no
#CpG_chrm   CpG_beg   CpG_end
#1     chr9 109780162 109780164
#2    chr19  45940542  45940544
#3    chr12 117860349 117860351
#4    chr10 132817110 132817112
#5    chr15 101469650 101469652
#6    chr17   5100887   5100889


verb_sges<- read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220722NOU_Output/resultats_nous_sges/QC_ResultsSGES/verbal/verbal_QCData.txt", header=TRUE)
head(verb_sges)
#probeID          BETA           SE        P_VAL  padj.fdr padj.bonf
#1 cg15480200 -1.691008e-03 3.322093e-04 3.577108e-07 0.1122911        no
#2 cg01566028 -2.958284e-04 5.902807e-05 5.396335e-07 0.1122911        no
#3 cg01572157 -2.257658e-03 4.529570e-04 6.219770e-07 0.1122911        no
#4 cg22965956 -1.070893e-03 2.150187e-04 6.343192e-07 0.1122911        no
#5 cg00856839  2.958225e-05 6.158409e-06 1.558736e-06 0.1223682        no
#6 cg18059802 -8.436252e-05 1.758898e-05 1.615998e-06 0.1223682        no
#CpG_chrm   CpG_beg   CpG_end
#1    chr10 132817110 132817112
#2     chr9 120578988 120578990
#3     chr7 125438150 125438152
#4    chr18  74212614  74212616
#5     chr9 109780162 109780164
#6    chr12  50953274  50953276


verb_all<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/verbal/verbal_QCData.txt", header=TRUE)

head(verb_all)
#probeID          BETA           SE        P_VAL  padj.fdr padj.bonf
#1 cg15480200 -0.0017058791 3.282710e-04 2.030122e-07 0.1025309        no
#2 cg01566028 -0.0002960664 5.772491e-05 2.914105e-07 0.1025309        no
#3 cg16897193 -0.0015097311 2.987653e-04 4.343885e-07 0.1025309        no
#4 cg01572157 -0.0022572901 4.576186e-04 8.110547e-07 0.1241466        no
#5 cg22965956 -0.0010666854 2.180143e-04 9.944656e-07 0.1241466        no
#6 cg01725682 -0.0002951650 6.077573e-05 1.194028e-06 0.1241466        no
#CpG_chrm   CpG_beg   CpG_end
#1    chr10 132817110 132817112
#2     chr9 120578988 120578990
#3    chr19  45940542  45940544
#4     chr7 125438150 125438152
#5    chr18  74212614  74212616
#6    chr17  50359485  50359487



probes<-verb_all$probeID #probes from the verbal model that are sig by pval*10^-4

verbmatord<-verb_mat[order(verb_mat$probeID),]
verbsgesord<-verb_sges[order(verb_sges$probeID),]
verballord<-verb_all[order(verb_all$probeID),]

verbmatordb<- verbmatord[,c("probeID","BETA")]
verbsgesordb<- verbsgesord[,c("probeID","BETA")]
verballordb<- verballord[,c("probeID","BETA")]

colnames(verbmatordb)<-c("probeID","BETA_verb_maternal_iq")
colnames(verbsgesordb)<-c("probeID","BETA_verb_wo_sges")
colnames(verballordb)<-c("probeID","BETA_main_verb_all")

corr_matrix<-merge(verbmatordb, verballordb,by="probeID") 
corr_matrix2<-merge(verbsgesordb, corr_matrix,by="probeID") 
head(corr_matrix2)
#probeID BETA_verb_wo_sges BETA_verb_maternal_iq BETA_main_verb_all
#1 cg00000029     -6.094347e-05         -4.243054e-05      -6.043817e-05
#2 cg00000103     -1.840856e-04         -1.745696e-04      -1.913580e-04
#3 cg00000109      1.530694e-05          3.889829e-06       1.477226e-05
#4 cg00000155     -7.594606e-06         -6.970487e-06      -7.754261e-06
#5 cg00000158     -2.651698e-06         -2.535844e-06      -2.657303e-06
#6 cg00000165     -1.182384e-04          3.994567e-05      -1.172157e-04



cor.test(corr_matrix2$BETA_verb_maternal_iq, corr_matrix2$BETA_main_verb_all)

#Pearson's product-moment correlation

#data:  corr_matrix2$BETA_verb_maternal_iq and corr_matrix2$BETA_main_verb_all
#t = 4035.1, df = 708103, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9788426 0.9790368
#sample estimates:
#      cor
#0.9789399


cor.test(corr_matrix2$BETA_verb_wo_sges, corr_matrix2$BETA_main_verb_all)
#Pearson's product-moment correlation

#data:  corr_matrix2$BETA_verb_wo_sges and corr_matrix2$BETA_main_verb_all
#t = 21344, df = 708103, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9992201 0.9992273
#sample estimates:
#      cor
#0.9992237



library(GGally)

corr_matrix2<-corr_matrix2[,2:4]

plot.cor.BETA<-ggpairs(corr_matrix2, title="Correlation between main verbal all vs verbal maternal iq all and wo gestational age model ")
file.create(paste0("QCCorrelation_verbal_matiq_sges_all.pdf"))
path1<-file.path(paste0("QCCorrelation_verbal_matiq_sges_all.pdf"))
pdf(path1)
print(plot.cor.BETA)
dev.off()

#  SAME CORRELATION given by the THE COR.TEST FUNCTION THAN IN THE PDF RESULTING FROM GGALLY package

###################################################
###verbal MAIN vs Mat_iq and model without sgest###
###################################################

### CpGs that are SIG with PVAL< 10-4 (IN THE MAIN MODEL)

setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/COMPARE_MODELS/net_NEW/")

verb_all<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/verbal/verbal_QCData.txt", header=TRUE)

head(verb_all)
#probeID          BETA           SE        P_VAL  padj.fdr padj.bonf
#1 cg15480200 -0.0017058791 3.282710e-04 2.030122e-07 0.1025309        no
#2 cg01566028 -0.0002960664 5.772491e-05 2.914105e-07 0.1025309        no
#3 cg16897193 -0.0015097311 2.987653e-04 4.343885e-07 0.1025309        no
#4 cg01572157 -0.0022572901 4.576186e-04 8.110547e-07 0.1241466        no
#5 cg22965956 -0.0010666854 2.180143e-04 9.944656e-07 0.1241466        no
#6 cg01725682 -0.0002951650 6.077573e-05 1.194028e-06 0.1241466        no
#CpG_chrm   CpG_beg   CpG_end
#1    chr10 132817110 132817112
#2     chr9 120578988 120578990
#3    chr19  45940542  45940544
#4     chr7 125438150 125438152
#5    chr18  74212614  74212616
#6    chr17  50359485  50359487


verb_all4<-verb_all[verb_all$P_VAL<0.0001,]

probes<-verb_all4$probeID #probes from the verbal model that are sig by pval*10^-4

verb_mat4<-verb_mat[verb_mat$probeID %in% probes,]
verb_sges4<-verb_sges[verb_sges$probeID %in% probes,]

dim(verb_sges4)
#103
dim(verb_mat4)
#103
dim(verb_all4)
#103

verbmatord<-verb_mat4[order(verb_mat4$probeID),]
verbsgesord<-verb_sges4[order(verb_sges4$probeID),]
verballord<-verb_all4[order(verb_all4$probeID),]

verbmatordb<- verbmatord[,c("probeID","BETA")]
verbsgesordb<- verbsgesord[,c("probeID","BETA")]
verballordb<- verballord[,c("probeID","BETA")]

colnames(verbmatordb)<-c("probeID","BETA_verb_maternal_iq_4")
colnames(verbsgesordb)<-c("probeID","BETA_verb_wo_sges_4")
colnames(verballordb)<-c("probeID","BETA_main_verb_4")

corr_matrix<-merge(verbmatordb, verballordb,by="probeID") 
corr_matrix2<-merge(verbsgesordb, corr_matrix,by="probeID") 
head(corr_matrix2)
#probeID BETA_verb_wo_sges_4 BETA_verb_maternal_iq_4 BETA_main_verb_4
#1 cg00295161       -1.584154e-04           -1.430049e-04    -1.584223e-04
#2 cg00431874        3.259013e-05            3.195708e-05     3.250708e-05
#3 cg00470352       -5.987081e-05           -5.704357e-05    -5.977362e-05
#4 cg00590134       -4.532358e-04           -4.194969e-04    -4.556521e-04
#5 cg00618746       -4.169330e-04           -3.860378e-04    -4.144913e-04
#6 cg00856839        2.958225e-05            3.247666e-05     2.917391e-05


cor.test(corr_matrix2$BETA_verb_maternal_iq_4, corr_matrix2$BETA_main_verb_4)

#Pearson's product-moment correlation

#data:  corr_matrix2$BETA_verb_maternal_iq_4 and corr_matrix2$BETA_main_verb_4
#t = 178.7, df = 170, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.9964176 0.9980384
#sample estimates:
#  cor
#0.997349


cor.test(corr_matrix2$BETA_verb_wo_sges_4, corr_matrix2$BETA_main_verb_4)
#Pearson's product-moment correlation

#data:  corr_matrix2$BETA_verb_wo_sges_4 and corr_matrix2$BETA_main_verb_4
#t = 1376.6, df = 170, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.9999394 0.9999668
#sample estimates:
#  cor
#0.9999551



library(GGally)

corr_matrix2<-corr_matrix2[,2:4]

plot.cor.BETA<-ggpairs(corr_matrix2, title="Correlation between main verbal sig -4 vs verbal maternal iq and wo gestational age model ")
file.create(paste0("QCCorrelation_verbal_matiq_sges_sig4.pdf"))
path1<-file.path(paste0("QCCorrelation_verbal_matiq_sges_sig4.pdf"))
pdf(path1)
print(plot.cor.BETA)
dev.off()

#  SAME CORRELATION given by the THE COR.TEST FUNCTION THAN IN THE PDF RESULTING FROM GGALLY package





###################################
###nonverbal MAIN vs SENSITIVITY###
###################################

### ALL

setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/COMPARE_MODELS/net_NEW/")

nonverb_sens<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220707_Output/SENSITIVITY_NOU/QC_ResultsSENSITIVITY/nonverbal/nonverbal_QCData.txt", header=TRUE)
head(nonverb_sens)
#probeID          BETA           SE        P_VAL   padj.fdr padj.bonf
#1 cg05742591 -1.566119e-03 3.027884e-04 2.312045e-07 0.09954241        no
#2 cg21142557  3.547450e-04 6.907483e-05 2.811515e-07 0.09954241        no
#3 cg01153132  1.137337e-04 2.275556e-05 5.790933e-07 0.10327737        no
#4 cg25096679 -1.493898e-04 3.021212e-05 7.626220e-07 0.10327737        no
#5 cg05409601 -6.534902e-05 1.322470e-05 7.754957e-07 0.10327737        no
#6 cg16495530  2.098862e-04 4.275209e-05 9.136493e-07 0.10327737        no
#CpG_chrm   CpG_beg   CpG_end
#1     chr5   1656679   1656681
#2     chr1 204401133 204401135
#3     chr6 119078948 119078950
#4    chr10 101841248 101841250
#5     chr4  77156758  77156760
#6    chr11 122726318 122726320

nonverb_all<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/nonverbal/nonverbal_QCData.txt", header=TRUE)

head(nonverb_all)
#probeID          BETA           SE        P_VAL  padj.fdr padj.bonf
#1 cg18892718 -1.934287e-03 3.761491e-04 2.713358e-07 0.1056628        no
#2 cg11995801 -4.945211e-05 9.650265e-06 2.984383e-07 0.1056628        no
#3 cg15110918 -2.337091e-03 4.731039e-04 7.815800e-07 0.1844802        no
#4 cg22034155 -9.322172e-05 1.998282e-05 3.084759e-06 0.2593288        no
#5 cg15856043  6.664539e-04 1.431903e-04 3.250469e-06 0.2593288        no
#6 cg00072720 -1.384769e-04 2.978712e-05 3.337357e-06 0.2593288        no
#CpG_chrm   CpG_beg   CpG_end
#1     chr7  87629447  87629449
#2     chr2 231961564 231961566
#3    chr11 126655602 126655604
#4     chr6  31835421  31835423
#5    chr12 118542075 118542077
#6    chr17   7262511   7262513



probes<-nonverb_all$probeID #probes from the nonverbal model that are sig by pval*10^-4

sensord<-nonverb_sens[order(nonverb_sens$probeID),]
nonverballord<-nonverb_all[order(nonverb_all$probeID),]

sensordb<- sensord[,c("probeID","BETA")]
nonverballordb<- nonverballord[,c("probeID","BETA")]

colnames(sensordb)<-c("probeID","BETA_nonverb_sensitivity")
colnames(nonverballordb)<-c("probeID","BETA_main_nonverb_all")

corr_matrix<-merge(sensordb, nonverballordb,by="probeID") 
head(corr_matrix)
#probeID BETA_nonverb_sensitivity BETA_main_nonverb_all
#1 cg00000029            -3.284833e-07         -3.942379e-08
#2 cg00000103             7.992723e-05          8.419662e-05
#3 cg00000109            -2.390160e-05         -1.185585e-05
#4 cg00000155            -2.178270e-05         -6.959237e-06
#5 cg00000158            -2.912062e-05         -3.234250e-05
#6 cg00000165            -6.865671e-04         -5.890364e-04


cor.test(corr_matrix$BETA_nonverb_sensitivity, corr_matrix$BETA_main_nonverb_all)
#Pearson's product-moment correlation

#data:  corr_matrix$BETA_nonverb_sensitivity and corr_matrix$BETA_main_nonverb_all
#t = 1742.5, df = 708103, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.9000516 0.9009325
#sample estimates:
#  cor
#0.900493


library(GGally)

corr_matrix2<-corr_matrix[,2:3]

plot.cor.BETA<-ggpairs(corr_matrix2, title="Correlation between main nonverbal all vs nonverbal sensitivity all ")
file.create(paste0("QCCorrelation_nonverbal_sens_all.pdf"))
path1<-file.path(paste0("QCCorrelation_nonverbal_sens_all.pdf"))
pdf(path1)
print(plot.cor.BETA)
dev.off()

#  SAME CORRELATION given by the THE COR.TEST FUNCTION THAN IN THE PDF RESULTING FROM GGALLY package

###################################
###nonverbal MAIN vs SENSITIVITY###
###################################

### CpGs that are SIG with PVAL< 10-4 (IN THE MAIN MODEL)

setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/COMPARE_MODELS/net_NEW/")

nonverb_sens<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220707_Output/SENSITIVITY_NOU/QC_ResultsSENSITIVITY/nonverbal/nonverbal_QCData.txt", header=TRUE)
head(nonverb_sens)
#probeID          BETA           SE        P_VAL   padj.fdr padj.bonf
#1 cg05742591 -1.566119e-03 3.027884e-04 2.312045e-07 0.09954241        no
#2 cg21142557  3.547450e-04 6.907483e-05 2.811515e-07 0.09954241        no
#3 cg01153132  1.137337e-04 2.275556e-05 5.790933e-07 0.10327737        no
#4 cg25096679 -1.493898e-04 3.021212e-05 7.626220e-07 0.10327737        no
#5 cg05409601 -6.534902e-05 1.322470e-05 7.754957e-07 0.10327737        no
#6 cg16495530  2.098862e-04 4.275209e-05 9.136493e-07 0.10327737        no
#CpG_chrm   CpG_beg   CpG_end
#1     chr5   1656679   1656681
#2     chr1 204401133 204401135
#3     chr6 119078948 119078950
#4    chr10 101841248 101841250
#5     chr4  77156758  77156760
#6    chr11 122726318 122726320



nonverb_all<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/nonverbal/nonverbal_QCData.txt", header=TRUE)

head(nonverb_all)
#probeID          BETA           SE        P_VAL  padj.fdr padj.bonf
#1 cg18892718 -1.934287e-03 3.761491e-04 2.713358e-07 0.1056628        no
#2 cg11995801 -4.945211e-05 9.650265e-06 2.984383e-07 0.1056628        no
#3 cg15110918 -2.337091e-03 4.731039e-04 7.815800e-07 0.1844802        no
#4 cg22034155 -9.322172e-05 1.998282e-05 3.084759e-06 0.2593288        no
#5 cg15856043  6.664539e-04 1.431903e-04 3.250469e-06 0.2593288        no
#6 cg00072720 -1.384769e-04 2.978712e-05 3.337357e-06 0.2593288        no
#CpG_chrm   CpG_beg   CpG_end
#1     chr7  87629447  87629449
#2     chr2 231961564 231961566
#3    chr11 126655602 126655604
#4     chr6  31835421  31835423
#5    chr12 118542075 118542077
#6    chr17   7262511   7262513



nonverb_all4<-nonverb_all[nonverb_all$P_VAL<0.0001,]

probes<-nonverb_all4$probeID #probes from the nonverbal model that are sig by pval*10^-4

nonverb_sens<-nonverb_sens[nonverb_sens$probeID %in% probes,]

dim(nonverb_sens)
#103

sensord<-nonverb_sens[order(nonverb_sens$probeID),]
nonverballord<-nonverb_all4[order(nonverb_all4$probeID),]

sensordb<- sensord[,c("probeID","BETA")]
nonverballordb<- nonverballord[,c("probeID","BETA")]

colnames(sensordb)<-c("probeID","BETA_nonverb_sensitivity")
colnames(nonverballordb)<-c("probeID","BETA_main_nonverb_all")

corr_matrix<-merge(sensordb, nonverballordb,by="probeID") 

cor.test(corr_matrix$BETA_nonverb_sensitivity, corr_matrix$BETA_main_nonverb_all)

#Pearson's product-moment correlation

#data:  corr_matrix$BETA_nonverb_sensitivity and corr_matrix$BETA_main_nonverb_all
#t = 95.558, df = 110, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9913225 0.9958954
#sample estimates:
#      cor
#0.9940307




library(GGally)

corr_matrix2<-corr_matrix[,2:3]

plot.cor.BETA<-ggpairs(corr_matrix2, title="Correlation between main nonverbal sig-4 vs nonverbal sensitivity sig-4 ")
file.create(paste0("QCCorrelation_nonverbal_sens_sig4.pdf"))
path1<-file.path(paste0("QCCorrelation_nonverbal_sens_sig4.pdf"))
pdf(path1)
print(plot.cor.BETA)
dev.off()

#  SAME CORRELATION given by the THE COR.TEST FUNCTION THAN IN THE PDF RESULTING FROM GGALLY package


######################################################
###nonverbal MAIN vs Mat_iq and model without sgest###
######################################################

### ALL

setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/COMPARE_MODELS/net_NEW/")

nonverb_mat<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220601NOUIQ_Output/Resultats_matIQ/QC_ResultsmaternalIQ/nonverbal/nonverbal_QCData.txt", header=TRUE)
head(nonverb_mat)
####ESTAVA POSAT EL VELL: /home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220627_Output/RESULTATS_sges_NOU/QC_ResultsSGES/nonverbal/nonverbal_QCData.txt'
### LI HE CANVIAT EL NOM PER NO. 
nonverb_sges<- read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_20220722NOU_Output/resultats_nous_sges/QC_ResultsSGES/nonverbal/nonverbal_QCData.txt", header=TRUE)
head(nonverb_sges)

nonverb_all<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/nonverbal/nonverbal_QCData.txt", header=TRUE)

head(nonverb_all)

probes<-nonverb_all$probeID #probes from the nonverbal model that are sig by pval*10^-4

nonverbmatord<-nonverb_mat[order(nonverb_mat$probeID),]
nonverbsgesord<-nonverb_sges[order(nonverb_sges$probeID),]
nonverballord<-nonverb_all[order(nonverb_all$probeID),]

nonverbmatordb<- nonverbmatord[,c("probeID","BETA")]
nonverbsgesordb<- nonverbsgesord[,c("probeID","BETA")]
nonverballordb<- nonverballord[,c("probeID","BETA")]

colnames(nonverbmatordb)<-c("probeID","BETA_nonverb_maternal_iq")
colnames(nonverbsgesordb)<-c("probeID","BETA_nonverb_wo_sges")
colnames(nonverballordb)<-c("probeID","BETA_main_nonverb_all")

corr_matrix<-merge(nonverbmatordb, nonverballordb,by="probeID") 
corr_matrix2<-merge(nonverbsgesordb, corr_matrix,by="probeID") 
head(corr_matrix2)
#probeID BETA_nonverb_wo_sges BETA_nonverb_maternal_iq
#1 cg00000029        -2.729840e-06             2.869676e-05
#2 cg00000103         1.754401e-04             1.589427e-04
#3 cg00000109        -1.343543e-05            -2.301129e-05
#4 cg00000155        -7.887482e-06            -3.039123e-06
#5 cg00000158        -3.196729e-05            -3.228232e-05
#6 cg00000165        -5.798512e-04            -3.857034e-04
#BETA_main_nonverb_all
#1         -3.942379e-08
#2          8.419662e-05
#3         -1.185585e-05
#4         -6.959237e-06
#5         -3.234250e-05
#6         -5.890364e-04



cor.test(corr_matrix2$BETA_nonverb_maternal_iq, corr_matrix2$BETA_main_nonverb_all)

#Pearson's product-moment correlation

#data:  corr_matrix2$BETA_nonverb_maternal_iq and corr_matrix2$BETA_main_nonverb_all
#t = 3341.4, df = 708103, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.9695821 0.9698600
#sample estimates:
#  cor
#0.9697214


cor.test(corr_matrix2$BETA_nonverb_wo_sges, corr_matrix2$BETA_main_nonverb_all)

#Pearson's product-moment correlation

#data:  corr_matrix2$BETA_nonverb_wo_sges and corr_matrix2$BETA_main_nonverb_all
#t = 6808.2, df = 708103, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.9924129 0.9924830
#sample estimates:
#  cor
#0.992448

library(GGally)

corr_matrix2<-corr_matrix2[,2:4]

plot.cor.BETA<-ggpairs(corr_matrix2, title="Correlation between main nonverbal all vs nonverbal maternal iq all and wo gestational age model ")
file.create(paste0("QCCorrelation_nonverbal_matiq_sges_all.pdf"))
path1<-file.path(paste0("QCCorrelation_nonverbal_matiq_sges_all.pdf"))
pdf(path1)
print(plot.cor.BETA)
dev.off()

#  SAME CORRELATION given by the THE COR.TEST FUNCTION THAN IN THE PDF RESULTING FROM GGALLY package

######################################################
###nonverbal MAIN vs Mat_iq and model without sgest###
######################################################

### CpGs that are SIG with PVAL< 10-4 (IN THE MAIN MODEL)

setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/COMPARE_MODELS/net_NEW/")

nonverb_all<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/nonverbal/nonverbal_QCData.txt", header=TRUE)

head(nonverb_all)

nonverb_all4<-nonverb_all[nonverb_all$P_VAL<0.0001,]

probes<-nonverb_all4$probeID #probes from the nonverbal model that are sig by pval*10^-4

nonverb_mat4<-nonverb_mat[nonverb_mat$probeID %in% probes,] 
nonverb_sges4<-nonverb_sges[nonverb_sges$probeID %in% probes,]

dim(nonverb_sges4)
#103
dim(nonverb_mat4)
#103
dim(nonverb_all4)
#103

nonverbmatord<-nonverb_mat4[order(nonverb_mat4$probeID),]
nonverbsgesord<-nonverb_sges4[order(nonverb_sges4$probeID),]
nonverballord<-nonverb_all4[order(nonverb_all4$probeID),]

nonverbmatordb<- nonverbmatord[,c("probeID","BETA")]
nonverbsgesordb<- nonverbsgesord[,c("probeID","BETA")]
nonverballordb<- nonverballord[,c("probeID","BETA")]

colnames(nonverbmatordb)<-c("probeID","BETA_nonverb_maternal_iq_4")
colnames(nonverbsgesordb)<-c("probeID","BETA_nonverb_wo_sges_4")
colnames(nonverballordb)<-c("probeID","BETA_main_nonverb_4")

corr_matrix<-merge(nonverbmatordb, nonverballordb,by="probeID") 
corr_matrix2<-merge(nonverbsgesordb, corr_matrix,by="probeID") 
head(corr_matrix2)
#probeID BETA_nonverb_wo_sges_4 BETA_nonverb_maternal_iq_4
#1 cg00072720          -0.0001374353              -0.0001453470
#2 cg00126261          -0.0014856187              -0.0016238756
#3 cg00184876          -0.0003076718              -0.0003049678
#4 cg00247269          -0.0004976955              -0.0005252350
#5 cg00734567          -0.0005592525              -0.0004959636
#6 cg01148113          -0.0015115625              -0.0016376028
#BETA_main_nonverb_4
#1       -0.0001384769
#2       -0.0015105157
#3       -0.0003050343
#4       -0.0005091527
#5       -0.0005824753
#6       -0.0016325517



cor.test(corr_matrix2$BETA_nonverb_maternal_iq_4, corr_matrix2$BETA_main_nonverb_4)

#Pearson's product-moment correlation

#data:  corr_matrix2$BETA_nonverb_maternal_iq_4 and corr_matrix2$BETA_main_nonverb_4
#t = 127.34, df = 110, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.9950916 0.9976806
#sample estimates:
#  cor
#0.9966255


cor.test(corr_matrix2$BETA_nonverb_wo_sges_4, corr_matrix2$BETA_main_nonverb_4)

#Pearson's product-moment correlation

#data:  corr_matrix2$BETA_nonverb_wo_sges_4 and corr_matrix2$BETA_main_nonverb_4
#t = 333.62, df = 110, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
 # 0.9992813 0.9996608
#sample estimates:
 # cor
#0.9995062

library(GGally)

corr_matrix2<-corr_matrix2[,2:4]

plot.cor.BETA<-ggpairs(corr_matrix2, title="Correlation between main nonverbal sig -4 vs nonverbal maternal iq and wo gestational age model ")
file.create(paste0("QCCorrelation_nonverbal_matiq_sges_sig4.pdf"))
path1<-file.path(paste0("QCCorrelation_nonverbal_matiq_sges_sig4.pdf"))
pdf(path1)
print(plot.cor.BETA)
dev.off()



