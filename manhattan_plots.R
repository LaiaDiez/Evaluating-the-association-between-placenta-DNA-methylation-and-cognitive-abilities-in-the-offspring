### Manhattan plots ###


setwd("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/")

# GENERAL 

general<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/general/general_QCData.txt", header=TRUE)

general$chr2<-gsub("[a-zA-Z ]", "", general$CpG_chrm)
general$chr<-as.numeric(general$chr2)

png("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/manhattan_GENERAL.png",height=600,width=800)

obspval <- (general$P_VAL)
chr <- (general$chr)
pos <- (general$CpG_beg)
obsmax <- trunc(max(-log10(obspval)))+1

sort.ind <- order(chr, pos) 
chr <- chr[sort.ind]
pos <- pos[sort.ind]
obspval <- obspval[sort.ind]

x <- 1:22
x2<- 1:22

for (i in 1:22)
{
  curchr=which(chr==i)
  x[i] <- trunc((max(pos[curchr]))/100) +100000
  x2[i] <- trunc((min(pos[curchr]))/100) -100000
}

x[1]=x[1]-x2[1]
x2[1]=0-x2[1]

for (i in 2:24)
{
  x[i] <- x[i-1]-x2[i]+x[i]
  x2[i] <- x[i-1]-x2[i]
  
}
locX = trunc(pos/100) + x2[chr]
locY = -log10(obspval)
col1=rgb(0,0,108,maxColorValue=255)
col2=rgb(100,149,237,maxColorValue=255)
col3=rgb(0,205,102,maxColorValue=255)
col4 <- ifelse (chr%%2==0, col1, col2)
curcol <- ifelse (obspval<5e-8, col3, col4) 
plot(locX,locY,pch=20,col=curcol,axes=F,ylab="-log10 p-value",xlab="",bty="n",ylim=c(0,8),cex=0.8, abline(h=c(7.0611,5), col=c("red","green"), lty = c(2,2)), yaxt = "none") 
# we write the abline in the bonferroni limit which is 0.05 divided by the total of cpg's = 0.05/708105 = 7.0611e-8, and because the y axis is -logpval then we write only 7.0611
axis(2, seq(0,8,1))

for (i in 1:22)
{
  labpos = (x[i] + x2[i]) / 2
  mtext(i,1,at=labpos,cex=0.8,line=0)
}
mtext("Chromosome",1,at=x[22]/2,cex=1,line=1)
dev.off()


###########

# VERBAL 
verbal<- read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/verbal/verbal_QCData.txt",header=TRUE)


verbal$chr2<-gsub("[a-zA-Z ]", "", verbal$CpG_chrm)
verbal$chr<-as.numeric(verbal$chr2)

png("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/manhattan_VERBAL.png",height=600,width=800)

obspval <- (verbal$P_VAL)
chr <- (verbal$chr)
pos <- (verbal$CpG_beg)
obsmax <- trunc(max(-log10(obspval)))+1

sort.ind <- order(chr, pos) 
chr <- chr[sort.ind]
pos <- pos[sort.ind]
obspval <- obspval[sort.ind]

x <- 1:22
x2<- 1:22

for (i in 1:22)
{
  curchr=which(chr==i)
  x[i] <- trunc((max(pos[curchr]))/100) +100000
  x2[i] <- trunc((min(pos[curchr]))/100) -100000
}

x[1]=x[1]-x2[1]
x2[1]=0-x2[1]

for (i in 2:24)
{
  x[i] <- x[i-1]-x2[i]+x[i]
  x2[i] <- x[i-1]-x2[i]
  
}
locX = trunc(pos/100) + x2[chr]
locY = -log10(obspval)
col1=rgb(0,0,108,maxColorValue=255)
col2=rgb(100,149,237,maxColorValue=255)
col3=rgb(0,205,102,maxColorValue=255)
col4 <- ifelse (chr%%2==0, col1, col2)
curcol <- ifelse (obspval<5e-8, col3, col4) 
plot(locX,locY,pch=20,col=curcol,axes=F,ylab="-log10 p-value",xlab="",bty="n",ylim=c(0,8),cex=0.8, abline(h=c(7.0611,5), col=c("red","green"), lty = c(2,2)), yaxt = "none") # we write the abline in the bonferroni limit which is 0.05 divided by the total of cpg's = 0.05/708105 = 7.0611e-8, and because the y axis is -logpval then we write only 7.0611
axis(2, seq(0,8,1))

for (i in 1:22)
{
  labpos = (x[i] + x2[i]) / 2
  mtext(i,1,at=labpos,cex=0.8,line=0)
}
mtext("Chromosome",1,at=x[22]/2,cex=1,line=1)
dev.off()


############
# NONVERBAL 
nonverbal<-read.table("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/QC_MAIN_NOU/nonverbal/nonverbal_QCData.txt", header=TRUE)

nonverbal$chr2<-gsub("[a-zA-Z ]", "", verbal$CpG_chrm)
nonverbal$chr<-as.numeric(verbal$chr2)

png("/home/isglobal.lan/ldiez/data/WS_INMA/Methylation_INMA/PACE/Pla_IQ_LD/results/INMA_202200601_Output/NOUsignificativeresultsMAIN/manhattan_NONVERBAL.png",height=600,width=800)

obspval <- (nonverbal$P_VAL)
chr <- (nonverbal$chr)
pos <- (nonverbal$CpG_beg)
obsmax <- trunc(max(-log10(obspval)))+1

sort.ind <- order(chr, pos) 
chr <- chr[sort.ind]
pos <- pos[sort.ind]
obspval <- obspval[sort.ind]

x <- 1:22
x2<- 1:22

for (i in 1:22)
{
  curchr=which(chr==i)
  x[i] <- trunc((max(pos[curchr]))/100) +100000
  x2[i] <- trunc((min(pos[curchr]))/100) -100000
}

x[1]=x[1]-x2[1]
x2[1]=0-x2[1]

for (i in 2:24)
{
  x[i] <- x[i-1]-x2[i]+x[i]
  x2[i] <- x[i-1]-x2[i]
  
}
locX = trunc(pos/100) + x2[chr]
locY = -log10(obspval)
col1=rgb(0,0,108,maxColorValue=255)
col2=rgb(100,149,237,maxColorValue=255)
col3=rgb(0,205,102,maxColorValue=255)
col4 <- ifelse (chr%%2==0, col1, col2)
curcol <- ifelse (obspval<5e-8, col3, col4) 
plot(locX,locY,pch=20,col=curcol,axes=F,ylab="-log10 p-value",xlab="",bty="n",ylim=c(0,8),cex=0.8, abline(h=c(7.0611,5), col=c("red","green"), lty = c(2,2)), yaxt = "none") # we write the abline in the bonferroni limit which is 0.05 divided by the total of cpg's = 0.05/708105 = 7.0611e-8, and because the y axis is -logpval then we write only 7.0611
axis(2, seq(0,8,1))

for (i in 1:22)
{
  labpos = (x[i] + x2[i]) / 2
  mtext(i,1,at=labpos,cex=0.8,line=0)
}
mtext("Chromosome",1,at=x[22]/2,cex=1,line=1)
dev.off()

