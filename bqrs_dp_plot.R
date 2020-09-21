library(ggplot2)
setwd("C:/Users/selee/Desktop/class/defaultCGP")
df <- read.table("random20000.genotypegvcfs.txt",header=T)
colnames(df)<-c('TRANSITION','REF','ALT','BaseQRankSum','DP')
dfSnv <- df[complete.cases(df),]
dfSnv <- dfSnv[dfSnv$DP<3000,]
nrow(dfSnv)
dfti <- dfSnv[dfSnv$TRANSITION==1,]
dftv <- dfSnv[dfSnv$TRANSITION==0,]
nrow(dfti)
nrow(dftv)
#ti
GtoA <- which(dfSnv$REF=='G'&dfSnv$ALT=='A')
AtoG <- which(dfSnv$REF=='A'&dfSnv$ALT=='G')
CtoT <- which(dfSnv$REF=='C'&dfSnv$ALT=='T')
TtoC <- which(dfSnv$REF=='T'&dfSnv$ALT=='C')

#tv
GtoC <- which(dfSnv$REF=='G'&dfSnv$ALT=='C')
GtoT <- which(dfSnv$REF=='G'&dfSnv$ALT=='T')
AtoT <- which(dfSnv$REF=='A'&dfSnv$ALT=='T')
AtoC <- which(dfSnv$REF=='A'&dfSnv$ALT=='C')

CtoA <- which(dfSnv$REF=='C'&dfSnv$ALT=='A')
CtoG <- which(dfSnv$REF=='C'&dfSnv$ALT=='G')
TtoA <- which(dfSnv$REF=='T'&dfSnv$ALT=='A')
TtoG <- which(dfSnv$REF=='T'&dfSnv$ALT=='G')
dfSnv$mutype<-dfSnv$TRANSITION
dfSnv$mutype[GtoA]<-'G>A'
dfSnv$mutype[AtoG]<-'A>G'
dfSnv$mutype[CtoT]<-'C>T'
dfSnv$mutype[TtoC]<-'T>C'
dfti <- dfSnv[dfSnv$TRANSITION==1,]
dftv <- dfSnv[dfSnv$TRANSITION==0,]
ggplot(dfti,aes(x=BaseQRankSum,y=DP,color=mutype))+geom_point(size=1)

#tv
dfSnv$mutype[GtoC]<-'G>C'
dfSnv$mutype[GtoT]<-'G>T'
dfSnv$mutype[AtoT]<-'A>T'
dfSnv$mutype[AtoC]<-'A>C'
dfSnv$mutype[CtoA]<-'C>A'
dfSnv$mutype[CtoG]<-'C>G'
dfSnv$mutype[TtoA]<-'T>A'
dfSnv$mutype[TtoG]<-'T>G'
dftv <- dfSnv[dfSnv$TRANSITION=='0',]
ggplot(dftv,aes(x=BaseQRankSum,y=DP,color=mutype))+geom_point(size=1)
