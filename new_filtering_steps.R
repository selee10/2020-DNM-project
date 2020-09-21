library(ggplot2)
setwd("C:/Users/selee/Desktop/class/defaultCGP")
df <- read.table("3_BQRS.all.dbsnp151.annotated.table.txt",header=T)
df <- df[complete.cases(df),]
dfSnv <- df[df$TYPE=='SNP',]
dfti <- dfSnv[dfSnv$TRANSITION=='1',]
dftv <- dfSnv[dfSnv$TRANSITION=='0',]
######separate ti######
attach(dfti)
lm1 <- lm(Chinese.AF~QD) #AF = 0.02736*QD+0.09290
xvec1 <- dfti$QD
yvec1 <- c()
for (val in xvec1){
  yval1 = 0.026*val+0.09290#0.02736 originally but the graph should have lower slope
  yvec1 <- c(yvec1,yval1)
}
dfti$lm.AF <- yvec1
ggplot(dfti,aes(x=QD))+
  geom_point(aes(y=Chinese.AF))+
  geom_line(aes(y=lm.AF))
#fp
d1 <- c()
for (i in seq(1,nrow(dfti))){
  a1 <- dfti$Chinese.AF[i]
  b1 <- dfti$lm.AF[i]
  d1 <- c(d1,a1<b1)
}
tifp <- dfti[d1,]
#tp
e1 <- c()
for (i in seq(1,nrow(dfti))){
  a1 <- dfti$Chinese.AF[i]
  b1 <- dfti$lm.AF[i]
  e1 <- c(e1,a1>=b1)
}
tiwithtp <- dfti[e1,]
nrow(tiwithtp)
######separate tv######
attach(dftv)
lm2 <- lm(Chinese.AF~QD) #AF = 0.02749 * QD + 0.09859
xvec <- dftv$QD
yvec <- c()
for (val in xvec){
  yval = 0.026*val+0.09859#0.02749 originally but the graph should have lower slope
  yvec <- c(yvec,yval)
}
dftv$lm.AF <- yvec
ggplot(dftv,aes(x=QD))+
  geom_point(aes(y=Chinese.AF))+
  geom_line(aes(y=lm.AF))
#fp
d <- c()
for (i in seq(1,nrow(dftv))){
  a <- dftv$Chinese.AF[i]
  b <- dftv$lm.AF[i]
  d <- c(d,a<b)
}
tvfp <- dftv[d,]
nrow(tvfp)#301
#tp
e <- c()
for (i in seq(1,nrow(dftv))){
  a <- dftv$Chinese.AF[i]
  b <- dftv$lm.AF[i]
  e <- c(e,a>=b)
}
tvwithtp <- dftv[e,]
nrow(tvwithtp)
######0.AF>=0.3(to filter out somatic variants)######
dfSnv <- dfSnv[dfSnv$Chinese.AF>=0.3,]
dfti <- dfSnv[dfSnv$TRANSITION=='1',]
dftv <- dfSnv[dfSnv$TRANSITION=='0',]
nrow(dfti)
nrow(dftv)#318/335
######1.af threshold to adjust ti/tv>=1######
#background information: 2011 conrad et al
dfSnv <- dfSnv[dfSnv$Chinese.AF>=0.4,]
dfti <- dfSnv[dfSnv$TRANSITION=='1',]
dftv <- dfSnv[dfSnv$TRANSITION=='0',]
nrow(dfti)
nrow(dftv)#113/95
######2. exclude variants in dbSNP######
dfti <- dfti[dfti$ID=='.',]
dftv <- dftv[dftv$ID=='.',]
nrow(dfti)#56 / dbsnp: 57
nrow(dftv)#73 / dbsnp: 22
######3.separate two lines in QD-AF######
#ti fp
dfti$lm.AF <- dfti$QD*0.026+0.09290
d <- c()
for (i in seq(1,nrow(n_dfti))){
  a <- dfti$Chinese.AF[i]
  b <- dfti$lm.AF[i]
  d <- c(d,a<b)
}
tifp <- dfti[d,]
nrow(tifp)#2 / dbsnp: 22
#ti tp
e <- c()
for (i in seq(1,nrow(dfti))){
  a <- dfti$Chinese.AF[i]
  b <- dfti$lm.AF[i]
  e <- c(e,a>=b)
}
tiwithtp <- dfti[e,]
nrow(tiwithtp)#54 / dbsnp: 35
#tv fp
dftv$lm.AF <- dftv$QD*0.026+0.09859
d <- c()
for (i in seq(1,nrow(dftv))){
  a <- dftv$Chinese.AF[i]
  b <- dftv$lm.AF[i]
  d <- c(d,a<b)
}
tvfp <- dftv[d,]
nrow(tvfp)#7 / dbsnp: 12
#tv tp
e <- c()
for (i in seq(1,nrow(dftv))){
  a <- dftv$Chinese.AF[i]
  b <- dftv$lm.AF[i]
  e <- c(e,a>=b)
}
tvwithtp <- dftv[e,]
nrow(tvwithtp)#66 / dbsnp: 10
#fp pos to file->observe in IGV
fp <- rbind(tifp,tvfp)
cat(fp$POS,file='C:/Users/selee/Desktop/fp.pos.txt')
######4.separate variants by BaseQRankSum(+/-)######
#ti
ti_posBQRS <- tiwithtp[tiwithtp$BaseQRankSum>=0,]
ti_negBQRS <- tiwithtp[tiwithtp$BaseQRankSum<0,]
nrow(ti_posBQRS)#23 / dbsnp: 28
nrow(ti_negBQRS)#31 / dbsnp: 7
#for (i in 1:nrow(tiwithtp)){
#  if (tiwithtp$BaseQRankSum[i]>=0){
#    tiwithtp$bqrs[i] <- 'pos'
#  }else{
#    tiwithtp$bqrs[i] <- 'neg'
#  }
#}
i<-1
while (i<=nrow(tiwithtp)){
  if (tiwithtp$BaseQRankSum[i]>=-5&tiwithtp$BaseQRankSum[i]<0){
    tiwithtp$bqrs[i] <- '-5<=BQRS<0'
  }else{
    tiwithtp$bqrs[i]<-'BQRS>=0'
  }
  i <- i + 1
}
tiwithtp$bqrs[which(tiwithtp$BaseQRankSum<(-5))]<-'BQRS<-5'
ggplot(tiwithtp,aes(x=BaseQRankSum,y=Chinese.AF,color=bqrs))+geom_point()
ggplot(tiwithtp,aes(x=QD,y=Chinese.AF,color=bqrs))+geom_point()
#tv
tv_posBQRS <- tvwithtp[tvwithtp$BaseQRankSum>=0,]
tv_negBQRS <- tvwithtp[tvwithtp$BaseQRankSum<0,]
nrow(tv_posBQRS)#32 -> 31 / dbsnp: 7
nrow(tv_negBQRS)#34 -> 35 / dbsnp: 3
i<-1
while (i<=nrow(tvwithtp)){
  if (tvwithtp$BaseQRankSum[i]>=-4&tvwithtp$BaseQRankSum[i]<0){
    tvwithtp$bqrs[i] <- 'neg1'
    }else{
      tvwithtp$bqrs[i]<-'pos'
    }
  i <- i + 1
}
tvwithtp$bqrs[which(tvwithtp$BaseQRankSum<(-4))]<-'neg2'
ggplot(tvwithtp,aes(x=BaseQRankSum,y=Chinese.AF,color=bqrs))+geom_point()
ggplot(tvwithtp,aes(x=QD,y=Chinese.AF,color=bqrs))+geom_point()
######4-1.dbSNP BaseQRankSum######
#change dfti<-dfti[dfti$ID!='.',]

#(test)filtering
tmp_tiwithtp <- tiwithtp[tiwithtp$QD>=10,]
tmp_tvwithtp <- tvwithtp[tvwithtp$QD>=10.2,]
nrow(tmp_tiwithtp)
nrow(tmp_tvwithtp)
ggplot(dfSnv1,aes(x=BaseQRankSum,y=Chinese.AF,color=TRANSITION))+geom_point()
dfSnv1 <- df[df$TYPE=='SNP',]
######annotated variants by mut spectra######
dfSnv <- df[df$TYPE=='SNP',]
dfSnv <- dfSnv[dfSnv$ID!='.',]
dfSnv <- dfSnv[dfSnv$Chinese.AF>=0.3,]
dfti <- dfSnv[dfSnv$TRANSITION=='1',]
dftv <- dfSnv[dfSnv$TRANSITION=='0',]
nrow(dfti)#136
nrow(dftv)#49
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
ggplot(dfSnv,aes(x=QD,y=Chinese.AF,color=mutype))+geom_point(size=1)
ggplot(dfSnv,aes(x=BaseQRankSum,y=Chinese.AF,color=mutype))+geom_point(size=1)
dfti <- dfSnv[dfSnv$TRANSITION==1,]
ggplot(dfti,aes(x=BaseQRankSum,y=Chinese.AF,color=mutype))+geom_point(size=1)
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
ggplot(dftv,aes(x=BaseQRankSum,y=Chinese.AF,color=mutype))+geom_point(size=1)
#T,A qual > C,G ?
ggplot(dfti,aes(x=BaseQRankSum,y=Chinese.AF,color=ALT))+geom_point()
