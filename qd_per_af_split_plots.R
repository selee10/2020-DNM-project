#install.packages("RColorBrewer")
#library(RColorBrewer)
library(ggplot2)
setwd("C:/Users/selee/Desktop/class/defaultCGP")
#df <- read.table("qual.table.Reannotated.hg19.common.SnpSift.masked.isHet.tranch99.9.vcf.txt",header=T)
df <- read.table("BQRS.all.dbsnp151.annotated.table.txt",header=T)
nrow(df)
df <- df[complete.cases(df),]
dfSnv <- df[df$TYPE=='SNP',]
#remove variants with AF<0.3
dfSnv <- dfSnv[dfSnv$Chinese.AF>=0.3,]

dfti <- dfSnv[dfSnv$TRANSITION=='1',]
dftv <- dfSnv[dfSnv$TRANSITION=='0',]
nrow(dfti)
nrow(dftv)
ggplot(dfti,aes(x=QD,y=Chinese.AF))+geom_point()
ggplot(dftv,aes(x=QD,y=Chinese.AF))+geom_point()

#####dbsnp build 151(all)#######
dfti$tmpID[which(dfti$ID!='.')]<-0
dfti$tmpID[which(dfti$ID=='.')]<-1
ggplot(dfti,aes(x=QD,y=Chinese.AF,color=tmpID))+geom_point(size=0.5)
ggplot(dfti,aes(x=BaseQRankSum,y=Chinese.AF,color=tmpID))+geom_point()

dftv$tmpID[which(dftv$ID!='.')]<-0
dftv$tmpID[which(dftv$ID=='.')]<-1
ggplot(dftv,aes(x=QD,y=Chinese.AF,color=tmpID))+geom_point(size=0.6)
ggplot(dftv,aes(x=BaseQRankSum,y=Chinese.AF,color=tmpID))+geom_point()

#exclude variants in dbSNP
dfti <- dfti[dfti$ID=='.',]
dftv <- dftv[dftv$ID=='.',]
nrow(dfti)
#(whole)774 to 281,(AF>=0.3)318 to 182
nrow(dftv)
#(whole)711 to 404,(AF>=0.3)335 to 286

#####################separate lines(ti)###########################
attach(dfti)
lm1 <- lm(Chinese.AF~QD) #AF = 0.02736*QD+0.09290
xvec1 <- dfti$QD
yvec1 <- c()

for (val in xvec1){
  yval1 = 0.026*val+0.09290#0.02736 originally but the graph should have lower slope
  yvec1 <- c(yvec1,yval1)
}
dfti$lm.AF <- yvec1
head(dfti)
ggplot(dfti,aes(x=QD))+
  geom_point(aes(y=Chinese.AF))+
  geom_line(aes(y=lm.AF))
d1 <- c()
for (i in seq(1,nrow(dfti))){
  a1 <- dfti$Chinese.AF[i]
  b1 <- dfti$lm.AF[i]
  d1 <- c(d1,a1<b1)
}
tifp <- dfti[d1,]
nrow(tifp)

#cat(tifp$POS,file='C:/Users/selee/Desktop/tifp_pos.txt',sep='\n')

ggplot(dfti,aes(x=BaseQRankSum,y=Chinese.AF,color=ID))+geom_point()

e1 <- c()
for (i in seq(1,nrow(dfti))){
  a1 <- dfti$Chinese.AF[i]
  b1 <- dfti$lm.AF[i]
  e1 <- c(e1,a1>=b1)
}
tiwithtp <- dfti[e1,]
nrow(tiwithtp) #409

#cat(tiwithtp$POS,file='C:/Users/selee/Desktop/tiwithtp_pos.txt',sep='\n')

ggplot(tiwithtp,aes(x=QD,y=Chinese.AF))+geom_point()
ggplot(tiwithtp,aes(x=SOR,y=Chinese.AF))+geom_point(size=.5)
ggplot(tiwithtp,aes(x=MQ,y=Chinese.AF))+geom_point()
ggplot(tiwithtp,aes(x=ReadPosRankSum,y=Chinese.AF))+geom_point()

nrow(dfti[dfti$tmpID==1,])#281
nrow(dftv[dftv$tmpID==1,])#404
#####################separate lines(tv)###########################
attach(dftv)
lm2 <- lm(Chinese.AF~QD) #AF = 0.02749 * QD + 0.09859

xvec <- dftv$QD
yvec <- c()
for (val in xvec){
  yval = 0.026*val+0.09859#0.02749 originally but the graph should have lower slope
  yvec <- c(yvec,yval)
}
dftv$lm.AF <- yvec
head(dftv)
ggplot(dftv,aes(x=QD))+
  geom_point(aes(y=Chinese.AF))+
  geom_line(aes(y=lm.AF))
d <- c()
for (i in seq(1,nrow(dftv))){
  a <- dftv$Chinese.AF[i]
  b <- dftv$lm.AF[i]
  d <- c(d,a<b)
}
tvfp <- dftv[d,]
nrow(tvfp)#301

ggplot(tvfp,aes(x=ReadPosRankSum,y=Chinese.AF))+geom_point()
#cat(tvfp$POS,file='C:/Users/selee/Desktop/tvfp_pos.txt',sep='\n')

ggplot(tvfp,aes(x=QD,y=Chinese.AF))+geom_point()

e <- c()
for (i in seq(1,nrow(dftv))){
  a <- dftv$Chinese.AF[i]
  b <- dftv$lm.AF[i]
  e <- c(e,a>=b)
}
tvwithtp <- dftv[e,]
nrow(tvwithtp) #410
#cat(tvwithtp$POS,file='C:/Users/selee/Desktop/tvwithtp_pos.txt',sep='\n')
ggplot(tvwithtp,aes(x=QD,y=Chinese.AF))+geom_point()
ggplot(tvwithtp,aes(x=SOR,y=Chinese.AF))+geom_point()
ggplot(tvwithtp,aes(x=FS,y=Chinese.AF))+geom_point()
ggplot(tvwithtp,aes(x=ReadPosRankSum,y=Chinese.AF))+geom_point()
ggplot(tvwithtp,aes(x=MQ,y=Chinese.AF))+geom_point()
ggplot(tvwithtp,aes(x=MQRankSum,y=Chinese.AF))+geom_point()

##########(tv)set apart possible somatic variants##########
tvwithtp$Chinese.AD<-round(tvwithtp$Chinese.AF*tvwithtp$Chinese.DP)
af_tvwithtp <- tvwithtp[tvwithtp$Chinese.AF>=0.3&
                          tvwithtp$Chinese.AF<=0.7,]
attach(af_tvwithtp)
lm(Chinese.AF~QD)#AF=0.02737*QD+0.14036
af_tvwithtp$lm.somatic.AF <- 0.0255*af_tvwithtp$QD+0.17
ggplot(af_tvwithtp,aes(x=QD))+
  geom_point(aes(y=Chinese.AF))+
  geom_line(aes(y=lm.somatic.AF,color='red'))
tmp <- which(af_tvwithtp$Chinese.AF>=
               af_tvwithtp$lm.somatic.AF)
length(tmp)#65
tvpossibletp <- af_tvwithtp
tvpossibletp$lm.somatic.AF[tmp]<-7
ggplot(tvpossibletp,aes(x=Chinese.DP))+
  geom_point(aes(y=Chinese.AF,color=lm.somatic.AF))
ggplot(tvpossibletp,aes(x=Chinese.GQ))+
  geom_point(aes(y=Chinese.AF,color=lm.somatic.AF))
ggplot(tvpossibletp,aes(x=Chinese.AD))+
  geom_point(aes(y=Chinese.AF,color=lm.somatic.AF))
candidate.tv <- tvpossibletp[tvpossibletp$lm.somatic.AF==7,]
hist(candidate.tv$Chinese.AF)

##########(ti)set apart possible somatic variants######
af_tiwithtp <- tiwithtp[tiwithtp$Chinese.AF>=0.47&
                          tiwithtp$Chinese.AF<=0.53,]
nrow(af_tiwithtp)
#0.4 to 0.6: 84
#0.41 to 0.59: 68
#0.42 to 0.58: 58
#0.43 to 0.57: 48
#0.44 to 0.56: 36
#0.45 to 0.55: 27
#0.46 to 0.54: 20
#0.47 to 0.53: 14
attach(af_tiwithtp)
lm(Chinese.AF~QD)#AF=0.0296*QD+0.1166
af_tiwithtp$lm.somatic.AF <- 0.0226*af_tiwithtp$QD+0.19
ggplot(af_tiwithtp,aes(x=QD))+
  geom_point(aes(y=Chinese.AF))+
  geom_line(aes(y=lm.somatic.AF,color='red'))
tmpti <- which(af_tiwithtp$Chinese.AF>=
               af_tiwithtp$lm.somatic.AF)
length(tmpti)#70
af_tiwithtp$Chinese.AD<-round(
  af_tiwithtp$Chinese.AF*
    af_tiwithtp$Chinese.DP)
tipossibletp <- af_tiwithtp
tipossibletp$lm.somatic.AF[tmpti]<-7
ggplot(tipossibletp,aes(x=ReadPosRankSum))+
  geom_point(aes(y=Chinese.AF,color=lm.somatic.AF))

ggplot(tipossibletp,aes(x=Chinese.DP))+
  geom_point(aes(y=Chinese.AF,color=lm.somatic.AF))
ggplot(tipossibletp,aes(x=Chinese.AD))+
  geom_point(aes(y=Chinese.AF,color=lm.somatic.AF))
candidate.ti <- tipossibletp[tipossibletp$lm.somatic.AF==7,]
hist(candidate.ti$Chinese.AF)

#ti plot(tiwithtp and tifp)
which(dfti$Chinese.AF>=dfti$lm.AF)
nadfti <- dfti
nadfti$lm.AF[which(nadfti$Chinese.AF>=nadfti$lm.AF)]<-NA
head(nadfti)
#nadfti$reciDP <- 1/nadfti$Chinese.DP
ggplot(nadfti,aes(x=QD,y=Chinese.AF,color=lm.AF))+geom_point()
#ggplot(nadfti,aes(x=reciDP,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadfti,aes(x=QUAL,y=Chinese.AF,color=lm.AF))+geom_point(size=.5)
ggplot(nadfti,aes(x=FS,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadfti,aes(x=MQ,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadfti,aes(x=MQRankSum,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadfti,aes(x=ReadPosRankSum,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadfti,aes(x=SOR,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadfti,aes(x=C_Father.GQ,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadfti,aes(x=C_Father.DP,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadfti,aes(x=C_Father.AF,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadfti,aes(x=Chinese.GQ,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadfti,aes(x=Chinese.DP,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadfti,aes(x=samplename.GQ,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadfti,aes(x=samplename.DP,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadfti,aes(x=samplename.AF,y=Chinese.AF,color=lm.AF))+geom_point()

colnames(nadfti)

#tp + somatic
nrow(tiwithtp)
ggplot(tiwithtp,aes(x=Chinese.DP,y=ReadPosRankSum,color=Chinese.AF))+geom_point()
ggplot(tiwithtp,aes(x=Chinese.DP,y=FS,color=Chinese.AF))+geom_point()
ggplot(tiwithtp,aes(x=ReadPosRankSum,y=FS,color=Chinese.AF))+geom_point()
#af0.4-0.6 to NA
siteti <- tiTPaf0.4_0.6NA[,seq(from=6,to=12)]#21
plot(siteti)
varti <- tiTPaf0.4_0.6NA[,seq(from=13,to=21)]
colnames(varti)
plot(varti)

x <- which(tiwithtp$Chinese.AF<0.4)
y <- which(tiwithtp$Chinese.AF>=0.4&tiwithtp$Chinese.AF<0.6)
z <- which(tiwithtp$Chinese.AF>0.6)
length(x)
length(y)
length(z)
tiTPaf0.4_0.6NA <- tiwithtp
tiTPaf0.4_0.6NA$Chinese.AF[x]<--99
tiTPaf0.4_0.6NA$Chinese.AF[y]<-30
tiTPaf0.4_0.6NA$Chinese.AF[z]<-99
summary(tiTPaf0.4_0.6NA$Chinese.AF)
ggplot(tiTPaf0.4_0.6NA,
       aes(x=MQ,y=MQRankSum,
           color=Chinese.AF))+geom_point()+
  scale_color_gradient2()
ggplot(tiTPaf0.4_0.6NA,
       aes(x=Chinese.DP,y=QUAL,
           color=Chinese.AF))+geom_point()+
  scale_color_gradient2()

hist(tiwithtp$Chinese.AF[y])
af0.4to0.6titp <- tiwithtp[tiwithtp$Chinese.AF>=0.4&tiwithtp$Chinese.AF<=0.6,]
nrow(af0.4to0.6titp)
hist(af0.4to0.6titp$Chinese.AF)
p <- which(af0.4to0.6titp$Chinese.AF<0.45)
q <- w
which(af0.4to0.6titp$Chinese.AF>=0.45&
             af0.4to0.6titp$Chinese.AF<0.5)
r <- which(af0.4to0.6titp$Chinese.AF>=0.5&
             af0.4to0.6titp$Chinese.AF<0.55)
s <- which(af0.4to0.6titp$Chinese.AF>=0.55)
Af0.4to0.6titp <- af0.4to0.6titp
Af0.4to0.6titp$Chinese.AF[p]<--99
Af0.4to0.6titp$Chinese.AF[q]<--20
Af0.4to0.6titp$Chinese.AF[r]<-20
Af0.4to0.6titp$Chinese.AF[s]<-99
ggplot(Af0.4to0.6titp[s,],
       aes(x=MQ,y=MQRankSum,
           color=Chinese.AF))+geom_point()+
  scale_color_gradient2()
mq60 <- Af0.4to0.6titp[Af0.4to0.6titp$MQ==60,]
mqrs0mq60 <- mq60[mq60$MQRankSum==0,]
nrow(mqrs0mq60)
nrow(Af0.4to0.6titp)
#features
mqtiWhole <- ggplot(dfti,aes(MQ))+geom_histogram()
ggplot(af0.4to0.6titp,
                 aes(x=Chinese.DP,y=Chinese.AF))+geom_point()
attach(af0.4to0.6titp)
cor(Chinese.DP, Chinese.AF,use='complete.obs',method='pearson')
ggplot(dfti,aes(C_Father.GQ))+geom_histogram()+
  geom_histogram(aes(samplename.GQ),fill='yellow',alpha=.5)
ggplot(tiwithtp,aes(C_Father.GQ))+geom_histogram()+
  geom_histogram(aes(samplename.GQ),fill='yellow',alpha=.5)

ggplot(tiwithtp,aes(x=samplename.GQ,y=Chinese.AF))+geom_point()
ggplot(tiwithtp,aes(x=C_Father.GQ,y=Chinese.AF))+geom_point()

############################tv plotting###############
nrow(dftv)
which(dftv$Chinese.AF>=dftv$lm.AF)
nadftv <- dftv
nadftv$lm.AF[which(dftv$Chinese.AF>=dftv$lm.AF)]<-7
head(nadfti)
ggplot(nadftv,aes(x=QD,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadftv,aes(x=QUAL,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadftv,aes(x=FS,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadftv,aes(x=MQ,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadftv,aes(x=MQRankSum,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadftv,aes(x=ReadPosRankSum,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadftv,aes(x=SOR,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadftv,aes(x=C_Father.GQ,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadftv,aes(x=C_Father.DP,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadftv,aes(x=C_Father.AF,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadftv,aes(x=Chinese.GQ,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadftv,aes(x=Chinese.DP,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadftv,aes(x=samplename.GQ,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadftv,aes(x=samplename.DP,y=Chinese.AF,color=lm.AF))+geom_point()
ggplot(nadftv,aes(x=samplename.AF,y=Chinese.AF,color=lm.AF))+geom_point()

af0.4to0.6tvtp <- dftv[dftv$Chinese.AF>=0.4&dftv$Chinese.AF<=0.6,]
nrow(af0.4to0.6tvtp)
hist(af0.4to0.6tvtp$Chinese.AF)
##ggplot(af0.4to0.6tvtp,aes(x=MQ,y=Chinese.AF))+geom_point()
#ggplot(af0.4to0.6tvtp,aes(x=MQRankSum,y=Chinese.AF))+geom_point()
ggplot(af0.4to0.6tvtp,aes(x=ReadPosRankSum,y=Chinese.AF))+geom_point()
afti <- tiwithtp[tiwithtp$Chinese.AF>=0.3&
                   tiwithtp$Chinese.AF<=0.7,]
aftv <- tvwithtp[tvwithtp$Chinese.AF>=0.47&
                   tvwithtp$Chinese.AF<=0.53,]
nrow(aftv)
#0.4 to 0.6: 74
#0.41 to 0.59: 63
#0.42 to 0.58: 51
#0.43 to 0.57: 42
#0.44 to 0.56: 31
#0.45 to 0.55: 23
#0.46 to 0.54: 19
#0.47 to 0.53: 13
ggplot(tvwithtp,aes(x=QD,y=Chinese.AF))+geom_point()
ggplot(tiwithtp,aes(x=QD,y=Chinese.AF))+geom_point()

which(dftv$POS==93699129)
dftv[110,]
dftv$TRANSITION[110]<--99
ggplot(dftv,aes(x=QD,y=Chinese.AF,color=TRANSITION))+geom_point()
aftitp <- tiwithtp[tiwithtp$Chinese.AF>=0.45&
                   tiwithtp$Chinese.AF<=0.55,]
nrow(aftitp)#27

aftifp <- tifp[tifp$Chinese.AF>=0.45&
                     tifp$Chinese.AF<=0.55,]
nrow(aftifp)#4

aftvtp <- tvwithtp[tvwithtp$Chinese.AF>=0.45&
                     tvwithtp$Chinese.AF<=0.55,]
nrow(aftvtp)#23

aftvfp <- tvfp[tvfp$Chinese.AF>=0.45&
                     tvfp$Chinese.AF<=0.55,]
nrow(aftvfp)#4

cat(aftitp$POS,file="C:/Users/selee/Desktop/af0.45-0.55.titp.pos.txt",sep='\n')
cat(aftifp$POS,file="C:/Users/selee/Desktop/af0.45-0.55.tifp.pos.txt",sep='\n')
cat(aftvtp$POS,file="C:/Users/selee/Desktop/af0.45-0.55.tvtp.pos.txt",sep='\n')
cat(aftvfp$POS,file="C:/Users/selee/Desktop/af0.45-0.55.tvfp.pos.txt",sep='\n')

#plotting previously (hard) filtered SNVs
matchpos <- match(qdti2$POS,dfti$POS)
dfti$matpos<-"filtered"
dfti$matpos[matchpos]<-"true"
ggplot(dfti,aes(x=QD,y=Chinese.AF,color=matpos))+geom_point()

dfti$POS<-as.numeric(dfti$POS)

ct<-0

for(i in 1:773){
  if(dfti$matpos[i]=="true"){ct=ct+1}
}