require(reshape2)
library(ggplot2)
library(rgl)
library(randomForest)
library(caret)
library(gridExtra)
library(scales)
#dataframe of z-scores
ceu.only <- ceu.only[complete.cases(ceu.only),]
yri.only <- yri.only[complete.cases(yri.only),]
compHcDfSnv3 <- compHcDfSnv3[complete.cases(compHcDfSnv3),]

z.ceu.only.new <- ceu.only#[,1:33]
z.yri.only.new <- yri.only#[,1:33]
z.hc.only.new <- compHcDfSnv3
colnames(ceu.only)
z.ceu.only.new$QUAL<-scale(ceu.only$QUAL)
z.ceu.only.new$BaseQRankSum<-scale(ceu.only$BaseQRankSum)
z.ceu.only.new$FS<-scale(ceu.only$FS)
z.ceu.only.new$MQ<-scale(ceu.only$MQ)
z.ceu.only.new$MQRankSum<-scale(ceu.only$MQRankSum)
z.ceu.only.new$QD<-scale(ceu.only$QD)
z.ceu.only.new$ReadPosRankSum<-scale(ceu.only$ReadPosRankSum)
z.ceu.only.new$VQSLOD<-scale(ceu.only$VQSLOD)
z.ceu.only.new$father.DP<-scale(ceu.only$father.DP)
z.ceu.only.new$mother.DP<-scale(ceu.only$mother.DP)
z.ceu.only.new$offspring.DP<-scale(ceu.only$offspring.DP)
z.ceu.only.new$father.PL<-scale(ceu.only$father.PL)
z.ceu.only.new$father.PP<-scale(ceu.only$father.PP)
z.ceu.only.new$mother.PL<-scale(ceu.only$mother.PL)
z.ceu.only.new$mother.PP<-scale(ceu.only$mother.PP)
z.ceu.only.new$offspring.PL<-scale(ceu.only$offspring.PL)
z.ceu.only.new$offspring.PP<-scale(ceu.only$offspring.PP)
z.ceu.only.new$offspring.REF<-scale(ceu.only$offspring.REF)
z.ceu.only.new$offspring.ALT<-scale(ceu.only$offspring.ALT)
z.ceu.only.new$father.REF<-scale(ceu.only$father.REF)
z.ceu.only.new$mother.REF<-scale(ceu.only$mother.REF)

z.yri.only.new$QUAL<-scale(yri.only$QUAL)
z.yri.only.new$BaseQRankSum<-scale(yri.only$BaseQRankSum)
z.yri.only.new$FS<-scale(yri.only$FS)
z.yri.only.new$MQ<-scale(yri.only$MQ)
z.yri.only.new$MQRankSum<-scale(yri.only$MQRankSum)
z.yri.only.new$QD<-scale(yri.only$QD)
z.yri.only.new$ReadPosRankSum<-scale(yri.only$ReadPosRankSum)
z.yri.only.new$VQSLOD<-scale(yri.only$VQSLOD)
z.yri.only.new$father.DP<-scale(yri.only$father.DP)
z.yri.only.new$mother.DP<-scale(yri.only$mother.DP)
z.yri.only.new$offspring.DP<-scale(yri.only$offspring.DP)
z.yri.only.new$father.PL<-scale(yri.only$father.PL)
z.yri.only.new$father.PP<-scale(yri.only$father.PP)
z.yri.only.new$mother.PL<-scale(yri.only$mother.PL)
z.yri.only.new$mother.PP<-scale(yri.only$mother.PP)
z.yri.only.new$offspring.PL<-scale(yri.only$offspring.PL)
z.yri.only.new$offspring.PP<-scale(yri.only$offspring.PP)
z.yri.only.new$offspring.REF<-scale(yri.only$offspring.REF)
z.yri.only.new$offspring.ALT<-scale(yri.only$offspring.ALT)
z.yri.only.new$father.REF<-scale(yri.only$father.REF)
z.yri.only.new$mother.REF<-scale(yri.only$mother.REF)

z.hc.only.new$QUAL<-scale(compHcDfSnv3$QUAL)
z.hc.only.new$BaseQRankSum<-scale(compHcDfSnv3$BaseQRankSum)
z.hc.only.new$FS<-scale(compHcDfSnv3$FS)
z.hc.only.new$MQ<-scale(compHcDfSnv3$MQ)
z.hc.only.new$MQRankSum<-scale(compHcDfSnv3$MQRankSum)
z.hc.only.new$QD<-scale(compHcDfSnv3$QD)
z.hc.only.new$ReadPosRankSum<-scale(compHcDfSnv3$ReadPosRankSum)
z.hc.only.new$VQSLOD<-scale(compHcDfSnv3$VQSLOD)
z.hc.only.new$father.DP<-scale(compHcDfSnv3$father.DP)
z.hc.only.new$mother.DP<-scale(compHcDfSnv3$mother.DP)
z.hc.only.new$offspring.DP<-scale(compHcDfSnv3$offspring.DP)
z.hc.only.new$father.PL<-scale(compHcDfSnv3$father.PL)
z.hc.only.new$father.PP<-scale(compHcDfSnv3$father.PP)
z.hc.only.new$mother.PL<-scale(compHcDfSnv3$mother.PL)
z.hc.only.new$mother.PP<-scale(compHcDfSnv3$mother.PP)
z.hc.only.new$offspring.PL<-scale(compHcDfSnv3$offspring.PL)
z.hc.only.new$offspring.PP<-scale(compHcDfSnv3$offspring.PP)
z.hc.only.new$offspring.REF<-scale(compHcDfSnv3$offspring.REF)
z.hc.only.new$offspring.ALT<-scale(compHcDfSnv3$offspring.ALT)
z.hc.only.new$father.REF<-scale(compHcDfSnv3$father.REF)
z.hc.only.new$mother.REF<-scale(compHcDfSnv3$mother.REF)

p.ceu <- ggplot(z.ceu.only.new,aes(x=mother.REF,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
p.yri <- ggplot(z.yri.only.new,aes(x=mother.PL,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
p.hc <- ggplot(z.hc.only.new,aes(x=mother.PL,y=offspring.AF))+geom_point(size=1)
grid.arrange(p.ceu,p.yri,p.hc,nrow=3)

ceu.only.dnm <- ceu.only.new1[ceu.only.new1$DNMSNV!='.',]
yri.only.dnm <- yri.only.new1[yri.only.new1$DNMSNV!='.',]
p.ceu.dnm <- ggplot(ceu.only.dnm,aes(x=FS,y=offspring.AF,color=DNMSNV))+geom_point(size=1)+coord_cartesian(xlim=c(0:1),ylim=c(0:1))
p.yri.dnm <- ggplot(yri.only.dnm,aes(x=FS,y=offspring.AF,color=DNMSNV))+geom_point(size=1)+coord_cartesian(xlim=c(0:1),ylim=c(0:1))
grid.arrange(p.ceu.dnm,p.yri.dnm,p.hc)

#####|z-score|<=5#####
ceu.z.inlier<-z.ceu.only.new[z.ceu.only.new$QUAL<=5&z.ceu.only.new$QUAL>=-5&
                           z.ceu.only.new$BaseQRankSum<=5&z.ceu.only.new$BaseQRankSum>=-5&
                           z.ceu.only.new$FS<=5&z.ceu.only.new$FS>=-5&
                           z.ceu.only.new$MQ<=5&z.ceu.only.new$MQ>=-5&
                           z.ceu.only.new$MQRankSum<=5&z.ceu.only.new$MQRankSum>=-5&
                           z.ceu.only.new$QD<=5&z.ceu.only.new$QD>=-5&
                           z.ceu.only.new$ReadPosRankSum<=5&z.ceu.only.new$ReadPosRankSum>=-5&
                           z.ceu.only.new$VQSLOD<=5&z.ceu.only.new$VQSLOD>=-5&
                           z.ceu.only.new$father.DP<=5&z.ceu.only.new$father.DP>=-5&
                           z.ceu.only.new$mother.DP<=5&z.ceu.only.new$mother.DP>=-5&
                           z.ceu.only.new$offspring.DP<=5&z.ceu.only.new$offspring.DP>=-5&
                           z.ceu.only.new$father.PL<=5&z.ceu.only.new$father.PL>=-5&
                           z.ceu.only.new$father.PP<=5&z.ceu.only.new$father.PP>=-5&
                           z.ceu.only.new$mother.PL<=5&z.ceu.only.new$mother.PL>=-5&
                           z.ceu.only.new$mother.PP<=5&z.ceu.only.new$mother.PP>=-5&
                           z.ceu.only.new$offspring.PL<=5&z.ceu.only.new$offspring.PL>=-5&
                           z.ceu.only.new$offspring.PP<=5&z.ceu.only.new$offspring.PP>=-5,]
nrow(ceu.z.inlier)#6405 to 6211

yri.z.inlier<-z.yri.only.new[z.yri.only.new$QUAL<=5&z.yri.only.new$QUAL>=-5&
                               z.yri.only.new$BaseQRankSum<=5&z.yri.only.new$BaseQRankSum>=-5&
                               z.yri.only.new$FS<=5&z.yri.only.new$FS>=-5&
                               z.yri.only.new$MQ<=5&z.yri.only.new$MQ>=-5&
                               z.yri.only.new$MQRankSum<=5&z.yri.only.new$MQRankSum>=-5&
                               z.yri.only.new$QD<=5&z.yri.only.new$QD>=-5&
                               z.yri.only.new$ReadPosRankSum<=5&z.yri.only.new$ReadPosRankSum>=-5&
                               z.yri.only.new$VQSLOD<=5&z.yri.only.new$VQSLOD>=-5&
                               z.yri.only.new$father.DP<=5&z.yri.only.new$father.DP>=-5&
                               z.yri.only.new$mother.DP<=5&z.yri.only.new$mother.DP>=-5&
                               z.yri.only.new$offspring.DP<=5&z.yri.only.new$offspring.DP>=-5&
                               z.yri.only.new$father.PL<=5&z.yri.only.new$father.PL>=-5&
                               z.yri.only.new$father.PP<=5&z.yri.only.new$father.PP>=-5&
                               z.yri.only.new$mother.PL<=5&z.yri.only.new$mother.PL>=-5&
                               z.yri.only.new$mother.PP<=5&z.yri.only.new$mother.PP>=-5&
                               z.yri.only.new$offspring.PL<=5&z.yri.only.new$offspring.PL>=-5&
                               z.yri.only.new$offspring.PP<=5&z.yri.only.new$offspring.PP>=-5,]
nrow(yri.z.inlier)#2599 to 2545

hc.z.inlier<-z.hc.only.new[z.hc.only.new$QUAL<=5&z.hc.only.new$QUAL>=-5&
                               z.hc.only.new$BaseQRankSum<=5&z.hc.only.new$BaseQRankSum>=-5&
                               z.hc.only.new$FS<=5&z.hc.only.new$FS>=-5&
                               z.hc.only.new$MQ<=5&z.hc.only.new$MQ>=-5&
                               z.hc.only.new$MQRankSum<=5&z.hc.only.new$MQRankSum>=-5&
                               z.hc.only.new$QD<=5&z.hc.only.new$QD>=-5&
                               z.hc.only.new$ReadPosRankSum<=5&z.hc.only.new$ReadPosRankSum>=-5&
                               z.hc.only.new$VQSLOD<=5&z.hc.only.new$VQSLOD>=-5&
                               z.hc.only.new$father.DP<=5&z.hc.only.new$father.DP>=-5&
                               z.hc.only.new$mother.DP<=5&z.hc.only.new$mother.DP>=-5&
                               z.hc.only.new$offspring.DP<=5&z.hc.only.new$offspring.DP>=-5&
                               z.hc.only.new$father.PL<=5&z.hc.only.new$father.PL>=-5&
                               z.hc.only.new$father.PP<=5&z.hc.only.new$father.PP>=-5&
                               z.hc.only.new$mother.PL<=5&z.hc.only.new$mother.PL>=-5&
                               z.hc.only.new$mother.PP<=5&z.hc.only.new$mother.PP>=-5&
                               z.hc.only.new$offspring.PL<=5&z.hc.only.new$offspring.PL>=-5&
                               z.hc.only.new$offspring.PP<=5&z.hc.only.new$offspring.PP>=-5,]
nrow(hc.z.inlier)#1485 to 1443

#####random forest#####
#merge ceu and yri to one dataframe
whole.ceu.yri.z.inlier <- rbind(z.ceu.only.new,z.yri.only.new)

nrow(whole.ceu.yri.z.inlier[whole.ceu.yri.z.inlier$DNMSNV!='.',])
#8756 / 75
colnames(whole.ceu.yri.z.inlier)
p1 <- ggplot(ceu.z.inlier,aes(x=QD,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
p2 <- ggplot(yri.z.inlier,aes(x=QD,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
p3 <- ggplot(hc.z.inlier,aes(x=QD,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
grid.arrange(p1,p2,p3,nrow=3)
test.w <- whole.ceu.yri.z.inlier
test.w$DNMSNV[test.w$DNMSNV!='.']<-'DNM'
whole.dnm <- whole.ceu.yri.z.inlier[whole.ceu.yri.z.inlier$DNMSNV!='.',]
whole.not.dnm <- whole.ceu.yri.z.inlier[whole.ceu.yri.z.inlier$DNMSNV=='.',]

set.seed(1000)
dnm.intrain <- createDataPartition(y=whole.ceu.yri.z.inlier$DNMSNV,p=0.7,list=FALSE)
ceu.and.yri.dnm$offspring.AD
#SMOTE
#install.packages('DMwR')
library(DMwR)
test.w$DNMSNV<-as.factor(test.w$DNMSNV)
rawdata <- SMOTE(DNMSNV~QUAL+BaseQRankSum+FS+QD+MQ+
                   MQRankSum+ReadPosRankSum+VQSLOD+offspring.AF+offspring.DP+
                   father.DP+mother.DP+offspring.PL+offspring.PP+
                   father.PL+father.PP+mother.PL+mother.PP,data=test.w,perc.over = 3500,perc.under = 100,k=3)
#over 500, under 200
#dnm 450, not dnm 750
table(rawdata$DNMSNV)
ggplot(rawdata,aes(x=QD,y=offspring.AF,color=DNMSNV))+geom_point()

test.w.2 <- test.w
v.index.w1 <- c()
for (i in 3:5){
  for (j in 7:10){
    model2 <- randomForest(DNMSNV~QUAL+BaseQRankSum+FS+QD+MQ+
                             MQRankSum+ReadPosRankSum+VQSLOD+offspring.AF+offspring.DP+
                             father.DP+mother.DP+offspring.PL+offspring.PP+
                             father.PL+father.PP+mother.PL+mother.PP,data=test.w,mtry=i,ntree=j,importance=T)
    test.w.2$predDNM <- predict(model2,test.w)
    indices <- which(test.w.2$predDNM=='DNM')
    v.index.w1 <- c(v.index.w1,indices)
  }
}
length(v.index.w1)#w0 : 759? ,w1: 767
length(unique(v.index.w1))#w0 : 78, w1 : 77
unique.v.index.w1 <- unique(v.index.w1)
test.w[unique.v.index.w1,]$DNMSNV# 75 / 78 , 75 / 77
length(which(test.w$DNMSNV=='DNM'))#75
####han chinese####
test.hc <- hc.z.inlier
w <- 1
total.v.index <- c()
while (w<3){
  v.index <- c()
  for (i in 3){
    for (j in 30:50){
      model2 <- randomForest(DNMSNV~QUAL+BaseQRankSum+FS+QD+MQ+
                               MQRankSum+ReadPosRankSum+VQSLOD+offspring.AF+offspring.DP+
                               father.DP+mother.DP+offspring.PL+offspring.PP+
                               father.PL+father.PP+mother.PL+mother.PP,data=test.w,mtry=i,ntree=j,importance=T)
      test.hc$predDNM <- predict(model2,hc.z.inlier)
      indices <- which(test.hc$predDNM=='DNM')
      v.index <- c(v.index,indices)
    }
  }
  unique.v.index <- unique(v.index)
  total.v.index <- c(total.v.index,unique.v.index)
  w <- w + 1
}
unique.total.v.index <- unique(total.v.index)
length(total.v.index)
length(unique.total.v.index)
hist(test.hc[unique.total.v.index,]$offspring.AF)
try30.hc.index <- unique.total.v.index #22
o.try30.hc.index <- try30.hc.index[order(try30.hc.index)]

test.hc$predDNM<-'.'
test.hc$predDNM[o.try30.hc.index]<-'DNM'
ggplot(test.hc,aes(x=QD,y=offspring.AF,color=predDNM))+geom_point()+coord_cartesian(ylim=c(0:1))
#ggplot(ceu.only,aes(x=QD,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
