#install.packages("ROCR")
library(ROCR)
library(randomForest)
library(dplyr)
library(DMwR)
require(reshape2)
#install.packages("yardstick")
library(yardstick)
#####SMOTE#####
#library(DMwR)
test.w$DNMSNV<-as.factor(test.w$DNMSNV)
rawdata <- SMOTE(DNMSNV~QUAL+BaseQRankSum+FS+QD+MQ+
                MQRankSum+ReadPosRankSum+VQSLOD+offspring.AF+offspring.DP+
                father.DP+mother.DP+offspring.PL+offspring.PP+
                father.PL+father.PP+mother.PL+mother.PP+mother.REF+mother.ALT+
                offspring.ALT+offspring.REF,data=test.w,
                perc.over = 3000,perc.under = 110,k=5)
head(test.w)
table(rawdata$DNMSNV)#non:2475,DNM:2325
ggplot(rawdata,aes(x=mother.REF,y=offspring.AF,color=DNMSNV))+geom_point()

#####SMOTE#####
dnm.intrain <- createDataPartition(y=rawdata$DNMSNV,p=0.7,list=FALSE)
rawdata.2 <- rawdata[dnm.intrain,]
length(dnm.intrain)#3361
table(rawdata.2$DNMSNV)

v1 <- c(1:4800)
v2 <- v1[! v1 %in% dnm.intrain]
length(v2)#1439
testdata <- rawdata[v2,]
table(testdata$DNMSNV)#DNM:697, else: 742
testdata$preDNMSNV<-vector(length=nrow(testdata))

model3 <- randomForest(DNMSNV~QUAL+BaseQRankSum+FS+QD+MQ+
                         MQRankSum+ReadPosRankSum+VQSLOD+offspring.AF+offspring.DP+
                         father.DP+mother.DP+offspring.PL+offspring.PP+
                         father.PL+father.PP+mother.PL+mother.PP+mother.REF+father.REF+
                         offspring.ALT+offspring.REF,data=rawdata.2,
                       mtry=4,ntree=500,nodesize=1,importance=T,norm.votes=T,proximity=T)
testdata$predDNM <- predict(model3,testdata)

indices <- which(testdata$predDNM=='DNM')
length(indices)#701 -> 705

nrow(testdata[indices,][testdata[indices,]$DNMSNV=='DNM',])
#roc curve
testdata$prob <- predict(model3,testdata,type='prob')
test.pr <- prediction(testdata$prob[,2],testdata$DNMSNV)
test.prf <- performance(test.pr, measure='tpr',x.measure = 'fpr')
win.graph(); plot(test.prf, main='ROC of SMOTEd CEU+YRI test data')

#####ceu#####
test.ceu <- z.ceu.only.new

test.ceu$predDNM <- predict(model3,z.ceu.only.new)
test.ceu$prob <- predict(model3,z.ceu.only.new,type='prob')
pred.ceu.index <- which(test.ceu$predDNM=='DNM')
ceu.pred <- test.ceu[pred.ceu.index,]
nrow(ceu.pred)#249 -> 282
length(which(ceu.pred$DNMSNV!='.'))#43 -> 42
hist(ceu.pred[ceu.pred$DNMSNV!='.',]$prob[,2])

ggplot(test.ceu,aes(x=QD,y=offspring.AF,color=predDNM))+geom_point()
ggplot(ceu.pred,aes(x=QD,y=offspring.AF,color=DNMSNV))+geom_point()
hist(ceu.pred$offspring.AF)
ceu.pred.af <- ceu.pred[ceu.pred$offspring.AF>=0.3&ceu.pred$offspring.AF<=0.7,]
nrow(ceu.pred.af)#248 -> 281

nrow(test.ceu[test.ceu$prob[,2]>=0.3&
                test.ceu$prob[,2]<0.5,])#347
table(test.ceu[test.ceu$prob[,2]>=0.3&
           test.ceu$prob[,2]<0.5,]$DNMSNV)#one DNM
table(test.ceu[test.ceu$prob[,2]<0.3,]$DNMSNV)
#PR curve
ceu.pr <- prediction(test.ceu$prob[,2],test.ceu$DNMSNV)
ceu.prf <- performance(ceu.pr,measure='ppv',x.measure='tpr')
win.graph(); plot(ceu.prf,main='PRC of CEU trio')

slotNames(ceu.prf)
rd <- data.frame(x=ceu.prf@x.values[[1]],y=ceu.prf@y.values[[1]],alpha.values=ceu.prf@alpha.values[[1]])
which(rd$alpha.values==0.5)
rd[168,]
which(rd$alpha.values==0.3)
rd[263,]
win.graph(); ggplot(rd,aes(x=x,y=y))+geom_line()+
  geom_point(x=rd[168,]$x,y=rd[168,]$y)+
  geom_point(x=rd[263,]$x,y=rd[263,]$y)+
  labs(x='Recall',y='Precision')
####yri####
test.yri <- z.yri.only.new

test.yri$predDNM <- predict(model3,z.yri.only.new)
test.yri$prob <- predict(model3,z.yri.only.new,type='prob')

pred.yri.index <- which(test.yri$predDNM=='DNM')
yri.pred <- test.yri[pred.yri.index,]
nrow(yri.pred)#45 -> 75
length(which(yri.pred$DNMSNV!='.'))#31 -> 31
ggplot(test.yri,aes(x=QD,y=offspring.AF,color=predDNM))+geom_point()
ggplot(yri.pred,aes(x=QD,y=offspring.AF,color=DNMSNV))+geom_point()
hist(yri.pred$offspring.AF)
hist(yri.pred[yri.pred$DNMSNV!='.',]$prob[,2])
hist(yri.pred$prob[,2])
nrow(test.yri[test.yri$prob[,2]>=0.3&
                test.yri$prob[,2]<0.5,])#59
table(test.yri[test.yri$prob[,2]>=0.3&
                 test.yri$prob[,2]<0.5,]$DNMSNV)#one
hist()
#PR curve
yri.pr <- prediction(test.yri$prob[,2],test.yri$DNMSNV)
yri.prf <- performance(yri.pr, measure='ppv',x.measure = 'tpr')
win.graph(); plot(yri.prf, main='PRC of YRI trio')
str(yri.prf)

slotNames(yri.prf)
yri.rd <- data.frame(x=yri.prf@x.values[[1]],y=yri.prf@y.values[[1]],alpha.values=yri.prf@alpha.values[[1]])
head(yri.rd)
which(yri.rd$alpha.values==0.5)
yri.rd[58,]
which(yri.rd$alpha.values==0.302)
yri.rd[99,]
win.graph(); ggplot(yri.rd,aes(x=x,y=y))+geom_line()+
  geom_point(x=yri.rd[58,]$x,y=yri.rd[58,]$y)+
  geom_point(x=yri.rd[99,]$x,y=yri.rd[99,]$y)+
  geom_point(x=yri.rd[29,]$x,y=yri.rd[29,]$y)+
  labs(x='Recall',y='Precision')
####han chinese#####
test.hc <- z.hc.only.new

test.hc$predDNM <- predict(model3,z.hc.only.new)
test.hc$prob <- predict(model3,z.hc.only.new,type='prob')
hc.index <- which(test.hc$predDNM=='DNM')

length(hc.index)#39 -> 61

test.hc$predDNM<-'.'
test.hc$predDNM[hc.index]<-'DNM'
ggplot(test.hc,aes(x=QD,y=offspring.AF,color=predDNM))+geom_point(size=1)+coord_cartesian(ylim=c(0:1))
hist(test.hc[test.hc$predDNM=='DNM',]$offspring.AF)

hist(test.hc[test.hc$predDNM!='.',]$prob[,2])
nrow(test.hc[test.hc$prob[,2]>=0.3&
               test.hc$prob[,2]<0.5,])#183
which(test.hc$predDNM=='DNM')
# 61   78  127  172  195  230  245  297  303
#345  395  417  419  509  663  706  717  722
#729  735  791  882  893  935  975 1043 1052
#1059 1071 1072 1074 1128 1129 1135 1154 1156
#1244 1278 1298 1314 1336 1343 1344
canpos43 <- c(61,78,127,172,195,230,245,297,303,
             345,395,417,419,509,663,706,717,722,
             729,735,791,882,893,935,975,1043,1052,
             1059,1071,1072,1074,1128,1129,1135,1154,1156,
             1244,1278,1298,1314,1336,1343,1344)

# 61  127  172  195  198  230  245  297  302
#303  345  395  419  474  509  579  663  705
#706  717  791  795  804  882  893  935  942
#975 1052 1059 1062 1065 1071 1072 1074 1086
#1129 1135 1139 1154 1156 1244 1269 1278 1298
#1314 1320 1336 1343
canpos49 <- c(61,127,172,195,198,230,245,297,302,303,345,395,419,474,509,579,663,705,
  706,717,791,795,804,882,893,935,942,975,1052,1059,1062,1065,1071,1072,1074,1086,
  1129,1135,1139,1154,1156,1244,1269,1278,1298,1314,1320,1336,1343)
#39 possible dnms
# 61  127  172  195  230  245  302  303  345  395  509  561  579
#663  706  717  791  795  804  882  893  935 1052 1057 1059 1065
#1071 1072 1074 1129 1135 1154 1244 1269 1278 1298 1336 1343 1395

#61 possible dnms
# 54   64   81  108  130  177  200  229  235  266
#269  309  360  411  433  435  526  528  530  533
#597  621  674  729  740  751  752  758  781  786
#787  814  851  916  935  960  967  983 1011 1025
#1072 1081 1088 1091 1095 1101 1102 1104 1165 1166
#1383
hc.possdnm <- test.hc[test.hc$predDNM=='DNM',]
head(hc.possdnm)
nrow(hc.possdnm[hc.possdnm$TRANSITION=='1',])#20 -> 25
nrow(hc.possdnm[hc.possdnm$TRANSITION=='0',])#19 -> 36

hc.possdnm.49 <- test.hc[canpos49,]
head(hc.possdnm.49)
nrow(hc.possdnm.49[hc.possdnm.49$TRANSITION=='1',])#27
nrow(hc.possdnm.49[hc.possdnm.49$TRANSITION=='0',])#22

hc.possdnm.43 <- test.hc[canpos43,]
nrow(hc.possdnm.43[hc.possdnm.43$TRANSITION=='1',])#22
nrow(hc.possdnm.43[hc.possdnm.43$TRANSITION=='0',])#21

#####ASJ trio#####
raw.asj <- read.table('ASJ.table.txt',header=T)
head(raw.asj)
colnames(raw.asj)
colnames(raw.asj)[25:29]<-c("offspring.AD","offspring.AF","offspring.DP","offspring.PL","offspring.PP")
onlysnp.asj <- raw.asj[raw.asj$TRANSITION!='-1',]
nrow(onlysnp.asj)
split.asj <- transform(onlysnp.asj,offspring.PL=colsplit(offspring.PL,pattern='\\,',names=c(1,2,3)))
split.asj <- transform(split.asj,offspring.PP=colsplit(offspring.PP,pattern='\\,',names=c(1,2,3)))
split.asj <- transform(split.asj,offspring.AD=colsplit(offspring.AD,pattern='\\,',names=c(1,2)))

split.asj <- transform(split.asj,father.PP=colsplit(father.PP,pattern='\\,',names=c(1,2,3)))
split.asj <- transform(split.asj,father.PL=colsplit(father.PL,pattern='\\,',names=c(1,2,3)))
split.asj <- transform(split.asj,father.AD=colsplit(father.AD,pattern='\\,',names=c(1,2)))

split.asj <- transform(split.asj,mother.PP=colsplit(mother.PP,pattern='\\,',names=c(1,2,3)))
split.asj <- transform(split.asj,mother.PL=colsplit(mother.PL,pattern='\\,',names=c(1,2,3)))
split.asj <- transform(split.asj,mother.AD=colsplit(mother.AD,pattern='\\,',names=c(1,2)))

select.asj <- split.asj[,c(-15,-19,-20,-24)]
colnames(select.asj)
select.asj$offspring.PL <- split.asj$POS
select.asj$offspring.PP <- split.asj$POS

select.asj$father.PL <- split.asj$POS
select.asj$father.PP <- split.asj$POS

select.asj$mother.PL <- split.asj$POS
select.asj$mother.PP <- split.asj$POS

for (i in 1:nrow(split.asj)){
  if (split.asj$offspring.PL[1][i,]>=split.asj$offspring.PL[3][i,]){
    select.asj$offspring.PL[i]<-split.asj$offspring.PL[3][i,]
  }else{
    select.asj$offspring.PL[i]<-split.asj$offspring.PL[1][i,]
  }
}
for (i in 1:nrow(split.asj)){
  if (split.asj$offspring.PP[1][i,]>=split.asj$offspring.PP[3][i,]){
    select.asj$offspring.PP[i]<-split.asj$offspring.PP[3][i,]
  }else{
    select.asj$offspring.PP[i]<-split.asj$offspring.PP[1][i,]
  }
}
for (i in 1:nrow(split.asj)){
  if (split.asj$father.PL[2][i,]>=split.asj$father.PL[3][i,]){
    select.asj$father.PL[i]<-split.asj$father.PL[3][i,]
  }else{
    select.asj$father.PL[i]<-split.asj$father.PL[2][i,]
  }
}
for (i in 1:nrow(split.asj)){
  if (split.asj$father.PP[2][i,]>=split.asj$father.PP[3][i,]){
    select.asj$father.PP[i]<-split.asj$father.PP[3][i,]
  }else{
    select.asj$father.PP[i]<-split.asj$father.PP[2][i,]
  }
}
for (i in 1:nrow(split.asj)){
  if (split.asj$mother.PL[2][i,]>=split.asj$mother.PL[3][i,]){
    select.asj$mother.PL[i]<-split.asj$mother.PL[3][i,]
  }else{
    select.asj$mother.PL[i]<-split.asj$mother.PL[2][i,]
  }
}
for (i in 1:nrow(split.asj)){
  if (split.asj$mother.PP[2][i,]>=split.asj$mother.PP[3][i,]){
    select.asj$mother.PP[i]<-split.asj$mother.PP[3][i,]
  }else{
    select.asj$mother.PP[i]<-split.asj$mother.PP[2][i,]
  }
}
for (i in 1:nrow(split.asj)){
  if (split.asj$mother.PP[2][i,]>=split.asj$mother.PP[3][i,]){
    select.asj$mother.PP[i]<-split.asj$mother.PP[3][i,]
  }else{
    select.asj$mother.PP[i]<-split.asj$mother.PP[2][i,]
  }
}
#AD
#AF calculation
for (i in 1:nrow(split.asj)){
  select.asj$mother.REF[i] <- split.asj$mother.AD[i,1]
  select.asj$mother.ALT[i] <- split.asj$mother.AD[i,2]
  #select.asj$mother.AF[i] <- round(split.asj$mother.AD[i,2] / split.asj$mother.DP[i],digits = 2)
  select.asj$father.REF[i] <- split.asj$father.AD[i,1]
  select.asj$father.ALT[i] <- split.asj$father.AD[i,2]
  #select.asj$father.AF[i] <- round(split.asj$father.AD[i,2] / split.asj$father.DP[i], digits=2)
  select.asj$offspring.REF[i] <- split.asj$offspring.AD[i,1]
  select.asj$offspring.ALT[i] <- split.asj$offspring.AD[i,2]
  #select.asj$offspring.AF[i] <- round(split.asj$offspring.AD[i,2] / split.asj$offspring.DP[i], digits=2)
}
head(select.asj)
nrow(select.asj)#19886
offspring.het.asj <- select.asj[select.asj$offspring.ALT!='0',]
nrow(offspring.het.asj)#11950
mom.hom.off.het.asj <- offspring.het.asj[offspring.het.asj$mother.ALT=='0'|
                                              offspring.het.asj$mother.ALT=='1',]
nrow(mom.hom.off.het.asj)#3862
parent.hom.off.het.asj <- mom.hom.off.het.asj[mom.hom.off.het.asj$father.ALT=='0'|
                                                mom.hom.off.het.asj$father.ALT=='1',]
nrow(parent.hom.off.het.asj)#3246
candidate.asj <- na.omit(parent.hom.off.het.asj)
nrow(candidate.asj)#3211
#z score
colnames(candidate.asj)

z.asj$QUAL<-scale(candidate.asj$QUAL)
z.asj$BaseQRankSum<-scale(candidate.asj$BaseQRankSum)
z.asj$FS<-scale(candidate.asj$FS)
z.asj$MQ<-scale(candidate.asj$MQ)
z.asj$MQRankSum<-scale(candidate.asj$MQRankSum)
z.asj$QD<-scale(candidate.asj$QD)
z.asj$ReadPosRankSum<-scale(candidate.asj$ReadPosRankSum)
z.asj$VQSLOD<-scale(candidate.asj$VQSLOD)
z.asj$father.DP<-scale(candidate.asj$father.DP)
z.asj$mother.DP<-scale(candidate.asj$mother.DP)
z.asj$offspring.DP<-scale(candidate.asj$offspring.DP)
z.asj$father.PL<-scale(candidate.asj$father.PL)
z.asj$father.PP<-scale(candidate.asj$father.PP)
z.asj$mother.PL<-scale(candidate.asj$mother.PL)
z.asj$mother.PP<-scale(candidate.asj$mother.PP)
z.asj$offspring.PL<-scale(candidate.asj$offspring.PL)
z.asj$offspring.PP<-scale(candidate.asj$offspring.PP)
head(z.asj)

z.asj$predDNM <- predict(model3,z.asj)
length(which(z.asj$predDNM=='DNM'))#90 -> 34
z.asj$prob <- predict(model3,z.asj,type='prob')

ggplot(z.asj,aes(x=QD,y=offspring.AF,color=predDNM))+geom_point()
win.graph(); hist(z.asj[z.asj$predDNM=='DNM',]$offspring.AF)
table(z.asj[z.asj$predDNM=='DNM',]$CHROM)
plot(z.asj[z.asj$predDNM=='DNM',]$CHROM)
table(z.asj$CHROM)

pdn.asj <- z.asj[z.asj$predDNM=='DNM',]
nrow(pdn.asj)
pdn.asj[pdn.asj$CHROM=='chr5',]
z.asj[z.asj$POS==15897669,]#prob 0.306
z.asj[z.asj$POS==29864622,]#80341085
hist(z.asj$prob[,2])
hist(pdn.asj$prob[,2])
nrow(z.asj[z.asj$prob[,2]>=0.306,])#334
hist(z.asj[z.asj$prob[,2]>=0.306,]$prob[,2])
whole.z.asj <- z.asj
whole.z.asj$denovoconf <- vector(length=nrow(z.asj))
whole.z.asj$denovoconf[whole.z.asj$prob[,2]>0.306&whole.z.asj$prob[,2]<0.5]<-'lowConf'
whole.z.asj$denovoconf[whole.z.asj$prob[,2]==0.306]<-'true'

whole.z.asj$denovoconf[whole.z.asj$prob[,2]>=0.5]<-'highConf'
whole.z.asj$denovoconf[whole.z.asj$prob[,2]<0.306]<-'.'
e <- whole.z.asj[whole.z.asj$denovoconf!='.',]
ggplot(e,aes(x=QD,y=offspring.AF,color=denovoconf))+geom_point()
ti.asj <- pdn.asj[pdn.asj$TRANSITION==1,]
tv.asj <- pdn.asj[pdn.asj$TRANSITION==0,]
nrow(ti.asj)/nrow(tv.asj)#1.9 -> 2.4
z.asj[z.asj$offspring.AF>=0.5&z.asj$offspring.AF<=0.6,c(1,2,34,35)]
pdn.asj$offspring.ALT

e[e$denovoconf=='lowConf',c(1,2,35)]
colnames(e)

hist(whole.z.asj[whole.z.asj$prob[,2]>=0.3&
              whole.z.asj$prob[,2]<0.5,]$offspring.AF)
af<-whole.z.asj[whole.z.asj$offspring.AF>=0.3&
  whole.z.asj$offspring.AF<=0.7,]
nrow(af)#645
nrow(whole.z.asj[whole.z.asj$prob[,2]>=0.3&
                   whole.z.asj$prob[,2]<0.5,])#310
af[250:300,c(1,2,22,36)]
colnames(af)
af[af$POS==43163632,]

split.asj[split.asj$POS==43163632,]

whole.z.asj[whole.z.asj$prob[,2]<0.3,c(1,2,22,34)]
split.asj[split.asj$POS==83958816,]#83958816,160572530
whole.z.asj[whole.z.asj$POS==160572530,]
