#install.packages('rgl')
#install.packages("randomForest")
require(reshape2)
library(ggplot2)
library(rgl)
library(randomForest)
DNMpos <- c(75884343,110583335,182974758,39556621,152899032,182693277,101454745,
            118900031,15227519,104624818,6466106,52638226,126385924,145107247,52120843,145808310,
            160334960,21568355,74680107,38096405,89775629,123350598,56256293,56256294,8598739,120325447,40851625,
            76051551,85631211,117264626,56763353,78809184,50953965,51248561,58669774,85715561,87791919,99175976,29601086,
            82540504,53502801,71375843,80712655,74117587,6661912,7195809,55356548,20165354,24262476)
getwd()
SMA.CEU.df <- read.table('table.SMA.includeFiltered.PL.PP.AF.dbSNP.annotated.qual.thres.0.MV.BI.Het.CEU.wgs.consensus.20131118.snps_indels.txt',header=T)
colnames(SMA.CEU.df)
SMA.CEU.df2 <- transform(SMA.CEU.df,NA12878.PL=colsplit(NA12878.PL,pattern='\\,',names=c(1,2,3)))
SMA.CEU.df2 <- transform(SMA.CEU.df2,NA12878.PP=colsplit(NA12878.PP,pattern='\\,',names=c(1,2,3)))
SMA.CEU.df2 <- transform(SMA.CEU.df2,NA12891.PP=colsplit(NA12891.PP,pattern='\\,',names=c(1,2,3)))
SMA.CEU.df2 <- transform(SMA.CEU.df2,NA12891.PL=colsplit(NA12891.PL,pattern='\\,',names=c(1,2,3)))
SMA.CEU.df2 <- transform(SMA.CEU.df2,NA12892.PP=colsplit(NA12892.PP,pattern='\\,',names=c(1,2,3)))
SMA.CEU.df2 <- transform(SMA.CEU.df2,NA12892.PL=colsplit(NA12892.PL,pattern='\\,',names=c(1,2,3)))
SMA.CEU.df3 <- SMA.CEU.df2[,c(1:17,20:23,26:29,32)]
SMA.CEU.df3$NA12878.PL <- SMA.CEU.df2$POS
SMA.CEU.df3$NA12878.PP <- SMA.CEU.df2$POS
SMA.CEU.df3$NA12891.PL <- SMA.CEU.df2$POS
SMA.CEU.df3$NA12891.PP <- SMA.CEU.df2$POS
SMA.CEU.df3$NA12892.PL <- SMA.CEU.df2$POS
SMA.CEU.df3$NA12892.PP <- SMA.CEU.df2$POS
colnames(SMA.CEU.df3)

for (i in 1:nrow(SMA.CEU.df2)){
  if (SMA.CEU.df2$NA12891.PL[2][i,]>=SMA.CEU.df2$NA12891.PL[3][i,]){
    SMA.CEU.df3$NA12891.PL[i]<-SMA.CEU.df2$NA12891.PL[3][i,]
  }else{
    SMA.CEU.df3$NA12891.PL[i]<-SMA.CEU.df2$NA12891.PL[2][i,]
  }
}
head(SMA.CEU.df3)

SMA.CEU.df$DNMSNV <- SMA.CEU.df$POS
for (pos in SMA.CEU.df$POS){
  if (pos %in% DNMpos){
    cs <- which(pos==SMA.CEU.df$POS)
    SMA.CEU.df$DNMSNV[cs]<-'red'
  }else{
    ws <- which(pos==SMA.CEU.df$POS)
    SMA.CEU.df$DNMSNV[ws]<-'blue'
  }
}
length(which(SMA.CEU.df$DNMSNV=='red'))#43
ValidatedDNMSNV<-SMA.CEU.df[SMA.CEU.df$DNMSNV=='red',]
nrow(ValidatedDNMSNV)
DNM_ti_tv <- nrow(ValidatedDNMSNV[ValidatedDNMSNV$TRANSITION==1,])/nrow(ValidatedDNMSNV[ValidatedDNMSNV$TRANSITION==0,])
DNM_ti_tv#3.3
dfSnv3<-SMA.CEU.df3[SMA.CEU.df3$TRANSITION!=-1,]
ggplot(dfSnv,aes(x=NA12878.DP,y=NA12878.AF,z=MQ,color=DNMSNV))+geom_point(size=1)
ggplot(dfSnv,aes(x=VQSLOD,y=NA12878.AF,z=MQ,color=DNMSNV))+geom_point(size=1)
ggplot(dfSnv,aes(x=BaseQRankSum,y=NA12878.AF,color=DNMSNV))+geom_point(size=.5)
ggplot(dfSnv,aes(x=ClippingRankSum,y=NA12878.AF,color=DNMSNV))+geom_point(size=.5)
ggplot(dfSnv,aes(x=MQRankSum,y=NA12878.AF,color=DNMSNV))+geom_point(size=.5)

dfSnv3$DNMSNV<-factor(dfSnv3$DNMSNV,levels = c('red','blue'))
str(dfSnv3$DNMSNV)
plot3d(dfSnv$VQSLOD,dfSnv$QD,dfSnv$NA12878.AF,col=as.integer(dfSnv$DNMSNV))

nrow(dfSnv)#6469
####random forest####
compDfSnv3 <- dfSnv3[complete.cases(dfSnv3),]
nrow(compDfSnv3)
head(compDfSnv3)
colnames(compDfSnv3)
compDfSnv3 <- compDfSnv3[,c(2,4,7:14,16:17,20:21,24:25,27:33)]
colnames(compDfSnv3)
colnames(compDfSnv3)[c(11:16,18:23)]<-c("offspring.AF","offspring.DP",
                                        "father.AF","father.DP",
                                        "mother.AF","mother.DP",
                                        "offspring.PL","offspring.PP","father.PL","father.PP","mother.PL","mother.PP")

model1 <- randomForest(DNMSNV~.,data=compDfSnv3,mtry=5,ntree=100,importance=T)
importance(model1)
#han chinese data
hcdf <- read.table('table.SMA.PL.PP.JL.JP.dbsnp151.grch37.All_20180423.masked.isHet.tranch99.9.txt',header=T)
nrow(hcdf)
hcDfSnv <- hcdf[hcdf$TRANSITION!=(-1),]
nrow(hcDfSnv)#1494
compHcDfSnv <- hcDfSnv[complete.cases(hcDfSnv),]
colnames(compHcDfSnv)
nrow(compHcDfSnv)#1485
compHcDfSnv<-compHcDfSnv[,c(2,4,7:14,16:17,19:20,24:25,27:28,32:33,35:36)]
colnames(compHcDfSnv)[11:22]<-c("father.AF","father.DP","father.PL","father.PP",
                                "offspring.AF","offspring.DP","offspring.PL","offspring.PP",
                                "mother.AF","mother.DP","mother.PL","mother.PP")

compHcDfSnv2 <- transform(compHcDfSnv,offspring.PL=colsplit(offspring.PL,pattern='\\,',names=c(1,2,3)))
compHcDfSnv2 <- transform(compHcDfSnv2,offspring.PP=colsplit(offspring.PP,pattern='\\,',names=c(1,2,3)))
compHcDfSnv2 <- transform(compHcDfSnv2,father.PP=colsplit(father.PP,pattern='\\,',names=c(1,2,3)))
compHcDfSnv2 <- transform(compHcDfSnv2,father.PL=colsplit(father.PL,pattern='\\,',names=c(1,2,3)))
compHcDfSnv2 <- transform(compHcDfSnv2,mother.PP=colsplit(mother.PP,pattern='\\,',names=c(1,2,3)))
compHcDfSnv2 <- transform(compHcDfSnv2,mother.PL=colsplit(mother.PL,pattern='\\,',names=c(1,2,3)))

compHcDfSnv3 <- compHcDfSnv
compHcDfSnv3$offspring.PL <- compHcDfSnv2$POS
compHcDfSnv3$offspring.PP <- compHcDfSnv2$POS
compHcDfSnv3$father.PL <- compHcDfSnv2$POS
compHcDfSnv3$father.PP <- compHcDfSnv2$POS
compHcDfSnv3$mother.PL <- compHcDfSnv2$POS
compHcDfSnv3$mother.PP <- compHcDfSnv2$POS
for (i in 1:nrow(compHcDfSnv2)){
  if (compHcDfSnv2$offspring.PL[1][i,]>=compHcDfSnv2$offspring.PL[3][i,]){
    compHcDfSnv3$offspring.PL[i]<-compHcDfSnv2$offspring.PL[3][i,]
  }else{
    compHcDfSnv3$offspring.PL[i]<-compHcDfSnv2$offspring.PL[1][i,]
  }
}
for (i in 1:nrow(compHcDfSnv2)){
  if (compHcDfSnv2$offspring.PP[1][i,]>=compHcDfSnv2$offspring.PP[3][i,]){
    compHcDfSnv3$offspring.PP[i]<-compHcDfSnv2$offspring.PP[3][i,]
  }else{
    compHcDfSnv3$offspring.PP[i]<-compHcDfSnv2$offspring.PP[1][i,]
  }
}
for (i in 1:nrow(compHcDfSnv2)){
  if (compHcDfSnv2$mother.PP[2][i,]>=compHcDfSnv2$mother.PP[3][i,]){
    compHcDfSnv3$mother.PP[i]<-compHcDfSnv2$mother.PP[3][i,]
  }else{
    compHcDfSnv3$mother.PP[i]<-compHcDfSnv2$mother.PP[2][i,]
  }
}
for (i in 1:nrow(compHcDfSnv2)){
  if (compHcDfSnv2$mother.PL[2][i,]>=compHcDfSnv2$mother.PL[3][i,]){
    compHcDfSnv3$mother.PL[i]<-compHcDfSnv2$mother.PL[3][i,]
  }else{
    compHcDfSnv3$mother.PL[i]<-compHcDfSnv2$mother.PL[2][i,]
  }
}
for (i in 1:nrow(compHcDfSnv2)){
  if (compHcDfSnv2$father.PP[2][i,]>=compHcDfSnv2$father.PP[3][i,]){
    compHcDfSnv3$father.PP[i]<-compHcDfSnv2$father.PP[3][i,]
  }else{
    compHcDfSnv3$father.PP[i]<-compHcDfSnv2$father.PP[2][i,]
  }
}
for (i in 1:nrow(compHcDfSnv2)){
  if (compHcDfSnv2$father.PL[2][i,]>=compHcDfSnv2$father.PL[3][i,]){
    compHcDfSnv3$father.PL[i]<-compHcDfSnv2$father.PL[3][i,]
  }else{
    compHcDfSnv3$father.PL[i]<-compHcDfSnv2$father.PL[2][i,]
  }
}
head(compHcDfSnv3)
length(compHcDfSnv3)
length(compDfSnv3)