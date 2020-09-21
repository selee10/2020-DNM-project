#install.packages('rgl')
#install.packages("randomForest")
#install.packages('caret')
require(reshape2)
library(ggplot2)
library(rgl)
library(randomForest)
library(caret)
library(gridExtra)
DNMpos <- c(75884343,110583335,182974758,39556621,152899032,182693277,101454745,
            118900031,15227519,104624818,6466106,52638226,126385924,145107247,52120843,145808310,
            160334960,21568355,74680107,38096405,89775629,123350598,56256293,56256294,8598739,120325447,40851625,
            76051551,85631211,117264626,56763353,78809184,50953965,51248561,58669774,85715561,87791919,99175976,29601086,
            82540504,53502801,71375843,80712655,74117587,6661912,7195809,55356548,20165354,24262476)
getwd()
SMA.CEU.df <- read.table('table.SMA.includeFiltered.PL.PP.AF.dbSNP.annotated.qual.thres.0.MV.BI.Het.CEU.wgs.consensus.20131118.snps_indels.txt',header=T)
SMA.CEU.df$DNMSNV <- SMA.CEU.df$POS
for (pos in SMA.CEU.df$POS){
  if (pos %in% DNMpos){
    cs <- which(pos==SMA.CEU.df$POS)
    SMA.CEU.df$DNMSNV[cs]<-'DNM'
  }else{
    ws <- which(pos==SMA.CEU.df$POS)
    SMA.CEU.df$DNMSNV[ws]<-'.'
  }
}
SMA.CEU.df2 <- transform(SMA.CEU.df,NA12878.PL=colsplit(NA12878.PL,pattern='\\,',names=c(1,2,3)))
SMA.CEU.df2 <- transform(SMA.CEU.df2,NA12878.PP=colsplit(NA12878.PP,pattern='\\,',names=c(1,2,3)))
SMA.CEU.df2 <- transform(SMA.CEU.df2,NA12891.PP=colsplit(NA12891.PP,pattern='\\,',names=c(1,2,3)))
SMA.CEU.df2 <- transform(SMA.CEU.df2,NA12891.PL=colsplit(NA12891.PL,pattern='\\,',names=c(1,2,3)))
SMA.CEU.df2 <- transform(SMA.CEU.df2,NA12892.PP=colsplit(NA12892.PP,pattern='\\,',names=c(1,2,3)))
SMA.CEU.df2 <- transform(SMA.CEU.df2,NA12892.PL=colsplit(NA12892.PL,pattern='\\,',names=c(1,2,3)))
SMA.CEU.df2 <- transform(SMA.CEU.df2,NA12878.AD=colsplit(NA12878.AD,pattern='\\,',names=c(1,2)))
SMA.CEU.df2 <- transform(SMA.CEU.df2,NA12891.AD=colsplit(NA12891.AD,pattern='\\,',names=c(1,2)))
SMA.CEU.df2 <- transform(SMA.CEU.df2,NA12892.AD=colsplit(NA12892.AD,pattern='\\,',names=c(1,2)))
SMA.CEU.df2$NA12892.AD
colnames(SMA.CEU.df2)
SMA.CEU.df3 <- SMA.CEU.df2[,c(1:14,16:17,19:20,22:23,25:26,28:29,31:33)]
SMA.CEU.df3$NA12878.PL <- SMA.CEU.df2$POS
SMA.CEU.df3$NA12878.PP <- SMA.CEU.df2$POS
SMA.CEU.df3$NA12891.PL <- SMA.CEU.df2$POS
SMA.CEU.df3$NA12891.PP <- SMA.CEU.df2$POS
SMA.CEU.df3$NA12892.PL <- SMA.CEU.df2$POS
SMA.CEU.df3$NA12892.PP <- SMA.CEU.df2$POS

colnames(SMA.CEU.df3)

for (i in 1:nrow(SMA.CEU.df2)){
  if (SMA.CEU.df2$NA12878.PL[1][i,]>=SMA.CEU.df2$NA12878.PL[3][i,]){
    SMA.CEU.df3$NA12878.PL[i]<-SMA.CEU.df2$NA12878.PL[3][i,]
  }else{
    SMA.CEU.df3$NA12878.PL[i]<-SMA.CEU.df2$NA12878.PL[1][i,]
  }
}

for (i in 1:nrow(SMA.CEU.df2)){
  if (SMA.CEU.df2$NA12878.PP[1][i,]>=SMA.CEU.df2$NA12878.PP[3][i,]){
    SMA.CEU.df3$NA12878.PP[i]<-SMA.CEU.df2$NA12878.PP[3][i,]
  }else{
    SMA.CEU.df3$NA12878.PP[i]<-SMA.CEU.df2$NA12878.PP[1][i,]
  }
}

for (i in 1:nrow(SMA.CEU.df2)){
  if (SMA.CEU.df2$NA12891.PL[2][i,]>=SMA.CEU.df2$NA12891.PL[3][i,]){
    SMA.CEU.df3$NA12891.PL[i]<-SMA.CEU.df2$NA12891.PL[3][i,]
  }else{
    SMA.CEU.df3$NA12891.PL[i]<-SMA.CEU.df2$NA12891.PL[2][i,]
  }
}

for (i in 1:nrow(SMA.CEU.df2)){
  if (SMA.CEU.df2$NA12891.PP[2][i,]>=SMA.CEU.df2$NA12891.PP[3][i,]){
    SMA.CEU.df3$NA12891.PP[i]<-SMA.CEU.df2$NA12891.PP[3][i,]
  }else{
    SMA.CEU.df3$NA12891.PP[i]<-SMA.CEU.df2$NA12891.PP[2][i,]
  }
}

for (i in 1:nrow(SMA.CEU.df2)){
  if (SMA.CEU.df2$NA12892.PL[2][i,]>=SMA.CEU.df2$NA12892.PL[3][i,]){
    SMA.CEU.df3$NA12892.PL[i]<-SMA.CEU.df2$NA12892.PL[3][i,]
  }else{
    SMA.CEU.df3$NA12892.PL[i]<-SMA.CEU.df2$NA12892.PL[2][i,]
  }
}

for (i in 1:nrow(SMA.CEU.df2)){
  if (SMA.CEU.df2$NA12892.PP[2][i,]>=SMA.CEU.df2$NA12892.PP[3][i,]){
    SMA.CEU.df3$NA12892.PP[i]<-SMA.CEU.df2$NA12892.PP[3][i,]
  }else{
    SMA.CEU.df3$NA12892.PP[i]<-SMA.CEU.df2$NA12892.PP[2][i,]
  }
}
#AD
for (i in 1:nrow(SMA.CEU.df2)){
  SMA.CEU.df3$NA12892.REF[i] <- SMA.CEU.df2$NA12892.AD[i,1]
  SMA.CEU.df3$NA12892.ALT[i] <- SMA.CEU.df2$NA12892.AD[i,2]
  #select.asj$mother.AF[i] <- round(split.asj$mother.AD[i,2] / split.asj$mother.DP[i],digits = 2)
  SMA.CEU.df3$NA12891.REF[i] <- SMA.CEU.df2$NA12891.AD[i,1]
  SMA.CEU.df3$NA12891.ALT[i] <- SMA.CEU.df2$NA12891.AD[i,2]
  #select.asj$father.AF[i] <- round(split.asj$father.AD[i,2] / split.asj$father.DP[i], digits=2)
  SMA.CEU.df3$NA12878.REF[i] <- SMA.CEU.df2$NA12878.AD[i,1]
  SMA.CEU.df3$NA12878.ALT[i] <- SMA.CEU.df2$NA12878.AD[i,2]
  #select.asj$offspring.AF[i] <- round(split.asj$offspring.AD[i,2] / split.asj$offspring.DP[i], digits=2)
}
head(SMA.CEU.df3)

length(which(SMA.CEU.df3$DNMSNV=='DNM'))#43
ValidatedDNMSNV<-SMA.CEU.df[SMA.CEU.df$DNMSNV=='DNM',]
nrow(ValidatedDNMSNV)
DNM_ti_tv <- nrow(ValidatedDNMSNV[ValidatedDNMSNV$TRANSITION==1,])/nrow(ValidatedDNMSNV[ValidatedDNMSNV$TRANSITION==0,])
DNM_ti_tv#3.3
ceu.df.snv3<-SMA.CEU.df3[SMA.CEU.df3$TRANSITION!=-1,]

ceu.df.snv3$DNMSNV<-factor(ceu.df.snv3$DNMSNV,levels = c('DNM','.'))
str(ceu.df.snv3$DNMSNV)

#dfSnv3<-SMA.CEU.df3[SMA.CEU.df3$TRANSITION!=-1,](See CEU_DNM_SNV.R)
#dfSnv3$DNMSNV<-factor(dfSnv3$DNMSNV,levels = c('red','blue'))
nrow(dfSnv3)#6469
####random forest####
compDfSnv3 <- ceu.df.snv3[complete.cases(ceu.df.snv3),]
nrow(compDfSnv3)#6405
head(compDfSnv3)
colnames(compDfSnv3)
compDfSnv3.1 <- compDfSnv3
colnames(compDfSnv3.1)

colnames(compDfSnv3.1)[c(15:26,28:33)]<-c("offspring.AF","offspring.DP","offspring.PL","offspring.PP",
                                 "father.AF","father.DP","father.PL","father.PP",
                                 "mother.AF","mother.DP","mother.PL","mother.PP",
                                 "mother.REF","mother.ALT",
                                 "father.REF","father.ALT",
                                 "offspring.REF","offspring.ALT")
compDfSnv3.2 <- compDfSnv3.1#[,c(27:33)]

#han chinese data
hcdf <- read.table('table.SMA.PL.PP.JL.JP.dbsnp151.grch37.All_20180423.masked.isHet.tranch99.9.txt',header=T)
nrow(hcdf)
hcDfSnv <- hcdf[hcdf$TRANSITION!=(-1),]
nrow(hcDfSnv)#1494
compHcDfSnv <- hcDfSnv[complete.cases(hcDfSnv),]
nrow(compHcDfSnv)#1485
colnames(compHcDfSnv)
compHcDfSnv<-compHcDfSnv[,c(1:17,19:20,23:25,27:28,31:33,35:36)]
colnames(compHcDfSnv)[15:29]<-c("father.AD","father.AF","father.DP","father.PL","father.PP",
                               "offspring.AD","offspring.AF","offspring.DP","offspring.PL","offspring.PP",
                               "mother.AD","mother.AF","mother.DP","mother.PL","mother.PP")

compHcDfSnv2 <- transform(compHcDfSnv,offspring.PL=colsplit(offspring.PL,pattern='\\,',names=c(1,2,3)))
compHcDfSnv2 <- transform(compHcDfSnv2,offspring.PP=colsplit(offspring.PP,pattern='\\,',names=c(1,2,3)))
compHcDfSnv2 <- transform(compHcDfSnv2,father.PP=colsplit(father.PP,pattern='\\,',names=c(1,2,3)))
compHcDfSnv2 <- transform(compHcDfSnv2,father.PL=colsplit(father.PL,pattern='\\,',names=c(1,2,3)))
compHcDfSnv2 <- transform(compHcDfSnv2,mother.PP=colsplit(mother.PP,pattern='\\,',names=c(1,2,3)))
compHcDfSnv2 <- transform(compHcDfSnv2,mother.PL=colsplit(mother.PL,pattern='\\,',names=c(1,2,3)))
compHcDfSnv2 <- transform(compHcDfSnv2,mother.AD=colsplit(mother.AD,pattern='\\,',names=c(1,2)))
compHcDfSnv2 <- transform(compHcDfSnv2,father.AD=colsplit(father.AD,pattern='\\,',names=c(1,2)))
compHcDfSnv2 <- transform(compHcDfSnv2,offspring.AD=colsplit(offspring.AD,pattern='\\,',names=c(1,2)))

compHcDfSnv3 <- compHcDfSnv2[,c(1:14,16:19,21:24,26:29)]
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
  if (compHcDfSnv2$father.PL[2][i,]>=compHcDfSnv2$father.PL[3][i,]){
    compHcDfSnv3$father.PL[i]<-compHcDfSnv2$father.PL[3][i,]
  }else{
    compHcDfSnv3$father.PL[i]<-compHcDfSnv2$father.PL[2][i,]
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
  if (compHcDfSnv2$mother.PL[2][i,]>=compHcDfSnv2$mother.PL[3][i,]){
    compHcDfSnv3$mother.PL[i]<-compHcDfSnv2$mother.PL[3][i,]
  }else{
    compHcDfSnv3$mother.PL[i]<-compHcDfSnv2$mother.PL[2][i,]
  }
}
for (i in 1:nrow(compHcDfSnv2)){
  if (compHcDfSnv2$mother.PP[2][i,]>=compHcDfSnv2$mother.PP[3][i,]){
    compHcDfSnv3$mother.PP[i]<-compHcDfSnv2$mother.PP[3][i,]
  }else{
    compHcDfSnv3$mother.PP[i]<-compHcDfSnv2$mother.PP[2][i,]
  }
}
#AD
for (i in 1:nrow(compHcDfSnv2)){
  compHcDfSnv3$mother.REF[i] <- compHcDfSnv2$mother.AD[i,1]
  compHcDfSnv3$mother.ALT[i] <- compHcDfSnv2$mother.AD[i,2]
  #select.asj$mother.AF[i] <- round(split.asj$mother.AD[i,2] / split.asj$mother.DP[i],digits = 2)
  compHcDfSnv3$father.REF[i] <- compHcDfSnv2$father.AD[i,1]
  compHcDfSnv3$father.ALT[i] <- compHcDfSnv2$father.AD[i,2]
  #select.asj$father.AF[i] <- round(split.asj$father.AD[i,2] / split.asj$father.DP[i], digits=2)
  compHcDfSnv3$offspring.REF[i] <- compHcDfSnv2$offspring.AD[i,1]
  compHcDfSnv3$offspring.ALT[i] <- compHcDfSnv2$offspring.AD[i,2]
  #select.asj$offspring.AF[i] <- round(split.asj$offspring.AD[i,2] / split.asj$offspring.DP[i], digits=2)
}
colnames(compHcDfSnv3)

##model apply on han chinese trio##
compHcDfSnv3.1 <- compHcDfSnv3
compDfSnv3.1$predDNMSNV<-vector(length=nrow(compDfSnv3.1))

###YRI trio###
yri.dnm.pos <- c(182996168,215782887,245546833,151668827,26093263,11818938,117201136,
                            72217598,41332942,160248391,3829079,74029053,20418371,75866623,125759035,
                            125760549,140510278,3733188,76747972,103741508,79433746,99106873,133574332,
                            126202977,131373565,24488204,25405510,60760728,34678317,41133189,98830142,54808415,46619702,1676341,59258550)
yri.df <- read.table('table.YRI.txt',header=T)
nrow(yri.df)#3553
yri.df$DNMSNV <- vector(length=nrow(yri.df))
for (pos in yri.df$POS){
  if (pos %in% yri.dnm.pos){
    vs <- which(pos==yri.df$POS)
    yri.df$DNMSNV[vs]<-'DNM'
  }else{
    ns <- which(pos==yri.df$POS)
    yri.df$DNMSNV[ns]<-'.'
  }
}

colnames(yri.df)
yri.df <- yri.df[,c(1:17,19:23,25:29,31:33)]
colnames(yri.df)[15:29]<-c("mother.AD","mother.AF","mother.DP","mother.PL","mother.PP",
                           "father.AD","father.AF","father.DP","father.PL","father.PP",
                           "offspring.AD","offspring.AF","offspring.DP","offspring.PL","offspring.PP")
yri.df.SNP <- yri.df[yri.df$TRANSITION!=-1,]
nrow(yri.df.SNP)#2645
yri.df.SNP2 <- transform(yri.df.SNP,mother.PL=colsplit(mother.PL,pattern='\\,',names=c(1,2,3)))
yri.df.SNP2 <- transform(yri.df.SNP2,mother.PP=colsplit(mother.PP,pattern='\\,',names=c(1,2,3)))
yri.df.SNP2 <- transform(yri.df.SNP2,father.PL=colsplit(father.PL,pattern='\\,',names=c(1,2,3)))
yri.df.SNP2 <- transform(yri.df.SNP2,father.PP=colsplit(father.PP,pattern='\\,',names=c(1,2,3)))
yri.df.SNP2 <- transform(yri.df.SNP2,offspring.PL=colsplit(offspring.PL,pattern='\\,',names=c(1,2,3)))
yri.df.SNP2 <- transform(yri.df.SNP2,offspring.PP=colsplit(offspring.PP,pattern='\\,',names=c(1,2,3)))
yri.df.SNP2 <- transform(yri.df.SNP2,offspring.AD=colsplit(offspring.AD,pattern='\\,',names=c(1,2)))
yri.df.SNP2 <- transform(yri.df.SNP2,father.AD=colsplit(father.AD,pattern='\\,',names=c(1,2)))
yri.df.SNP2 <- transform(yri.df.SNP2,mother.AD=colsplit(mother.AD,pattern='\\,',names=c(1,2)))
yri.df.SNP2 <- transform(yri.df.SNP2,offspring.AD=colsplit(offspring.AD,pattern='\\,',names=c(1,2)))

yri.df.SNP3 <- yri.df.SNP2[,c(1:14,16:19,21:24,26:30)]
colnames(yri.df.SNP3)
yri.df.SNP3$father.PP<-yri.df.SNP2$POS
yri.df.SNP3$father.PL<-yri.df.SNP2$POS
yri.df.SNP3$mother.PP<-yri.df.SNP2$POS
yri.df.SNP3$mother.PL<-yri.df.SNP2$POS
yri.df.SNP3$offspring.PP<-yri.df.SNP2$POS
yri.df.SNP3$offspring.PL<-yri.df.SNP2$POS

for (i in 1:nrow(yri.df.SNP2)){
  if (yri.df.SNP2$father.PP[2][i,]>=yri.df.SNP2$father.PP[3][i,]){
    yri.df.SNP3$father.PP[i]<-yri.df.SNP2$father.PP[3][i,]
  }else{
    yri.df.SNP3$father.PP[i]<-yri.df.SNP2$father.PP[2][i,]
  }
}
for (i in 1:nrow(yri.df.SNP2)){
  if (yri.df.SNP2$father.PL[2][i,]>=yri.df.SNP2$father.PL[3][i,]){
    yri.df.SNP3$father.PL[i]<-yri.df.SNP2$father.PL[3][i,]
  }else{
    yri.df.SNP3$father.PL[i]<-yri.df.SNP2$father.PL[2][i,]
  }
}

for (i in 1:nrow(yri.df.SNP2)){
  if (yri.df.SNP2$mother.PP[2][i,]>=yri.df.SNP2$mother.PP[3][i,]){
    yri.df.SNP3$mother.PP[i]<-yri.df.SNP2$mother.PP[3][i,]
  }else{
    yri.df.SNP3$mother.PP[i]<-yri.df.SNP2$mother.PP[2][i,]
  }
}

for (i in 1:nrow(yri.df.SNP2)){
  if (yri.df.SNP2$mother.PL[2][i,]>=yri.df.SNP2$mother.PL[3][i,]){
    yri.df.SNP3$mother.PL[i]<-yri.df.SNP2$mother.PL[3][i,]
  }else{
    yri.df.SNP3$mother.PL[i]<-yri.df.SNP2$mother.PL[2][i,]
  }
}

for (i in 1:nrow(yri.df.SNP2)){
  if (yri.df.SNP2$offspring.PP[1][i,]>=yri.df.SNP2$offspring.PP[3][i,]){
    yri.df.SNP3$offspring.PP[i]<-yri.df.SNP2$offspring.PP[3][i,]
  }else{
    yri.df.SNP3$offspring.PP[i]<-yri.df.SNP2$offspring.PP[1][i,]
  }
}

for (i in 1:nrow(yri.df.SNP2)){
  if (yri.df.SNP2$offspring.PL[1][i,]>=yri.df.SNP2$offspring.PL[3][i,]){
    yri.df.SNP3$offspring.PL[i]<-yri.df.SNP2$offspring.PL[3][i,]
  }else{
    yri.df.SNP3$offspring.PL[i]<-yri.df.SNP2$offspring.PL[1][i,]
  }
}
#AD
for (i in 1:nrow(yri.df.SNP2)){
  yri.df.SNP3$mother.REF[i] <- yri.df.SNP2$mother.AD[i,1]
  yri.df.SNP3$mother.ALT[i] <- yri.df.SNP2$mother.AD[i,2]
  #select.asj$mother.AF[i] <- round(split.asj$mother.AD[i,2] / split.asj$mother.DP[i],digits = 2)
  yri.df.SNP3$father.REF[i] <- yri.df.SNP2$father.AD[i,1]
  yri.df.SNP3$father.ALT[i] <- yri.df.SNP2$father.AD[i,2]
  #select.asj$father.AF[i] <- round(split.asj$father.AD[i,2] / split.asj$father.DP[i], digits=2)
  yri.df.SNP3$offspring.REF[i] <- yri.df.SNP2$offspring.AD[i,1]
  yri.df.SNP3$offspring.ALT[i] <- yri.df.SNP2$offspring.AD[i,2]
  #select.asj$offspring.AF[i] <- round(split.asj$offspring.AD[i,2] / split.asj$offspring.DP[i], digits=2)
}
head(yri.df.SNP3)

###concatenate two trio data###
colnames(yri.df.SNP3)
colnames(compDfSnv3)
ceu.only <- compDfSnv3.2
yri.only <- yri.df.SNP3
ceu.only$DNMSNV<-as.character(ceu.only$DNMSNV)
yri.only$DNMSNV<-as.character(yri.only$DNMSNV)
for (i in c(which(ceu.only$DNMSNV=='DNM'))){
  ceu.only$DNMSNV[i]<-'CEU.DNM'
}
for (i in c(which(yri.only$DNMSNV=='DNM'))){
  yri.only$DNMSNV[i]<-'YRI.DNM'
}

ceu.and.yri.whole <- rbind(ceu.only,yri.only)
ceu.and.yri.dnm <- ceu.and.yri.whole[ceu.and.yri.whole$DNMSNV!='.',]
ceu.and.yri.no <- ceu.and.yri.whole[ceu.and.yri.whole$DNMSNV=='.',]

compHcDfSnv3$n.offspring.PL<-compHcDfSnv3$offspring.PL/max(compHcDfSnv3$offspring.PL)
yri.only$n.offspring.PL<-yri.only$offspring.PL/max(yri.only$offspring.PL)
ceu.only$n.offspring.PL<-ceu.only$offspring.PL/max(ceu.only$offspring.PL)
ceu.only1<-ceu.only[ceu.only$n.offspring.PL<=0.75,]

compHcDfSnv3.2 <- compHcDfSnv3
compHcDfSnv3.2$DNMSNV<-'.'
p1 <- ggplot(ceu.only1,aes(x=n.offspring.PL,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
p2 <- ggplot(yri.only,aes(x=n.offspring.PL,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
p3 <- ggplot(compHcDfSnv3.2,aes(x=n.offspring.PL,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
grid.arrange(p1,p2,p3,nrow=3)

ggplot(ceu.and.yri.dnm,aes(x=QD,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.no,aes(x=QD,y=offspring.AF,color=DNMSNV))+geom_point(size=1)

ggplot(ceu.and.yri.whole,aes(x=BaseQRankSum,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.dnm,aes(x=BaseQRankSum,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.no[ceu.and.yri.no$BaseQRankSum<5&ceu.and.yri.no$BaseQRankSum>-4&ceu.and.yri.no$offspring.AF>=0.3&ceu.and.yri.no$offspring.AF<=0.7,],aes(x=BaseQRankSum,y=offspring.AF,color=DNMSNV))+geom_point(size=1)

ggplot(ceu.and.yri.dnm,aes(x=MQ,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.no,aes(x=MQ,y=offspring.AF,color=DNMSNV))+geom_point(size=1)

ggplot(ceu.and.yri.whole,aes(x=MQRankSum,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.dnm,aes(x=MQRankSum,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.no[ceu.and.yri.no$MQRankSum>=-3&ceu.and.yri.no$MQRankSum<=2&ceu.and.yri.no$offspring.AF>=0.3&ceu.and.yri.no$offspring.AF<=0.7,],aes(x=MQRankSum,y=offspring.AF,color=DNMSNV))+geom_point(size=1)

ggplot(ceu.and.yri.whole,aes(x=VQSLOD,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.dnm,aes(x=VQSLOD,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.no,aes(x=VQSLOD,y=offspring.AF,color=DNMSNV))+geom_point(size=1)

ggplot(ceu.and.yri.whole,aes(x=offspring.PL,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.dnm,aes(x=offspring.PL,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.no,aes(x=offspring.PL,y=offspring.AF,color=DNMSNV))+geom_point(size=1)

ggplot(ceu.and.yri.whole,aes(x=offspring.PP,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.dnm,aes(x=offspring.PP,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.no,aes(x=offspring.PP,y=offspring.AF,color=DNMSNV))+geom_point(size=1)

ggplot(ceu.and.yri.whole,aes(x=father.PL,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.dnm,aes(x=father.PL,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.no,aes(x=father.PL,y=offspring.AF,color=DNMSNV))+geom_point(size=1)

ggplot(ceu.and.yri.whole,aes(x=father.PP,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.dnm,aes(x=father.PP,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.no,aes(x=father.PP,y=offspring.AF,color=DNMSNV))+geom_point(size=1)

ggplot(ceu.and.yri.whole,aes(x=mother.PL,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.dnm,aes(x=mother.PL,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.no,aes(x=mother.PL,y=offspring.AF,color=DNMSNV))+geom_point(size=1)

ggplot(ceu.and.yri.whole,aes(x=mother.PP,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.dnm,aes(x=mother.PP,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.no[ceu.and.yri.no$mother.PP<=290&ceu.and.yri.no$offspring.AF>=0.3&ceu.and.yri.no$offspring.AF<=0.7,],aes(x=mother.PP,y=offspring.AF,color=DNMSNV))+geom_point(size=1)

ggplot(ceu.and.yri.whole,aes(x=offspring.PL,y=father.PL,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.dnm,aes(x=FS,y=offspring.PP,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.whole,aes(x=MQ,y=offspring.PL,color=DNMSNV))+geom_point(size=1)

ggplot(ceu.and.yri.whole,aes(x=mother.DP,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.dnm,aes(x=mother.DP,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
ggplot(ceu.and.yri.no,aes(x=mother.PL,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
