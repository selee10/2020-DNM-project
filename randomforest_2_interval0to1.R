library(scales)
#turn every feature value in range of [0,1]
#CEU trio: ceu.only
#YRI trio: yri.only
#Han Chinese trio: compHcDfSnv3.2
ceu.only <- ceu.only[complete.cases(ceu.only),]
yri.only <- yri.only[complete.cases(yri.only),]
#yri: 2645 -> 2599
compHcDfSnv3.2 <- compHcDfSnv3.2[complete.cases(compHcDfSnv3.2),]

ceu.only.new <- ceu.only[,1:33]
yri.only.new <- yri.only[,1:33]
hc.only.new <- compHcDfSnv3.2

#values to adjust:QUAL,BaseQRankSum,FS,MQ,MQRankSum,
#QD,ReadPosRankSum,VQSLOD,father.DP,offspring.DP,
#mother.DP,offspring.PL,offspring.PP,father.PL,
#father.PP,mother.PL,mother.PP
#####ceu#####
ceu.only.new$QUAL<-ceu.only$QUAL/max(ceu.only$QUAL)
ceu.only.new$BaseQRankSum<-
  (ceu.only$BaseQRankSum-min(ceu.only$BaseQRankSum))/(max(ceu.only$BaseQRankSum)-min(ceu.only$BaseQRankSum))
ceu.only.new$FS<-ceu.only$FS/max(ceu.only$FS)
ceu.only.new$MQ<-ceu.only$MQ/max(ceu.only$MQ)
ceu.only.new$MQRankSum<-
  (ceu.only$MQRankSum-min(ceu.only$MQRankSum))/(max(ceu.only$MQRankSum)-min(ceu.only$MQRankSum))
ceu.only.new$QD<-ceu.only$QD/max(ceu.only$QD)
ceu.only.new$ReadPosRankSum<-
  (ceu.only$ReadPosRankSum-min(ceu.only$ReadPosRankSum))/(max(ceu.only$ReadPosRankSum)-min(ceu.only$ReadPosRankSum))
ceu.only.new$VQSLOD<-
  (ceu.only$VQSLOD-min(ceu.only$VQSLOD))/(max(ceu.only$VQSLOD)-min(ceu.only$VQSLOD))
ceu.only.new$father.DP<-ceu.only$father.DP/max(ceu.only$father.DP)
ceu.only.new$mother.DP<-ceu.only$mother.DP/max(ceu.only$mother.DP)
ceu.only.new$offspring.DP<-ceu.only$offspring.DP/max(ceu.only$offspring.DP)
ceu.only.new$father.PL<-ceu.only$father.PL/max(ceu.only$father.PL)
ceu.only.new$father.PP<-ceu.only$father.PP/max(ceu.only$father.PP)
ceu.only.new$mother.PL<-ceu.only$mother.PL/max(ceu.only$mother.PL)
ceu.only.new$mother.PP<-ceu.only$mother.PP/max(ceu.only$mother.PP)
ceu.only.new$offspring.PL<-ceu.only$offspring.PL/max(ceu.only$offspring.PL)
ceu.only.new$offspring.PP<-ceu.only$offspring.PP/max(ceu.only$offspring.PP)

####yri####
yri.only.new$QUAL<-yri.only$QUAL/max(yri.only$QUAL)
yri.only.new$BaseQRankSum<-
  (yri.only$BaseQRankSum-min(yri.only$BaseQRankSum))/(max(yri.only$BaseQRankSum)-min(yri.only$BaseQRankSum))
yri.only.new$FS<-yri.only$FS/max(yri.only$FS)
yri.only.new$MQ<-yri.only$MQ/max(yri.only$MQ)
yri.only.new$MQRankSum<-
  (yri.only$MQRankSum-min(yri.only$MQRankSum))/(max(yri.only$MQRankSum)-min(yri.only$MQRankSum))
yri.only.new$QD<-yri.only$QD/max(yri.only$QD)
yri.only.new$ReadPosRankSum<-
  (yri.only$ReadPosRankSum-min(yri.only$ReadPosRankSum))/(max(yri.only$ReadPosRankSum)-min(yri.only$ReadPosRankSum))
yri.only.new$VQSLOD<-
  (yri.only$VQSLOD-min(yri.only$VQSLOD))/(max(yri.only$VQSLOD)-min(yri.only$VQSLOD))
yri.only.new$father.DP<-yri.only$father.DP/max(yri.only$father.DP)
yri.only.new$mother.DP<-yri.only$mother.DP/max(yri.only$mother.DP)
yri.only.new$offspring.DP<-yri.only$offspring.DP/max(yri.only$offspring.DP)
yri.only.new$father.PL<-yri.only$father.PL/max(yri.only$father.PL)
yri.only.new$father.PP<-yri.only$father.PP/max(yri.only$father.PP)
yri.only.new$mother.PL<-yri.only$mother.PL/max(yri.only$mother.PL)
yri.only.new$mother.PP<-yri.only$mother.PP/max(yri.only$mother.PP)
yri.only.new$offspring.PL<-yri.only$offspring.PL/max(yri.only$offspring.PL)
yri.only.new$offspring.PP<-yri.only$offspring.PP/max(yri.only$offspring.PP)

#####han chinese#####
hc.only.new$QUAL<-compHcDfSnv3.2$QUAL/max(compHcDfSnv3.2$QUAL)
hc.only.new$BaseQRankSum<-
  (compHcDfSnv3.2$BaseQRankSum-min(compHcDfSnv3.2$BaseQRankSum))/(max(compHcDfSnv3.2$BaseQRankSum)-min(compHcDfSnv3.2$BaseQRankSum))
hc.only.new$FS<-compHcDfSnv3.2$FS/max(compHcDfSnv3.2$FS)
hc.only.new$MQ<-compHcDfSnv3.2$MQ/max(compHcDfSnv3.2$MQ)
hc.only.new$MQRankSum<-
  (compHcDfSnv3.2$MQRankSum-min(compHcDfSnv3.2$MQRankSum))/(max(compHcDfSnv3.2$MQRankSum)-min(compHcDfSnv3.2$MQRankSum))
hc.only.new$QD<-compHcDfSnv3.2$QD/max(compHcDfSnv3.2$QD)
hc.only.new$ReadPosRankSum<-
  (compHcDfSnv3.2$ReadPosRankSum-min(compHcDfSnv3.2$ReadPosRankSum))/(max(compHcDfSnv3.2$ReadPosRankSum)-min(compHcDfSnv3.2$ReadPosRankSum))
hc.only.new$VQSLOD<-
  (compHcDfSnv3.2$VQSLOD-min(compHcDfSnv3.2$VQSLOD))/(max(compHcDfSnv3.2$VQSLOD)-min(compHcDfSnv3.2$VQSLOD))
hc.only.new$father.DP<-compHcDfSnv3.2$father.DP/max(compHcDfSnv3.2$father.DP)
hc.only.new$mother.DP<-compHcDfSnv3.2$mother.DP/max(compHcDfSnv3.2$mother.DP)
hc.only.new$offspring.DP<-compHcDfSnv3.2$offspring.DP/max(compHcDfSnv3.2$offspring.DP)
hc.only.new$father.PL<-compHcDfSnv3.2$father.PL/max(compHcDfSnv3.2$father.PL)
hc.only.new$father.PP<-compHcDfSnv3.2$father.PP/max(compHcDfSnv3.2$father.PP)
hc.only.new$mother.PL<-compHcDfSnv3.2$mother.PL/max(compHcDfSnv3.2$mother.PL)
hc.only.new$mother.PP<-compHcDfSnv3.2$mother.PP/max(compHcDfSnv3.2$mother.PP)
hc.only.new$offspring.PL<-compHcDfSnv3.2$offspring.PL/max(compHcDfSnv3.2$offspring.PL)
hc.only.new$offspring.PP<-compHcDfSnv3.2$offspring.PP/max(compHcDfSnv3.2$offspring.PP)


p.ceu <- ggplot(ceu.only.new,aes(x=QUAL,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
p.yri <- ggplot(yri.only.new,aes(x=QUAL,y=offspring.AF,color=DNMSNV))+geom_point(size=1)
p.hc <- ggplot(hc.only.new,aes(x=QUAL,y=offspring.AF))+geom_point(size=1)+coord_cartesian(xlim=c(0:1),ylim=c(0,1))
grid.arrange(p.ceu,p.yri,p.hc,nrow=3)

ceu.only.dnm <- ceu.only.new[ceu.only.new$DNMSNV!='.',]
yri.only.dnm <- yri.only.new[yri.only.new$DNMSNV!='.',]
p.ceu.dnm <- ggplot(ceu.only.dnm,aes(x=QUAL,y=offspring.AF,color=DNMSNV))+geom_point(size=1)+coord_cartesian(xlim=c(0:1),ylim=c(0,1))
p.yri.dnm <- ggplot(yri.only.dnm,aes(x=QUAL,y=offspring.AF,color=DNMSNV))+geom_point(size=1)+coord_cartesian(xlim=c(0:1),ylim=c(0:1))
grid.arrange(p.ceu.dnm,p.yri.dnm,p.hc)
