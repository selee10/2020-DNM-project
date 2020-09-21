v.auc <- c()
for (nodesize in 1:10){
  t.model <- randomForest(DNMSNV~QUAL+BaseQRankSum+FS+QD+MQ+
                           MQRankSum+ReadPosRankSum+VQSLOD+offspring.AF+offspring.DP+
                           father.DP+mother.DP+offspring.PL+offspring.PP+
                           father.PL+father.PP+mother.PL+mother.PP,data=rawdata.2,
                         mtry=4,ntree=500,nodesize=nodesize,
                         importance=T,norm.votes=T,proximity=T)
  prob <- predict(t.model,testdata,type='prob')
  pr <- prediction(prob[,2],testdata$DNMSNV)
  auc.tmp <- performance(pr,"auc"); auc <- as.numeric(auc.tmp@y.values)
  v.auc <- c(v.auc,auc)
}
nodesize <- c(1:10)
AUC <- v.auc
win.graph(); plot(x=nodesize,y=AUC,type='b')

c.auc <- c()
for (cutoff in seq(0.5,1.0,0.01)){
  t.model <- randomForest(DNMSNV~QUAL+BaseQRankSum+FS+QD+MQ+
                            MQRankSum+ReadPosRankSum+VQSLOD+offspring.AF+offspring.DP+
                            father.DP+mother.DP+offspring.PL+offspring.PP+
                            father.PL+father.PP+mother.PL+mother.PP,data=rawdata.2,
                          mtry=4,ntree=500,nodesize=1,cutoff=c(cutoff,1-cutoff),
                          importance=T,norm.votes=T,proximity=T)
  prob <- predict(t.model,testdata,type='prob')
  pr <- prediction(prob[,2],testdata$DNMSNV)
  auc.tmp <- performance(pr,"auc"); auc <- as.numeric(auc.tmp@y.values)
  c.auc <- c(c.auc,auc)
}
AUC <- c.auc
cutoff <- seq(0.5,0.99,0.01)
plot(x=cutoff,y=AUC,type='b')
