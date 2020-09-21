test.yri$predDNM <- predict(model3,yri.z.inlier)
test.yri$prob <- predict(model3,yri.z.inlier,type='prob')

yri.pr <- prediction(test.yri$prob,test.yri$DNMSNV)
yri.prf <- performance(yri.pr, measure='ppv',x.measure = 'tpr')
win.graph; plot(yri.prf, main='PR curve of YRI trio')

#which cutoff RF used to determine a site as dnm
range(test.yri[test.yri$predDNM!='.',]$prob[,2])
range(test.ceu[test.ceu$predDNM!='.',]$prob[,2])
colnames(test.yri)
test.yri.whprob <- test.yri[,1:34]
test.yri.whprob$prob <- vector(length=nrow(test.yri))
test.yri.whprob$prob <- test.yri$prob[,2]
win.graph(); ggplot(test.yri.whprob[test.yri.whprob$prob>=0.4,],
       aes(x=prob,fill=DNMSNV))+geom_histogram()
##rf parameter tuning##
data.in.use <- rawdata[,c(3:12,14,15,17,18,20:25)]
x <- rawdata[,c(3:12,14,15,17,18,21:25)]
y <- rawdata[,20]
#trainClasses <- rawdata[,20]
seed<-1000
metric <- 'Accuracy'
set.seed(seed)
mtry<-sqrt(ncol(x))
control <- trainControl(method='boot',search='random')
set.seed(seed)
#rf_random <- train(DNMSNV~.,data=x,method='rf',tunelength=15,
#                   metric=metric,trControl=control)
#######install.packages('e1071',dependencies=T)#####
library(e1071)
library(caret)
rf_random_2 <- train(DNMSNV~.,data=x,method='rf',
                   metric=metric,trControl=control)
print(rf_random_2)
win.graph(); plot(rf_random_2)
#####OOB error#####
set.seed(seed)
bestmtry <- tuneRF(x,y,stepFactor = 1.5,improve = 1e-5,ntree=500)
#mtry = 4  OOB error = 1.62% 
#Searching left ...
#mtry = 3 	OOB error = 1.85% 
#-0.1410256 1e-05 
#Searching right ...
#mtry = 6 	OOB error = 1.83% 
#-0.1282051 1e-05
win.graph(); plot(bestmtry,type='b')

bestmtry2 <- tuneRF(x,y,stepFactor = 1.5,improve = 0.01,ntree=500)
#mtry = 4  OOB error = 1.67% 
#Searching left ...
#mtry = 3 	OOB error = 1.73% 
#-0.0375 0.01 
#Searching right ...
#mtry = 6 	OOB error = 1.9% 
#-0.1375 0.01

####ntree tuning####
customRF <- list(type="Classification",library="randomForest",loop=NULL)
customRF$parameters <- data.frame(parameter=c("mtry","ntree"),class=rep("numeric",2),
                                  label=c("mtry","ntree"))
customRF$grid <- function(x,y,len=NULL,search='grid'){}
customRF$fit <- function(x,y,wts,param,lev,last,weights,classProbs, ...){
  randomForest(x,y,mtry = param$mtry,ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit,newdata,preProc=NULL,submodels=NULL)
predict(modelFit,newdata)  
customRF$prob <- function(modelFit,newdata,preProc=NULL,submodels=NULL)
predict(modelFit,newdata,type='prob')  
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes
control <- trainControl(method='boot',search='grid')
tunegrid <- expand.grid(.mtry=c(1:7),.ntree=c(100,300,500,700))
set.seed(seed)
custom.re <- train(DNMSNV~.,data=data.in.use,method=customRF,metric=metric,
                tuneGrid=tunegrid, trControl = control)
win.graph(); plot(custom.re)
