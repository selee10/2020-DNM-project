library(ggplot2)
TP <- c(55,26,113,53,96,73)
FP <- c(36,37,107,32,73,38)
real <- c(4,4,9,7,3,11,3)#400:c(40,40,96,68,24,108,24)
aftp <- c(3,2,7,7,6,12,2)#39
rnames <- c('C>A','C>G','C>T','CpG>TpG','T>A','T>C','T>G')
dframe <- data.frame(real,aftp,row.names = rnames)
chisq.test(dframe)
a <- ggplot(dframe,aes(x=rnames))+
  geom_bar(aes(y=real),stat='identity',position='dodge')
a+geom_bar(aes(y=aftp),stat='identity',
           fill='yellow',alpha=.5,
           position=position_dodge(10))

#C>A
v1 <- c('ACA','ACC','ACG','ACT','CCA','CCC','CCG','CCT','GCA','GCC','GCG','GCT','TCA','TCC','TCG','TCT')
N1 <- c(0,0,0,0,1/39*100,0,0,0,0,1/39*100,0,0,0,0,0,1/39*100)
D1 <- data.frame(N1,row.names = v1)
ggplot(D1,aes(x=v1,y=N1))+geom_bar(stat='identity',fill='blue',alpha=.5)
#C>G
v2<-c('ACA1','ACC1','ACG1','ACT1','CCA1','CCC1','CCG1','CCT1','GCA1','GCC1','GCG1','GCT1','TCA1','TCC1',
      'TCG1','TCT1')
N2 <- c(0,0,0,1/39*100,0,0,0,0,0,0,0,0,0,1/39*100,0,0) 
D2 <- data.frame(N2,row.names = v2)
ggplot(D2,aes(x=v2,y=N2))+geom_bar(stat='identity')
#C>T
v3<-c('ACA2','ACC2','ACG2','ACT2','CCA2','CCC2','CCG2','CCT2','GCA2','GCC2','GCG2','GCT2',
      'TCA2','TCC2','TCG2','TCT2')
N3<-c(0,0,4,1,0,0,0,3,0,0,1,0,0,1,2,1)/39*100
D3 <- data.frame(N3,row.names = v3)
ggplot(D3,aes(x=v3,y=N3))+geom_bar(stat='identity',fill='red',alpha=.5)
#T>A
v4 <- c('ATA','ATC','ATG','ATT','CTA','CTC','CTG','CTT','GTA','GTC','GTG','GTT','TTA','TTC','TTG','TTT')
N4 <-c(0,0,0,2,0,0,0,1,0,1,0,0,1,0,0,1)/39*100
#T>C
v5 <- c('ATA1','ATC1','ATG1','ATT1','CTA1','CTC1','CTG1','CTT1','GTA1','GTC1',
        'GTG1','GTT1','TTA1','TTC1','TTG1','TTT1')
N5 <- c(2,1,0,2,0,2,1,0,1,0,1,0,0,0,1,1)/39*100
#T>G
v6 <- c('ATA2','ATC2','ATG2','ATT2','CTA2','CTC2','CTG2','CTT2','GTA2',
        'GTC2','GTG2','GTT2','TTA2','TTC2','TTG2','TTT2')
N6 <- c(0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0)/39*100

name <- c(v1,v2,v3,v4,v5,v6)
num <- c(N1,N2,N3,N4,N5,N6)
tdframe <- data.frame(num,row.names = name)
name1 <- seq(from=1,to=96)
ggplot(tdframe,aes(x=name1,y=num))+geom_bar(stat='identity')
