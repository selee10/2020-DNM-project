setwd("C:/Users/selee/Desktop/class/defaultCGP")
library("readxl")
#####CEU preprocessing#####
ceu.v <- read_excel("CEU_and_YRI_variants.xls",sheet='CEU',col_names=TRUE)
colnames(ceu.v)
ceu.va <- ceu.v[,c(2:5,50)]
#tibble to dataframe
ceu.var <- as.data.frame(ceu.va)
nrow(ceu.var)

#####YRI preprocessing#####
yri.v <- read_excel("CEU_and_YRI_variants.xls",sheet='YRI',col_names=TRUE)
colnames(yri.v)
yri.va <- yri.v[,c(2:5,61)]
yri.var <- as.data.frame(yri.va)
head(yri.var)

##need to liftover hg18 to hg19...
colnames(ceu.var)
nrow(ceu.var)
for (i in 1:3236){
  line <- paste('chr',as.character(ceu.var$chr[i]),':',as.character(ceu.var$pos[i]),
                '-',as.character(ceu.var$pos[i]),sep='')
  cat(line,file='C:/Users/selee/Desktop/test.ceu.hg18.txt',sep='\n',append=T)
}
for (i in 1:2750){
  line <- paste('chr',as.character(yri.var$chr[i]),':',as.character(yri.var$pos[i]),
                '-',as.character(yri.var$pos[i]),sep='')
  cat(line,file='C:/Users/selee/Desktop/yri.hg18.txt',sep='\n',append=T)
}

for (i in 1:2750){
  line1 <- paste('chr',as.character(yri.var$chr[i]),':',as.character(yri.var$pos[i]),
                '-',as.character(yri.var$pos[i]),sep='')
  cat(line1,file='C:/Users/selee/Desktop/test.yri.hg18.txt',sep='\n',append=T)
}
c <- read.table('CEU.hg19.txt',sep=':')
y <- read.table('YRI.hg19.txt',sep=':')
head(c)
colnames(c) <- c('chr','pos')
colnames(y) <- c('chr','pos')
split.c <- data.frame(do.call('rbind',strsplit(as.character(c$pos),'-',fixed=TRUE)))
split.y <- data.frame(do.call('rbind',strsplit(as.character(y$pos),'-',fixed=TRUE)))
c.1 <- within(c,pos <- as.vector(split.c[1]))
y.1 <- within(y,pos<-as.vector(split.y[1]))

ceu.var[,c(1,2)]<-c.1
head(ceu.var)
yri.var[,c(1,2)]<-y.1
colnames(ceu.var)[2] <- 'pos'
colnames(yri.var)[2] <- 'pos'

ceu.var$chr <- gsub('chr([0-9A-Z]+)','\\1',ceu.var$chr)
yri.var$chr <- gsub('chr([0-9A-Z]+)','\\1',yri.var$chr)

head(ceu.var)
tail(ceu.var)

ceu.var$pos <- as.integer(as.vector(unlist(ceu.var$pos)))
yri.var$pos <- as.integer(as.vector(unlist(yri.var$pos)))

#####^^^^^^^- liftover done(hg18 to hg19) -^^^^^^^######

#####vvvvvvvv- classify by variant type -vvvvvv#####
###merged_vald.code
## 0 - no call, data not informative enough
## 1 - germline de novo
## 2 - somatic/cell line de novo
## 3 - inherited polymorphism
## 4 - false positive call, no variation in any sample
## 5 - cell line mutation/revertant in parents (only in YRI)

ceu.nocall <- ceu.var[ceu.var$merged_vald.code==0,]
ceu.dnm <- ceu.var[ceu.var$merged_vald.code==1,]
ceu.soma <- ceu.var[ceu.var$merged_vald.code==2,]
ceu.inherited <- ceu.var[ceu.var$merged_vald.code==3,]
ceu.fp <- ceu.var[ceu.var$merged_vald.code==4,]

test.ceu$type <- vector(length=nrow(test.ceu))

for (i in 1:nrow(ceu.dnm)){
  for (j in 1:nrow(test.ceu)){
    if (ceu.dnm$chr[i]==test.ceu$CHROM[j]){
      if (ceu.dnm$pos[i]==test.ceu$POS[j]){
        test.ceu$type[j] <- 'DNM'
      }
    }
  }
}
table(test.ceu$type)

for (i in 1:nrow(ceu.soma)){
  for (j in 1:nrow(test.ceu)){
    if (ceu.soma$chr[i]==test.ceu$CHROM[j]){
      if (ceu.soma$pos[i]==test.ceu$POS[j]){
        test.ceu$type[j] <- 'somatic'
      }
    }
  }
}
table(test.ceu$type)

for (i in 1:nrow(ceu.inherited)){
  for (j in 1:nrow(test.ceu)){
    if (ceu.inherited$chr[i]==test.ceu$CHROM[j]){
      if (ceu.inherited$pos[i]==test.ceu$POS[j]){
        test.ceu$type[j] <- 'germline'
      }
    }
  }
}
table(test.ceu$type)

for (i in 1:nrow(ceu.fp)){
  for (j in 1:nrow(test.ceu)){
    if (ceu.fp$chr[i]==test.ceu$CHROM[j]){
      if (ceu.fp$pos[i]==test.ceu$POS[j]){
        test.ceu$type[j] <- 'FP'
      }
    }
  }
}
table(test.ceu$type)

for (i in 1:nrow(ceu.nocall)){
  for (j in 1:nrow(test.ceu)){
    if (ceu.nocall$chr[i]==test.ceu$CHROM[j]){
      if (ceu.nocall$pos[i]==test.ceu$POS[j]){
        test.ceu$type[j] <- 'no call'
      }
    }
  }
}
#########yri - classify by variant type######
yri.nocall <- yri.var[yri.var$merged_vald.code==0,]
yri.dnm <- yri.var[yri.var$merged_vald.code==1,]
yri.soma <- yri.var[yri.var$merged_vald.code==2,]
yri.inherited <- yri.var[yri.var$merged_vald.code==3,]
yri.fp <- yri.var[yri.var$merged_vald.code==4,]
yri.type5 <- yri.var[yri.var$merged_vald.code==5,]
nrow(yri.type5)

test.yri$type <- vector(length=nrow(test.yri))

for (i in 1:nrow(yri.dnm)){
  for (j in 1:nrow(test.yri)){
    if (yri.dnm$chr[i]==test.yri$CHROM[j]){
      if (yri.dnm$pos[i]==test.yri$POS[j]){
        test.yri$type[j] <- 'DNM'
      }
    }
  }
}

for (i in 1:nrow(yri.soma)){
  for (j in 1:nrow(test.yri)){
    if (yri.soma$chr[i]==test.yri$CHROM[j]){
      if (yri.soma$pos[i]==test.yri$POS[j]){
        test.yri$type[j] <- 'somatic'
      }
    }
  }
}

for (i in 1:nrow(yri.inherited)){
  for (j in 1:nrow(test.yri)){
    if (yri.inherited$chr[i]==test.yri$CHROM[j]){
      if (yri.inherited$pos[i]==test.yri$POS[j]){
        test.yri$type[j] <- 'germline'
      }
    }
  }
}

for (i in 1:nrow(yri.fp)){
  for (j in 1:nrow(test.yri)){
    if (yri.fp$chr[i]==test.yri$CHROM[j]){
      if (yri.fp$pos[i]==test.yri$POS[j]){
        test.yri$type[j] <- 'FP'
      }
    }
  }
}
for (i in 1:nrow(yri.type5)){
  for (j in 1:nrow(test.yri)){
    if (yri.type5$chr[i]==test.yri$CHROM[j]){
      if (yri.type5$pos[i]==test.yri$POS[j]){
        test.yri$type[j] <- 'cell line mut/rev in parents'
      }
    }
  }
}
for (i in 1:nrow(yri.nocall)){
  for (j in 1:nrow(test.yri)){
    if (yri.nocall$chr[i]==test.yri$CHROM[j]){
      if (yri.nocall$pos[i]==test.yri$POS[j]){
        test.yri$type[j] <- 'no call'
      }
    }
  }
}
table(test.yri$type)
table(test.ceu$type)
#####plot#####
library(ggplot2)
fp.test.yri <- test.yri[test.yri$type=='FP',]
soma.test.yri <- test.yri[test.yri$type=='somatic',]
dnm.test.yri <- test.yri[test.yri$type=='DNM',]
ggplot(soma.test.yri,
       aes(x=QD,y=offspring.AF,color=type))+geom_point()
ggplot(soma.test.yri,
       aes(x=offspring.PL,y=offspring.AF,color=type))+geom_point(size=.5)
win.graph(); ggplot(test.yri,
       aes(x=prob[,2],y=offspring.AF,color=type))+geom_point()
hist(soma.test.yri$prob[,2])
hist(dnm.test.yri$prob[,2])

win.graph(); ggplot(test.ceu,
                    aes(x=prob[,2],y=offspring.AF,color=type))+geom_point()

win.graph(); ggplot(test.hc,
                    aes(x=prob[,2],y=offspring.AF))+geom_point()
win.graph(); ggplot(whole.z.asj,
                    aes(x=prob[,2],y=offspring.AF))+geom_point()

soma.test.ceu <- test.ceu[test.ceu$type=='somatic',]
win.graph(); ggplot(soma.test.ceu,
                    aes(x=prob[,2],y=offspring.AF,color=type))+geom_point()
table(test.ceu$type)
table(test.yri$type)
