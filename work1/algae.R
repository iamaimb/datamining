###################################################
### Loading the Data into R
###################################################
library(DMwR)
setwd("D:\\R\\datamining1\\")
algae <- read.table('Analysis.txt',
          header=F,
          dec='.',
          col.names=c('season','size','speed','mxPH','mnO2','Cl',
          'NO3','NH4','oPO4','PO4','Chla','a1','a2','a3','a4',
          'a5','a6','a7'),
          na.strings=c('XXXXXXX'))

colum <- list('season','size','speed','mxPH','mnO2','Cl',
            'NO3','NH4','oPO4','PO4','Chla','a1','a2','a3','a4',
            'a5','a6','a7')
###################################################
### Data Visualization and Summarization
###################################################
sink("Summary.txt")   #把结果保存到summary.txt
summary(algae)
sink()

library(car)
for (i in 4:18)
{
a <- colum[[i]]

jpeg(file=paste('Histogram of ',a,'.jpg') )
#par(mfrow=c(1,2))
hist(algae[,i], prob=T, xlab='',main=paste('Histogram of ',a))
lines(density(algae[,i],na.rm=T))
rug(jitter(algae[,i]))
dev.off( )
jpeg(file=paste('QQ of ',a,'.jpg') )
qq.plot(algae[,i],main=paste('Normal QQ plot of ',a),ylab= as.character(a))
dev.off( )


jpeg(file=paste('box of ',a,'.jpg') )
#par(mfrow=c(1,1))
boxplot(algae[,i],ylab= as.character(a))
rug(jitter(algae[,i]),side=2)
abline(h=mean(algae[,i],na.rm=T),lty=2)
dev.off( )
#identify(algae[,i])
}

library(lattice)
for (i in 12:18)
{
a <- colum[[i]]
jpeg(file=paste('Season and ',a,'.jpg') )
print(bwplot(season ~ algae[,i], data=algae,ylab='Season',xlab=paste('Algal ',a)))
dev.off( )
jpeg(file=paste('River Size and ',a,'.jpg') )
print(bwplot(size ~ algae[,i], data=algae,ylab='River Size',xlab=paste('Algal ',a)))
dev.off( )
jpeg(file=paste('River Speed and ',a,'.jpg') )
print(bwplot(speed ~ algae[,i], data=algae,ylab='River Speed',xlab=paste('Algal ',a)))
dev.off( )
}

###################################################
### Unkwnon Values
###################################################
#algae <- algae[-manyNAs(algae),]
###0,orginal
write.csv(algae,file = "Analysis_Orginal.csv",na = "XXXXXXX")

###1,Delete
algae <- na.omit(algae)
write.csv(algae,file = "Analysis_Delete.csv",na = "XXXXXXX")

###2,most Frequency
algae <- read.table('Analysis.txt',
                    header=F,
                    dec='.',
                    col.names=c('season','size','speed','mxPH','mnO2','Cl',
                                'NO3','NH4','oPO4','PO4','Chla','a1','a2','a3','a4',
                                'a5','a6','a7'),
                    na.strings=c('XXXXXXX'))
for (i in 4:11)
{
a <- colum[[i]]
as.numeric(names(table(algae[,as.character(a)])))[which.max(table(algae[,as.character(a)]))]
algae[is.na(algae[,i]),as.character(a)] <- as.numeric(names(table(algae[,as.character(a)])))[which.max(table(algae[,as.character(a)]))]

}
write.csv(algae,file = "Analysis_Frequency.csv",na = "XXXXXXX") 
###3,Correlation
algae <- read.table('Analysis.txt',
                    header=F,
                    dec='.',
                    col.names=c('season','size','speed','mxPH','mnO2','Cl',
                                'NO3','NH4','oPO4','PO4','Chla','a1','a2','a3','a4',
                                'a5','a6','a7'),
                    na.strings=c('XXXXXXX'))
algae <- algae[-manyNAs(algae),]
sink("CorrelationMatrix.txt")   #把结果保存到CorrelationMatrix.txt
symnum(cor(algae[,4:18],use="complete.obs"))
sink()  

sink("CorrelationCoefficient.txt")   #把结果保存到CorrelationCoefficient.txt
lm(PO4 ~ oPO4,data=algae)
lm(NH4 ~ NO3,data=algae)
lm(Chla ~ mxPH,data=algae)
lm(mxPH ~Chla ,data=algae)
lm(mnO2 ~oPO4 ,data=algae)
lm(Cl ~oPO4 ,data=algae)
sink()

jpeg(file=paste('season-speed-size of mnO2','.jpg') )
stripplot(~mnO2|season *speed *size,data=algae,jitter=T)
dev.off()
jpeg(file=paste('season-speed-size of Cl','.jpg') )
stripplot(~Cl|season *speed *size,data=algae,jitter=T)
dev.off()
jpeg(file=paste('season-speed-size of Chla','.jpg') )
stripplot(~Chla|season *speed *size,data=algae,jitter=T)
dev.off()
jpeg(file=paste('season-speed-size of mxPH','.jpg') )
stripplot(~mxPH|season *speed *size,data=algae,jitter=T)
dev.off()


fillPO4 <- function(oP) {
  if (is.na(oP)) return(NA)
  else return(42.897 + 1.293 * oP)
}
algae[is.na(algae$PO4),'PO4'] <- sapply(algae[is.na(algae$PO4),'oPO4'],fillPO4)

fillChla <- function(mP) {
  if (is.na(mP)) return(NA)
  else return(abs(-139.4 + 19 * mP))
}
algae[is.na(algae$Chla),'Chla'] <- sapply(algae[is.na(algae$Chla),'mxPH'],fillChla)

fillmxPH <- function(Ch) {
  if (is.na(Ch)) return(NA)
  else return(7.92896 + 0.01047 * Ch)
}
algae[is.na(algae$mxPH),'mxPH'] <- sapply(algae[is.na(algae$mxPH),'Chla'],fillmxPH)

fillmnO2 <- function(oP) {
  if (is.na(oP)) return(NA)
  else return(9.93341 -0.01093 * oP)
}
algae[is.na(algae$mnO2),'mnO2'] <- sapply(algae[is.na(algae$mnO2),'oPO4'],fillmnO2)

fillCl <- function(oP) {
  if (is.na(oP)) return(NA)
  else return(28.4771 + 0.1992 * oP)
}
algae[is.na(algae$Cl),'Cl'] <- sapply(algae[is.na(algae$Cl),'oPO4'],fillCl)

#algae[is.na(algae$mnO2),"mnO2"] <- mean(algae$mnO2,na.rm = T)

write.csv(algae,file = "Analysis_Correlation.csv",na = "XXXXXXX") 

###4,Similarity
algae <- read.table('Analysis.txt',
                    header=F,
                    dec='.',
                    col.names=c('season','size','speed','mxPH','mnO2','Cl',
                                'NO3','NH4','oPO4','PO4','Chla','a1','a2','a3','a4',
                                'a5','a6','a7'),
                    na.strings=c('XXXXXXX'))
algae <- knnImputation(algae,k=10)
write.csv(algae,file = "Analysis_Similarity.csv",na = "XXXXXXX") 
