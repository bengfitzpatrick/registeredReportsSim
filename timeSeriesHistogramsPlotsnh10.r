library(ggplot2)
library(ggpubr)
#
# Case: valRR = 1.2
#
load("D:\\zdrive\\tempest\\projects\\integrity\\registeredReports\\nh10\\pH0p8dp5r10_50ka1.Rdata")

sPowerList <- truePosList
sAlfaList  <- falsePosList
sValueList <- valueList
sEffortList<- effortList
sRetryList <- rrList
sFPList    <- falsePosList[,2]/(truePosList[,2]+falsePosList[,2])

print(30*mean(labs$value/labs$age))
print(nLabs)

finalVal1  <- colMeans(finalVal) 
finalEff1  <- colMeans(finalEff) 
finalRR1   <- colMeans(finalRR) 
finalDelt1 <- colMeans(finalDelt) 
finalAge1  <- colMeans(finalAge) 
finalFP1   <- colSums(finalcFP)
finalTP1   <- colSums(finalcTP)
finalFPP1  <- finalFP1/(finalFP1+finalTP1)
finalcFP1   <- matrix(finalcFP,ncol=1)
finalcTP1   <- matrix(finalcTP,ncol=1)
finalcFPP1  <- finalcFP1/(finalcFP1+finalcTP1)

valRR1      <- valRR

fppDF               <- data.frame(fppC=finalcFPP1,fpp=finalFPP1)
fppDF$SimCondition  <- paste("RR value = ",valRR)

histDF              <- data.frame(value=finalVal1,effort=finalEff1,rr=finalRR1)
histDF$SimCondition <- paste("RR value = ",valRR)
#
# Case: valRR = 1.6
#
load("D:\\zdrive\\tempest\\projects\\integrity\\registeredReports\\nh10\\pH0p8dp5r10_50ka5.Rdata")

rPowerList <- truePosList
rAlfaList  <- falsePosList
rValueList <- valueList
rEffortList<- effortList
rRetryList <- rrList
rFPList    <- falsePosList[,2]/(truePosList[,2]+falsePosList[,2])

print(30*mean(labs$value/labs$age))

print(nLabs)

finalVal2  <- colMeans(finalVal) 
finalEff2  <- colMeans(finalEff) 
finalRR2   <- colMeans(finalRR) 
finalDelt2 <- colMeans(finalDelt) 
finalAge2  <- colMeans(finalAge) 
finalFP2   <- colSums(finalcFP)
finalTP2   <- colSums(finalcTP)
finalFPP2  <- finalFP2/(finalFP2+finalTP2)
finalcFP2   <- matrix(finalcFP,ncol=1)
finalcTP2   <- matrix(finalcTP,ncol=1)
finalcFPP2  <- finalcFP2/(finalcFP2+finalcTP2)
valRR2      <- valRR

hist2DF     <- data.frame(value=finalVal2,effort=finalEff2,rr=finalRR2)
hist2DF$SimCondition <- paste("RR value = ",valRR)
fpp2DF      <- data.frame(fppC=finalcFPP2,fpp=finalFPP2)
fpp2DF$SimCondition  <- paste("RR value = ",valRR)

#
# Case: valRR = 2.0
#
load("D:\\zdrive\\tempest\\projects\\integrity\\registeredReports\\nh10\\pH0p8dp5r10_50ka9.Rdata")

pPowerList <- truePosList
pAlfaList  <- falsePosList
pValueList <- valueList
pEffortList<- effortList
pRetryList <- rrList
pFPList    <- falsePosList[,2]/(truePosList[,2]+falsePosList[,2])

yrs        <- seq(1,nSimSteps)


print(30*mean(labs$value/labs$age))
print(nLabs)
fE3  <- finalEff
fE3a <- matrix(fE3,ncol=1)
finalVal3  <- colMeans(finalVal) 
finalEff3  <- colMeans(finalEff) 
finalRR3   <- colMeans(finalRR) 
finalDelt3 <- colMeans(finalDelt) 
finalAge3  <- colMeans(finalAge) 
finalFP3   <-colSums(finalcFP)
finalTP3   <-colSums(finalcTP)
finalFPP3  <- finalFP3/(finalFP3+finalTP3)
finalcFP3   <- matrix(finalcFP,ncol=1)
finalcTP3   <- matrix(finalcTP,ncol=1)
finalcFPP3  <- finalcFP3/(finalcFP3+finalcTP3)

valRR3       <- valRR
simCondition <- c(rep(paste("RR value = ",valRR1),nSimSteps),rep(paste("RR value = ",valRR2),nSimSteps),rep(paste("RR value = ",valRR3),nSimSteps))
valRRL       <- as.factor(c(rep(valRR1,nSimSteps),rep(valRR2,nSimSteps),rep(valRR3,nSimSteps)))
timeL        <- c(seq(1,nSimSteps),seq(1,nSimSteps),seq(1,nSimSteps))
effL         <- c(sEffortList[,2],rEffortList[,2],pEffortList[,2])
valL         <- c(sValueList[,2],rValueList[,2],pValueList[,2])
retryL       <- c(sRetryList[,2],rRetryList[,2],pRetryList[,2])
FPL          <- c(sFPList,rFPList,pFPList)

tsDF         <- data.frame(x=timeL,valueRR=valRRL,eff=effL,val=valL,nh=retryL,fp=FPL)

hist3DF     <- data.frame(value=finalVal3,effort=finalEff3,rr=finalRR3)
hist3DF$SimCondition <- paste("RR value = ",valRR)
fpp3DF      <- data.frame(fppC=finalcFPP3,fpp=finalFPP3)
fpp3DF$SimCondition  <- paste("RR value = ",valRR)

histDF <- rbind(histDF,hist2DF,hist3DF)
fppDF <- rbind(fppDF,fpp2DF,fpp3DF)

#
# generate time series plot objects
#
effPlot <- ggplot(data=tsDF,aes(x=x,y=eff))+geom_line(aes(color = simCondition),linewidth=2) + 
  scale_color_manual(values = c("black","blue","red")) + #ylim(c(0,1)) + #xlim(c(1.5,3.5)) + 
  xlab("Time Step") + ylab("Effort") + ggtitle("Median Effort vs Time Step")
 valPlot <- ggplot(data=tsDF,aes(x=x,y=val))+geom_line(aes(color = simCondition),linewidth=2) + 
  scale_color_manual(values = c("black","blue","red")) + #ylim(c(0,1)) + #xlim(c(1.5,3.5)) + 
  xlab("Time Step") + ylab("Value") + ggtitle("Median Value vs Time Step") 
 nhPlot <- ggplot(data=tsDF,aes(x=x,y=nh))+geom_line(aes(color = simCondition),linewidth=2) + 
  scale_color_manual(values = c("black","blue","red")) + ylim(c(0,10)) + #xlim(c(1.5,3.5)) + 
  xlab("Time Step") + ylab("Number of Hypotheses") + ggtitle("Median Number of Hypotheses vs Time Step")
 fpPlot <- ggplot(data=tsDF,aes(x=x,y=fp))+geom_line(aes(color = simCondition),linewidth=2) + 
  scale_color_manual(values = c("black","blue","red")) + ylim(c(0,1)) + #xlim(c(1.5,3.5)) + 
  xlab("Time Step") + ylab("False Positive Rate") + ggtitle("Fraction of False Positive Pubs vs Time Step") 
#
# generate time series plot grid
#
fPlot <- ggarrange(effPlot,valPlot,nhPlot,fpPlot)
dev.new(width=12,height=8)
print(fPlot)
#
# generate histogram plots
#
dev.new(width=12,height=8)
p1 <- ggplot(histDF, aes(x=effort, color=SimCondition, fill=SimCondition)) + geom_histogram(position = "identity", alpha = 0.5, binwidth = 5, boundary = 0)  +
scale_color_manual(values=c("#000000", "#090de0","#de0a02" ))+scale_fill_manual(values=c("#404040", "#090de0", "#de0a02")) +
xlim(c(0,100))+xlab("Effort")
print(p1)

dev.new(width=12,height=8)
p2<-ggplot(fppDF, aes(x=fpp, color=SimCondition, fill=SimCondition)) + geom_histogram(position = "identity", alpha = 0.5, binwidth = 0.02, boundary = 0) +
scale_color_manual(values=c("#000000", "#090de0","#de0a02" ))+scale_fill_manual(values=c("#404040", "#090de0", "#de0a02")) +
xlim(c(0,1)) + xlab("Fraction of False Positive Publications")
print(p2)
dev.new(width=12,height=8)
p3<-ggplot(histDF, aes(x=rr, color=SimCondition, fill=SimCondition)) + geom_histogram(position = "identity", alpha = 0.5, binwidth = 1, boundary = 0) +
scale_color_manual(values=c("#000000", "#090de0","#de0a02" ))+scale_fill_manual(values=c("#404040", "#090de0", "#de0a02"))+
xlim(c(0,10)) + xlab("Number of Hypotheses") + scale_x_continuous(breaks=seq(0,10,by=1))
print(p3)
