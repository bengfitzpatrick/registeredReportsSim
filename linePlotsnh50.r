library("ggplot2")
library("tidyverse")

nFiles <- 12 
fName  <- "pH0p8dp5r50_50k_r_"
gTitle <- "Pr[H0] = 0.8, effects below 0.5 considered null, max number of hypotheses = 50"


pRRa <- rep(0,nFiles)
pFP <- rep(0,nFiles)
vRR <- rep(0,nFiles)
pPP <- rep(0,nFiles)


for (ii in 1:nFiles){


fileName <- paste("D:\\zdrive\\tempest\\projects\\integrity\\registeredReports\\nH50\\",fName,ii,".Rdata",sep="")

load(fileName)

vRR[ii] <- valRR
pRRa[ii] <- mean(finalfRRP)       



finalFP <- as.vector(finalNegP)
finalTP <- as.vector(finalPosP)

s       <- finalFP+finalTP
frTP    <- finalTP/s
pPP[ii] <- sum(finalTP)/sum(s)

finalP   <- as.vector(finalcTP) + as.vector(finalcFP)
finalFPa <- as.vector(finalcFP);
s        <- finalP
frFP     <- finalFPa/s
pFP[ii]  <- sum(finalFPa)/sum(s)

}

myDF         <- data.frame(vRR = vRR,pRRa = pRRa,pPP = pPP,pFP = pFP)
names(myDF)  <- c("Value of RR publication","Proportion of pubs that are RR","Proportion of pos pubs","Proportion of false pos pubs")

nRR <- rep(names(myDF)[2],length(vRR))
nPP <- rep(names(myDF)[3],length(vRR))
nFP <- rep(names(myDF)[4],length(vRR))
newDF <- data.frame(x=c(vRR,vRR,vRR),y=c(pRRa,pPP,pFP),z=c(nRR,nPP,nFP))
myLinePlot <- ggplot(data=newDF,aes(x=x,y=y))+geom_line(aes(color = z, linetype = z),linewidth=2) + 
  scale_color_manual(values = c("darkred","steelblue","black")) + ylim(c(0,1)) + xlim(c(1.5,3.5)) + 
  xlab("Value of RR Publication") + ylab("Proportions") + ggtitle(gTitle)
 print(myLinePlot)
