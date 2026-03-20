rm(list=ls())
library(dplyr)
library("tictoc")
tic()
#
# This code is designed to run many cases, looping over the value of an RR publication.
#     Warning: if looping over 12 cases, it will several days to run
#
#     To speed up, run a single case, or reduce the number of Monte Carlo replications, 
#                  or reduce the number of sim steps
#
#     Warning: changing maxHypotheses to 50 will drastically increase the run time.

myGrid <- seq(1.2,2.4,by=0.2)
nGrid <- length(myGrid)

for (iRun in 1:nGrid){  

#
# Specify simulation parameters
#

myFileName         <- paste("pH0p8dp5r10_50ka",iRun,".Rdata",sep="")         # file for saving output
nSimSteps          <- 50000                  # number of simulation time steps
nMC                <- 100                    # number of Monte Carlo replication runs
nLabs              <- 2000                   # number of researchers in the sim
propRR             <- 0.05
nRR                <- ceiling(propRR*nLabs)
valRR              <- myGrid[iRun]  

maxHypotheses      <- 10;  # upper limit on number of hypotheses 
nHypotheses        <- sample(seq(1,maxHypotheses), nLabs, replace = TRUE)    # initial number of hypotheses (p-hacking model)
alfa               <- 0.05                   # p-value threshold for testing hypotheses  
                

pH0True   		   <- 0.8                   # prior probability of true null
deltMin            <- 0.5                   # minimum true effect size for false null


                     # initial effort
maxEff      	   <- 200                    # maximum allowable effort
minEff             <- 5                      # minimum allowable effort
effRR              <- 64;

effort             <- sample(seq(minEff,maxEff), nLabs, replace = TRUE)

LM                 <- 0.1                    # fraction of additional labs to draw from in evolution
deathRate          <- 0.002                  # probability of lab death
delayTime          <- 100;                   # time step at which evolution begins

sigma_RR           <- 1                      # standard deviation for nHypotheses evolution
sigma_e            <- 10                     # standard deviation for effort evolution
sigma_prRR         <- 0.05;                  # std dev prob RR evolution

eta                <- 0.01*1.5;                  # shape parameter for effort-to-experiment-attempt function

#
# Initialize arrays for the simulation
#
lab_id             <- 1:nLabs
testMat            <- matrix(0,ncol=nLabs,nrow=maxHypotheses)  # testMat is an indicator for hypothesis tests for each lab
deltMat            <- testMat;                                 # deltMat is a matrix for the true effect sizes, labs X nHypotheses
#
# labs is a data frame containing the records for all living researchers
#
labs               <- data.frame(matrix(ncol=20, nrow= nLabs))
colnames(labs)     <- c('age','id','effort','truePos','nHypotheses','value','falsePos','runExp','cumTruePos','cumFalsePos','pub_rate','cumPosPubs','cumNegPubs','cumRRTruePos','cumRRTrueNeg','cumRRFalsePos','cumRRFalseNeg','prRR','nhSave','fpRR')


#
# save quartiles, across MC runs, for various quantities of interest
#
falsePosList   <- matrix(0,ncol=3,nrow=nSimSteps)
truePosList    <- matrix(0,ncol=3,nrow=nSimSteps)
falseNegList   <- matrix(0,ncol=3,nrow=nSimSteps)
trueNegList    <- matrix(0,ncol=3,nrow=nSimSteps)
effortList     <- matrix(0,ncol=3,nrow=nSimSteps)
valueList      <- matrix(0,ncol=3,nrow=nSimSteps)
rrList         <- matrix(0,ncol=3,nrow=nSimSteps)
mutList        <- matrix(0,ncol=3,nrow=nSimSteps)
rrAdopters     <- matrix(0,ncol=3,nrow=nSimSteps)
rrPapers       <- matrix(0,ncol=3,nrow=nSimSteps)

#
# save final results for each living lab and each MC run
#
finalRR   <- matrix(0,nrow=nLabs,ncol=nMC)
finalEff  <- matrix(0,nrow=nLabs,ncol=nMC)
finalVal  <- matrix(0,nrow=nLabs,ncol=nMC)
finalFP   <- matrix(0,nrow=nLabs,ncol=nMC)
finalTP   <- matrix(0,nrow=nLabs,ncol=nMC)
finalFN   <- matrix(0,nrow=nLabs,ncol=nMC)
finalTN   <- matrix(0,nrow=nLabs,ncol=nMC)
finalAge  <- matrix(0,nrow=nLabs,ncol=nMC)
finalDelt <- matrix(0,nrow=nLabs,ncol=nMC)
finalcFP  <- matrix(0,nrow=nLabs,ncol=nMC)
finalcTP  <- matrix(0,nrow=nLabs,ncol=nMC)
finalcFN  <- matrix(0,nrow=nLabs,ncol=nMC)
finalcTN  <- matrix(0,nrow=nLabs,ncol=nMC)
finalcExp <- matrix(0,nrow=nLabs,ncol=nMC)
finalYNRR <- matrix(0,nrow=nLabs,ncol=nMC)
finalPosP <- matrix(0,nrow=nLabs,ncol=nMC)
finalNegP <- matrix(0,nrow=nLabs,ncol=nMC)
finalfRRP <- matrix(0,nrow=nLabs,ncol=nMC)


#
#  indxMat is used for indexing the researchers' number of hypotheses over time
#
indxMat    <- rep(seq(1,maxHypotheses),times=nLabs);
indxMat    <- t(matrix(indxMat,ncol=nLabs))
#
# monte carlo replication loop
#
for (iMC in 1:nMC){

#
#  initialize the lab dataframe
#
    print(paste("iMC = ",iMC,sep=""))
	labs$age               <- 0
	labs$id                <- seq(1,nLabs)
	labs$effort            <- effort
	labs$truePos           <- 0
	labs$nHypotheses       <- nHypotheses
	labs$value             <- 0
	labs$falsePos          <- 0
	labs$runExp            <- 0
	labs$cumTruePos        <- 0
	labs$cumFalsePos       <- 0
	labs$cumTrueNeg        <- 0
	labs$cumFalseNeg       <- 0
	labs$cumPosPubs        <- 0
	labs$cumNegPubs        <- 0
	labs$pub_rate          <- 0
	labs$cumRRTruePos      <- 0
	labs$cumRRTrueNeg      <- 0
	labs$cumRRFalsePos     <- 0
	labs$cumRRFalseNeg     <- 0
	labs$prRR              <- runif(nLabs,0,1)*ifelse(labs$effort<effRR,0,1) #ifelse(labs$effort<effRR,0,1) #runif(nLabs,0,1)*ifelse(labs$effort<effRR,0,1)
	labs$nhSave            <- nHypotheses
	labs$fpRR              <- 0
#
#  simulate the base effect sizes for the labs 
#
	delta1                 <- rgamma(nLabs,1,rate = (-log(1-pH0True)/deltMin)); #*ifelse(labs$YNRR==1,2,1)
	myDelta                <- matrix(rep(delta1,times=maxHypotheses),ncol=maxHypotheses)		


	
#
# time stepping loop
#
	for (i in 1:nSimSteps) {

#
# labs doing RRs
# 

	labs$nHypotheses <- labs$nhSave
	
	pRR         <- runif(nLabs,0,1)
#	indxPossRR  <- which((labs$effort>=effRR)&(delta1>=deltMin)&(pRR<=labs$prRR));
	indxPossRR  <- which((labs$effort>=effRR)&(pRR<=labs$prRR));
	numRR       <- length(indxPossRR)#min(length(indxPossRR),nRR)
	indxRR      <- indxPossRR[1:numRR]
	labs$nHypotheses[indxRR] <- 1
	
	labs$YNRR   <- 0
	labs$YNRR[indxRR] <- 1
	indxNot   <- which(labs$YNRR==0)
#
#  run an experiment?
#	
	  h       <- exp(-(labs$effort-minEff)*eta)  # probability of conducting experiment this time step 
      hRR     <- h; #/2;
	  w       <- runif(nLabs,0,1)

	  
	  runExp       <- ifelse(h>=w,1,0);          # simulate whether or not to run experiment
	  labs$runExp  <- labs$runExp+runExp;
#
#     Set the t-critical values for alfa and Sidak correction
#	  
	  alfSidak   <- 1-(1-alfa)^(1/labs$nHypotheses)     
	  t_alfSidak <- qt(1-alfSidak/2, 2*labs$effort-2)   #critical T for multiple testing
	  t_alfa     <- qt(1-alfa/2, 2*labs$effort-2)       #critical T for one-test 

#
#  initialize arrays for hypothesis testing
#
	  detct    <- rep(0,nLabs);
	  truth    <- rep(0,nLabs);
	  detCorr  <- rep(0,nLabs);	  
	  rrMat    <- matrix(rep(labs$nHypotheses,times=maxHypotheses),ncol=maxHypotheses);
	  testMat  <- ifelse(indxMat<=rrMat,1,0)	  
#
#    Simulate effect sizes using the lab effect
#	  
	  myDelta4 <- matrix(rgamma(nLabs,1,scale = myDelta),ncol=1)	   #matrix(abs(rnorm(nLabs,myDelta,0.1*myDelta)),ncol=1) #matrix(rgamma(nLabs,1,scale = myDelta),ncol=1)	  
	  myDelta3 <- matrix(rep(myDelta4,times=maxHypotheses),ncol=maxHypotheses);
	  myDelta2 <- (myDelta3)*testMat
#
#     find the true positives: myTruth is an indicator matrix of false nulls
#
	  myTruth  <- ifelse(myDelta2<=deltMin,0,1)
#
#     generate t-statistics for each lab, each hypothesis
#     non-central t's used for true effect sizes
#	  
	  effMat   <- matrix(rep(labs$effort,times=maxHypotheses),ncol=maxHypotheses); # effect sizes
	  lamda    <- sqrt(effMat)*myDelta2/sqrt(2)                                    # noncentrality
	  tVals    <- abs(rt(rrMat,2*effMat-2,lamda)*testMat)                          # t statistics
	  tMax     <- apply(tVals,1,max)	                                           # maximize over hypotheses  
	  tMaxMat  <- matrix(rep(tMax,times=maxHypotheses),ncol=maxHypotheses);        
	  chkMax   <- ifelse(tMaxMat == tVals,1,0);                                    # which hypothesis meets the max
	  detct    <- tMax>t_alfa                                                      # max beats standard T
	  detCorr  <- tMax>t_alfSidak                                                  # max beats Sidak T
	  truth    <- ifelse(rowSums(chkMax*myTruth)>0,1,0)                            # max was at a false null?
	  
	  truePos <- ifelse(((truth==1)&(detCorr==1)&(h>=w)),1,0)                      # false null, beat Sidak, experiment ran
	  falsPos <- ifelse(((truth==0)&(detct==1)&(h>=w)),1,0)                        # true null, beat T, experiment ran
	  trueNeg <- ifelse(((truth==0)&(detCorr==0)&(h>=w)),1,0)                      # false null, beat Sidak, experiment ran
	  falsNeg <- ifelse(((truth==1)&(detct==0)&(h>=w)),1,0)                        # true null, beat T, experiment ran
#
#  record experiment outcomes in labs data frame
#	  
	  labs$falsePos      <- falsPos
	  labs$truePos       <- truePos
	  #labs$value         <- labs$value + (truePos + falsPos)
	  labs$age           <- labs$age + 1
	  labs$cumTruePos    <- labs$cumTruePos  + truePos
	  labs$cumFalsePos   <- labs$cumFalsePos + falsPos 
	  labs$cumTrueNeg    <- labs$cumTrueNeg  + trueNeg
	  labs$cumFalseNeg   <- labs$cumFalseNeg + falsNeg
	  labs$cumRRTruePos    <- labs$cumTruePos  + truePos*ifelse(labs$YNRR==1,1,0)
	  labs$cumRRFalsePos   <- labs$cumFalsePos + falsPos*ifelse(labs$YNRR==1,1,0)
	  labs$cumRRTrueNeg    <- labs$cumTrueNeg  + trueNeg*ifelse(labs$YNRR==1,1,0)
	  labs$cumRRFalseNeg   <- labs$cumFalseNeg + falsNeg*ifelse(labs$YNRR==1,1,0)
#	  labs$pub_rate      <- labs$value/labs$age
	  
	  

	  valRRNow           <- 0
	  pubRRNow           <- 0
	  
	  if (numRR>0){
		valRRNow                 <- ifelse(hRR[indxRR] >= w[indxRR],valRR,0)
		pubRRNow                 <- ifelse(hRR[indxRR] >= w[indxRR],1,0)
		labs$cumPosPubs[indxRR]  <- labs$cumPosPubs[indxRR] + pubRRNow*(truePos[indxRR] + falsPos[indxRR])
		labs$cumNegPubs[indxRR]  <- labs$cumNegPubs[indxRR] + pubRRNow*(1-(truePos[indxRR] + falsPos[indxRR]))
	  
		labs$pub_rate[indxRR]      <- labs$pub_rate[indxRR] + pubRRNow 	  
		labs$value[indxRR]         <- labs$value[indxRR] + valRRNow
		
	  }
	  
	  nindxNot           <- length(indxNot)
	  rrPapersN          <- 0

	  if(numRR<nLabs){
		labs$value[indxNot]  <- labs$value[indxNot] + (truePos[indxNot] + falsPos[indxNot])
	  
		labs$cumPosPubs[indxNot] <- labs$cumPosPubs[indxNot] + (truePos[indxNot] + falsPos[indxNot])

		labs$pub_rate[indxNot]     <- labs$value[indxNot]
	  	rrPapersN <- sum(truePos[indxNot] + falsPos[indxNot])
		
	  }
	  rrPapersY <- sum(pubRRNow)

	  if (is.na(sum(pubRRNow))){browser()}
	  
	  denom         <- rrPapersY+ rrPapersN
	  denom         <- max(1,denom)
	  rrPapersY <- rrPapersY/denom
	  rrPapersN <- rrPapersN/denom
#
# if time is right, perform evolution
#	  

	  if(i>delayTime){
#
#  who dies 
#	  
	  	deaths      <- runif(nLabs,0,1)
		indxDeath   <- which(deaths<=deathRate);
		
		#indxDeath    <- which(labs$value == min(labs$value))
#
#       sort the labs
#	  
	    L           <- length(indxDeath);   		# number of deaths
		
		if(L>0){
			L2          <- min(floor(L*(1+LM)),nLabs);                             # number of "best" labs to draw from
			bestLabs    <- sort.int(labs$value,index.return=TRUE,decreasing=TRUE)  # sort labs in terms of value
			#indxDeath   <- bestLabs$ix[(nLabs-L+1):nLabs]
			indxRep     <- sample(bestLabs$ix[1:L2],L,replace=FALSE)               # indices of L labs sampled from the best L*(1+LM) 
								   
			labs[indxDeath,]       <- labs[indxRep,]   # initialize new labs with data from "best"
	#
	#  reset parameters not associated with evolution inheritance
	#
			labs$age[indxDeath]            <-1
			labs$runExp[indxDeath]         <-0
			labs$value[indxDeath]          <-0
			labs$cumFalsePos[indxDeath]    <-0
			labs$cumTruePos[indxDeath]     <-0
			labs$cumFalseNeg[indxDeath]    <-0
			labs$cumTrueNeg[indxDeath]     <-0
			labs$pub_rate[indxDeath]       <-0
			labs$cumNegPubs[indxDeath]     <-0
			labs$cumPosPubs[indxDeath]     <-0
			labs$id[indxDeath]             <- seq(max(labs$id)+1,max(labs$id)+L)
			
			indxRR  <- which(labs$YNRR==1)
			indxNot <- which(labs$YNRR==0)

            newRR   <- intersect(indxRR,indxDeath)
	#
	#  perform evolution in terms of effort and nHypotheses
	#	  
			labs$effort[indxDeath] <- round(labs$effort[indxDeath] + rnorm(L,0,sigma_e))
			labs$effort            <- ifelse(labs$effort>maxEff,maxEff,labs$effort)
			labs$effort            <- ifelse(labs$effort<minEff,minEff,labs$effort)
			
			labs$nHypotheses[indxDeath] <- round(labs$nHypotheses[indxDeath] + rnorm(L,0,sigma_RR))
			labs$nHypotheses            <- ifelse(labs$nHypotheses>maxHypotheses,maxHypotheses,labs$nHypotheses)
			labs$nHypotheses            <- ifelse(labs$nHypotheses<1,1,labs$nHypotheses)  
			
			labs$nhSave[indxDeath] <- labs$nHypotheses[indxDeath]
			labs$prRR[indxDeath]   <- labs$prRR[indxDeath]+rnorm(L,0,sigma_prRR)
			labs$prRR[indxDeath]   <- ifelse(labs$prRR[indxDeath]>=1,1,labs$prRR[indxDeath])
			labs$prRR[indxDeath]   <- ifelse(labs$prRR[indxDeath]<=0,0,labs$prRR[indxDeath])			
	#
	#   set effect size of new labs from old labs
	#

			deltaD                 <- rgamma(L,1,rate = (-log(1-pH0True)/deltMin));
			myDelta[indxDeath,]    <- matrix(rep(deltaD,times=maxHypotheses),ncol=maxHypotheses)
	        }
	  }
	  
	  truePosList[i,]  <- truePosList[i,] + c(sum(labs$truePos)/nLabs-1/sqrt(nLabs),sum(labs$truePos)/nLabs,sum(labs$truePos)/nLabs+1/sqrt(nLabs))/nMC 
	  
	  falsePosList[i,]   <- falsePosList[i,] + c(sum(labs$falsePos)/nLabs-1/sqrt(nLabs),sum(labs$falsePos)/nLabs,sum(labs$falsePos)/nLabs+1/sqrt(nLabs))/nMC  
	  
	  trueNegList[i,]  <- trueNegList[i,] + c(sum(runExp*(1-labs$falsePos))/nLabs-1/sqrt(nLabs),sum(runExp*(1-labs$falsePos))/nLabs,sum(runExp*(1-labs$falsePos))/nLabs+1/sqrt(nLabs))/nMC 
	  
	  falseNegList[i,]   <- falseNegList[i,] + c(sum(runExp*(1-labs$truePos))/nLabs-1/sqrt(nLabs),sum(runExp*(1-labs$truePos))/nLabs,sum(runExp*(1-labs$truePos))/nLabs+1/sqrt(nLabs))/nMC   
	  
	  effortList[i,]  <- effortList[i,] + quantile(labs$effort,probs = c(0.25,0.5,0.75))/nMC
	  valueList[i,]   <- valueList[i,] + quantile(labs$value,probs = c(0.25,0.5,0.75))/nMC
	  rrList[i,]      <- rrList[i,] + quantile(labs$nHypotheses,probs = c(0.25,0.5,0.75))/nMC
	  
	  rrAdopters[i,1] <- rrAdopters[i,1] + (length(indxRR)/nLabs)/nMC
	  rrPapers[i,1]   <- rrPapers[i,1]  + rrPapersY/nMC
      rrPapers[i,2]   <- rrPapers[i,2]  + rrPapersN/nMC
	}
	medDelt         <- 0*labs$effort
	for (iii in 1:nLabs){
		 medDelt[iii]<- median(myDelta[iii,(1:labs$nHypotheses[iii])])
	}
	labs$delt <- medDelt
	
	finalRR[,iMC]   <- labs$nHypotheses
	finalEff[,iMC]  <- labs$effort
	finalVal[,iMC]  <- labs$value
	finalFP[,iMC]   <- labs$falsePos
	finalTP[,iMC]   <- labs$truePos
	finalAge[,iMC]  <- labs$age
	finalDelt[,iMC] <- medDelt
	finalcFP[,iMC]  <- labs$cumFalsePos
	finalcTP[,iMC]  <- labs$cumTruePos
	finalcFN[,iMC]  <- labs$cumFalseNeg
	finalcTN[,iMC]  <- labs$cumTrueNeg
	finalcExp[,iMC] <- labs$runExp
	finalYNRR[,iMC] <- labs$YNRR
	finalPosP[,iMC] <- labs$cumPosPubs
	finalNegP[,iMC] <- labs$cumNegPubs
	finalfRRP[,iMC] <- labs$prRR
}

save.image(file=myFileName)
print(iRun)


}

toc()

# print(labs)

# dev.new()
# plot(truePosList[,2],ylim=c(0,1),xlab="time step",ylab="True Discovery Rate")
# lines(truePosList[,1],col="red")
# lines(truePosList[,3],col="red")
# dev.new()
# plot(falsePosList[,2],ylim=c(0,1),xlab="time step",ylab="False Discovery Rate")
# lines(falsePosList[,1],col="red")
# lines(falsePosList[,3],col="red")
# dev.new()
# plot(effortList[,2],xlab="time step",ylab="Effort",ylim=c(0,maxEff))
# lines(effortList[,1],col="red")
# lines(effortList[,3],col="red")
# dev.new()
# plot(valueList[,2],xlab="time step",ylab="Value",ylim=c(0,max(valueList)))
# lines(valueList[,1],col="red")
# lines(valueList[,3],col="red")
# dev.new()
# plot(rrList[,2],ylim=c(0,maxHypotheses),xlab="time step",ylab="re-try rate")
# lines(rrList[,1],col="red")
# lines(rrList[,3],col="red")
# 
# save.image(file=myFileName)
