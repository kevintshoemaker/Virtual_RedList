##############################
###### R script for analyzing the NASA redlist data in terms of the robustness
######   of the red list categories. Specifically, this script generates histograms
######    of time listed to visualize the results. 


KEVIN = TRUE # FALSE # TRUE
KEVINLAB = FALSE # TRUE
JESSIE = FALSE


########################
### load predictor variables and response variables from all 'realized species' from the NASA project

if(KEVIN) setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\RData")
 # if(KEVINLAB) setwd("")
load("FinalNASAData.RData")

##########################
###   load model-specific redlist data

if(KEVIN) setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\RData")
if(KEVINLAB) setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\RData")
load("ModelSpecificRLData_2014-04-22.RData")      # ModelSpecificRLData_2014-02-14.RData

## ls()


###################################
#######  DEVELOP MASTER DATA FRAME

keepers <- c("mp.run","clim.chg","model","spec")    # "spec",

RLMaster <- subset(dfmain, subset=(match.stat=="Match" & clim.chg=='WRE'), select=keepers)   # subset only non-mismatches with climate scenario WRE
RLMaster$mp.run <- as.character(RLMaster$mp.run)

RLMaster_old <- RLMaster

cattimes.ex$spec <- as.character(cattimes.ex$spec)
cattimes.noex$spec <- as.character(cattimes.noex$spec)

cattimes.ex$mp.run <- cattimes.ex$spec
cattimes.noex$mp.run <- cattimes.noex$spec

newspec <- substr(cattimes.ex$spec,1,4)
drop.cols <- c("model","spec")
cattimes2.ex <- subset(cattimes.ex,select = -which(names(cattimes.ex)%in%drop.cols))
RLMaster <- merge(cattimes2.ex, RLMaster, by='mp.run')     ##  head(merge(cattimes.ex, RLMaster, by='mp.run'),100)

newspec <- substr(cattimes.noex$spec,1,4)
drop.cols <- c("model","spec")
cattimes2.noex <- subset(cattimes.noex,select = -which(names(cattimes.noex)%in%drop.cols))
RLMaster.noex <- merge(cattimes2.noex, RLMaster_old, by='mp.run')     ##  head(merge(cattimes.ex, RLMaster, by='mp.run'),100)

RLMaster$model<-as.character(RLMaster$model)
RLMaster$spec <- as.character(RLMaster$spec)

RLMaster.noex$model<-as.character(RLMaster.noex$model)
RLMaster.noex$spec <- as.character(RLMaster.noex$spec)

 # RLMaster.noex$vu.lcens

#####################################
          # determine the extinction risk under the no climate change scenario
dfmain4 <- subset(dfmain, subset=(match.stat=="Match" & clim.chg=='NoCC'), select=c(keepers,"ext.risk"))
  #head(dfmain4)

dfmain4$mp.run <- as.character(dfmain4$mp.run)

dfmain4 <- merge(cattimes2.ex, dfmain4, by='mp.run')     ##  head(merge(cattimes.ex, RLMaster, by='mp.run'),100)
  #nrow(dfmain4)
  #nrow(cattimes2.ex)

 #head(RLMaster)
RLMaster$ext.risk.nocc <- dfmain4$ext.risk


       ## And for non-extinction runs
          # determine the extinction risk under the no climate change scenario
dfmain4 <- subset(dfmain, subset=(match.stat=="Match" & clim.chg=='NoCC'), select=c(keepers,"ext.risk"))
   #head(dfmain4)

dfmain4$mp.run <- as.character(dfmain4$mp.run)

dfmain4 <- merge(cattimes2.noex, dfmain4, by='mp.run')     ##  head(merge(cattimes.ex, RLMaster, by='mp.run'),100)
  #nrow(dfmain4)
  #nrow(cattimes2.ex)

  #head(RLMaster.noex)
RLMaster.noex$ext.risk.nocc <- dfmain4$ext.risk


#####################################
          # determine the extinction risk under the WRE climate change scenario
dfmain2 <- subset(dfmain, subset=(match.stat=="Match" & clim.chg=='WRE'), select=c(keepers,"ext.risk"))
  #head(dfmain2)

dfmain2$mp.run <- as.character(dfmain2$mp.run)

dfmain2 <- merge(cattimes2.ex, dfmain2, by='mp.run')     ##  head(merge(cattimes.ex, RLMaster, by='mp.run'),100)
  #nrow(dfmain2)
  #nrow(cattimes2.ex)


RLMaster$ext.risk.wre <- dfmain2$ext.risk
#head(RLMaster)
#RLMaster$ext.risk.wre

     ## and for non-extinct runs

dfmain2 <- subset(dfmain, subset=(match.stat=="Match" & clim.chg=='WRE'), select=c(keepers,"ext.risk"))
 #head(dfmain2)

dfmain2$mp.run <- as.character(dfmain2$mp.run)

dfmain2 <- merge(cattimes2.noex, dfmain2, by='mp.run')     ##  head(merge(cattimes.ex, RLMaster, by='mp.run'),100)
  #nrow(dfmain2)
  #nrow(cattimes2.ex)


RLMaster.noex$ext.risk.wre <- dfmain2$ext.risk
 #head(RLMaster)
 # max(RLMaster.noex$ext.risk.wre)



   ## first combine the two data sets...

RLMaster_comb <- rbind(RLMaster,RLMaster.noex)

  ## add extinction variable!

RLMaster_comb$is.ext <- c(rep(TRUE,times=nrow(RLMaster)),rep(FALSE,times=nrow(RLMaster.noex)))

    ## filter out the runs that went extinct during the last 20 years of the simulation   
 #RLMaster_comb <- subset(RLMaster_comb,ex.year<=90|is.ext==F,na.rm=F)
   #RLMaster_comb$ex.year



###############################################################
###############################################################
###############################################################
###############################################################
###################################################
######################################        HISTOGRAMS




###############################################
#########   Histograms of years continuously listed as VU+ prior to extinction:

pnl1Var <- RLMaster$vu.timeto[-which(RLMaster$vu.lcens==1)]
pnl2Var <- RLMaster$vu.timeto5[-which(RLMaster$vu.lcens5==1)]
pnl3Var <- RLMaster$vu.timeto10[-which(RLMaster$vu.lcens10==1)]

pnl1Var[which(pnl1Var<=0)] <- 0
pnl2Var[which(pnl2Var<=0)] <- 0
pnl3Var[which(pnl3Var<=0)] <- 0

pnl1Var[which(pnl1Var>=100)] <- 100
pnl2Var[which(pnl2Var>=100)] <- 100
pnl3Var[which(pnl3Var>=100)] <- 100

pnl1_n <- length(pnl1Var)
pnl2_n <- length(pnl2Var)
pnl3_n <- length(pnl3Var)

pnl1_n
pnl2_n
pnl3_n


  #head(RLMaster,100)

breaks <- c(0,seq(10,90,10),100)

graphics.off()

setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\RawSVG")

filename=paste("YearsToExt_Histograms_",Sys.Date(),".svg",sep="") 
svg(filename=filename,width=4,height=7,onefile=TRUE)

par(mfrow=c(3,1))

hist(pnl1Var,freq=F,main="Evaluation interval = annual",xlab="",xlim=c(0,100),col=gray(0.9),breaks=breaks)
#lines(density(pnl1Var,bw=10),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl1Var),lwd=2)
points(median(pnl1Var),max(hist(pnl1Var,plot=F)$density),pch='|',cex=2)     # max(hist(pnl1Var,freq=F)$density)
#abline(v=20,lwd=2,lty=2)

hist(pnl2Var,freq=F,main="Evaluation interval = 5 yrs",xlab="",xlim=c(0,100),col=gray(0.9),breaks=breaks)
#lines(density(pnl2Var,bw=10),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl2Var),lwd=2)
points(median(pnl2Var),max(hist(pnl2Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

hist(pnl3Var,freq=F,main="Evaluation interval = 10 yrs",
         xlab="Years continously listed as threatened before extinction",xlim=c(0,100),col=gray(0.9),breaks=breaks)
#lines(density(pnl3Var,bw=10),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl3Var),lwd=2)
points(median(pnl3Var),max(hist(pnl3Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

dev.off()

graphics.off()


#### compute the proportion below 20 years

length(which(pnl1Var<=20))/length(pnl1Var)
length(which(pnl2Var<=20))/length(pnl1Var)
length(which(pnl3Var<=20))/length(pnl1Var)


#### for supplemental material: compute proportion of runs meeting or exceeding a range of thresholds

wt_thresholds <- c(5,10,15,20,25,30)       # thresholds for adequate warning times, in years
temp1 <- numeric(length(wt_thresholds))
temp2 <- numeric(length(wt_thresholds))
temp3 <- numeric(length(wt_thresholds))
for(i in 1:length(wt_thresholds)){
   temp1[i] <- length(which(pnl1Var<=wt_thresholds[i]))/length(pnl1Var)
   temp2[i] <- length(which(pnl2Var<=wt_thresholds[i]))/length(pnl2Var)
   temp3[i] <- length(which(pnl3Var<=wt_thresholds[i]))/length(pnl3Var)
}
tempdf <- data.frame(Eval1=temp1,Eval5=temp2,Eval10=temp3)
rownames(tempdf) <- paste(wt_thresholds," years",sep="")

setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\Results")
filename=paste("Sens_to_WTthreshold_reev_int",Sys.Date(),".csv",sep="")
write.csv(tempdf,row.names=T,file=filename)   # col.names=T,  ,sep=","

###########################################
########   Histograms for 10 year reevaluation interval: separate for VU, EN, and CR

response <- "timeto"
criteria <- "" # "A" # 
category <- "vu"  # "cr"    #  critically endangered 
eval.time <- "10" # ""   # "5" #      
resp.name <- "Time continuously listed as Vulnerable \nbefore extinction"  #  "Time continuously listed as critical \nbefore extinction"  
dframe <-  "RLMaster"   # "dfmain" #  "df_rel" # dfmain or df_rel for absolute or relative metrics       

full.response <- paste(category,".",response,eval.time,criteria,sep="")
full.respname <- ifelse(nchar(criteria)>0 ,paste(resp.name," under criterion ",criteria,sep=""),resp.name)

left.censored <- paste(category,".","lcens",eval.time,sep="")

eval(parse(text=paste("subset2 <- subset(",dframe,",",left.censored,"== 0",")",sep="") ))

nrow(subset2)


################################
######  Plot histograms

pnl1Var <- subset2$vu.timeto10
pnl2Var <- subset2$en.timeto10
pnl3Var <- subset2$cr.timeto10

pnl1Var[which(pnl1Var<=0)] <- 0
pnl2Var[which(pnl2Var<=0)] <- 0
pnl3Var[which(pnl3Var<=0)] <- 0

pnl1Var[which(pnl1Var>=100)] <- 100
pnl2Var[which(pnl2Var>=100)] <- 100
pnl3Var[which(pnl3Var>=100)] <- 100


pnl1_n <- length(pnl1Var)
pnl2_n <- length(pnl2Var)
pnl3_n <- length(pnl3Var)

pnl1_n
pnl2_n
pnl3_n


 # head(RLMaster,100)


graphics.off()

setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\RawSVG")

filename=paste("VUvsENvsCR_Histograms_",Sys.Date(),".svg",sep="") 

svg(filename=filename,width=4,height=7,onefile=TRUE)

par(mfrow=c(3,1))

hist(pnl1Var,freq=F,main="Years listed as Vulnerable",xlab="",
                 xlim=c(0,100),ylim=c(0,0.035),col=gray(0.9),breaks=breaks)
#lines(density(pnl1Var,bw=10),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl1Var),lwd=2)
points(median(pnl1Var),max(hist(pnl1Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

hist(pnl2Var,freq=F,main="Years listed as Endangered",xlab="",
                  xlim=c(0,100),ylim=c(0,0.035),col=gray(0.9),breaks=breaks)
#lines(density(pnl2Var,bw=10),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl2Var),lwd=2)
points(median(pnl2Var),max(hist(pnl2Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

hist(pnl3Var,freq=F,main="Years listed as Critically Endangered",
         xlab="Years continously listed before extinction",xlim=c(0,100),ylim=c(0,0.035),col=gray(0.9),breaks=breaks)
#lines(density(pnl3Var,bw=10),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl3Var),lwd=2)
points(median(pnl3Var),max(hist(pnl3Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

dev.off()

graphics.off()

#### compute the proportion below 20 years

length(which(pnl1Var<=20))/length(pnl1Var)
length(which(pnl2Var<=20))/length(pnl1Var)
length(which(pnl3Var<=20))/length(pnl1Var)


#### for supplemental material: compute proportion of runs meeting or exceeding a range of thresholds

wt_thresholds <- c(5,10,15,20,25,30)       # thresholds for adequate warning times, in years
temp1 <- numeric(length(wt_thresholds))
temp2 <- numeric(length(wt_thresholds))
temp3 <- numeric(length(wt_thresholds))
for(i in 1:length(wt_thresholds)){
   temp1[i] <- length(which(pnl1Var<=wt_thresholds[i]))/length(pnl1Var)
   temp2[i] <- length(which(pnl2Var<=wt_thresholds[i]))/length(pnl2Var)
   temp3[i] <- length(which(pnl3Var<=wt_thresholds[i]))/length(pnl3Var)
}
tempdf <- data.frame(VU=temp1,EN=temp2,CR=temp3)
rownames(tempdf) <- paste(wt_thresholds," years",sep="")

setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\Results")
filename=paste("Sens_to_WTthreshold_listStatus",Sys.Date(),".csv",sep="")
write.csv(tempdf,row.names=T,file=filename)   # col.names=T,  ,sep=","



################################
###########   Histograms- break out by criteria

################################
######  Plot histograms

pnl1Var <- subset2$vu.timeto10A   # subset2$vu.timeto10A   # subset2$vu.timeto10A   #    
pnl2Var <- subset2$vu.timeto10B   # subset2$vu.timeto10B   # 
pnl3Var <- subset2$vu.timeto10C   # subset2$vu.timeto10C   # 
pnl4Var <- subset2$vu.timeto10D   # subset2$vu.timeto10D   # 
backVar <- subset2$vu.timeto10   # subset2$vu.timeto10   # 

pnl1Var[which(pnl1Var<=0)] <- 0
pnl2Var[which(pnl2Var<=0)] <- 0
pnl3Var[which(pnl3Var<=0)] <- 0
pnl4Var[which(pnl4Var<=0)] <- 0
backVar[which(backVar<=0)] <- 0

pnl1Var[which(pnl1Var>=100)] <- 100
pnl2Var[which(pnl2Var>=100)] <- 100
pnl3Var[which(pnl3Var>=100)] <- 100
pnl4Var[which(pnl4Var>=100)] <- 100
backVar[which(backVar>=100)] <- 100

pnl1_n <- length(pnl1Var)
pnl2_n <- length(pnl2Var)
pnl3_n <- length(pnl3Var)
pnl4_n <- length(pnl4Var)

pnl1_n
pnl2_n
pnl3_n
pnl4_n


   #head(RLMaster,100)


graphics.off()

setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\RawSVG")

filename=paste("YearsToExt_Histograms_byCriteria_annual",Sys.Date(),".svg",sep="") 

svg(filename=filename,width=4,height=9,onefile=TRUE)

par(mfrow=c(4,1))

hist(pnl1Var,freq=F,main="Years listed as threatened under Criterion A",xlab="",
                 xlim=c(0,100),ylim=c(0,0.04),col=gray(0.9),breaks=breaks)
#lines(density(backVar,bw=5),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl1Var),lwd=2)
points(median(pnl1Var),max(hist(pnl1Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

hist(pnl2Var,freq=F,main="Years listed as threatened under Criterion B",xlab="",
                  xlim=c(0,100),ylim=c(0,0.04),col=gray(0.9),breaks=breaks)
#lines(density(pnl2Var,bw=10),lwd=1.5,col=gray(0.2))
#lines(density(backVar,bw=5),lwd=1.5,col=gray(0.2))
polygon(x=c(density(backVar,bw=5)$x,rev(density(backVar,bw=5)$x)),
		y=c(density(backVar,bw=5)$y,rep(0,times=length(density(backVar,bw=5)$y))),
		lwd=1.5,col=gray(0.2))
#abline(v=median(pnl2Var),lwd=2)
points(median(pnl2Var),max(hist(pnl2Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

hist(pnl3Var,freq=F,main="Years listed as threatened under Criterion C",
         xlab="",xlim=c(0,100),ylim=c(0,0.04),col=gray(0.9),breaks=breaks)
#lines(density(pnl3Var,bw=10),lwd=1.5,col=gray(0.2))
#lines(density(backVar,bw=5),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl3Var),lwd=2)
points(median(pnl3Var),max(hist(pnl3Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

hist(pnl4Var,freq=F,main="Years listed as threatened under Criterion D",
         xlab="Years continously listed before extinction",xlim=c(0,100),ylim=c(0,0.04),col=gray(0.9),breaks=breaks)
#lines(density(pnl4Var,bw=10),lwd=1.5,col=gray(0.2))
#lines(density(backVar,bw=10),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl4Var),lwd=2)
points(median(pnl4Var),max(hist(pnl4Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)


dev.off()

graphics.off()

#### compute the proportion below 20 years

length(which(pnl1Var<=20))/length(pnl1Var)
length(which(pnl2Var<=20))/length(pnl2Var)
length(which(pnl3Var<=20))/length(pnl3Var)
length(which(pnl4Var<=20))/length(pnl4Var)


###############################
###    Criteria A,B, A,C, A,D, B,C, B,D, and C,D    (re-craft the process trajectories set first...)

################################
###########   Histograms- break out by criteria

################################
######  Plot histograms

pnl1Var <- subset2$vu.timeto10AB       
pnl2Var <- subset2$vu.timeto10AC   
pnl3Var <- subset2$vu.timeto10AD    
pnl4Var <- subset2$vu.timeto10BC   
pnl5Var <- subset2$vu.timeto10BD  
pnl6Var <- subset2$vu.timeto10CD   
backVar <- subset2$vu.timeto10  

pnl1Var[which(pnl1Var<=0)] <- 0
pnl2Var[which(pnl2Var<=0)] <- 0
pnl3Var[which(pnl3Var<=0)] <- 0
pnl4Var[which(pnl4Var<=0)] <- 0
pnl5Var[which(pnl5Var<=0)] <- 0
pnl6Var[which(pnl6Var<=0)] <- 0
backVar[which(backVar<=0)] <- 0

pnl1Var[which(pnl1Var>=100)] <- 100
pnl2Var[which(pnl2Var>=100)] <- 100
pnl3Var[which(pnl3Var>=100)] <- 100
pnl4Var[which(pnl4Var>=100)] <- 100
pnl5Var[which(pnl5Var>=100)] <- 100
pnl6Var[which(pnl6Var>=100)] <- 100
backVar[which(backVar>=100)] <- 100

pnl1_n <- length(pnl1Var)
pnl2_n <- length(pnl2Var)
pnl3_n <- length(pnl3Var)
pnl4_n <- length(pnl4Var)
pnl5_n <- length(pnl5Var)
pnl6_n <- length(pnl6Var)

pnl1_n
pnl2_n
pnl3_n
pnl4_n
pnl5_n
pnl6_n

   #head(RLMaster,100)


graphics.off()

setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\RawSVG")

filename=paste("YearsToExt_Histograms_byTwoCriteria_10year",Sys.Date(),".svg",sep="") 

svg(filename=filename,width=8,height=9,onefile=TRUE)

par(mfrow=c(3,2))

hist(pnl1Var,freq=F,main="Years listed as threatened \nunder Criteria A,B",xlab="",
                 xlim=c(0,100),ylim=c(0,0.04),col=gray(0.9),breaks=breaks)
#lines(density(backVar,bw=5),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl1Var),lwd=2)
points(median(pnl1Var),max(hist(pnl1Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

hist(pnl2Var,freq=F,main="Years listed as threatened \nunder Criteria A,C",xlab="",
                  xlim=c(0,100),ylim=c(0,0.04),col=gray(0.9),breaks=breaks)
#lines(density(pnl2Var,bw=10),lwd=1.5,col=gray(0.2))
#lines(density(backVar,bw=5),lwd=1.5,col=gray(0.2))
polygon(x=c(density(backVar,bw=5)$x,rev(density(backVar,bw=5)$x)),
		y=c(density(backVar,bw=5)$y,rep(0,times=length(density(backVar,bw=5)$y))),
		lwd=1.5,col=gray(0.2))
#abline(v=median(pnl2Var),lwd=2)
points(median(pnl2Var),max(hist(pnl2Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

hist(pnl3Var,freq=F,main="Years listed as threatened \nunder Criteria A,D",
         xlab="",xlim=c(0,100),ylim=c(0,0.04),col=gray(0.9),breaks=breaks)
#lines(density(pnl3Var,bw=10),lwd=1.5,col=gray(0.2))
#lines(density(backVar,bw=5),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl3Var),lwd=2)
points(median(pnl3Var),max(hist(pnl3Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

hist(pnl4Var,freq=F,main="Years listed as threatened \nunder Criteria B,C",
         xlab="Years continously listed before extinction",xlim=c(0,100),ylim=c(0,0.04),col=gray(0.9),breaks=breaks)
#lines(density(pnl4Var,bw=10),lwd=1.5,col=gray(0.2))
#lines(density(backVar,bw=10),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl4Var),lwd=2)
points(median(pnl4Var),max(hist(pnl4Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

hist(pnl5Var,freq=F,main="Years listed as threatened \nunder Criteria B,D",
         xlab="Years continously listed before extinction",xlim=c(0,100),ylim=c(0,0.04),col=gray(0.9),breaks=breaks)
#lines(density(pnl4Var,bw=10),lwd=1.5,col=gray(0.2))
#lines(density(backVar,bw=10),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl5Var),lwd=2)
points(median(pnl5Var),max(hist(pnl5Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

hist(pnl6Var,freq=F,main="Years listed as threatened \nunder Criteria C,D",
         xlab="Years continously listed before extinction",xlim=c(0,100),ylim=c(0,0.04),col=gray(0.9),breaks=breaks)
#lines(density(pnl6Var,bw=10),lwd=1.5,col=gray(0.2))
#lines(density(backVar,bw=10),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl6Var),lwd=2)
points(median(pnl6Var),max(hist(pnl6Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

dev.off()

graphics.off()

#### compute the proportion below 20 years

length(which(pnl1Var<=20))/length(pnl1Var)
length(which(pnl2Var<=20))/length(pnl2Var)
length(which(pnl3Var<=20))/length(pnl3Var)
length(which(pnl4Var<=20))/length(pnl4Var)
length(which(pnl5Var<=20))/length(pnl4Var)
length(which(pnl6Var<=20))/length(pnl4Var)





######################################
#########   Q: can we develop a way to use RF importance measures that are comparable across RL scenarios?
#                 that is, can we (e.g.) compare the RF importance of a variable (say, occ.area) across
#                 scenarios such as "only criterion A" vs: "all criteria considered"



###################################################################
#####################################################
#################################   INVESTIGATE FALSE ALARMS

###########################################
########   Histograms for 10 year reevaluation interval: separate for VU, EN, and CR

response <- "timein"
criteria <- ""  #  "A"  # 
category <- "cr"  #    "vu"  #   critically endangered 
eval.time <- "10" # """  # 5" #       
resp.name <- "Time continuously listed as Vulnerable \nbefore extinction"  #  "Time continuously listed as critical \nbefore extinction"  
dframe <-  "RLMaster.noex"   # "dfmain" #  "df_rel" # dfmain or df_rel for absolute or relative metrics       

full.response <- paste(category,".",response,eval.time,criteria,sep="")
full.respname <- ifelse(nchar(criteria)>0 ,paste(resp.name," under criterion ",criteria,sep=""),resp.name)

#left.censored <- paste(category,".","lcens",eval.time,sep="")

#eval(parse(text=paste("subset2 <- subset(",dframe,",",left.censored,"== 0",")",sep="") ))   ### don't left censor for false alarms!
#eval(parse(text=paste("subset2 <- subset(",dframe,",",left.censored,"== 0",")",sep="") ))   ### don't left censor for false alarms!

subset2 <- RLMaster.noex 

nrow(subset2)

###########
################################
######  Plot histograms

#pnl1Var <- eval(parse(text=paste("subset2$vu.timein",eval.time,criteria,sep="") ))
#pnl2Var <- eval(parse(text=paste("subset2$en.timein",eval.time,criteria,sep="") ))
#pnl3Var <- eval(parse(text=paste("subset2$cr.timein",eval.time,criteria,sep="") ))

pnl1Var <- subset2$vu.timein10 #- ifelse(subset2$vu.timeto10>=20,20,subset2$vu.timeto10)      # subtract the years continuously listed in the last 20 years- these may have gone extinct right afterwards... 
pnl2Var <- subset2$en.timein10 #- ifelse(subset2$en.timeto10>=20,20,subset2$en.timeto10)
pnl3Var <- subset2$cr.timein10 #- ifelse(subset2$cr.timeto10>=20,20,subset2$cr.timeto10)

pnl1Var[which(pnl1Var<=0)] <- 0
pnl2Var[which(pnl2Var<=0)] <- 0
pnl3Var[which(pnl3Var<=0)] <- 0

pnl1_n <- length(pnl1Var)
pnl2_n <- length(pnl2Var)
pnl3_n <- length(pnl3Var)

pnl1_n
pnl2_n
pnl3_n


##  cbind(subset2$vu.timein10,subset2$vu.timeto10)

  #head(RLMaster,100)


graphics.off()

setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\RawSVG")

filename=paste("VUvsENvsCR_Histograms_noex",eval.time,criteria,"_",Sys.Date(),".svg",sep="") 


svg(filename=filename,width=4,height=7,onefile=TRUE)

par(mfrow=c(3,1))

hist(pnl1Var,freq=F,main="Years listed as Vulnerable",xlab="",
                 xlim=c(0,110),ylim=c(0,0.04),col=gray(0.9),breaks=breaks)
#lines(density(pnl1Var,bw=10),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl1Var),lwd=2)
points(median(pnl1Var),max(hist(pnl1Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

hist(pnl2Var,freq=F,main="Years listed as Endangered",xlab="",
                  xlim=c(0,110),ylim=c(0,0.04),col=gray(0.9),breaks=breaks)
#lines(density(pnl2Var,bw=10),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl2Var),lwd=2)
points(median(pnl2Var),max(hist(pnl2Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

hist(pnl3Var,freq=F,main="Years listed as Critically Endangered",
         xlab="Total years listed",xlim=c(0,110),ylim=c(0,0.07),col=gray(0.9),breaks=breaks)
#lines(density(pnl3Var,bw=10),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl3Var),lwd=2)
points(median(pnl3Var),max(hist(pnl3Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

dev.off()

graphics.off()

#### compute the proportion above 50 years

length(which(pnl1Var>=50))/length(pnl1Var)
length(which(pnl2Var>=50))/length(pnl2Var)
length(which(pnl3Var>=50))/length(pnl3Var)


              ### across all runs, both extinct and not extinct, most species were listed for long periods as vulnerable... 

			  
#### for supplemental material: compute proportion of runs meeting or exceeding a range of false alarm thresholds

fa_thresholds <- c(35,40,45,50,55,60)       # thresholds for adequate warning times, in years
temp1 <- numeric(length(fa_thresholds))
temp2 <- numeric(length(fa_thresholds))
temp3 <- numeric(length(fa_thresholds))
for(i in 1:length(wt_thresholds)){
   temp1[i] <- length(which(pnl1Var>=fa_thresholds[i]))/length(pnl1Var)
   temp2[i] <- length(which(pnl2Var>=fa_thresholds[i]))/length(pnl2Var)
   temp3[i] <- length(which(pnl3Var>=fa_thresholds[i]))/length(pnl3Var)
}
tempdf <- data.frame(VU=temp1,EN=temp2,CR=temp3)
rownames(tempdf) <- paste(fa_thresholds," years",sep="")

setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\Results")
filename=paste("Sens_to_FAthreshold_listStatus",Sys.Date(),".csv",sep="")
write.csv(tempdf,row.names=T,file=filename)   # col.names=T,  ,sep=","




################################
###########   Histograms- break out by criteria

################################
######  Plot histograms

pnl1Var <- subset2$vu.timein10A
pnl2Var <- subset2$vu.timein10B
pnl3Var <- subset2$vu.timein10C
pnl4Var <- subset2$vu.timein10D
backVar <- subset2$vu.timein10   #subset2$vu.timeto10   #

pnl1Var[which(pnl1Var<=0)] <- 0
pnl2Var[which(pnl2Var<=0)] <- 0
pnl3Var[which(pnl3Var<=0)] <- 0
pnl4Var[which(pnl4Var<=0)] <- 0
backVar[which(backVar<=0)] <- 0


pnl1_n <- length(pnl1Var)
pnl2_n <- length(pnl2Var)
pnl3_n <- length(pnl3Var)
pnl4_n <- length(pnl4Var)

pnl1_n
pnl2_n
pnl3_n
pnl4_n


  # cbind(subset2$vu.timein10, pnl1Var,pnl2Var,pnl3Var,pnl4Var)   # checks out
  # cbind(subset2$en.timein10, pnl1Var,pnl2Var,pnl3Var,pnl4Var)   # checks out

   #head(RLMaster,100)


graphics.off()

setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\RawSVG")

filename=paste("YearsListedVU_noExtinct_Histograms_byCriteria",Sys.Date(),".svg",sep="") 

svg(filename=filename,width=4,height=9,onefile=TRUE)

par(mfrow=c(4,1))

hist(pnl1Var,freq=F,main="Years listed as endangered under Criterion A",xlab="",
                 xlim=c(0,100),ylim=c(0,max(hist(pnl1Var,plot=F)$density)*1.1),col=gray(0.9),breaks=breaks)
#lines(density(pnl1Var,bw=10),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl1Var),lwd=2)
polygon(x=c(density(backVar,bw=5)$x,rev(density(backVar,bw=5)$x)),
		y=c(density(backVar,bw=5)$y,rep(0,times=length(density(backVar,bw=5)$y))),
		lwd=1.5,col=gray(0.2))
points(median(pnl1Var),max(hist(pnl1Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

hist(pnl2Var,freq=F,main="Years listed as endangered under Criterion B",xlab="",
                  xlim=c(0,100),ylim=c(0,max(hist(pnl2Var,plot=F)$density)*1.1),col=gray(0.9),breaks=breaks)
#lines(density(pnl2Var,bw=10),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl2Var),lwd=2)
polygon(x=c(density(backVar,bw=5)$x,rev(density(backVar,bw=5)$x)),
		y=c(density(backVar,bw=5)$y,rep(0,times=length(density(backVar,bw=5)$y))),
		lwd=1.5,col=gray(0.2))
points(median(pnl2Var),max(hist(pnl2Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

hist(pnl3Var,freq=F,main="Years listed as endangered under Criterion C",xlab="",
         xlim=c(0,100),ylim=c(0,max(hist(pnl3Var,plot=F)$density)*1.1),col=gray(0.9),breaks=breaks)
#lines(density(pnl3Var,bw=10),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl3Var),lwd=2)
polygon(x=c(density(backVar,bw=5)$x,rev(density(backVar,bw=5)$x)),
		y=c(density(backVar,bw=5)$y,rep(0,times=length(density(backVar,bw=5)$y))),
		lwd=1.5,col=gray(0.2))
points(median(pnl3Var),max(hist(pnl3Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)

hist(pnl4Var,freq=F,main="Years listed as endangered under Criterion D",
         xlab="Total years listed",xlim=c(0,100),ylim=c(0,max(hist(pnl4Var,plot=F)$density)*1.1),col=gray(0.9),breaks=breaks)
#lines(density(pnl4Var,bw=10),lwd=1.5,col=gray(0.2))
#abline(v=median(pnl4Var),lwd=2)
polygon(x=c(density(backVar,bw=5)$x,rev(density(backVar,bw=5)$x)),
		y=c(density(backVar,bw=5)$y,rep(0,times=length(density(backVar,bw=5)$y))),
		lwd=1.5,col=gray(0.2))
points(median(pnl4Var),max(hist(pnl4Var,plot=F)$density),pch='|',cex=2)
#abline(v=20,lwd=2,lty=2)


dev.off()

graphics.off()

#### compute the proportion above 50 years

length(which(pnl1Var>=50))/length(pnl1Var)
length(which(pnl2Var>=50))/length(pnl2Var)
length(which(pnl3Var>=50))/length(pnl3Var)
length(which(pnl4Var>=50))/length(pnl4Var)

# 

#############################################
######  Barplot: x-axis: time listed    y-axis: time listed    bars: extinct vs. not extinct

sequence <- c(0,seq(10,90,10),110)
names <- as.character(sequence[-length(sequence)]+5)

RLMaster_comb$cuttimeVU <- cut(RLMaster_comb$vu.timein10,breaks=sequence)
RLMaster_comb$cuttimeEN <- cut(RLMaster_comb$en.timein10,breaks=sequence)
RLMaster_comb$cuttimeCR <- cut(RLMaster_comb$cr.timein10,breaks=sequence)
levels(RLMaster_comb$cuttimeVU) <- names
levels(RLMaster_comb$cuttimeEN) <- names
levels(RLMaster_comb$cuttimeCR) <- names

ext_df <- subset(RLMaster_comb,is.ext==T)
noext_df <- subset(RLMaster_comb,is.ext==F)

VUbarsEx <- table(ext_df$cuttimeVU)
VUbarsNoEx <- table(noext_df$cuttimeVU)
ENbarsEx <- table(ext_df$cuttimeEN)
ENbarsNoEx <- table(noext_df$cuttimeEN)
CRbarsEx <- table(ext_df$cuttimeCR)
CRbarsNoEx <- table(noext_df$cuttimeCR)

VUmatr <- rbind(VUbarsEx,VUbarsNoEx)
ENmatr <- rbind(ENbarsEx,ENbarsNoEx)
CRmatr <- rbind(CRbarsEx,CRbarsNoEx)

row.names(VUmatr) <- c("EX","NoEX")
row.names(ENmatr) <- c("EX","NoEX")
row.names(CRmatr) <- c("EX","NoEX")

graphics.off()
par(mfrow=c(3,1))


barplot(VUmatr,beside=T,legend.text=T,ylab="Number of reps",xlab="years listed, binned",main="VU")

barplot(ENmatr,beside=T,legend.text=T,ylab="Number of reps",xlab="years listed, binned",main="EN")

barplot(CRmatr,beside=T,legend.text=T,ylab="Number of reps",xlab="years listed, binned",main="CR")

  ## TODO:

#   investigate the reasons why some species were listed for long times but didn't go extinct, vs. species that were listed for 
#          long times and did go extinct.   [[what predicts the time wasted on conservation since those were not ultimately going to go extinct]]

#       [how many years did a non-extinct run spend listed as CR, VU, and EN.  How did this frequency shift over time?
#         KTS: were the non-extinct species listed as VU likely to move to EN- were there any trends?
#                 


#         e.g., filter only those species that were listed for 50 years or more, and see which variables predict
#                 whether the species went extinct or didn't go extinct. E.g., those that were listed as EN for 50 years,
#                 are expected to have a high risk of extinction- runs that didn't go extinct went against expectation- what explains that?



### why were so many "species" listed as vulnerable even though they didn't go extinct? What part of the Red List was malfunctioning? 




#####################################################
###############  DEVELOP TABLE FOR MANUSCRIPT

### for each red list category, compute the percentage of total listings that met each of the four criteria...

### use 10-year interval...


table <- matrix(nrow=4,ncol=8) 

names(RLMaster)

colnames(table) <- c("A","B","C","D","A","B","C","D")
rownames(table) <- c("CR","EN","VU","NT")

categories <- c("cr","en","vu","nt")
criteria <- c("A","B","C","D")

###  for extinct runs... 

subset2 <- subset(RLMaster,en.lcens10==0)
nrow(subset2)

i=1;j=2
for(i in 1:length(categories)){
  for(j in 1:length(criteria)){
    temp1 <- paste(categories[i],".timein10",sep="")
    temp2 <- paste(categories[i],".timein10",criteria[j],sep="")
    temp3 <- sum(eval(parse(text=paste("subset2$",temp1,sep=""))))     # denom
    temp4 <- sum(eval(parse(text=paste("subset2$",temp2,sep=""))))
    table[i,j] <- temp4/temp3
  }
}

### for non-extinct runs...

subset2 <- RLMaster.noex
nrow(subset2)

i=1;j=2
for(i in 1:length(categories)){
  for(j in 1:length(criteria)){
    temp1 <- paste(categories[i],".timein10",sep="")
    temp2 <- paste(categories[i],".timein10",criteria[j],sep="")
    temp3 <- sum(eval(parse(text=paste("subset2$",temp1,sep=""))))     # denom
    temp4 <- sum(eval(parse(text=paste("subset2$",temp2,sep=""))))
    table[i,(j+4)] <- temp4/temp3
  }
}


### Note: something seems strange with the "EN" category here... look into this.


############# WRITE OUT TO EXCEL

setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\Results") 

filename1 <- "Table1_raw"
filename <- paste(filename1,"_",Sys.Date(),".csv",sep="")

write.table(table,file=filename,row.names=T,col.names=T,sep=",")






