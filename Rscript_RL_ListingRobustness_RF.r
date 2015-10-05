##############################
###### R script for analyzing the NASA redlist data in terms of the robustness
######   of the red list categories. Specifically, this script evaluates which 
######   species traits might affect the number of years a species is listed as 
######   critical (e.g.) before it goes extinct, using Random Forest analy 


KEVIN = TRUE # FALSE # TRUE
KEVINLAB = FALSE # TRUE
JESSIE = FALSE

##########################
### load random forest functions

             # version controlled using GIT...
if(KEVIN) setwd("C:\\Users\\Kevin\\GIT\\Random-Forest-Functions")
if(KEVIN) source("RF_Extensions.R")


########################
### load predictor variables and response variables from all 'realized species' from the NASA project

if(KEVIN) setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\RData")
 # if(KEVINLAB) setwd("")
load("FinalNASAData.RData")

##########################
###   load model-specific redlist data

if(KEVIN) setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\RData")
if(KEVINLAB) setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\RData")
load("ModelSpecificRLData_2014-04-22.RData")   #("ModelSpecificRLData_2014-02-14.RData")

## ls()

###################################################
 ######            DEFINE RESPONSE AND PREDICTOR VARIABLES

predictors <- c(	"log.metapop.initab",
			"Gen.Tgen",
			"GrowthRt",
			"cv.avg",
			"NewCorrelationIndex",
			"log.patch.n.10",
			"log.tot.patch.area.10",
			"overall.frac.dim.10",
			"LargestPatchFrac",
			"ConnectivityIndex",
			"FragIndex",
			"Dispersability",
			"NicheBreadthT",
			"NicheBreadthP",
			"RecentAreaChange",
			"RecentNChange",
			"RecentSubpopChange",
			"RecentFracDimChange",
			"RecentConnectivity",
			"RecentFragmentation",
			"RecentLPFChange",
                  #"spec",                   #### Add species to list of predictors
                  "init.cat"                    #### add initial listing category to list of predictors...
)     

predictorNames <- c(	"Initial Abundance",
				"Generation Length",
				"Maximum Growth Rate",
				"Variability in Vital Rates",
				"Spatial Correlation of Variability",
				"Number of patches",
				"Total Occupied Area",
				"Overall Fractal Dimension",
				"Largest Patch Fraction",
                   	"Connectivity Index",
				"Fragmentation Index",
				"Innate Dispersal Ability",
				"Climate Breadth at t0, temp",
				"Climate Breadth at t0, precip",
				"Recent change in occupied area",
				"Recent change in abundance",
				"Recent change in number of subpops",
				"Recent change in fractal dimension",
				"Recent Change in Connectivity",
				"Recent fragmentation", 
				"Recent Change in Patch Dominance" ,
                        #"Species"  ,
                        "Initial Listing Category"
)  

cbind(predictors,predictorNames)


###################################
#######  DEVELOP MASTER DATA FRAME

keepers <- c("mp.run","clim.chg","model", predictors[-length(predictors)] ,"spec") # note: init.cat is not in this data frame

RLMaster <- subset(dfmain, subset=(match.stat=="Match" & clim.chg=='WRE'), select=keepers)   # subset only non-mismatches with climate scenario WRE
RLMaster$mp.run <- as.character(RLMaster$mp.run)

RLMaster_old <- RLMaster

#  table(RLMaster$mp.run)

#  head(RLMaster)
#  RLMaster$mp.run
#  nrow(RLMaster)
#  nrow(cattimes.ex)

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

  # nrow(RLMaster.noex)

  #predictor_cols <- which(names(RLMaster)%in%predictors)    # is this used at all?

RLMaster$model<-as.character(RLMaster$model)
RLMaster$spec <- as.character(RLMaster$spec)

RLMaster.noex$model<-as.character(RLMaster.noex$model)
RLMaster.noex$spec <- as.character(RLMaster.noex$spec)

 # RLMaster.noex$vu.lcens

#####################################
          # determine the extinction risk under the no climate change scenario
dfmain4 <- subset(dfmain, subset=(match.stat=="Match" & clim.chg=='NoCC'), select=c(keepers,"ext.risk"))
head(dfmain4)

dfmain4$mp.run <- as.character(dfmain4$mp.run)

dfmain4 <- merge(cattimes2.ex, dfmain4, by='mp.run')     ##  head(merge(cattimes.ex, RLMaster, by='mp.run'),100)
  #nrow(dfmain4)
  #nrow(cattimes2.ex)

head(RLMaster)
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
head(dfmain2)

dfmain2$mp.run <- as.character(dfmain2$mp.run)

dfmain2 <- merge(cattimes2.ex, dfmain2, by='mp.run')     ##  head(merge(cattimes.ex, RLMaster, by='mp.run'),100)
  #nrow(dfmain2)
  #nrow(cattimes2.ex)


RLMaster$ext.risk.wre <- dfmain2$ext.risk
#head(RLMaster)
#RLMaster$ext.risk.wre

     ## and for non-extinct runs

dfmain2 <- subset(dfmain, subset=(match.stat=="Match" & clim.chg=='WRE'), select=c(keepers,"ext.risk"))
head(dfmain2)

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



########## write test file for Resit

 # setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\RData")
 #  write.table(RLMaster,file="RLMaster.csv",sep=",",row.names=FALSE,col.names=TRUE)




########################
###  START ANALYSIS




########################
### EXPLORE DATA: Visualize predictor variables

graphics.off()
par(mfrow=c(2,2))
hist(subset1$Gen.Tgen)
hist(subset1$log.metapop.initab)
hist(subset1$GrowthRt)
hist(subset1$NewCorrelationIndex)

graphics.off()
par(mfrow=c(2,2))
hist(subset1$log.patch.n.10)
hist(subset1$log.tot.patch.area.10)
hist(subset1$LargestPatchFrac)
hist(subset1$ConnectivityIndex)

 
graphics.off()
par(mfrow=c(2,2))
hist(subset1$FragIndex)
hist(subset1$Dispersability)
hist(subset1$NicheBreadthT)
hist(subset1$NicheBreadthP)

graphics.off()
par(mfrow=c(2,2))
hist(subset1$RecentAreaChange)
hist(subset1$RecentNChange)
hist(subset1$RecentSubpopChange)
hist(subset1$RecentFracDimChange)

graphics.off()
par(mfrow=c(2,2))
hist(subset1$RecentAreaChange)
hist(subset1$RecentNChange)
hist(subset1$RecentSubpopChange)
hist(subset1$RecentFracDimChange)

#[a few not visualized]




######################################################
#######################################
##################  DEFINE MODEL AND DEVELOP FINAL DATA SET (subset1) 

#  possible generic models: Snake, Turtle, Small_Salamander, Large_Salamander, Tortoise




###############################################################
###############################################################
###############################################################
#####  What predicts the amount of time listed before extinction?


response <- "timeto" 
with.vars <-  "with_GLH but not spec"  #  "with_GLH_and_spec_only"  #"no_spec"  #"no_corr" #"spec_and_initCat"
criteria <- ""  #  "A"  # use "" for all criterial
category <- "vu"  # "cr"    #  critically endangered 
eval.time <- "10" # "5" # ""        
resp.name <- "Time continuously listed as Vulnerable \nbefore extinction"  #  "Time continuously listed as critical \nbefore extinction"  
binaryresponse <- FALSE  
sp.group <-  "All species"  #   "Small_Salamander"  # "Snake"    #  "Turtle" #   "Tortoise"  #  
species <- "all"  #  "FPSN"   #   "BLTU"   #   "AWSN"  #  "BGTU" #    "CTSA"  #  "DNSA" #  "JPSA" #  "OSSA"  #          
dframe <-  "RLMaster"   #  RLMaster.noex  #  
climate <- "WRE"  #  "LEV"  #   "NoCC"  #     
rtype <- "absolute" #  "relative"  #        

  # predictors <- c(predictors,"model")    #  predictors <- c("spec","model")
  # predictorNames <- c(predictorNames,"GLH Model") # c("Species","GLH Model")

predictors2 <- paste(c(predictors),collapse="+")  #  "IsolationIndex+log.metapop.initab+Gen.Tgen+GrowthRt+log.N.CV.10+avg.corr.dist.b+log.patch.n.10+log.tot.patch.area.10+overall.frac.dim.10+LargestPatchFrac"   #     # spec+   log.mean.t0.disp.rate+

full.response <- paste(category,".",response,eval.time,criteria,sep="")
left.censored <- paste(category,".","lcens",eval.time,sep="")
full.respname <- ifelse(nchar(criteria)>0 ,paste(resp.name," under criterion ",criteria,sep=""),resp.name)

if(sp.group=="All species"){
  sp.group <- c("Snake", "Turtle", "Small_Salamander", "Large_Salamander", "Tortoise", "Lizard")
  eval(parse(text=paste("subset1 <- subset(",dframe,",model%in%","sp.group","&",left.censored,"== 0",")",sep="") ))
} else{
  eval(parse(text=paste("subset1 <- subset(",dframe,",(model==\"",sp.group,"\"&",left.censored,"== 0))",sep="") ))
}

nrow(subset1)

#if(species!="all"){
#  eval(parse(text=paste("subset1 <- subset1[which(subset1$spec==\"",species,"\"&match.stat==\"Match\"),]",sep="") ))
#}

sp.group2 <- "all species"
if(length(sp.group)==1&species=="all"){
  title <- paste(resp.name,", all ",sp.group," species",sep="")
} else{
 title <- paste(resp.name,", ",rtype,", ",sp.group2,", ",climate,sep="")
}

if(species=="all"&length(sp.group)>1) title <- paste(full.respname,", ",rtype,", ",sp.group2,", ",climate,sep="")

formula <- eval(parse(text=paste("as.formula(",full.response,"~",predictors2,")",sep="")))



###############################################################
###############################################################
###############################################################
#####  ALTERNATIVE: What predicts the time listed for non-extinct runs?

response <- "timein" 
criteria <- ""  #  "A"  # use "" for all criterial
category <- "en"  # "vu"  # "cr"    #  critically endangered 
eval.time <- "10" # "5" # ""        
resp.name <- "Total time listed as Endangered \nbefore extinction"  #  "Time continuously listed as critical \nbefore extinction"  
binaryresponse <- FALSE  
sp.group <-  "All species"  #   "Small_Salamander"  # "Snake"    #  "Turtle" #   "Tortoise"  #  
species <- "all"  #  "FPSN"   #   "BLTU"   #   "AWSN"  #  "BGTU" #    "CTSA"  #  "DNSA" #  "JPSA" #  "OSSA"  #          
dframe <-  "RLMaster.noex"  # "RLMaster"   #   
climate <- "WRE"  #  "LEV"  #   "NoCC"  #     
rtype <- "absolute" #  "relative"  #        

predictors2 <- paste(c(predictors),collapse="+")  #  "IsolationIndex+log.metapop.initab+Gen.Tgen+GrowthRt+log.N.CV.10+avg.corr.dist.b+log.patch.n.10+log.tot.patch.area.10+overall.frac.dim.10+LargestPatchFrac"   #     # spec+   log.mean.t0.disp.rate+

full.response <- paste(category,".",response,eval.time,criteria,sep="")
left.censored <- paste(category,".","lcens",eval.time,sep="")
full.respname <- ifelse(nchar(criteria)>0 ,paste(resp.name," under criterion ",criteria,sep=""),resp.name)

if(sp.group=="All species"){
  sp.group <- c("Snake", "Turtle", "Small_Salamander", "Large_Salamander", "Tortoise", "Lizard")
  eval(parse(text=paste("subset1 <- subset(",dframe,",model%in%","sp.group","&",left.censored,"== 0",")",sep="") ))
} else{
  eval(parse(text=paste("subset1 <- subset(",dframe,",(model==\"",sp.group,"\"&",left.censored,"== 0))",sep="") ))
}

nrow(subset1)

#if(species!="all"){
#  eval(parse(text=paste("subset1 <- subset1[which(subset1$spec==\"",species,"\"&match.stat==\"Match\"),]",sep="") ))
#}

sp.group2 <- "all species"
if(length(sp.group)==1&species=="all"){
  title <- paste(resp.name,", all ",sp.group," species",sep="")
} else{
 title <- paste(resp.name,", ",rtype,", ",sp.group2,", ",climate,sep="")
}

if(species=="all"&length(sp.group)>1) title <- paste(full.respname,", ",rtype,", ",sp.group2,", ",climate,sep="")

formula <- eval(parse(text=paste(full.response,"~",predictors2,sep="")))




#############################################################
#######   RANDOM FOREST ANALYSIS 



###################################
########  CONDITIONAL INFERENCE TREE ANALYSIS

subset1$spec <- as.factor(subset1$spec)
subset1$model <- as.factor(subset1$model)
con <- ctree_control(maxdepth = 4)
treeobj <- ctree(formula, data=subset1,control=con)
graphics.off()

filename=paste("RL_CITree_",with.vars,"_",dframe,"_",full.response,"_",Sys.Date(),".svg",sep="")
setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\RawSVG")
svg(file=filename,width=8,height=6,onefile=TRUE)

plot(treeobj,main=full.respname,type="simple" , inner_panel=node_inner(treeobj, digits = 2, abbreviate = FALSE, 
  fill = "white", pval = FALSE, id = FALSE),terminal_panel=node_terminal(treeobj, digits = 2, abbreviate = FALSE, 
  fill = c("lightgray", "white"), id = FALSE))    #  #   # type="extended"

dev.off()


######################################
###############  RANDOM FOREST: use conditional inference trees   (should try subsampling etc...)
##note that these RFs use the specifications in the EXPLORATORY REGRESSION TREES block, above

  #set.seed(47) 
data.controls <- cforest_unbiased(ntree=1000, mtry=2)   # ntree=5000, mtry=3 #SELECT RF PARAMETERS HERE: mtry is the number of randomly 
									#selected variables; ntree is the number of trees.

               # find a fraction that generally yields independent observations
 #data.controls@fraction <- find_fraction(subset1)
data.controls@fraction <- 0.3

data.cforest <- cforest(formula, data = subset1, controls=data.controls)  # controls=data.controls

   ### NOTE: use "getAnywhere(RandomForest)" to find the source code
      ### type 'party:::RandomForest' to access the source code (but also note that much of the source code
          ### for these analyses is in C++ and R just calls pre-compiled routines.)


data.cforest.varimp <- varimp(data.cforest, conditional = FALSE)   # conditional = TRUE produces memory limitation error



  ### plot importance values...   [[if possible probably best to use conditional- takes intercorrelations into account ]]
graphics.off()

filename=paste("RL_","VarImp_",with.vars,"_",dframe,"_",full.response,"_",Sys.Date(),".svg",sep="")
setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\RawSVG")
svg(file=filename,width=6,height=5,onefile=TRUE)

lengthndx <- length(data.cforest.varimp)
par(mai=c(0.95,3.1,0.6,0.4))
col <- rainbow(lengthndx, start = 3/6, end = 4/6)      # rep(brewer.pal(6,"Blues"),each=2)
barplot(height=data.cforest.varimp[order(data.cforest.varimp,decreasing = FALSE)],
            horiz=T,las=1,main=paste(resp.name, ", ",climate, sep=""),
            xlab="Index of overall importance",col=col,           
            names.arg=predictorNames[match(names(data.cforest.varimp),predictors)][order(data.cforest.varimp,decreasing = FALSE)])


dev.off()

graphics.off()  ##run this to reset graphics devise AFTER running the above plot##


###################################
#################### RANDOM FOREST: develop performance metrics
            ##################
                      #  EXTRACT PERFORMANCE METRIC(S)...

   ### use cross validation algorithm...

n.folds = length(unique(subset1$spec))  #  10 # 10  #   #  to make leave-one-out, choose number of species: n.folds = length(unique(df_rel$spec))

foldVector <- make_foldvec(n.folds=n.folds,foldvar=as.character(subset1$spec)) #to run validtion by species (best)
  #foldVector <- rep(c(1:10),length=(nrow(subset1)))[sample(nrow(subset1))]   # to run random validation

counter = 1
CVprediction <- numeric(nrow(subset1))
CVobserved <- numeric(nrow(subset1))
realprediction <- numeric(nrow(subset1))
#test <- numeric(nrow(subset1))
i=0
i=i+1
for(i in 1:n.folds){
   #data.controls@fraction <- 0.64  # find_fraction(subset1[which(foldVector!=i),])
  model <- cforest(formula, data = subset1[which(foldVector!=i),], controls=data.controls)
  CVprediction[counter:(counter+(length(which(foldVector==i))-1))] <- predict(model,newdata=subset1[which(foldVector==i),])
  CVobserved[counter:(counter+(length(which(foldVector==i))-1))] <- eval(parse(text=paste("subset1$",full.response,"[which(foldVector==i)]",sep="")))  #subset1$ext.risk[which(foldVector==i)]
  #test[counter:(counter+(length(which(foldVector==i))-1))] <- i
  realprediction[counter:(counter+(length(which(foldVector==i))-1))] <- predict(data.cforest,newdata=subset1[which(foldVector==i),]) 
  counter = counter + length(which(foldVector==i))
}

   #realdata <- eval(parse(text=paste("subset1$",response,sep="")))
   #realprediction <- data.cforest@predict_response()
CV_RMSE = sqrt(mean((CVobserved-CVprediction)^2))       # root mean squared error for holdout samples in 10-fold cross-validation ...
real_RMSE = sqrt(mean((CVobserved-realprediction)^2))  # root mean squared error for residuals from final model

       # print RMSE statistics
CV_RMSE 
real_RMSE

graphics.off()
par(mfrow=c(1,2))
plot(CVprediction~CVobserved,xlab="Observed",ylab="Predicted, cross-validation")
plot(realprediction~CVobserved,xlab="Observed",ylab="Predicted, full model")


fit_deviance_CV <- mean((CVobserved-CVprediction)^2)
if(binaryresponse) fit_deviance_CV <- mean(-2*(dbinom(CVobserved,1,CVprediction,log=T)-dbinom(realdata,1,realdata,log=T)))
fit_deviance_real <- mean((CVobserved-realprediction)^2)
if(binaryresponse) fit_deviance_real <- mean(-2*(dbinom(CVobserved,1,realprediction,log=T)-dbinom(realdata,1,realdata,log=T)))
null_deviance <- mean((CVobserved-mean(CVobserved))^2)
if(binaryresponse) null_deviance <- mean(-2*(dbinom(CVobserved,1,mean(CVobserved),log=T)-dbinom(realdata,1,realdata,log=T)))
deviance_explained_CV <- (null_deviance-fit_deviance_CV)/null_deviance   # based on holdout samples
deviance_explained_real <- (null_deviance-fit_deviance_real)/null_deviance   # based on full model...

deviance_explained_CV
deviance_explained_real

##################################
##################### PRINT PERFORMANCE METRICS

filename=paste("Performance_stats_leaveSpecOut_",with.vars,"_",full.response,Sys.Date(),".csv",sep="")
setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\Results")
temp <- data.frame(Rsquared_CV = deviance_explained_CV,Rsquared_real = deviance_explained_real,
                     RMSE_CV = CV_RMSE, RMSE_real = real_RMSE )
write.table(temp,file=filename,sep=",",row.names=F)




			### RESIDUAL PLOTS 
                   ###########   eval(parse(text=paste(workingModel,"$fitted",sep="")))     # gbm.fracChange.tc4.lr01$fitted

           # compare to original response variable
  #cbind(eval(,eval(parse(text=paste("subset1$",response,sep=""))) )

           # residuals
residuals <- CVobserved - realprediction 
graphics.off()
hist(residuals,freq=F,breaks=25)
curve(dnorm(x,mean(residuals),sd(residuals)),add=T)

qqnorm(residuals)



#############################  RANDOM FOREST
#############################  Display univariate plots
  
 	 #The call should now follow this syntax:
#RF_UnivariatePlots(object=<<random forest object>>, varimp=<<importance object from random forest>>, data=<<main data frame>>,
	#predictors=<<names of predictor variables>>, labels=<<readable predictor names>>, allpredictors=<<names of predictor variables>>,plot.layout=<<number of rows and columns to plot in each graphics window>>)
	 
graphics.off()

setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\RawSVG")

filename=paste("UnivariatePlots_top9_",with.vars,"_",full.response,Sys.Date(),".svg",sep="")
svg(filename=filename,width=7,height=6,onefile=TRUE)
 #pdf(file=filename,width=3,height=3,onefile=TRUE)


predictors_varimp <- predictors[order(data.cforest.varimp,decreasing=T)]
predictorNames_varimp <- predictorNames[order(data.cforest.varimp,decreasing=T)]
RF_UnivariatePlots(object=data.cforest, varimp=data.cforest.varimp, data=subset1, 
                   predictors=predictors_varimp[1:9],labels=predictorNames_varimp[1:9],allpredictors=predictors_varimp,plot.layout=c(3,3),plot=T)  # predictors




dev.off()
 #dev.off()  
graphics.off()



#### look at species effect more closely

graphics.off()
par(las=3)
plot(en.timein10~spec,data=subset1)     # might be some interesting interactions here?


graphics.off()
par(las=0)
plot(en.timein10~as.factor(model),data=subset1)     # might be some interesting interactions here?




####################################
#######################   RANDOM FOREST FIND AND PLOT INTERACTIONS

          # NOTE: this one can take a very long time- maybe up to 2 hours...
# rf_findint <- RF_FindInteractions(object=data.cforest,data=subset1,predictors=predictors)

   # display and plot out interactions...
rf_findint$interactions1

rf_findint$rank.list1

  ### plot interaction strength
graphics.off()
lengthndx <- min(7,nrow(rf_findint$rank.list))
par(mai=c(0.95,3.1,0.6,0.4))
barplot(height=(rf_findint$rank.list1[c(1:min(9,nrow(rf_findint$rank.list))),5][c(lengthndx:1)]),
            horiz=T,las=1,main=paste(full.response, sep=""),
            xlab="Index of interaction strength",col=brewer.pal(lengthndx,"Blues"),           
            names.arg=paste("",predictorNames[match(rf_findint$rank.list1[,2][c(lengthndx:1)],predictors)],"\n",predictorNames[match(rf_findint$rank.list1[,4][c(lengthndx:1)],predictors)],sep="") )

graphics.off()

rf_findint$rank.list1


      ########### New code (Feb 2013) to output 3d plots as vector, for publication quality:

graphics.off()

#svg(filename = "Fig3_RF_reference_rank4.svg",
#    width = 6, height = 6, pointsize = 10,
#    onefile = TRUE, family = "sans", bg = "white")

filename="InteractionPlot_part2.pdf" # 
 #svg(filename=filename,width=3,height=3,onefile=TRUE)
pdf(file=filename,width=6,height=6,onefile=TRUE,pointsize=8)

par(mai=c(1,1.5,0,1))

fam <- "gaussian"  #  "bernoulli"  #

#RF_InteractionPlots(6,13,object=data.cforest,data=subset1,predictors=predictors,family=fam) #version that does not set z axis scale
RF_InteractionPlots(20,2,object=data.cforest,data=subset1,predictors=predictors,family=fam)  # zlim=c(0,0.5) #version for setting z axis scale (zlim parameters)

dev.off()

graphics.off()

###END OF FIGURE 3 CODE






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






