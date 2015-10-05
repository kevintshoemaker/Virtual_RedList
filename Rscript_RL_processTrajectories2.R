#######################################
####  Red List analysis for the NASA project
####     This script reads raw time-specific RL criteria/status data (in turn compiled from MP files in a separate script)
####     to compile, calculate and collate the data in a form usable for futher visualization and analysis

####  J. Stanton and K. Shoemaker Jan 2014
####
#### KTS 29 Jan: code now accommodates 5 and 10 year survey intervals
#### KTS 11 April: main processes are now flexible functions- code is more streamlined t
####       
#######################################

#######################################
#########################
#### SET GLOBAL VARIABLES

KEVIN = TRUE
JESSIE = FALSE


## Create numerical column for assessment category- for reclassification
cat.lookup= data.frame( 
  cat.num=c(5,4,3,2,1,0), 
  AssessmentCategory=c("LC","NT","VU","EN","CR","EX"))


eval.intervals = c(1,5,10)  # species are reevaluated after 1, 5, or 10 year intervals.

last.year = 111

           # columns to keep for each raw spreadsheet for red-listing
RLColNames <- c("Model", "FileName", "Replicate", "Year", "AssessmentCategory", 
"Crit.A2", "Crit.B1", "Crit.B2", "Crit.C", "Crit.D", "CRA", "CRB", 
"CRC", "CRD", "ENA", "ENB", "ENC", "END", "VUA", "VUB", "VUC", 
"VUD", "NTA", "NTB", "NTC", "NTD")


burnin.year <- 10

final.year_noex <- 90   # last year to be considered for non-extinct scenarios...

          # different categories of interest
categories <- c("cr","en","vu","nt")

          # criteria of interest
criteria <- c("","A","B","C","D","AB","AC","AD","BC","BD","CD","ABD","ACD","BCD")    # note: "" signified all criteria


intervals <- c("","5","10")

   ### listing will be performed for all intervals, categories, and criteria


#####################################################
#####    INITIALIZE STORAGE STRUCTURES 

catcols <- c('cat','tmstp','count')
cattable.risk <- data.frame()
cattable.lowrisk <- data.frame()
cattimes.ex <- data.frame()     	# times in each category (at least) before extinction
cattimes.noex <- data.frame()   	# times in each category (at least) before sim end
ex.reps<-data.frame()           	# all info for extinction runs
non.ex.reps <- data.frame()     	# all info for non-extinct runs

TimeInCategories <- data.frame(temp = rep(0,times=100000))     # Analysis Resit wants to run: less time spent in more severe categories??

CategoryMatrices <- list()  # new analysis: look at the transition rate between and among categories

##########################
### SET DIRECTORY WHERE RAW SPREADSHEETS ARE

     # if(JESSIEH) setwd('G:/NASA_REDLIST') #working at home
if(JESSIE) setwd('F:/NASA_REDLIST') #working at school
if(KEVIN) setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres")

 # getwd()


#########################
###  read in names of files that store raw information for redlisting

        ## populate list of CSV filesrzl.files <-list.files(,pattern='*new_newRL.csv')
rl.files <-list.files(,pattern='*new_newRL.csv')


#########################
###  PROCESS GLOBAL VARIABLES

tempdf <- read.csv("ADTO_new_newRL.csv",h=T)
keep.cols <- which(names(tempdf)%in%RLColNames)    #c(1,2,4,5,6,7,8,9,10,11,12:27)  # columns of redlist files to keep 



##########################
### LOAD DEPENDENCIES

suppressWarnings(suppressMessages(library(ggplot2)))   # suppressMessages(library(foo))
suppressWarnings(suppressMessages(library(zoo)))
suppressWarnings(suppressMessages(library(plyr)))
suppressWarnings(suppressMessages(library(reshape)))
suppressWarnings(suppressMessages(library(dismo)))


##########################
###  LOAD FUNCTIONS

 ## dataFile=rl.mod; replicate=2; model=models[3]

processTrajectory <- function(dataFile,replicate,model){

       ### KTS: this function performs all red-listing tasks for a given replicate.

  rl.rep <- subset(dataFile, Replicate==replicate)      # data frame: only data for this replicate for this SRPM for this species
  rl.rep[,11:26] <- apply(rl.rep[,11:26],2,as.logical)           # process data: make sure the single criteria columns are all logical


	      ##### set up the re-evaluation intervals
  i=1
  for(i in 1:length(intervals)){
    tempname <- paste("evalyears",eval.intervals[i],sep="")
      #eval(parse(text=paste(tempname,"<-numeric(0)",sep="")))
    samp <- sample(x=rep(c(burnin.year:(burnin.year+eval.intervals[i]-1)),each=2),size=1)
    sequence <- seq(from=samp,to=last.year,by=eval.intervals[i])
    eval(parse(text=paste(tempname,"<- sequence",sep="")))
    if(eval.intervals[i]==10)   evalyears10_alt <- seq(from=(samp+5),to=last.year,by=eval.intervals[i])     
  }

      	#########  process data: make sure all years are represented in the data frame even if all reps went extinct      
  if (max(rl.rep$Year)<last.year){ 				# If years end because all replicates go extinct before end
    n.missing <- last.year-max(rl.rep$Year)
    rep.missing <- c(max(rl.rep$Year)+1:n.missing)
    add.ex <- data.frame(Model= rep(model,n.missing), FileName = rep(rl.rep[1,2],n.missing),
                           Replicate = rep(n,n.missing),Year = rep.missing, AssessmentCategory = rep('EX',n.missing), 
                           Crit.A2 = rep(FALSE,n.missing), Crit.B1 = rep(FALSE,n.missing),
                           Crit.B2 = rep(FALSE,n.missing), Crit.C = rep(FALSE,n.missing),
                           Crit.D = rep(FALSE,n.missing), 
                           CRA = rep(FALSE,n.missing),CRB = rep(FALSE,n.missing),CRC = rep(FALSE,n.missing),CRD = rep(FALSE,n.missing),
                           ENA = rep(FALSE,n.missing),ENB = rep(FALSE,n.missing),ENC = rep(FALSE,n.missing),END = rep(FALSE,n.missing),
                           VUA = rep(FALSE,n.missing),VUB = rep(FALSE,n.missing),VUC = rep(FALSE,n.missing),VUD = rep(FALSE,n.missing),
                           NTA = rep(FALSE,n.missing),NTB = rep(FALSE,n.missing),NTC = rep(FALSE,n.missing),NTD = rep(FALSE,n.missing))
    rl.rep <- rbind(rl.rep,add.ex) 
  } 	#end if too few years
     
      	# Add numerical column for assessment cat 
  rl.rep <- join(rl.rep, cat.lookup, by='AssessmentCategory') 	# a left join, in sql terminology: add numeric representation of category

      
      	# head(rl.rep)								# Drop first 9 replicate years (burn-in from NASA project)
  rl.rep <- subset(rl.rep, Year>=burnin.year) 

##################################################################
##############   IF EXTINCT
##################################################################      
  if(rl.rep$AssessmentCategory[101]=='EX'){         
    ex.status <- TRUE    # flag for extinction...
    ex.reps <<- rbind(ex.reps,rl.rep) #compile extinction replicate runs       # add to extinction club
    ex.year <- min(subset(rl.rep, select=Year, subset=(cat.num<1))) # Extinction year (year of simulation, not years after burnin)
           # set up the data frame for storing all information from this replicate...this is what will be returned from this function.
    cat.rep.ex <- data.frame(spec=rl.rep[1,2],init.cat=rl.rep$cat.num[1],model = models[mod], rep=n, ex.year = ex.year)
    masterFileName <- "cat.rep.ex"
    yearsin_years <- c(burnin.year:last.year)

                  # SET UP CATEGORY TRANSITION MATRIX
        ###  KTS: add new transition matrix structure:
    CategoryMatrices[[CMCounter]] <<- list()       ########<<-#########
    CategoryMatrices[[CMCounter]]$spec <<- rl.rep[1,2] 
    CategoryMatrices[[CMCounter]]$model <<- models[mod]
    CategoryMatrices[[CMCounter]]$ex.year <<- ex.year             # model = models[mod], rep=n, ex.year = ex.year
    CategoryMatrices[[CMCounter]]$MAT <<- matrix(NA,nrow=6,ncol=6)
    from <- 6-rl.rep$cat.num[1:(length(rl.rep$cat.num)-1)]
    to <- 6-rl.rep$cat.num[2:(length(rl.rep$cat.num))]
    for(cgy1 in min(6-rl.rep$cat.num):6){    # loop through categories and compute transitions
      for(cgy2 in 1:6){
        tmp <- length(which((from==cgy1)&(to==cgy2)))
        CategoryMatrices[[CMCounter]]$MAT[cgy2,cgy1] <<- tmp
      }
    }
 
  } else{       # IF NOT EXTINCT
    ex.status <- FALSE
    non.ex.reps <<- rbind(non.ex.reps,rl.rep) #compile extinction replicate runs
    ex.year <- last.year
    cat.rep.noex <- data.frame(spec=rl.rep[1,2],init.cat=rl.rep$cat.num[1],model = models[mod], rep=n, ex.year = NA)
    masterFileName <- "cat.rep.noex"
    yearsin_years <- c(burnin.year:final.year_noex)
  }        
########################################
             # LOOP THROUGH CATEGORIES
   #  c=1
    for(c in 1:length(categories)){
      cat.num2 <- cat.lookup$cat.num[which(cat.lookup$AssessmentCategory==toupper(categories[c]))]

      if(any(rl.rep$cat.num<(cat.num2+1)&rl.rep$cat.num>=1)){          # check whether there are any listings at this level...
##################################### 
             # LOOP THROUGH INTERVALS
         # i=1
        for(i in 1:length(intervals)){
          eval.yrs2 <- eval(parse(text=paste("evalyears",eval.intervals[i],sep="")))  # select evaluation years...
          if(eval.intervals[i]==10) eval.yrs2_alt <- evalyears10_alt


#####################################
             # LOOP THROUGH LISTING CRITERIA
          l=1
          for(l in 1:length(criteria)){
            if(nchar(criteria[l])==0){
              listed.yrs2 <- subset(rl.rep,select=Year, subset=((cat.num<(cat.num2+1))&(cat.num>=1)))[,1]    # store years listed
            } else if(nchar(criteria[l])==1){
              temp1 <- paste(toupper(categories[c]),criteria[l],sep="")
              cond1 <- paste("(",temp1,"==T)",sep="")
              listed.yrs2 <- eval(parse(text=paste("subset(rl.rep,select=Year, subset=",cond1,")[,1]",sep="")))
            } else if(nchar(criteria[l])==2){
              temp1 <- paste(toupper(categories[c]),substr(criteria[l],1,1),sep="")
              temp2 <- paste(toupper(categories[c]),substr(criteria[l],2,2),sep="")
              cond1 <- paste("(",temp1,"==T)",sep="")
              cond2 <- paste("(",temp2,"==T)",sep="")
              listed.yrs2 <- eval(parse(text=paste("subset(rl.rep,select=Year, subset=(",cond1,"|",cond2,"))[,1]",sep="")))
            } else if(nchar(criteria[l])==3){
              temp1 <- paste(toupper(categories[c]),substr(criteria[l],1,1),sep="")
              temp2 <- paste(toupper(categories[c]),substr(criteria[l],2,2),sep="")
              temp3 <- paste(toupper(categories[c]),substr(criteria[l],3,3),sep="")
              cond1 <- paste("(",temp1,"==T)",sep="")
              cond2 <- paste("(",temp2,"==T)",sep="")
              cond3 <- paste("(",temp3,"==T)",sep="")
              listed.yrs2 <- eval(parse(text=paste("subset(rl.rep,select=Year, subset=(",cond1,"|",cond2,"|",cond3,"))[,1]",sep="")))
            }
            listed.yrs3 <- intersect(listed.yrs2,eval.yrs2)

                      ## FOR 10-YEAR INTERVAL
            if(eval.intervals[i]==10){

             			# To handle the complexities of the 10-year re-evaluation interval
 					# KTS: use a while loop- as long as there are one or more "YN" signals, do a re-evaluation...
              listed.yrs_alt <- intersect(listed.yrs2,eval.yrs2_alt)   		# correct listing status for all potential re-evaluations to confirm downlisting status
              listed.yrs_doubled <- sort(c(listed.yrs3,(listed.yrs3+5)),decreasing=F)
              all.years <- sort(c(eval.yrs2,eval.yrs2_alt),decreasing=F)          
              logicalvec <- (all.years%in%listed.yrs_doubled)|(all.years>=ex.year)      # boolean vector: naive, considering only original evaluations, no re-evaluations
              reeval.ndx <- as.numeric(gregexpr("10",paste(as.numeric(logicalvec),collapse=""))[[1]])+1  # re-evaluate if this pattern is matched
              reeval.yrs <- all.years[reeval.ndx]+5       		# years for which a re-evaluation is needed
              while(any(reeval.ndx>0)&any(reeval.yrs%in%listed.yrs_alt)){
                logicalvec[reeval.ndx[which(reeval.yrs%in%listed.yrs_alt)]-1] <- TRUE     # set the listing status to TRUE
                logicalvec[reeval.ndx[which(reeval.yrs%in%listed.yrs_alt)]] <- TRUE 
                reeval.ndx <- as.numeric(gregexpr("10",paste(as.numeric(logicalvec),collapse=""))[[1]])+1  # re-evaluate if this pattern is matched
                reeval.yrs <- all.years[reeval.ndx]+5       		# years for which a re-evaluation is needed
              }
              listed.yrs <- all.years[logicalvec]             
            }

                        ## FOR 5-year interval
            if(eval.intervals[i]==5){
                            				 # Handle the same issue for 5-year interval (much easier!)
              logicalvec <- (eval.yrs2%in%listed.yrs3)|(eval.yrs2>=ex.year)      # boolean vector: naive, considering only original evaluations, no re-evaluations          
              reeval.ndx <- as.numeric(gregexpr("101",paste(as.numeric(logicalvec),collapse=""))[[1]])+1  # re-evaluate if this pattern is matched
              while(any(reeval.ndx>0)){
                logicalvec[reeval.ndx] <- TRUE    # fill in "gaps"
                reeval.ndx <- as.numeric(gregexpr("101",paste(as.numeric(logicalvec),collapse=""))[[1]])+1  # re-evaluate if this pattern is matched
              } 
              listed.yrs <- eval.yrs2[logicalvec]
            }

            if(eval.intervals[i]==1){
              listed.yrs <- listed.yrs2
            }

                      # determine time length between listings...
            tempdiff <- diff(listed.yrs)

            timeInVarName <- paste(categories[c],".timein",intervals[i],criteria[l],sep="")
            multiplier <- ifelse(eval.intervals[i]>1,5,1)               # multiplier for total time listed
            eval(parse(text=paste(timeInVarName,"<-length(which(listed.yrs%in%yearsin_years))*multiplier",sep="")))
            
            eval(parse(text=paste(masterFileName,"$",timeInVarName,"<-",timeInVarName,sep=""))) 
                ###[[masterFileName]]$[[timeInVarName]] <- [[timeInVarName]]   

     
                       # determine time continuously listed...

            firstVarName <- paste("first.",categories[c],intervals[i],criteria[l],sep="")
            timeToVarName <- paste(categories[c],".timeto",intervals[i],criteria[l],sep="")

            if(length(listed.yrs)>0){         # if there are any listings 
              eval(parse(text=paste(firstVarName,"<-min(listed.yrs)",sep="")))
              if (all(tempdiff<=5) ){ #Test if no more than five years between listings.. 
                eval(parse(text=paste(timeToVarName,"<-ex.year - min(listed.yrs)",sep="")))   
              } else {                        # if not "continuous" listing (if species was downlisted after initially being listed)
                last.run <- c((max(which(tempdiff>5))+1):length(listed.yrs)) #index vector of last constant run before extinction
                listed.yrs <- listed.yrs[last.run]
                eval(parse(text=paste(timeToVarName,"<-ex.year - min(listed.yrs)",sep="")))
              }
            } else {
              eval(parse(text=paste(firstVarName,"<-NA",sep="")))
              eval(parse(text=paste(timeInVarName,"<-0",sep="")))
              eval(parse(text=paste(timeToVarName,"<-0",sep="")))
            }    # if no listings

            eval(parse(text=paste(masterFileName,"$",timeToVarName,"<-",timeToVarName,sep="")))
            eval(parse(text=paste(masterFileName,"$",firstVarName,"<-",firstVarName,sep="")))

            ############################
            ####### DETERMINE TIME SPENT IN EACH CATEGORY PER RESIT'S REQUEST
            ####### build into new data frame

            if((categories[c]=="vu")&(intervals[i]=="10")&(criteria[l]=="") ){ # if "vu" and no special conditions
              TimeInCategories$spec[counter] <<- rl.rep[1,2]   
              TimeInCategories$model[counter] <<- models[mod]
              TimeInCategories$rep[counter] <<- n 
              TimeInCategories$ex.year[counter] <<- ex.year
              TimeInCategories$is.ex[counter] <<- ex.status 
              if(length(listed.yrs)>0){
                vu.years <- listed.yrs
                yearsOfInterest <- c(min(vu.years):(ex.year-1))
                ndx <- which((rl.rep$cat.num<3)&(rl.rep$cat.num>=1))
                en.years_t1 <- rl.rep$Year[ndx]
                ndx <- which((rl.rep$cat.num<2)&(rl.rep$cat.num>=1))
                cr.years_t1 <- rl.rep$Year[ndx]

                TimeInCategories$TotalTimeVU[counter] <<- ex.year - min(vu.years)
                TimeInCategories$TotalTimeCR[counter] <<- length(intersect(cr.years_t1,yearsOfInterest))
                TimeInCategories$TotalTimeEN[counter] <<- length(intersect(en.years_t1,yearsOfInterest))
              
              } else{
                TimeInCategories$TotalTimeVU[counter] <<- 0   #ex.year - min(listed.yrs) 
                TimeInCategories$TotalTimeCR[counter] <<- NA
                TimeInCategories$TotalTimeEN[counter] <<- NA
              }
              counter <<- counter+1
            }

            ##############
            ##############

          }  # END listing criteria
        }  # END eval intervals
      }else{  # END if any annual listings

            # loop through evaluation intervals
        i=1
        for(i in 1:length(intervals)){

             # loop through listing criteria
          l=1
          for(l in 1:length(criteria)){
            timeInVarName <- paste(categories[c],".timein",intervals[i],criteria[l],sep="")
            timeToVarName <- paste(categories[c],".timeto",intervals[i],criteria[l],sep="")
            firstVarName <- paste("first.",categories[c],intervals[i],criteria[l],sep="")
            eval(parse(text=paste(firstVarName,"<-NA",sep="")))
            eval(parse(text=paste(timeToVarName,"<-0",sep="")))
            eval(parse(text=paste(timeInVarName,"<-0",sep="")))
            eval(parse(text=paste(masterFileName,"$",timeInVarName,"<-",timeInVarName,sep="")))
            eval(parse(text=paste(masterFileName,"$",timeToVarName,"<-",timeToVarName,sep="")))
            eval(parse(text=paste(masterFileName,"$",firstVarName,"<-",firstVarName,sep="")))
          }
        } 

      } # END if no annual listings
    } # END categories
  #} # END if extinct??  [try to get by without this loop...]

  tempdf <- eval(parse(text=masterFileName))
  
  #return(tempdf)  # return all information about this model run

  if(ex.status){
    cattimes.ex <<- rbind(cattimes.ex,tempdf)
  }else{
    cattimes.noex <<- rbind(cattimes.noex,tempdf)
  }

  if(ex.status) CMCounter <<- CMCounter + 1

}    # END FUNCTION



       

#######################
### END LOADING FUNCTIONS



########################
### LOAD ANCILLARY DATA and response variables from all 'specific realized population models' 
       # from the NASA project (stores all information for each SRPM)

if(KEVIN) setwd("C:\\Users\\Kevin\\Desktop\\NASA Analysis")
load("FinalNASAData.RData")
if(KEVIN) setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres")



############################################################################################
# find models that show either increased/or no increased risk of extinction risk due to 
# climate change. 
# Flag models that show ex. risk >10% under no climate change?

mod.pex <- c("mp.run","ext.risk" )   # names of columns to select (to extract extinction probability)  

WREpex <- subset(dfmain, subset=(match.stat=="Match" & clim.chg=='WRE'), select=mod.pex)

NOCCpex <- subset(dfmain, subset=(match.stat=="Match" & clim.chg=='NoCC'), select=mod.pex)

p.ex <- merge(x=WREpex, y=NOCCpex, by='mp.run')
p.ex['ch.pex'] <- NA
p.ex$ch.pex <- p.ex$ext.risk.x - p.ex$ext.risk.y

p.ex <- subset(p.ex, ext.risk.y<0.1)  # only models with low prob extinction without climate change
p.ex.risk <- subset(p.ex, ch.pex>0)  # Only models with increase in risk due to climate change
p.ex.lowrisk <- subset(p.ex, ch.pex<=0) # Only models with no increase in risk due to climate change
mp.risk <- as.character(unique(p.ex.risk$mp.run))  # identify those models of interest: those MP files that exhibited low risk of extinction without climate change but higher risk with climate change.
mp.lowrisk <- as.character(unique(p.ex.lowrisk$mp.run))

#######################
##############################################################################################

# KTS: new downlisting model (e.g., for 10 year evaluation interval): if a downlisting evaluation occurs after 10 years then 
      #  re-evaluate after 5 years. If still meeting the downlisting criteria, then downlist. For each 5-year reevaluation interval,
      #  store whether a downlisting occurred, and what category it was listed as... 

###All Models##################
#Assemble the category count by year results from each .mp into single data frame 


counter <<- 1   # set counter variable for the TimeInCategories analysis    

CMCounter <<- 1  # set counter variable for the stage transition matrix 

###############################
####    LOOP THROUGH RED-LISTING RAW SPREADSHEETS


###################################
######   LOOP THROUGH SPECIES
     #cnt = 32 #18 #  12 # 1 # 7 #   1 # 20  #7 #18     # 1 #  
for (cnt in 1:length(rl.files)){   				# loop through species
###################################

  rlres <- read.csv(rl.files[cnt]) 				# Read in model results file (all SRPMs for that species)
  rlres$FileName <- as.character(rlres$FileName)
  mp <- unlist(strsplit(rlres$FileName,split='/')) 	# get file name of .mp
  rlres$FileName <- mp[seq(5, length(mp), 5)]    	# unique for each realized species (unique mp file)
  rlres <- rlres[,keep.cols]
  models <- unique(rlres$Model)

###################################
######   LOOP THROUGH SRPMs
         #mod=which(models==3); n=2    # which(models==11);n=10
  for (mod in 1:length(models)){  				#separate by distinct MP models (realized species)
###################################

    rl.mod <- subset(rlres, Model==models[mod])    # data frame: only data for this SRPM

###################################
######  LOOP THROUGH REPLICATES
    for (n in 1:max(rl.mod$Replicate)) {  		#separate by replicate (10 per MP file)
###################################
      ### CALL NEW FUNCTION HERE!!! ######
  
      processTrajectory(rl.mod,n,mod)      
      #cattimes.ex <-  rbind(cattimes.ex, cat.rep.ex)

    }# end loop for replicates
  } # end loop for separating models
}# end loop results file


#######  remove excess rows from the TimeInCategories data frame

head(TimeInCategories)
tail(TimeInCategories)

counter
CMCounter

TimeInCategories2 <- TimeInCategories[1:counter,]

###############################
##########  reminder of the key files in the workspace at this point

####  ex.reps ===> redlist model runs that went extinct, each year of each simulation, starting at year 10 
####  non.ex.reps ===> redlist model runs that weren't extinct by the end of the simulation
####  cattimes.ex ===> for each redlist replicate that went extinct (up to 10 per model), key summary statistics: especially time spent at or above a particular classification before going extinct...
####  cattimes.noex ===> for each redlist replicate that didn't go extinct (up to 10 per model), key summary statistics: especially time spent at or above a particular classification before the end of the simulation...


#####  Remove models that don't show increased extinction risk due to climate change
ex.reps.risk <- subset(ex.reps, FileName %in% mp.risk)  
non.ex.reps.risk <- subset(non.ex.reps, FileName %in% mp.risk)

risk.reps <- rbind(ex.reps.risk, non.ex.reps.risk) #re-combine extinction and non-extinction reps: now filtered so that only those models with higher extinction risk due to climate change are included 
#  head(risk.reps)

ex.reps.lowrisk <- subset(ex.reps, FileName %in% mp.lowrisk)  
non.ex.reps.lowrisk <- subset(non.ex.reps, FileName %in% mp.lowrisk)

lowrisk.reps <- rbind(ex.reps.lowrisk, non.ex.reps.lowrisk) #re-combine extinction and non-extinction reps    # now filtered so that only those models with little climate related risk are included

##################################################################################
                #   KTS: I'm not sure this is what I want: need to keep only cattimes.ex, right? if simulations don't go extinct before the end, then this is censored data... 

cattimes <- rbind(cattimes.noex, cattimes.ex)       # conflates time to extinction with time to end of simulation, right?
cattimes.low <- subset(cattimes, spec %in% mp.lowrisk)    

cattimes.risk <- subset(cattimes, spec %in% mp.risk)      # now only for those models that have higher risk due to CC


#################################################################################
                 # Now determine which replicates are left censored for subsetting purposes...
			# KTS: now, the filter is only "CR"

      # KTS: note: filtering will need to change for 5 and 10 year runs...

  #  head(cattimes.ex)   # names(cattimes.ex)

cattimes.ex$cr.lcens <- ifelse((cattimes.ex$ex.year-cattimes.ex$cr.timeto)==10,1,0)
cattimes.ex$en.lcens <- ifelse((cattimes.ex$ex.year-cattimes.ex$cr.timeto)==10,1,0)
cattimes.ex$vu.lcens <- ifelse((cattimes.ex$ex.year-cattimes.ex$cr.timeto)==10,1,0)
cattimes.ex$nt.lcens <- ifelse((cattimes.ex$ex.year-cattimes.ex$cr.timeto)==10,1,0)

cattimes.ex$cr.lcens5 <- ifelse((cattimes.ex$ex.year-cattimes.ex$cr.timeto5)%in%c(10:14),1,0)
cattimes.ex$en.lcens5 <- ifelse((cattimes.ex$ex.year-cattimes.ex$cr.timeto5)%in%c(10:14),1,0)
cattimes.ex$vu.lcens5 <- ifelse((cattimes.ex$ex.year-cattimes.ex$cr.timeto5)%in%c(10:14),1,0)
cattimes.ex$nt.lcens5 <- ifelse((cattimes.ex$ex.year-cattimes.ex$cr.timeto5)%in%c(10:14),1,0)

cattimes.ex$cr.lcens10 <- ifelse((cattimes.ex$ex.year-cattimes.ex$cr.timeto10)%in%c(10:19),1,0)
cattimes.ex$en.lcens10 <- ifelse((cattimes.ex$ex.year-cattimes.ex$cr.timeto10)%in%c(10:19),1,0)
cattimes.ex$vu.lcens10 <- ifelse((cattimes.ex$ex.year-cattimes.ex$cr.timeto10)%in%c(10:19),1,0)
cattimes.ex$nt.lcens10 <- ifelse((cattimes.ex$ex.year-cattimes.ex$cr.timeto10)%in%c(10:19),1,0)



cattimes.noex$cr.lcens <- ifelse((111-cattimes.noex$cr.timeto)==10,1,0)
cattimes.noex$en.lcens <- ifelse((111-cattimes.noex$cr.timeto)==10,1,0)
cattimes.noex$vu.lcens <- ifelse((111-cattimes.noex$cr.timeto)==10,1,0)
cattimes.noex$nt.lcens <- ifelse((111-cattimes.noex$cr.timeto)==10,1,0)

cattimes.noex$cr.lcens5 <- ifelse((111-cattimes.noex$cr.timeto5)%in%c(10:14),1,0)
cattimes.noex$en.lcens5 <- ifelse((111-cattimes.noex$cr.timeto5)%in%c(10:14),1,0)
cattimes.noex$vu.lcens5 <- ifelse((111-cattimes.noex$cr.timeto5)%in%c(10:14),1,0)
cattimes.noex$nt.lcens5 <- ifelse((111-cattimes.noex$cr.timeto5)%in%c(10:14),1,0)

cattimes.noex$cr.lcens10 <- ifelse((111-cattimes.noex$cr.timeto10)%in%c(10:19),1,0)
cattimes.noex$en.lcens10 <- ifelse((111-cattimes.noex$cr.timeto10)%in%c(10:19),1,0)
cattimes.noex$vu.lcens10 <- ifelse((111-cattimes.noex$cr.timeto10)%in%c(10:19),1,0)
cattimes.noex$nt.lcens10 <- ifelse((111-cattimes.noex$cr.timeto10)%in%c(10:19),1,0)



################################################################################
                # Save key time-specific variables for further use (e.g., plotting)

 #if(KEVIN) setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\RData")
if(KEVIN) setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlist_overflow")
filename <- paste("TimeSpecificRLData_",Sys.Date(),".RData",sep="")
save(ex.reps,non.ex.reps,risk.reps,lowrisk.reps,file=filename)


################################################################################
                # Save key model-specific variables for further use (e.g., plotting)

if(KEVIN) setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\RData")
filename <- paste("ModelSpecificRLData_",Sys.Date(),".RData",sep="")
save(cattimes.ex,cattimes.noex,file=filename)


###############################################################################
                # Save extra storage variables
if(KEVIN) setwd("C:\\Users\\Kevin\\Dropbox\\NASA_redlistNewres\\RData")
filename <- paste("AncillaryRLData_",Sys.Date(),".RData",sep="")
save(CategoryMatrices,TimeInCategories2,file=filename)




