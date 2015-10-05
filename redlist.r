redlist <- function(mp.file.name, raw.file.name, mac=FALSE, ptc=FALSE, ptcFiles='no file', 
                     habdyn=FALSE, hdhFile='no file') {
# RAMAS Redlist Continuous Analysis Program
# Authors: Jessica Stanton, Matt Aiello-Lammens, and Ben Weinstein
# Date: October 2012
#
# Args:
#	mp.file.name: the name and path of the *.mp file to be read.  Path can be either from the 
#	current working directory or from the root directory
#
#	raw.file.name: the name and path of the rawfile.txt that corresponds to the results in the mpFile to be read.
#
# 	mac: TRUE/FALSE: Is this being run on a mac or *nix machine vs. a Windows machine? If 
#     TRUE the program Wine must be installed http://www.winehq.org/
# 
# 	ptc: TRUE/FALSE - Are there ptc files to extract area data? If not, critera based on area
#     measures will not be calculated.
#  
# 	ptcFiles: A single path to a ptc file or a vector of ptc files with which to extract
#  	  data from.
#  
#	habdyn: TRUE/FALSE - Does the model include a dynamic habitat with a habitat dynamics history file? If true,
#		the HabDyn history file must be specified.
#
# 	hdhFile: The full path name to the HabDyn history (.txt) file
#
# Returns:
# redlist.pars - A data.frame of Redlist parameters and category designations for each timestep of each replicate
###################################################################################################################
# This script is written to calculate threat categories of extinction risk according to the IUCN Red List criteria
# continuously from population simulations in RAMAS.
#
#
#
# THRESHOLDS
MinYears <- c(3, 5, 10)
SmallPopSize <- 100  # for severely fragmented, defn of "small" 

# Read in basic variables for calculations from the .mp file
mp.file <- mp.read(mpFile=mp.file.name)
mp.file <- mp.file$mp.file
Stages <- mp.file$Stages

mp.raw <- mp.read.raw(rawFile=raw.file.name)
nPop <- mp.raw$npop
RepRun <- mp.raw$RepRun
tsteps <- mp.raw$tsteps
rawOut <- mp.raw$rawOut

#calculate the mature stages  ###This will only read the first matrix if there are more than one###
Mtrx<-mp.file$StMatr[[1]]$Matr
const.Mtrx<-mp.file$ConstraintsMatr
Mtrx.const<-(1-const.Mtrx)*Mtrx
max.fec <- apply(Mtrx.const,2,max)
mature.stages <- as.numeric(max.fec>0)
#Add column to raw output matrix to for total mature individuals in each population
which.mature <- which(mature.stages > 0)
totNmat <- apply(as.matrix(rawOut[,which.mature+4]),1,sum)
rawOut <- cbind(rawOut,totNmat)
#Total individuals in all populations
Ntotal <- tapply(rawOut$totN, INDEX = list(rawOut$rep, rawOut$ts), sum)
Ntotal <- t(Ntotal)
Ntotal <- replace(Ntotal, is.na(Ntotal), 0)
#Total mature individuals in all populations 
Nmature <- tapply(rawOut$totNmat, INDEX = list(rawOut$rep, rawOut$ts), sum)
Nmature <- t(Nmature)
Nmature <- replace(Nmature, is.na(Nmature), 0)
#Count number of occupied (with mature individuals) populations
popOcc <- as.numeric(rawOut$totNmat>0)
rawOut <- cbind(rawOut, popOcc)
subpopOcc <- tapply(rawOut$popOcc, INDEX = list(rawOut$rep, rawOut$ts), sum) 
subpopOcc <- t(subpopOcc)
subpopOcc <- replace(subpopOcc, is.na(subpopOcc), 0)
singleLoc<- subpopOcc==1
fiveLoc<- subpopOcc<=5
tenLoc<-subpopOcc<=10
fifteenLoc <- subpopOcc<=15

#calculate proportion of total population in largest subpopulation at a given timestep (mature indv.)
####TO DO:   DEAL WITH ZEROS 
LargestPopSize <- tapply(rawOut$totNmat, INDEX = list(rawOut$rep, rawOut$ts), max)
LargestPopSize <- t(LargestPopSize)
LargestPopSize <- replace(LargestPopSize, is.na(LargestPopSize), 0 )
PropInLargestPop <- LargestPopSize / Nmature
PropInLargestPop <- replace(PropInLargestPop, is.nan(PropInLargestPop), 1)
#calculate proportion of total (mature) population in small populations
popSm <- as.numeric(rawOut$totNmat<=SmallPopSize)
NpopSm <- rawOut$totNmat*popSm
rawOut <- cbind(rawOut, NpopSm)
Nsmallpop <- tapply(rawOut$NpopSm, INDEX = list(rawOut$rep, rawOut$ts), sum)
Nsmallpop <-t(Nsmallpop)
#Determine if criteria for severe fragmentation is met
#####TO DO: Add Isolation criteria
SevereFrag <- (Nsmallpop/Nmature) > 0.5
SevereFrag <- replace(SevereFrag, is.na(SevereFrag),FALSE)

#Determine if the population goes extinct
Extinct <- Ntotal == 0

#_________________________________________________________________________#
#Calculate generation timeCall GenTime.exe to determine various calculations of generation time and eigen values
if( mac ) {
  # Development Notes: figured out a way to run GenTime on Mac or Linux if Wine is installed
  wine <- '/Applications/Wine.app/Contents/Resources/bin/wine' # Works for Matt's Mac
  # Assuming GenTime.exe is in same directory as scripts
  gen.exe <- paste( redlist.base.dir, 'GenTime.exe ', sep='')
  gen.call <- paste( wine, gen.exe , '"', mp.file.name , '"' )
  gen <- system( gen.call, intern=TRUE, ignore.stderr = TRUE )
  gen <- strsplit( gen, " ")
  gen <- unlist( lapply( gen, as.numeric ) )
  Gen.Abar <- gen[1]
  Gen.Mu1 <- gen[2]
  Gen.Tgen <- gen[3]
  
} else {
  gen.call <- paste( redlist.base.dir,'GenTime ', '"', mp.file.name , '"',sep="" )
  #print( paste('GenTime.exe called as: ', gen.call) ) ### DEBUG LINE
  gen <- system( gen.call, intern = TRUE )
  gen <- strsplit( gen, " ")
  gen <- unlist( lapply( gen, as.numeric ) )
  Gen.Abar <- gen[1]
  Gen.Mu1 <- gen[2]
  Gen.Tgen <- gen[3]
}

all.gen.times <- sort(c(Gen.Abar,Gen.Mu1,Gen.Tgen))
GTime <- all.gen.times[2]
#Determine interval lengths for percent reduction
gen <- round(c(GTime, GTime*2, GTime*3))
intlength <- pmax(gen, MinYears)


#calculate slope and percent reduction for the ln(total mature) for each interval
#percentreduction gives a list of 3 elements corresponding to the 3 interval lengths
Nmature.nozero <- replace(Nmature, Nmature ==0, 0.001)
lnNmat <- log(Nmature.nozero)
slope <- function(data){
  coef(lm(data ~ as.numeric(rownames(data)), data = as.data.frame(data)))[2,]
}
percentreduction <- vector(mode='list',length = 0)
for(i in 1:length(intlength)){
  loopcalc <-data.frame()
  loopcalc <- rbind(loopcalc, rep(NA, ncol(lnNmat)), rep(NA, ncol(lnNmat)))
  timeperiod <- intlength[i]
  if(timeperiod>3){ #Extrapolate slope for timesteps less than the first interval
    for (x in 3:(timeperiod-1)){
      data.rowx <- lnNmat[1:x,] 
      rate <- slope(data.rowx)
      calc.rowx <- 100-exp(rate*timeperiod)*100
      loopcalc <- rbind(loopcalc,calc.rowx)
    }
  }
  rate <- rollapply(lnNmat, width = timeperiod, FUN = slope, by.column = FALSE)
  rollcalc <- as.data.frame(100-exp(rate*timeperiod)*100)
  names(loopcalc)<-names(rollcalc)
  loopcalc <- rbind(loopcalc, rollcalc)
  percentreduction[[i]] <- loopcalc
}
for(i in 1:3){
percentreduction[[i]] <- replace(percentreduction[[i]], is.na(percentreduction[[i]]),0) #remove NA values
percentreduction[[i]] <- replace(percentreduction[[i]], percentreduction[[i]]<0,0)  #Set populaton increases (neg percent reduction) to zero
percentreduction[[i]] <- round(percentreduction[[i]], 1)  #round to one decimal place
}

#-----------------------------------------------------------------------------------#
##### Begin section concerning spatial parameters ################################
occupied.area <- data.frame()
occupied.extent <- data.frame()

if (ptc) {
  ##list patch map files - RAMAS output files .P_S and .RDCS
  ## are equivelent to .RST and .RDS raster files
  spatial.info <- unlist(strsplit(ptcFiles[1], split='/'))
  pth.lgth <- length(spatial.info)
  spatial.path <- paste(spatial.info[1:pth.lgth-1], collapse='/')
  
  ######  !!!!!!! NOTE: If running in Parallel, code below should be done outside of function
  ######   Raster maps generated from RAMAS are in .P_S format and must be converted to .RST
  ######   to be read-in. Beginning at line 338 files are converted back    !!!!!!!!!!
  
  #patchmaps <- list.files(path=spatial.path,pattern='*LO.....P_S', full.names=TRUE)
  #patch_rdcs<- list.files(path=spatial.path,pattern='*LO.....RDCS', full.names=TRUE)
  #sapply(patchmaps, FUN=function(eachfile){
  #  file.rename(from = eachfile, to = sub(pattern=".P_S",replacement=".RST",eachfile))
  #})
  #sapply(patch_rdcs, FUN=function(eachfile){
  #  file.rename(from = eachfile, to = sub(pattern=".RDCS",replacement=".RDC",eachfile))
  #})
  patchrast <- list.files(path=spatial.path, pattern='*.RST', full.names=TRUE)
  
  ## Check that the ptcFiles exist
  ptcFilesExist <- lapply(ptcFiles, file.exists)
  ptcFilesExist <- unlist(ptcFilesExist)
  # If any file does not exist, then stop the program
  if( !all(ptcFilesExist) ){
    stop('Did not find all *.ptc file(s).')
  }
  ## Check that there are an equal number of ptcFiles and ptcFileIter
  ptcFileIter <- c(1:mp.raw$tsteps)
  if( !(length(ptcFiles)==length(ptcFileIter)) ) {
    stop('Uneven number of ptcFiles and ptcFileIter.')
  }
  ## Check that there are an equal number of raster files and ptcfiles
  if( !(length(ptcFiles)==length(patchrast)) )  {
    nopatch<-length(patchrast)
	warning('Uneven number of ptcFiles and Rasters. 
          Timesteps with no habitat maps will be zeros for EOO & AOO')
  } else {
  nopatch <- length(ptcFiles)+1
  }
  ## Read the ptcFile(s)
  ptc.files <- lapply( ptcFiles, ptc.read )
  
  ## Begin a for loop that goes through each of the *.ptc files included
  for ( ptc.cnt in 1:length(ptcFiles) ){
    if( ptc.cnt<=nopatch) {  # Check if there are any patches at all
	# Check that the program is reading in the raster that matches the ptc file
	ptc.year <- as.numeric(gsub("\\D","",ptcFiles[ptc.cnt]))
	rast.year <- as.numeric(gsub("\\D","",patchrast[ptc.cnt]))
	if( !(ptc.year==rast.year)) {
	stop('ptc file does not match raster')
	}

    # Load in the first raster and convert to points
    print("Reading in raster")
	  allpatch.rst <- raster(patchrast[ptc.cnt])
    allpatch.pts<- rasterToPoints(allpatch.rst)
    colnames(allpatch.pts)<- c("x","y","patch")
    
	## Begin another loop that goes through each replicate in the raw output file
    area.ts.rep <- vector()
    EOO.ts.rep <- vector()
	pop.ab.iter <- subset(rawOut, ts == ptc.cnt, select= c(rep,pop,totNmat))
    for (rep.cnt in 1:mp.raw$MaxRep){
      # Look at the ptc file and determine the number of patches that exist
      # in the current file
      PatchN <- ptc.files[[ ptc.cnt ]]$PopN
      pop.ab.iter <- subset(pop.ab.iter, rep == rep.cnt)
      
      if (habdyn) {
        ## First, let's consider the scenario in which habdyn = TRUE
        ## In this scenario the HabDyn Module in RAMAS GIS was used to incorporate 
        ## habitat dynamics into the *.mp file.
        
        ## Check that hdhFile exists
        if( !file.exists( hdhFile ) ){
          stop( paste( 'Did not find HabDyn History file: ', hdhFile ) )
        }  
        ## Read the HabDyn History file
        hdhist <-hdhist.read( hdhFile )
        
        # Look at hdhist file and determine which *.mp populations exist at 
        # ptcFileIter[ptc.cnt] iteration (time step).  There can be a greater number of
        # mp populations than ptc populations, since some populations in the mp file will
        # have K=0 for some time steps
        patch.mat.sub <- subset( hdhist$patch.info.mat, 
                                 subset=hdhist$patch.info.mat$iter==ptcFileIter[ptc.cnt])
        print(paste('Determining occupied patches at timestep ',ptc.cnt, ' replicate ', rep.cnt))
		# Limit this subset to only the patches that exist in the ptc file for this 
        # iteration step
        patch.mat.sub <- patch.mat.sub[1:PatchN,]
        # Address the case where there is a zero in the new2old column (usually 
        # only for the first iteration step)
        iter.mp.pops <- ifelse( test=patch.mat.sub$new2old==0,
                                yes=patch.mat.sub$patch, no=patch.mat.sub$new2old)
        
        # Determine the matching mp simulation start year for a particular iteration
        iter.start.yr <- unique( patch.mat.sub$Fyear )
        if ( length(iter.start.yr) > 1 ){
          stop('mp.results: HabDyn History file is incorrectly formated.
  			  There should only be one Fyear per iter step.')
        } 
        # Now check which populations have non-zero abundance.  If iter.start.yr==0, then use
        # initial abundance from PopData_df.  If iter.start.yr>1, then use populations with non-zero adult
        # abundance from the rawfile
        if ( iter.start.yr==0 ){
          pop.ab <- mp.file$PopData_df[ iter.mp.pops, 4 ] # InitAb is column 4
        } else {
          pop.ab <- pop.ab.iter[ iter.mp.pops, 3]
        }
        
        # Determine the ptc populations with non-zero abundance
        patch.nonzero <- patch.mat.sub$patch[ which( pop.ab > 0 ) ]
        #browser()  #REMOVE AFTER DEBUG
      } else {
        ## If habdyn=FALSE, then it is assumed the that there is only one *.ptc file associated
        ## with the *.mp file. That is, no habitat dynamics were calculated using the HabDyn
        ## module. Using HabDyn is NOT the only way to incorporate habitat dynamics (e.g. simple
        ## declines in K can be added to the *.mp file directly), but for the sake of this 
        ## calculation, we are assuming that if there is no HabDyn History file included, then
        ## the *.ptc file sets an unchanging spatial structure for the duration of the simulation.
        ## Occupied area may change in this instance due to colonization and extinction of existing patches.
    
        patch.nonzero <- which( pop.ab.iter[3] > 0 )
        
      } #END IF HABDYN
      
      ## Calculate Area and Shape Indices for non-zero abundance populations
      #
      # Number of non-zero patches
      patch.n.nonzero <- length( patch.nonzero )
      
	  if ( patch.n.nonzero > 0 ) {  #determine if there are any occupied patches
      # Extract non-zero population (patch) information from the 
      # ptc.file$PopLandInd_df data frame and from ptc.file$PopPatchChar_df
      PopLandInd_nonzero <- ptc.files[[ptc.cnt]]$PopLandInd_df[ patch.nonzero, ]
      PopPatchChar_nonzero <- ptc.files[[ptc.cnt]]$PopPatchChar_df[ patch.nonzero, ]
      
      # Extract non-zero patch locations
      occpatch.pts <- allpatch.pts[which(allpatch.pts[,3] %in% patch.nonzero),]
      
      # Total area of patches (km^2)
      # NOTE: must convert from area in cells to km^2 by multiplying the number
      # of cells by the cell length squared
      tot.patch.area.cells <- sum(PopPatchChar_nonzero$Area)
      tot.patch.area.km2 <- tot.patch.area.cells * ptc.files[[ptc.cnt]]$cell.length^2
      area.ts.rep <- c(area.ts.rep, tot.patch.area.km2)  
      
      # Total extent of occupied patches (km^2) - assuming original file is in meters
      z<- rep(1,nrow(occpatch.pts))
      xyz<-as.data.frame(cbind(occpatch.pts[,1:2],z))
      tot.extent.area.km2<- mcp.area(xyz[,1:2],xyz[,3], percent = 100,
                                     unin= 'm', unout = 'km2', plotit=FALSE )
      EOO.ts.rep <- c(EOO.ts.rep, tot.extent.area.km2[,1])
	  } else {  # If there are no occupied patches insert zeros for both area measures
	  area.ts.rep <- c(area.ts.rep, patch.n.nonzero)
	  EOO.ts.rep <- c(EOO.ts.rep, patch.n.nonzero)
	  }
    } # End for rep.cnt loop
	} else {   # If there are no available patches occupied or not insert Zeros for both measures
	area.ts.rep <- rep(0,mp.raw$MaxRep)
	EOO.ts.rep <- rep(0,mp.raw$MaxRep)
	}
    occupied.area <- rbind(occupied.area, area.ts.rep)
    occupied.extent <- rbind(occupied.extent, EOO.ts.rep )

	
  } # End for ptc.cnt loop
  ###!!!!!!!  NOTE: RST files are returned to RAMAS specific file extentions 
  ##Return patch map files to RAMAS output format .P_S and .RDCS
  #patchmaps <- list.files(path=spatial.path,pattern='*LO.....RST', full.names=TRUE)
  #patch_rdc<- list.files(path=spatial.path,pattern='*LO.....RDC$', full.names=TRUE)
  #sapply(patchmaps, FUN=function(eachfile){
  #  file.rename(from = eachfile, to = sub(pattern=".RST",replacement=".P_S",eachfile))
  #})
  #sapply(patch_rdc, FUN=function(eachfile){
  #  file.rename(from = eachfile, to = sub(pattern=".RDC$",replacement=".RDCS",eachfile))
  #})  
  
} else { #End if PTC
  occupied.area <- (Nmature*4) 
  occupied.extent <- replace(Nmature,Nmature>=0,9999999) } # Marker value to indicate that EOO has not been calculated

AOO <- pmin(as.matrix(occupied.area), (Nmature*4))
AOO <- AOO[1:nrow(Nmature),]
EOO <- as.matrix(occupied.extent)
EOO <- EOO[1:nrow(Nmature),]

########################################  END SPATIAL CALCULATIOnS ##################


#Calculate Continuing decline for B1 or B2 part b. "Continuing" is defined as the slope of the function 
	# observed over the last 1 generation or 3 years (whichever is longer). To meet the criteria any of the 
	# following may be met:
		# .EOO  decline in Extent of Occurrence with slope less than -.1
		# .AOO  decline in Area of Occupancy with slope less than -.1
		# .subpopOcc decline in the number of occupied locations with slope less than -.1
		# .lnNmat  decline in ln(N mature) of 1% or more
 
 	###Continuing Decline in EOO
  loopcalc <-data.frame()
  loopcalc <- rbind(loopcalc, rep(NA, ncol(EOO)), rep(NA, ncol(EOO)))
  timeperiod <- intlength[1]
  if(timeperiod>3){ #Extrapolate slope for timesteps less than the first interval
    for (x in 3:(timeperiod-1)){
      data.rowx <- EOO[1:x,] 
      rate <- slope(data.rowx)
      calc.rowx <- rate<(-0.1)
      loopcalc <- rbind(loopcalc,calc.rowx)
    }
  }
  rate <- rollapply(EOO, width = timeperiod, FUN = slope, by.column = FALSE)
  rollcalc <- as.data.frame(rate<(-.1))
  names(loopcalc)<-names(rollcalc)
  loopcalc <- rbind(loopcalc, rollcalc)
  ContDecline.EOO <- loopcalc
  ContDecline.EOO <- replace(ContDecline.EOO, is.na(ContDecline.EOO), FALSE)
  names(ContDecline.EOO)<-c(1:RepRun)
 
	###Continuing Decline in AOO
  loopcalc <-data.frame()
  loopcalc <- rbind(loopcalc, rep(NA, ncol(AOO)), rep(NA, ncol(AOO)))
  timeperiod <- intlength[1]
  if(timeperiod>3){ #Extrapolate slope for timesteps less than the first interval
    for (x in 3:(timeperiod-1)){
      data.rowx <- AOO[1:x,] 
      rate <- slope(data.rowx)
      calc.rowx <- rate<(-0.1)
      loopcalc <- rbind(loopcalc,calc.rowx)
    }
  }
  rate <- rollapply(AOO, width = timeperiod, FUN = slope, by.column = FALSE)
  rollcalc <- as.data.frame(rate<(-.1))
  names(loopcalc)<-names(rollcalc)
  loopcalc <- rbind(loopcalc, rollcalc)
  ContDecline.AOO <- loopcalc
  ContDecline.AOO <- replace(ContDecline.AOO, is.na(ContDecline.AOO), FALSE)
  names(ContDecline.AOO)<-c(1:RepRun)

  ###Continuing Decline in subpopOcc
  loopcalc <-data.frame()
  loopcalc <- rbind(loopcalc, rep(NA, ncol(subpopOcc)), rep(NA, ncol(subpopOcc)))
  timeperiod <- intlength[1]
  if(timeperiod>3){ #Extrapolate slope for timesteps less than the first interval
    for (x in 3:(timeperiod-1)){
      data.rowx <- subpopOcc[1:x,] 
      rate <- slope(data.rowx)
      calc.rowx <- rate<(-0.1)
      loopcalc <- rbind(loopcalc,calc.rowx)
    }
  }
  rate <- rollapply(subpopOcc, width = timeperiod, FUN = slope, by.column = FALSE)
  rollcalc <- as.data.frame(rate<(-.1))
  names(loopcalc)<-names(rollcalc)
  loopcalc <- rbind(loopcalc, rollcalc)
  ContDecline.subpopOcc <- loopcalc
  ContDecline.subpopOcc <- replace(ContDecline.subpopOcc, is.na(ContDecline.subpopOcc), FALSE)

  ###Continuing Decline in lnNmat
ContDecline.lnNmat <- percentreduction[[1]]>=1
ContDecline.lnNmat <- replace(ContDecline.lnNmat, is.na(ContDecline.lnNmat), FALSE)

	##Continuing Decline of any type
ContDecline <- (ContDecline.lnNmat + ContDecline.subpopOcc + ContDecline.AOO + ContDecline.EOO) >= 1

#___________________________________________________________________________________#
#####################Extreme fluctuations: Nmax/Nmin >= 10 calculated for 2 or 3 generations 
################### Parameter ExFluc can be true for extreme fluctuations if any of the 
	# following may be met:
		# .EOO  Extreme fluctuations in Extent of Occurrence 
		# .AOO  Extreme fluctuations in Area of Occupancy 
		# .subpopOcc Extreme fluctuations in number of subpopulations
		# .Nmature  Extreme fluctuations in number of mature individuals

#Find minium values 
colmin <- function(data){
  apply(data,2,min)
}
#Find Maximums
colmax <- function(data){
  apply(data,2,max)
}

#calculate Extreme Fluctuations for EOO ExFluc.EOO  
MinN <- vector(mode='list',length = 0)
for(i in 2:3){
  mincalc <-data.frame()
  mincalc <- rbind(mincalc,EOO[1,])
  timeperiod <- gen[i]
    for (x in 2:(timeperiod-1)){
      data.rowx <- colmin(EOO[1:x,]) 
      mincalc <- rbind(mincalc,data.rowx)
    }
  rollmin <- rollapply(EOO, width = timeperiod, FUN = colmin, by.column = FALSE)
  rollmin <- as.data.frame(rollmin)
  names(mincalc)<-names(rollmin)
  mincalc <- rbind(mincalc, rollmin)
  MinN[[i]] <- mincalc
}
MaxN <- vector(mode='list',length = 0)
for(i in 2:3){
  maxcalc <-data.frame()
  maxcalc <- rbind(maxcalc,EOO[1,])
  timeperiod <- gen[i]
  for (x in 2:(timeperiod-1)){
    data.rowx <- colmax(EOO[1:x,]) 
    maxcalc <- rbind(maxcalc,data.rowx)
  }
  rollmax <- rollapply(EOO, width = timeperiod, FUN = colmax, by.column = FALSE)
  rollmax <- as.data.frame(rollmax)
  names(maxcalc)<-names(rollmax)
  maxcalc <- rbind(maxcalc, rollmax)
  MaxN[[i]] <- maxcalc
}
#Calculate Extreme Fluctuations over 2 generations
max2gen <- as.data.frame(MaxN[[2]])
min2gen <- as.data.frame(MinN[[2]])
min2gen <- replace(min2gen, min2gen==0, 1)
ratio2gen <- max2gen/min2gen
ExFluc2gen.EOO <- ratio2gen>=10
#Calculate Extreme Fluctuations over 3 generations
max3gen <- as.data.frame(MaxN[[3]])
min3gen <- as.data.frame(MinN[[3]])
min3gen <- replace(min3gen, min3gen==0, 1)
ratio3gen <- max3gen/min3gen
ExFluc3gen.EOO <- ratio3gen>=10


#calculate Extreme Fluctuations for AOO ExFluc.AOO  
MinN <- vector(mode='list',length = 0)
for(i in 2:3){
  mincalc <-data.frame()
  mincalc <- rbind(mincalc,AOO[1,])
  timeperiod <- gen[i]
    for (x in 2:(timeperiod-1)){
      data.rowx <- colmin(AOO[1:x,]) 
      mincalc <- rbind(mincalc,data.rowx)
    }
  rollmin <- rollapply(AOO, width = timeperiod, FUN = colmin, by.column = FALSE)
  rollmin <- as.data.frame(rollmin)
  names(mincalc)<-names(rollmin)
  mincalc <- rbind(mincalc, rollmin)
  MinN[[i]] <- mincalc
}
MaxN <- vector(mode='list',length = 0)
for(i in 2:3){
  maxcalc <-data.frame()
  maxcalc <- rbind(maxcalc,AOO[1,])
  timeperiod <- gen[i]
  for (x in 2:(timeperiod-1)){
    data.rowx <- colmax(AOO[1:x,]) 
    maxcalc <- rbind(maxcalc,data.rowx)
  }
  rollmax <- rollapply(AOO, width = timeperiod, FUN = colmax, by.column = FALSE)
  rollmax <- as.data.frame(rollmax)
  names(maxcalc)<-names(rollmax)
  maxcalc <- rbind(maxcalc, rollmax)
  MaxN[[i]] <- maxcalc
}
#Calculate Extreme Fluctuations over 2 generations
max2gen <- as.data.frame(MaxN[[2]])
min2gen <- as.data.frame(MinN[[2]])
min2gen <- replace(min2gen, min2gen==0, 1)
ratio2gen <- max2gen/min2gen
ExFluc2gen.AOO <- ratio2gen>=10
#Calculate Extreme Fluctuations over 3 generations
max3gen <- as.data.frame(MaxN[[3]])
min3gen <- as.data.frame(MinN[[3]])
min3gen <- replace(min3gen, min3gen==0, 1)
ratio3gen <- max3gen/min3gen
ExFluc3gen.AOO <- ratio3gen>=10

#calculate Extreme Fluctuations for number ccupied populalations ExFluc.subpopOcc  
MinN <- vector(mode='list',length = 0)
for(i in 2:3){
  mincalc <-data.frame()
  mincalc <- rbind(mincalc,subpopOcc[1,])
  timeperiod <- gen[i]
    for (x in 2:(timeperiod-1)){
      data.rowx <- colmin(subpopOcc[1:x,]) 
      mincalc <- rbind(mincalc,data.rowx)
    }
  rollmin <- rollapply(subpopOcc, width = timeperiod, FUN = colmin, by.column = FALSE)
  rollmin <- as.data.frame(rollmin)
  names(mincalc)<-names(rollmin)
  mincalc <- rbind(mincalc, rollmin)
  MinN[[i]] <- mincalc
}
MaxN <- vector(mode='list',length = 0)
for(i in 2:3){
  maxcalc <-data.frame()
  maxcalc <- rbind(maxcalc,subpopOcc[1,])
  timeperiod <- gen[i]
  for (x in 2:(timeperiod-1)){
    data.rowx <- colmax(subpopOcc[1:x,]) 
    maxcalc <- rbind(maxcalc,data.rowx)
  }
  rollmax <- rollapply(subpopOcc, width = timeperiod, FUN = colmax, by.column = FALSE)
  rollmax <- as.data.frame(rollmax)
  names(maxcalc)<-names(rollmax)
  maxcalc <- rbind(maxcalc, rollmax)
  MaxN[[i]] <- maxcalc
}
#Calculate Extreme Fluctuations over 2 generations
max2gen <- as.data.frame(MaxN[[2]])
min2gen <- as.data.frame(MinN[[2]])
min2gen <- replace(min2gen, min2gen==0, 1)
ratio2gen <- max2gen/min2gen
ExFluc2gen.subpopOcc <- ratio2gen>=10
#Calculate Extreme Fluctuations over 3 generations
max3gen <- as.data.frame(MaxN[[3]])
min3gen <- as.data.frame(MinN[[3]])
min3gen <- replace(min3gen, min3gen==0, 1)
ratio3gen <- max3gen/min3gen
ExFluc3gen.subpopOcc <- ratio3gen>=10

#calculate Extreme Fluctuations for number of mature individuals  ExFluc.Nmature  
MinN <- vector(mode='list',length = 0)
for(i in 2:3){
  mincalc <-data.frame()
  mincalc <- rbind(mincalc,Nmature[1,])
  timeperiod <- gen[i]
    for (x in 2:(timeperiod-1)){
      data.rowx <- colmin(Nmature[1:x,]) 
      mincalc <- rbind(mincalc,data.rowx)
    }
  rollmin <- rollapply(Nmature, width = timeperiod, FUN = colmin, by.column = FALSE)
  rollmin <- as.data.frame(rollmin)
  names(mincalc)<-names(rollmin)
  mincalc <- rbind(mincalc, rollmin)
  MinN[[i]] <- mincalc
}
MaxN <- vector(mode='list',length = 0)
for(i in 2:3){
  maxcalc <-data.frame()
  maxcalc <- rbind(maxcalc,Nmature[1,])
  timeperiod <- gen[i]
  for (x in 2:(timeperiod-1)){
    data.rowx <- colmax(Nmature[1:x,]) 
    maxcalc <- rbind(maxcalc,data.rowx)
  }
  rollmax <- rollapply(Nmature, width = timeperiod, FUN = colmax, by.column = FALSE)
  rollmax <- as.data.frame(rollmax)
  names(maxcalc)<-names(rollmax)
  maxcalc <- rbind(maxcalc, rollmax)
  MaxN[[i]] <- maxcalc
}
#Calculate Extreme Fluctuations over 2 generations
max2gen <- as.data.frame(MaxN[[2]])
min2gen <- as.data.frame(MinN[[2]])
min2gen <- replace(min2gen, min2gen==0, 1)
ratio2gen <- max2gen/min2gen
ExFluc2gen.Nmature <- ratio2gen>=10
#Calculate Extreme Fluctuations over 3 generations
max3gen <- as.data.frame(MaxN[[3]])
min3gen <- as.data.frame(MinN[[3]])
min3gen <- replace(min3gen, min3gen==0, 1)
ratio3gen <- max3gen/min3gen
ExFluc3gen.Nmature <- ratio3gen>=10

#Combine to find Extreme Fluctuations of any type

ExFluc2gen <- (ExFluc2gen.Nmature + ExFluc2gen.subpopOcc + ExFluc2gen.AOO + ExFluc2gen.EOO)>= 1
ExFluc3gen <- (ExFluc3gen.Nmature + ExFluc3gen.subpopOcc + ExFluc3gen.AOO + ExFluc3gen.EOO)>= 1


######################################
##################################
##########################Task 2: Create conditional statements
##IN FINAL TABLE
##   1 Extinct (EX)
##   2 Critically Endagered (CR)
##   3 Endangered (EN)
##   4 Vulnerable (VU)
##   5 Near Threatened (NT)


######################################Extinct##############################
# population is zero

EX<-Extinct*1
EX[EX==0]<-NA

##############################Critically Endangered#######################

## Criteria A ####
# Criteria A here is applied as if the observed trajectory from the beginning of the simulation
# up to the year for which the criteria is being applied is regarded as direct observation.
# Potential declines (future timesteps) are not considered. The criterion of threats being understood, reversible, 
# and ceased are assumed to not be met for simulations. Criteria A1, A3, and A4 are not considered. 
#
# CR.A.2 80% reduction over 3 generations or 10 years wichever is longer

CRA2<-apply(percentreduction[[3]],2,function(x){
  x >= 80
})*1
CRA2<-replace(CRA2, is.na(CRA2),0)
CRA<-CRA2

##Criteria B ##
### Criteria B is met if B1 OR B2 is met
# CR.B.1 extent of occurrence (EOO) is less than 100km AND at least two a-c criteria

CRB1.part1 <- apply(EOO,2,function(x) {
	x < 100
	})*1
#### a-c criteria calculated below (same for B2) 
   
# CR.B.2 area of occupancy (AOO), is less than 10km (part1) AND two of a-c criteria (part2)

#Part 1:AOO is less than 10km
CRB2.part1<-apply(AOO,2,function(x){
  x < 10
})*1

#Part 2: At least two of a-c:
#CR.B Part 2 a: Severely Fragmented OR known to exist at only a single location
CRB.part2.a<-(SevereFrag+singleLoc >=1)
#CR.B Part 2 b: Continuing decline (Calculated above combined to single var)
#CR.B Part 2 c: Extreme fluctuations (Calculated above combined to single var)
#NOTE: This is for Extreme fluctions over 2 generations - 3 gen is option
#Determine if 2 of a-c is met
CRB.part2<-(CRB.part2.a+ContDecline+ExFluc2gen >= 2)*1

# CR.B1 and CR.B2 (combine parts 1 and 2)
CRB1<-CRB1.part1*CRB.part2
CRB2<-CRB2.part1*CRB.part2

CRB<-CRB2 + CRB1 >=1

###Criteria C ######
### To meet this critera Part 1 AND Part 2 most both be met
### Critera C - Part 1: Population size estimated to number fewer than 250 mature 
CRC.part1<-apply(Nmature,2,function(x){
  x < 250
})*1

#Crtiteria C - Part 2:  Contdition 1 OR condition 2
#Cond 1. Continuing decline of at least 25% within 3 years
#or one generation whichever is longer OR 
#Cond 2. Continuing decline in N mature AND sub condition a or b: 
#a: popualtion structure such that condition (i) OR (ii)
# (i) there are no subpopulation with at least 50 mature 
# (ii)  at least 90% of mature indv. in single pop
#b: Extreme fluctuations in N mat (calculated over 2 generations change for 3)
CRC.part2.cond1<-apply(percentreduction[[1]],2,function(x){
  x >= 25
})*1
CRC.part2.cond1 <- replace(CRC.part2.cond1, is.na(CRC.part2.cond1),0)

CRC.part2.cond2.ai<-apply(LargestPopSize,2,function(x){
  x<=50
})*1

CRC.part2.cond2.aii<-apply(PropInLargestPop,2,function(x){
  x >= .90
})*1

CRC.part2.cond2.a <- (CRC.part2.cond2.ai + CRC.part2.cond2.aii >= 1)*1
CRC.part2.cond2.b <- ExFluc2gen.Nmature *1 
CRC.part2.cond2 <- (CRC.part2.cond2.a + CRC.part2.cond2.b >=1)
CRC.part2.cond2 <- (CRC.part2.cond2*ContDecline.lnNmat)*1
CRC.part2 <- (CRC.part2.cond1 + CRC.part2.cond2 >=1)*1

CRC<-CRC.part1*CRC.part2 

#D Population size estimated fewer 50 mature

CRD<-Nmature < 50
CRD <-CRD*1
#E SKIP!

CR<-(CRA+CRB+CRC+CRD)>=1
CR <- CR*1
#change 0's to NA's
CR[CR==0]<-NA
CR<-CR+1
#In the next step (i.e for Endangered" change 1 to 3, to create hierarchy)

##############################Endangered#######################

## Criteria A ####
# EN.A.1 70% reduction where cause is understood and ceased - not modeled
# EN.A.2 50% reduction over 3 generations or 10 years (and possibly ongoing)

ENA2<-apply(percentreduction[[3]],2,function(x){
  x >= 50
})*1
ENA2<-replace(ENA2, is.na(ENA2),0)
ENA<-ENA2

##Criteria B ##
### Criteria B is met if B1 OR B2 is met
#EN.B.1 extent of occurrence (EOO) is less than 5000km AND at least two a-c criteria
ENB1.part1<-apply(EOO,2,function(x){
  x < 5000
})*1    

#EN.B.2 area of occupancy (AOO), is less than 500km (part1) AND *two* of a-c criteria (part2)

#ENB2 Part 1:AOO is less than 500km
ENB2.part1<-apply(AOO,2,function(x){
  x < 500
})*1

#Part 2: At least two of a-c:
#ENB Part 2 a: Severely Fragmented OR known to exist at no more than 5 populations
ENB.part2.a<-(SevereFrag+fiveLoc >=1)
#ENB Part 2 b: Continuing decline (Calculated above combined to single var)
#ENB Part 2 c: Extreme fluctuations (Calculated above combined to single var)
#NOTE: This is for Extreme fluctions over 2 generations - 3 gen is option
#Determine if 2 of a-c is met
ENB.part2<-(ENB.part2.a+ContDecline+ExFluc2gen >= 2)*1

# ENB1 and ENB2 (including parts 1 and 2)
ENB1 <- ENB1.part1*ENB.part2
ENB2 <- ENB2.part1*ENB.part2

ENB<-ENB2 + ENB1 >= 1

###Criteria C ######
### To meet this critera Part 1 AND Part 2 most both be met
### Critera C - Part 1: Population size estimated to number fewer than 2500 mature 
ENC.part1<-apply(Nmature,2,function(x){
  x < 2500
})*1

#Crtiteria C - Part 2:  Contdition 1 OR condition 2
#Cond 1. Continuing decline of at least 20% within 5 years
#or two generations whichever is longer OR 
#Cond 2. Continuing decline in N mature AND sub condition a or b: 
#a: popualtion structure such that condition (i) OR (ii)
# (i) there are no subpopulation with at least 250 mature 
# (ii)  at least 95% of mature indv. in single pop
#b: Extreme fluctuations in N mat (calculated over 2 generations change for 3)
ENC.part2.cond1<-apply(percentreduction[[2]],2,function(x){
  x >= 20
})*1
ENC.part2.cond1 <- replace(ENC.part2.cond1, is.na(ENC.part2.cond1),0)

ENC.part2.cond2.ai<-apply(LargestPopSize,2,function(x){
  x<=250
})*1

ENC.part2.cond2.aii<-apply(PropInLargestPop,2,function(x){
  x >= .95
})*1

ENC.part2.cond2.a <- (ENC.part2.cond2.ai + ENC.part2.cond2.aii >= 1)*1
ENC.part2.cond2.b <- ExFluc2gen.Nmature *1 
ENC.part2.cond2 <- (ENC.part2.cond2.a + ENC.part2.cond2.b >=1)
ENC.part2.cond2 <- (ENC.part2.cond2*ContDecline.lnNmat)*1
ENC.part2 <- (ENC.part2.cond1 + ENC.part2.cond2 >=1)*1

ENC<-ENC.part1*ENC.part2 

#D Population size estimated fewer 250 mature

END<-Nmature < 250
END <-END*1
#E SKIP!

EN<-(ENA+ENB+ENC+END)>=1
EN <- EN*1
#change 0's to NA's
EN[EN==0]<-NA
EN<-EN+2

##############################Vulnerable#######################

## Criteria A ####
# VU.A.1 50% reduction where cause is understood and ceased - not modeled
# VU.A.2 30% reduction over 3 generations or 10 years (and possibly ongoing)

VUA2<-apply(percentreduction[[3]],2,function(x){
  x >= 30
})*1
VUA2<-replace(VUA2, is.na(VUA2),0)
VUA<-VUA2

##Criteria B ##
### Criteria B is met if B1 OR B2 is met
#VU.B.1 extent of occurrence (EOO) is less than 20,000km (part1) AND *two* of a-c criteria (part2)

VUB1.part1<-apply(EOO,2,function(x){
  x < 20000
})*1

#VU.B.2 area of occupancy (AOO), is less than 2,000km (part1) AND *two* of a-c criteria (part2)

#VUB2 Part 1:AOO is less than 2,000km
VUB2.part1<-apply(AOO,2,function(x){
  x < 2000
})*1

#Part 2: At least two of a-c:
#VUB Part 2 a: Severely Fragmented OR known to exist at no more than 10 locations
VUB.part2.a<-(SevereFrag+tenLoc >=1)
#VUB Part 2 b: Continuing decline (Calculated above combined to single var)
#VUB Part 2 c: Extreme fluctuations (Calculated above combined to single var)
#NOTE: This is for Extreme fluctions over 2 generations - 3 gen is option
#Determine if 2 of a-c is met
VUB.part2<-(VUB.part2.a+ContDecline+ExFluc2gen >= 2)*1

# VUA1 and VUB2 (including parts 1 and 2)
VUB1<-VUB1.part1*VUB.part2
VUB2<-VUB2.part1*VUB.part2

VUB<-VUB2 + VUB1 >= 1

###Criteria C ######
### To meet this critera Part 1 AND Part 2 most both be met
### Critera C - Part 1: Population size estimated to number fewer than 10,000 mature 
VUC.part1<-apply(Nmature,2,function(x){
  x < 10000
})*1

#Crtiteria C - Part 2:  Contdition 1 OR condition 2
#Cond 1. Continuing decline of at least 10% within 10 years
#or three generations whichever is longer OR 
#Cond 2. Continuing decline in N mature AND sub condition a or b: 
#a: popualtion structure such that condition (i) OR (ii)
# (i) there are no subpopulation with at least 1000 mature 
# (ii)  100% of mature indv. in single pop
#b: Extreme fluctuations in N mat (calculated over 2 generations change for 3)
VUC.part2.cond1<-apply(percentreduction[[3]],2,function(x){
  x >= 10
})*1
VUC.part2.cond1 <- replace(VUC.part2.cond1, is.na(VUC.part2.cond1),0)

VUC.part2.cond2.ai<-apply(LargestPopSize,2,function(x){
  x<=1000
})*1

VUC.part2.cond2.aii<-apply(PropInLargestPop,2,function(x){
  x >= 1
})*1

VUC.part2.cond2.a <- (VUC.part2.cond2.ai + VUC.part2.cond2.aii >= 1)*1
VUC.part2.cond2.b <- ExFluc2gen.Nmature *1 
VUC.part2.cond2 <- (VUC.part2.cond2.a + VUC.part2.cond2.b >=1)
VUC.part2.cond2 <- (VUC.part2.cond2*ContDecline.lnNmat)*1
VUC.part2 <- (VUC.part2.cond1 + VUC.part2.cond2 >=1)*1

VUC<-VUC.part1*VUC.part2 

#D Either (1) Population size estimated fewer 1000 mature  OR
# (2) AOO < 20km or Number of locations <= 5

VUD.1<-Nmature < 1000
if (ptc) { # Criteria 2 may not make sense to evaluate without spatial information
  VUD.2<-apply(AOO,2,function(x){
    x < 20
  })*1
  VUD.2<-(VUD.2 + fiveLoc)>=1
  VUD <-(VUD.1 + VUD.2 >=1)*1
} else {
  VUD <- VUD.1
  VUD.2 <- replace(VUD, VUD==TRUE, 'Not calculated') 
  VUD.2 <- replace(VUD.2, VUD.2==FALSE, 'Not calculated')
}
#E SKIP!

VU<-(VUA+VUB+VUC+VUD)>=1
VU <- VU*1
#change 0's to NA's
VU[VU==0]<-NA
VU<-VU+3

##############################Near Threatened#######################

## Criteria A ####
# NT - Decline over 3 generations of 20%

NTA2<-apply(percentreduction[[3]],2,function(x){
  x >= 20
})*1
NTA2<-replace(NTA2, is.na(NTA2),0)
NTA<-NTA2

##Criteria B ##
### Criteria B is met if B1 OR B2 is met
#Near threatened category is less precisely defined. Here it is defined as meeting the 
#area requirement for threatned with only one (rather than 2) of a-c
#OR
#meeting a NT threshold of 30,000km2 for EOO or 3,000km2 for AOO and 2 of a-c

#NTB1 Part 1:EOO is less than 30,000km
NTB1.part1<-apply(EOO,2,function(x){
  x < 30000
})*1

#NTB2 Part 1:AOO is less than 3,000km
NTB2.part1<-apply(AOO,2,function(x){
  x < 3000
})*1

#Part 2: At least one of a-c:
#NTB2 Part 2 a: Severely Fragmented OR known to exist at no more than 15 locations
NTB.part2.a<-(SevereFrag+fifteenLoc >=1)
#NTB Part 2 b: Continuing decline (Calculated above combined to single var)
#NTB Part 2 c: Extreme fluctuations (Calculated above combined to single var)
#NOTE: This is for Extreme fluctions over 2 generations - 3 gen is option
#Determine if any 1 or at least 2 of a-c is met
NTB.part2.any<-(NTB.part2.a+ContDecline+ExFluc2gen >= 1)*1
NTB.part2.two<-(NTB.part2.a+ContDecline+ExFluc2gen >= 2)*1

# NTB1 (including parts 1 and 2)
NTB1a<-VUB1.part1*NTB.part2.any
NTB1b<-NTB1.part1*NTB.part2.two

NTB1<-(NTB1a+NTB1b>=1)*1

# NTB2 (including parts 1 and 2)
NTB2a<-VUB2.part1*NTB.part2.any
NTB2b<-NTB2.part1*NTB.part2.two

NTB2<-(NTB2a+NTB2b>=1)*1


NTB<-NTB2 + VUB1 >= 1

###Criteria C ######
### To meet this critera Part 1 AND Part 2 most both be met
### Critera C - Part 1: Population size estimated to number fewer than 15,000 mature 
NTC.part1<-apply(Nmature,2,function(x){
  x < 15000
})*1

#Crtiteria C - Part 2:  Contdition 1 OR condition 2
#Cond 1. Continuing decline of at least 10% within 10 years
#or three generations whichever is longer (note: this is the same as VU) OR 
#Cond 2. Continuing decline in N mature AND sub condition a or b: 
#a: popualtion structure such that condition (i) OR (ii)
# (i) there are no subpopulation with at least 1500 mature 
# (ii)  100% of mature indv. in single pop (note: this is the same as VU)
#b: Extreme fluctuations in N mat (calculated over 2 generations change for 3)


NTC.part2.cond2.ai<-apply(LargestPopSize,2,function(x){
  x<=1500
})*1


NTC.part2.cond2.a <- (NTC.part2.cond2.ai + VUC.part2.cond2.aii >= 1)*1
NTC.part2.cond2 <- (NTC.part2.cond2.a + VUC.part2.cond2.b >=1)
NTC.part2.cond2 <- (NTC.part2.cond2*ContDecline.lnNmat)*1
NTC.part2 <- (VUC.part2.cond1 + NTC.part2.cond2 >=1)*1

NTC<-NTC.part1*NTC.part2 

#D Either (1) Population size estimated fewer 1500 mature  OR
# (2) AOO < 30km or Number of locations <= 5

NTD.1<-Nmature < 1500
if (ptc) { # Criteria 2 may not make sense to evaluate without spatial information
  NTD.2<-apply(AOO,2,function(x){
    x <= 30
  })*1
  NTD.2<-(NTD.2 + fiveLoc)>=1
  NTD <-(NTD.1 + NTD.2 >=1)*1
} else {
  NTD <- NTD.1
}
#E SKIP!

NT<-(NTA+NTB+NTC+NTD)>=1
NT <- NT*1
#change 0's to NA's
NT[NT==0]<-NA
NT<-NT+4
#############################################################################
redlist.cat <- pmin(EX,CR,EN,VU,NT, na.rm =TRUE) #combine categories into single table
redlist.cat <- replace(redlist.cat, is.na(redlist.cat), 6)

#melt all the output variable tables so they can be assembled 
rlcat.mltd<-melt(redlist.cat)
CRA2.mltd<-melt(CRA2)
CRB1.mltd <- melt(CRB1)
CRB2.mltd <- melt(CRB2)
CRC.mltd <- melt(CRC)
CRD.mltd <- melt(CRD)
ENA2.mltd<-melt(ENA2)
ENB1.mltd <- melt(ENB1)
ENB2.mltd <- melt(ENB2)
ENC.mltd <- melt(ENC)
END.mltd <- melt(END)
VUA2.mltd<-melt(VUA2)
VUB1.mltd <- melt(VUB1)
VUB2.mltd <- melt(VUB2)
VUC.mltd <- melt(VUC)
VUD.mltd <- melt(VUD)
NTA2.mltd<-melt(NTA2)
NTB1.mltd <- melt(NTB1)
NTB2.mltd <- melt(NTB2)
NTC.mltd <- melt(NTC)
NTD.mltd <- melt(NTD)

ExFluc.mltd <- melt(ExFluc2gen.Nmature)
percred1gen <- as.data.frame(percentreduction[1])
percred2gen <- as.data.frame(percentreduction[2])
percred3gen <- as.data.frame(percentreduction[3])
colnames(percred1gen) <- c(1:RepRun)
colnames(percred2gen) <- c(1:RepRun)
colnames(percred3gen) <- c(1:RepRun)
percred1gen.mltd <- melt(percred1gen)
percred2gen.mltd <- melt(percred2gen)
percred3gen.mltd <- melt(percred3gen)

ContDec.mltd <- melt(ContDecline.lnNmat)
EOO.mltd <- melt(EOO)
ContDecline.EOO.mltd <- melt(ContDecline.EOO)
ExFluc2gen.EOO.mltd <- melt(ExFluc2gen.EOO)
AOO.mltd <- melt(AOO)
ContDecline.AOO.mltd <- melt(ContDecline.AOO)
ExFluc2gen.AOO.mltd <- melt(ExFluc2gen.AOO)
VUD.2.mltd <- melt(VUD.2)

SubpopCount <- replace(subpopOcc, subpopOcc==0, 1)
SubpopCount.mltd <- melt(SubpopCount)

ContDecline.subpopOcc.mltd <- melt(ContDecline.subpopOcc)
ExFluc2gen.subpopOcc.mltd <- melt(ExFluc2gen.subpopOcc)
LargestPopSize.mltd <- melt(LargestPopSize)
SevereFrag.mltd <- melt(SevereFrag)
singleLoc.mltd <- melt(singleLoc)

# create a lookup table to match RL category to name
cat.lookup <- cbind(rep('EX',nrow(rlcat.mltd)), rep('CR',nrow(rlcat.mltd)), 
					rep('EN',nrow(rlcat.mltd)), rep('VU',nrow(rlcat.mltd)), rep('NT',nrow(rlcat.mltd)), rep('LC', nrow(rlcat.mltd)))
					
# create a lookup table for each criteria to match the qualifying criteria to the RL category
Crit.A2.lookup <- cbind( rep(0,nrow(rlcat.mltd)), CRA2.mltd[,3], ENA2.mltd[,3],
						VUA2.mltd[,3], NTA2.mltd[,3], rep(0,nrow(rlcat.mltd)))
Crit.A2.lookup <- Crit.A2.lookup==1  # Convert to TRUE/FALSE
Crit.B1.lookup <-cbind( rep(0,nrow(rlcat.mltd)), CRB1.mltd[,3], ENB1.mltd[,3],
						VUB1.mltd[,3], NTB1.mltd[,3], rep(0,nrow(rlcat.mltd)))
Crit.B1.lookup <- Crit.B1.lookup==1  # Convert to TRUE/FALSE
Crit.B2.lookup <-cbind( rep(0,nrow(rlcat.mltd)), CRB2.mltd[,3], ENB2.mltd[,3],
						VUB2.mltd[,3], NTB2.mltd[,3], rep(0,nrow(rlcat.mltd)))
Crit.B2.lookup <- Crit.B2.lookup==1  # Convert to TRUE/FALSE
Crit.C.lookup <- cbind( rep(0,nrow(rlcat.mltd)), CRC.mltd[,3], ENC.mltd[,3],
						VUC.mltd[,3], NTC.mltd[,3], rep(0,nrow(rlcat.mltd)))
Crit.C.lookup <- Crit.C.lookup==1  # Convert to TRUE/FALSE
Crit.D.lookup <- cbind( rep(0,nrow(rlcat.mltd)), CRD.mltd[,3], END.mltd[,3],
						VUD.mltd[,3], NTD.mltd[,3], rep(0,nrow(rlcat.mltd)))
Crit.D.lookup <- Crit.D.lookup==1  # Convert to TRUE/FALSE
					
## use lookup table to match criteria to highest risk category 
rl.category <- vector()
Crit.A2 <- vector()
Crit.B1 <- vector()
Crit.B2 <- vector()
Crit.C <- vector()
Crit.D <- vector()

for (i in 1:nrow(rlcat.mltd)){
  cat.i <- cat.lookup[i, rlcat.mltd[i,3]]
  critA2.i <- Crit.A2.lookup[i, rlcat.mltd[i,3]]
  critB1.i <- Crit.B1.lookup[i, rlcat.mltd[i,3]]
  critB2.i <- Crit.B2.lookup[i, rlcat.mltd[i,3]]
  critC.i <- Crit.C.lookup[i, rlcat.mltd[i,3]]
  critD.i <- Crit.D.lookup[i, rlcat.mltd[i,3]]
  rl.category <- c(rl.category,cat.i)
  Crit.A2 <- c(Crit.A2, critA2.i)
  Crit.B1 <- c(Crit.B1, critB1.i)
  Crit.B2 <- c(Crit.B2, critB2.i)
  Crit.C <- c(Crit.C, critC.i)
  Crit.D <- c(Crit.D, critD.i)
}

redlist.pars <- data.frame( Model = rep(as.numeric(gsub("\\D","",mp.file.name)),length(rl.category)),
							FileName = rep(mp.file.name,length(rl.category)),
							Species = rep(substr(mp.file.name,1,4),length(rl.category)),
							Replicate = rlcat.mltd[,2], Year = rlcat.mltd[,1],
							AssessmentCategory = rl.category, Crit.A2 = Crit.A2,
							Crit.B1 = Crit.B1, Crit.B2 = Crit.B2, Crit.C = Crit.C,
							Crit.D = Crit.D, GenerationTime = rep(GTime, length(rl.category)),
							ExtremePopFluctuations = ExFluc.mltd[,3], PastReduction = percred3gen.mltd[,2],
							ContinuingDecline = ContDec.mltd[,3], ContinuingDeclineCR = percred1gen.mltd[,2],
							ContinuingDeclineEN = percred2gen.mltd[,2], ContinuingDeclineVU = percred3gen.mltd[,2],
							OccurrenceExtent = round(EOO.mltd[,3]), OccurrenceContinuingDecline = ContDecline.EOO.mltd[,2],
							OccurrenceExtremeFluctuations = ExFluc2gen.EOO.mltd[,3], 
							OccupancyArea = round(AOO.mltd[,3]), OccupancyContinuingDecline = ContDecline.AOO.mltd[,2],
							OccupancyExtremeFluctuations = ExFluc2gen.AOO.mltd[,3], VeryRestricted = VUD.2.mltd[,3],
							SubpopCount = SubpopCount.mltd, SubpopContinuingDecline = ContDecline.subpopOcc.mltd[,2],
							SubpopExtremeFluctuations = ExFluc2gen.subpopOcc.mltd[,3], LargestSubpopSize = LargestPopSize.mltd[,3],
							SeverelyFragmented = SevereFrag.mltd[,3], AllInOnePop = singleLoc.mltd[,3])
return(redlist.pars)               
} #End function
 










