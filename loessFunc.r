library("ncdf")

	

period <- 120
x <- 1:120
y <- sin(2*pi*x/period) + 10   # runif(length(x),-1,1)

x<-1:20
y<-c(0,0,0,0,0,0,0,0,0,10,9,8,7,6,5,4,3,2,1,0)

playwith({
			plot(x,y, main="Sine Curve + 'Uniform' Noise")
			mtext("showing loess smoothing (local regression smoothing)")
		}, new=TRUE) 

y.loess <- loess(y ~ x, span=0.75, data.frame(x=x, y=y))


playwith({
			plot(x,curVecInt, main="loess_fitted")
			mtext("showing loess smoothing fitted(local regression smoothing)")
		}, new=TRUE) 

y.predict <- predict(y.loess, data.frame(x=x))

###############
	#read EIC data
	################
	#inFilePath <- DataFilelist[[fileindex]]
	fileName<-parseFileName(inFilePath)
	EICFile<-paste(WorkDir,"/output/EIC/",fileName,"EIC.cdf",sep="")
	cat(fileName,"reading EIC data...\n")
	
	ncid <- open.ncdf(EICFile)
	vecInt <- get.var.ncdf(ncid, varid="intVec")
	mzVec<-get.var.ncdf(ncid, varid="mzVec")
	totalscan<-length(vecInt)/length(mzVec)
	
	##############
	#smoothing
	################
	cat(fileName,"smoothing EIC data...\n")
	
	DenoisedTIC <- 0
	
	time1<-Sys.time()
	for(mzInd in 1:length(mzVec))
	{
		mz<-mzVec[mzInd]	
		startInd<-(mzInd-1)*totalscan+1
		endInd<-startInd+totalscan-1
		curVecInt <- vecInt[startInd:endInd]
		
		if(isSm==T)
		{
			smoothingWindow<-20
			ma = rep(1, smoothingWindow)/smoothingWindow
			sn.ma=filter(curVecInt,ma)
			sn.ma[which(is.na(sn.ma))]<-0
	#
#				playwith({
#					plot((1:length(curVecInt))*params$ScanInterval+params$delaytime,curVecInt,type="l",main=paste(smoothingWindow,"-Point Moving Average Filter",sep=""))
#					points((1:length(curVecInt))*params$ScanInterval+params$delaytime,sn.ma,type="l",col="red")
#					})
			curVecInt<-sn.ma
		}	
		##########
		#baseline correction
		###########
		if(isBslnCrt==T)
		{   
			#nbinSize<-240
		    nbinSize<-5
			curVecInt <- y
			totalscan <- 20
			
			isMin <- peaks(-curVecInt, span=nbinSize)
			minInd<-which(isMin==TRUE)	
			
			if (length(minInd)==0) {bsln=0} # avoiding no minInd found - yan
			else {
				intervalOfScans<-c(minInd[1],minInd[-1]-minInd[-length(minInd)])
				
				bgs<-c(rep(curVecInt[minInd],intervalOfScans),rep(0,totalscan-minInd[length(minInd)]))
				f.lo <- loess(bgs ~ c(1:totalscan), span =0.05, degree = 2)
				bsln <- f.lo$fitted
			}
			
#			playwith(
#					{ 
#					
#						plot((1:length(curVecInt))*params$ScanInterval+params$delaytime,curVecInt, type = "l", col = "blue")
#						points((1:length(curVecInt))*params$ScanInterval+params$delaytime,bgs,type="l")
#						points((1:length(curVecInt))*params$ScanInterval+params$delaytime,bsln, type="l",col = "red")
#					})
			###baseline subtraction
			curVecInt<-curVecInt-bsln
			curVecInt[curVecInt<0]<-0
		}
		
		#####update vector of intensity
		vecInt[startInd:endInd]<-curVecInt
		DenoisedTIC <- DenoisedTIC + curVecInt
	}
	
	

peaks<-function (x, span = 3) 
{
	z <- embed(as.vector(x), span)
	s <- span%/%2
	result <- max.col(z,ties.method="first")== s+1
	c(rep(FALSE, s), result, rep(FALSE, s))      # (c(rep(FALSE, s), result, rep(FALSE, s))    
}