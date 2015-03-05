#Estimate and Remove the baseline according to the local valley(local minimum value)
peaks<-function (x, span = 3) 
{			
	z <- embed(as.vector(x), span)
	s <- span%/%2
	if(span%%2 ==0) #span is even number, such as 4, 6, 8,... 100, etc
	   result <- max.col(z,ties.method="first") == s
    else  # span is odd number, such as 3,5,9,....
	   result <- max.col(z,ties.method="first") == 1+s	
 	c(rep(FALSE, s), result, rep(FALSE, (s-1)))         # c(rep(FALSE, s), result, rep(FALSE, s))
}

BaseLine_Remove <- function(EIC_Smooth, nbinSize =6)
{
	totalpoint <- length(EIC_Smooth)
	
	isMin <- peaks(-EIC_Smooth, span=nbinSize)
	minInd<-which(isMin==TRUE)
	
	if (length(minInd)==0) 
	   {bsln=0} # avoiding no minInd found - yan
	else {
		intervalOfScans<-c(minInd[1],minInd[-1]-minInd[-length(minInd)])
		
		bgs<-c(rep(EIC_Smooth[minInd],intervalOfScans),rep(0,totalpoint-minInd[length(minInd)]))
		f.lo <- loess(bgs ~ c(1:totalpoint), span =0.05, degree = 2)
		bsln <- f.lo$fitted
	}	

	###baseline subtraction
	EIC_Baseline_Remove <-EIC_Smooth-bsln
	Neg_Ind<-which(EIC_Baseline_Remove< 0)
	EIC_Baseline_Remove[Neg_Ind]<- 0 #After removing baseline, some data point can be negtative, need to keep them as 0.  
	EIC_Baseline_Remove
}