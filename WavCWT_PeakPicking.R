# The following function are for CentWave based peak picking, which include three modules: 
#1) WavCWT    : CWT Transform 1D signal data into 2D coefficient plane {scale, time}
#2) WavCWTTree: Get the CWT widget Tree (include two submodules: (a) WTMM: try to find the 2D map of local maximum/minimum values point in the CWT plane, 1: local maximum/minum value, 0: not the maximum/minimum value 
#                                                                (b) wtmmBranches: Try to find the CWT_Tree branch according to the 2D CWT Coefficients and the local maximum/minimum map)
#3) WavCWTPeak: Identify the peak according to the CWT widget Tree.
wavCWT_Modify<-function (x, scale.range = deltat(x) * c(1, length(x)), n.scale = 100, 
		wavelet = "gaussian2", shift = 5, variance = 1) 
{	
	
#	cat("Entering wavCWT_Modify \n")
	
	checkVectorType(scale.range, "numeric")
	checkScalarType(n.scale, "integer")
	checkScalarType(wavelet, "character")
	checkScalarType(shift, "numeric")
	checkScalarType(variance, "numeric")
	checkRange(n.scale, c(1, Inf))
	series.name <- deparse(substitute(x))
	if (length(scale.range) != 2) 
		stop("scale.range must be a two-element numeric vector")
	if (variance <= 0) 
		stop("variance input must be positive")
	sampling.interval <- deltat(x)
	octave <- logb(scale.range, 2)
	scale <- ifelse1(n.scale > 1, 2^c(octave[1] + seq(0, n.scale - 
									2) * diff(octave)/(floor(n.scale) - 1), octave[2]), scale.range[1])
	scale <- unique(round(scale/sampling.interval) * sampling.interval)
	n.scale <- length(scale)
	if (abs(min(scale) - sampling.interval) > .Machine$double.eps) 
		stop("Minimum scale must be greater than or equal to sampling interval ", 
				"of the time series")
	if (inherits(x, "signalSeries")) 
		times <- as(x@positions, "numeric")
	else times <- time(x)
	x <- as.vector(x)
	storage.mode(x) <- "double"
	gauss1 <- c("gaussian1", "gauss1")
	gauss2 <- c("gaussian2", "gauss2", "mexican hat", "sombrero")
	supported.wavelets <- c("haar", gauss1, gauss2, "morlet")
	wavelet <- match.arg(lowerCase(wavelet), supported.wavelets)
	filter <- mutilsFilterTypeContinuous(wavelet)
	if (filter == 4) {
		filter.arg <- sqrt(variance)
		wavelet <- "gaussian1"
	}
	else if (filter == 5) {
		filter.arg <- sqrt(variance)
		wavelet <- "gaussian2"
	}
	else if (filter == 6) {
		filter.arg <- shift
		wavelet <- "morlet"
	}
	else if (filter == 7) {
		filter.arg <- 0
		wavelet <- "haar"
		scale <- sampling.interval * unique(round(scale/sampling.interval))
	}
	else stop("Unsupported filter type")
	z <- .Call("RS_wavelets_transform_continuous_wavelet", as.numeric(x), 
			as.numeric(sampling.interval), as.integer(filter), as.numeric(filter.arg), 
			as.numeric(scale), COPY = rep(FALSE, 5), CLASSES = c("numeric", 
					"numeric", "integer", "numeric", "numeric"), PACKAGE = "ifultools")
	if (wavelet != "morlet") 
		z <- Re(z)
	attr(z, "scale") <- scale
	attr(z, "time") <- as.vector(times)
	attr(z, "wavelet") <- wavelet
	attr(z, "series") <- x
	attr(z, "sampling.interval") <- sampling.interval
	attr(z, "series.name") <- series.name
	attr(z, "n.sample") <- length(x)
	attr(z, "n.scale") <- n.scale
	attr(z, "filter.arg") <- filter.arg
	oldClass(z) <- "wavCWT"
	z
}

wavCWTTree_Modify<-function (x, n.octave.min = 1, tolerance = 0, type = "maxima") 
{		
	
#	cat("In wavCWTTree_Modify \n")
#	browser()
	
	"WTMM" <- function(x, tolerance = NULL, type = "maxima") {
		
		if (!is(x, "wavCWT")) 
			stop("Input object must be of class wavCWT")
		x.attr <- attributes(x)
		times <- x.attr$time
		scales <- x.attr$scale
		n.sample <- x.attr$n.sample
		series <- x.attr$series
		if (is.null(tolerance)) {
			tolerance <- mad(Mod(x[, 1]))/scales
		}
		if (length(tolerance) < length(scales)) 
			tolerance <- tolerance[1]/sqrt(scales)
		wtmmz <- .Call("RS_wavelets_transform_continuous_wavelet_modulus_maxima", 
				as.matrix(x) + (0+0i), tolerance, mutilsTransformPeakType(type), 
				CLASSES = c("matrix", "numeric", "integer"), COPY = rep(FALSE, 
						3), PACKAGE = "ifultools")
		z <- matrix(0, nrow = nrow(x), ncol = ncol(x))
		z[matrix(unlist(wtmmz), ncol = 2) + 1] <- 1
		z
	}
	"wtmmBranches" <- function(wtmm, extrema.mask, times, scales, 
			span.min = 5, gap.max = 3, skip = NULL, sampling.interval = 1) {   
				
		scales <- as.integer(scales/sampling.interval)
		n.scale <- ncol(extrema.mask)
		n.sample <- nrow(extrema.mask)
		if (is.null(scales)) 
			scales <- 1:n.scale
		iwtmm <- which(extrema.mask[, n.scale] > 0)
		iscale <- seq(n.scale - 1, 1, -1)
		tree <- as.list(iwtmm)
		names(tree) <- iwtmm
		peakStatus <- as.list(rep(0, length(iwtmm)))
		names(peakStatus) <- iwtmm
		orphanRidgeList <- NULL
		orphanRidgeName <- NULL
		n.level <- length(iscale)
		for (j in seq(n.level)) {
			iscale.j <- iscale[j]
			scale.j <- scales[iscale.j]
			if (length(iwtmm) == 0) {
				iwtmm <- which(extrema.mask[, iscale.j] > 0)
				next
			}			
			span <- scale.j * 2 + 1
			if (span < span.min) 
				span <- span.min
			remove.j <- selPeak.j <- NULL
			for (k in seq(along = iwtmm)) {
				itime <- iwtmm[k]
				itime.start <- itime - span
				if (itime.start < 1) 
					itime.start <- 1
				itime.end <- itime + span
				if (itime.end > n.sample) 
					itime.end <- n.sample
				itime.candidates <- which(extrema.mask[itime.start:itime.end, 
								iscale.j] > 0) + itime.start - 1
				if (length(itime.candidates) == 0) {
					status.k <- peakStatus[[as.character(itime)]]
					if (status.k > gap.max & scale.j >= 2) {
						temp <- tree[[as.character(itime)]]
						orphanRidgeList <- c(orphanRidgeList, list(temp[1:(length(temp) - 
															status.k)]))
						orphanRidgeName <- c(orphanRidgeName, paste(iscale.j + 
												status.k + 1, itime, sep = "_"))
						remove.j <- c(remove.j, as.character(itime))
						next
					}
					else {
						itime.candidates <- itime
						peakStatus[[as.character(itime)]] <- status.k + 
								1
					}
				}
				else {
					peakStatus[[as.character(itime)]] <- 0
					if (length(itime.candidates) >= 2) 
						itime.candidates <- itime.candidates[which.min(abs(itime.candidates - 
														itime))]
				}
				tree[[as.character(itime)]] <- c(tree[[as.character(itime)]], 
						itime.candidates)
				selPeak.j <- c(selPeak.j, itime.candidates)
			}
			if (length(remove.j) > 0) {
				bad.tree <- which(is.element(names(tree), remove.j))
				tree <- tree[-bad.tree]
				peakStatus <- peakStatus[-bad.tree]
			}
			dupPeak.j <- unique(selPeak.j[duplicated(selPeak.j)])
			if (length(dupPeak.j) > 0) {
				bad.tree <- NULL
				for (dupPeak.jk in dupPeak.j) {
					selInd <- which(selPeak.j == dupPeak.jk)
					selLen <- sapply(tree[selInd], length)
					bad.tree.jk <- which.max(selLen)
					bad.tree <- c(bad.tree, selInd[-bad.tree.jk])
					orphanRidgeList <- c(orphanRidgeList, tree[bad.tree.jk])
					orphanRidgeName <- c(orphanRidgeName, paste(iscale.j, 
									selPeak.j[bad.tree.jk], sep = "_"))
				}
				selPeak.j <- selPeak.j[-bad.tree]
				tree <- tree[-bad.tree]
				peakStatus <- peakStatus[-bad.tree]
			}
			names(tree) <- selPeak.j
			names(peakStatus) <- selPeak.j
			if (scale.j >= 2) {
				maxInd.next <- which(extrema.mask[, iscale.j] > 
								0)
				unSelPeak.j <- maxInd.next[!is.element(maxInd.next, 
								selPeak.j)]
				newPeak.j <- as.list(unSelPeak.j)
				names(newPeak.j) <- unSelPeak.j
				tree <- c(tree, newPeak.j)
				iwtmm <- c(selPeak.j, unSelPeak.j)
				newPeakStatus <- as.list(rep(0, length(newPeak.j)))
				names(newPeakStatus) <- newPeak.j
				peakStatus <- c(peakStatus, newPeakStatus)
			}
			else {
				iwtmm <- selPeak.j
			}
		}
		names(tree) <- paste(1, names(tree), sep = "_")
		names(orphanRidgeList) <- orphanRidgeName
		tree <- c(tree, orphanRidgeList)
		tree <- lapply(tree, rev)
		tree <- tree[unique(names(tree))]
		tree <- lapply(seq(along = tree), function(i, tree, iscale.min, 
						times, scales, wtmm) {
					itime <- tree[[i]]
					iscale <- seq(iscale.min[i], length = length(itime))
					list(itime = itime, iscale = iscale, time = times[itime], 
							scale = scales[iscale], extrema = wtmm[cbind(itime, 
											iscale)])
				}, tree = tree, iscale.min = as.integer(gsub("_.*", "", 
								names(tree))), times = times, scales = scales * sampling.interval, 
				wtmm = wtmm)
		iflat <- lapply(tree, function(x, nr) (x$iscale - 1) * 
							nr + x$itime, nr = nrow(wtmm))
		flatset <- iflat[[1]]
		bad <- NULL
		if(length(iflat)>1)
		{
			for (i in seq(2, length(iflat))) {
				if (any(is.element(iflat[[i]], flatset))) 
					bad <- c(bad, i)
				else flatset <- c(flatset, iflat[[i]])
			}
		}
		
		if (length(bad) > 0) 
			tree <- tree[-bad]
		tree
	}
	x.attr <- attributes(x)
	times <- x.attr$time
	scales <- x.attr$scale
	n.sample <- x.attr$n.sample
	sampling.interval <- x.attr$sampling.interval
	border.times <- range(times) + sampling.interval * c(1, -1)
	extrema.mask <- WTMM(x, tolerance = tolerance, type = type)	
	
	if (!identical(dim(x), dim(extrema.mask))) 
		stop("Input WTMM dimenions do not match those of the input CWT matrix")
		 
	z <- wtmmBranches(ifelse1(is.complex(x), Mod(as.matrix(x)), 
					as.matrix(x)), extrema.mask, times, scales, sampling.interval = sampling.interval)
	n.scale <- length(scales)
	n.octave <- log2(max(scales)/min(scales))
	n.voice <- (n.scale - 1)/n.octave
	n.scale.min <- as.integer(n.voice * n.octave.min)
	good <- which(unlist(lapply(z, function(x, n.scale.min) length(x[[1]]) > 
										n.scale.min, n.scale.min = n.scale.min)))
	z <- z[good]
	endtime <- unlist(lapply(z, function(x, iscale) x$itime[iscale], 
					iscale = which.min(scales)))
	isort <- order(endtime)
	z <- z[isort]
	names(z) <- seq(z)
	attr(z, "iendtime") <- endtime[isort]
	attr(z, "endtime") <- times[endtime[isort]]
	attr(z, "time") <- times
	attr(z, "scale") <- scales
	attr(z, "extrema.mask") <- extrema.mask
	attr(z, "noise") <- x[, 1]
	attr(z, "branch.hist") <- colSums(extrema.mask * abs(x))
	attr(z, "wavelet") <- attr(x, "wavelet")
	attr(z, "filter.arg") <- attr(x, "filter.arg")
	attr(z, "series.name") <- attr(x, "series.name")
	attr(z, "series") <- attr(x, "series")
	attr(z, "sampling.interval") <- attr(x, "sampling.interval")
	oldClass(z) <- "wavCWTTree"
	z
}

wavCWTPeaks_Modify <- function( x, snr.min = 3, scale.range = NULL, length.min = 3, 
		noise.span = NULL, noise.fun = "quantile", noise.min = NULL)
{	
#	cat("In wavCWTPeaks_Modify \n")
#	browser()
	
	if (!is(x, "wavCWTTree")) 
		stop("Input must be an object of class wavCWTTree")
	xatt <- attributes(x)
	endtimes <- attr(x, "endtime")
	times <- attr(x, "time")
	scale <- attr(x, "scale")
	noise <- attr(x, "noise")
	wavelet <- attr(x, "wavelet")
	series <- attr(x, "series")
	branch.hist <- attr(x, "branch.hist")
	sampling.interval <- abs(diff(times[1:2]))
	if (!is.element(wavelet, "gaussian2")) 
		stop("Only CWT developed using the Mexican hat (gaussian2) filter are supported")
	if (is.null(noise.min)) 
		noise.min <- quantile(abs(attr(x, "noise")), prob = 0.05)
	
#	#############################################################################
#	Modified on 6/1/2012 by Xiuxia: Decrease prob to accommodate narrow peaks
	if (is.null(scale.range)) 
		scale.range <- scale[range(which(branch.hist > quantile(branch.hist, 
										prob = 0.3)))]
#		scale.range <- scale[range(which(branch.hist > quantile(branch.hist, 
#										prob = 0.5)))]
#	#############################################################################
	
	if (is.null(noise.span)) 
		noise.span <- max(0.01 * diff(range(times)), 5 * sampling.interval)
	noise.levels <- unlist(lapply(endtimes, function(x, noise.fun, 
							times, times.range, noise, noise.min, noise.span) {
						time.start <- x - noise.span
						if (time.start < times.range[1]) 
							time.start <- times.range[1]
						time.end <- x + noise.span
						if (time.end < times.range[2]) 
							time.end <- times.range[2]
						ix <- which(times >= time.start & times <= time.end)
						noise.local <- noise.fun(abs(noise[ix]))
						if (noise.local < noise.min) 
							noise.local <- noise.min
						noise.local
					}, noise.fun = switch(noise.fun, quantile = function(x) {
								quantile(x, probs = 0.95)
							}, sd = sd, mad = function(x) {
								mad(x, center = 0)
							}), times = times, times.range = range(times), noise = noise, 
					noise.min = noise.min, noise.span = noise.span))
	tmpargs <- lapply(x, function(x) unlist(lapply(x, function(x, 
										imax) x[imax], imax = which.max(x$extrema))))
	peaks <- data.frame(do.call("rbind", tmpargs))
	peaks <- cbind(data.frame(branch = row.names(peaks)), peaks, 
			data.frame(iendtime = attr(x, "iendtime")))
	peak.snr <- peaks[["extrema"]]/noise.levels
		
#	cat(paste("extrema =", peaks[["extrema"]], "\n"))
#	
#	cat(paste("noise =", noise.levels, "\n"))
	
	peak.scale <- peaks[["scale"]]
	branch.lengths <- unlist(lapply(x, function(x, scale.range) length(which(x$scale >= 
												scale.range[1] & x$scale <= scale.range[2])), scale.range = scale.range))
	good.snr <- peak.snr >= snr.min
	
#	6/1/2012: Xiuxia changed upper scale to accommodate wide peaks
	good.scale <- peak.scale >= scale.range[1]  
#	good.scale <- ((peak.scale >= max(scale.range[1],3))& (peak.scale <= min(scale.range[2], 64)))
#	good.scale <- ((peak.scale >= max(scale.range[1],3))& (peak.scale <= min(scale.range[2], 25)))   #Modifed by Wenchao Zhang, Use lower scale boundary and Upper scale boundary to filter  
	


	good.length <- branch.lengths >= length.min
	iendtime.min <- max(as.integer(noise.span/sampling.interval/4), 3)
	iendtime.max <- length(times) - iendtime.min + 1
	good.end <- peaks[["iendtime"]] > iendtime.min & peaks[["iendtime"]] < iendtime.max
	
	# Added some codes to filter out some very very small peaks. Added by wenchao zhang
	peak_extrema <- peaks[["extrema"]]
	
	good.peak_extrema <- peak_extrema> max(peak_extrema)/8.0
#	good.peak_extrema <- peak_extrema> max(peak_extrema)/40.0	
	peaks <- peaks[which(good.snr & good.scale & good.length & good.end & good.peak_extrema), ]	
	
	#Save the snr of peaks, added in 10.07.2011
    peaks_snr<-peak.snr[which(good.snr & good.scale & good.length & good.end & good.peak_extrema)] 
	
	#Debug for none good peak cases, 
	if(nrow(peaks)>0) row.names(peaks) <- as.character(seq(nrow(peaks)))
	#z <- list(x = times[peaks$iendtime], y = series[peaks$iendtime])
	z <- list(x = times[peaks$iendtime], y = series[peaks$iendtime], snr = peaks_snr)
	attr(z, "peaks") <- peaks
	attr(z, "snr.min") <- snr.min
	attr(z, "scale.range") <- scale.range
	attr(z, "length.min") <- length.min
	attr(z, "noise.span") <- noise.span
	attr(z, "noise.fun") <- noise.fun
	attr(z, "noise.min") <- noise.min
		
	att<- attr(z, "peaks")
	z
}