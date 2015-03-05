#Extract the corresponding EIC by Bining the points between mass range of (mass_upper, mass_bottom)
EIC_Binning<- function(mass_values, intensity_values, scan_index, scan_number, mass_bottom, mass_upper)
{
	Start_Scan <- min(scan_number)
	End_Scan <- max(scan_number-1)
	EIC <- numeric(length = End_Scan-Start_Scan+1)
	
	for(i in Start_Scan:End_Scan)
	{
		ind_start<- scan_index[i]+1
		ind_end  <- scan_index[i+1]
		mass_values_tempscan <- mass_values[ind_start:ind_end]
		intensity_values_tempscan <- intensity_values[ind_start:ind_end]
		ind <- which((mass_values_tempscan > mass_bottom) & (mass_values_tempscan < mass_upper)) # For each scan, find all of the points that fall between the mass range 	 		
		EIC[i] <- sum(intensity_values_tempscan[ind]) #Binning summing the corresponding point intensity 		
	}	
	EIC
}

Binning_Mass<- function(mass_values, intensity_values, scan_index, scan_number, mass_bottom, mass_upper)
{
	Start_Scan <- min(scan_number)
	End_Scan <- max(scan_number-1)
	Binning_mass <- numeric(length = End_Scan-Start_Scan+1)
	
	for(i in Start_Scan:End_Scan)
	{
		ind_start<- scan_index[i]+1
		ind_end  <- scan_index[i+1]
		mass_values_tempscan <- mass_values[ind_start:ind_end]
		ind <- which((mass_values_tempscan > mass_bottom) & (mass_values_tempscan < mass_upper)) # For each scan, find all of the points that fall between the mass range 	 		
		Binning_mass[i] <- median(mass_values_tempscan[ind]) #mean Binning summing the corresponding point intensity 		
	}	
	Binning_mass
}

EIC_Binning_PeakDetect_Plot <-function(mass_values, intensity_values, scan_index, scan_number, scan_acquisition_time, mass_bottom, mass_upper)
{
	EIC_Mass_Bining <- EIC_Binning(mass_values, intensity_values, scan_index, scan_number, mass_bottom, mass_upper)
	Binningmass<- Binning_Mass(mass_values, intensity_values, scan_index, scan_number, mass_bottom, mass_upper)	
	
	Peaklist <- PeakDetect_CentWave(EIC_Mass_Bining,Binningmass)
	
	peak_sim<-numeric(length(scan_number))
	for (i in 1:length(Peaklist))
	{
		for (j in (Peaklist[[i]]$Left_Boundary):(Peaklist[[i]]$Right_Boundary))
		{   
			peak_scale <- ((Peaklist[[i]]$Right_Boundary) - (Peaklist[[i]]$Left_Boundary))/2.0 
			peak_sim[j]<- (Peaklist[[i]]$Peak_Intensity) *exp(-((j-(Peaklist[[i]]$Scan_Pos))*1.0/peak_scale)^2)
		}
	}
	
	bound_pos <-numeric(2*length(Peaklist))
	bound_int <-numeric(2*length(Peaklist))
	for (i in 1:length(Peaklist))
	{
		bound_pos[i*2-1]= c(Peaklist[[i]]$Left_Boundary)	
		bound_pos[i*2]  = c(Peaklist[[i]]$Right_Boundary)	
		bound_int[i*2-1]= c(Peaklist[[i]]$Peak_Intensity)
		bound_int[i*2]  = c(Peaklist[[i]]$Peak_Intensity)
	}
	
	playwith({
				lower_scan <- min(scan_number)
				upper_scan <- max(scan_number)
				lower_time <- scan_acquisition_time[lower_scan]
				upper_time <- scan_acquisition_time[upper_scan]				
				
				plot(EIC_Mass_Bining, type="l", xlab="time (min)", ylab="intensity")
				axis(1, at=lower_scan:upper_scan-1, lab= scan_acquisition_time[lower_scan:upper_scan]/60.0)
				legend(50,max(EIC_Mass_Bining),paste("EIC_Binning_Mass:", mass_bottom, "~", mass_upper, sep=""))
				lines(bound_pos, bound_int, type="h", col ="green", cex =1)
				lines(peak_sim, type="l", col="blue")
				Point_x <- numeric()
				Point_y <- numeric()
				for (i in 1: length(Peaklist))
				{
					Point_x[i] <-Peaklist[[i]]$Scan_Pos
					Point_y[i] <-Peaklist[[i]]$Peak_Intensity
				}
				points(Point_x,Point_y,pch=16,col="red",cex=1.2)	
				
			}, new=TRUE)
	
	playwith({
				plot(EIC_Mass_Bining, type="l", xlab="scan", ylab="intensity") 
				legend(50,max(EIC_Mass_Bining),paste("EIC_Binning_Mass:", mass_bottom, "~", mass_upper, sep=""))
			}, new=TRUE)
}
