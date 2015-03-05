# Calculate Peak area by summing all of the intensity from left boundary to right boundary
Peak_Area_Calculation <- function(EIC_MassTrace, Left_Boundary, Right_Boundary)
{
	Peak_Area <- 0.0
	if(Left_Boundary<1) Left_Boundary <-1
	if(Right_Boundary> length(EIC_MassTrace)) Right_Boundary <- length(EIC_MassTrace) 
	for (i in Left_Boundary:Right_Boundary)
	{
		Peak_Area <- Peak_Area + EIC_MassTrace[i]
	}	
	
	Base_Area <- 0.5*(Right_Boundary-Left_Boundary+1)*(EIC_MassTrace[Right_Boundary]+EIC_MassTrace[Left_Boundary])
	Peak_Area <- Peak_Area-Base_Area 
}
# Calculate the Peak Significance Value by dividing the sum of the peak_apex and it left and right neighboring point to the mean of other points 
Peak_Significance_Calculation<- function(EIC_MassTrace, Peak_Apex_Pos,Left_Boundary, Right_Boundary)
{
	Peak_Significance <- 0
	if(((Right_Boundary- Left_Boundary)>=5) &((Left_Boundary+2)<Peak_Apex_Pos) &((Peak_Apex_Pos+2) < Right_Boundary ))
	{
		Part1<- mean(EIC_MassTrace[(Peak_Apex_Pos-1):(Peak_Apex_Pos+1)])
		Part2<- mean(c(EIC_MassTrace[Left_Boundary:(Peak_Apex_Pos-2)], EIC_MassTrace[(Peak_Apex_Pos+2):Right_Boundary]))
		
		if((EIC_MassTrace[Left_Boundary]>0)&&(EIC_MassTrace[Right_Boundary]>0)) Peak_Significance <- Part1/Part2
	}	
	Peak_Significance 	
}




PeakDetect_Process<- function(EIC_Masstrace_Array, scan_num, PeakFileName)
{
	cat("In PeakDetect_Process \n")
#	browser()
	
	fid<- file(PeakFileName, "w")	
	cat(PeakFileName)	
	cat(c("PeakID,EIC_ID,Apex_Pos,acute_mass,Peak_Intensity,Peak_SNR,Left_Boundary,Right_Boundary,Peak_Area,Peak_Significance"),file=fid)
	cat("\n",file = fid)
	
	
	
	row_index <-1
	Sel_BaseLine_Remove <- 0
	
	for (m in 1419:length(EIC_Masstrace_Array))
#	for (m in 89:89)
	{
		cat(paste("current EIC Index is", m, "\n"))
		#Use Centwave method to do peak detection  
		
		
		if (length(EIC_Masstrace_Array[[m]]$intensity)>4) # Only consider the EIC/mass trace with more than 5 valid points.
		{
			EIC_Smooth <- EIC_FillMissing_Smoothing(EIC_Masstrace_Array[[m]]$scan_index, EIC_Masstrace_Array[[m]]$intensity,scan_num)
			
#			6/1/2012: Xiuxia changed smooth times from 5 to 1 to prevent over_smoothing
			EIC_Smooth <- Chromatogram_WindowFilter_Smooth(EIC_Smooth, 1)
#			EIC_Smooth <- Chromatogram_WindowFilter_Smooth(EIC_Smooth, 5)
			EIC_BaseLine_Remove <- BaseLine_Remove(EIC_Smooth, 100)
			Acute_MassList <- EIC_FillMissing_Smoothing(EIC_Masstrace_Array[[m]]$scan_index, EIC_Masstrace_Array[[m]]$mass,scan_num)
			
			if(Sel_BaseLine_Remove == 0)
				EIC_Peak_Detect<- PeakDetect_CentWave(EIC_Smooth, Acute_MassList)
			else
				EIC_Peak_Detect<- PeakDetect_CentWave(EIC_BaseLine_Remove, Acute_MassList)   # Use the EIC after removing baseline to do centwave based peak detection
			
			if(length(EIC_Peak_Detect)>0)
			{
				# Save the corresponding peak information to csv file
				for (k in 1:length(EIC_Peak_Detect))
				{	
					cat(c(row_index,m,EIC_Peak_Detect[[k]]$Scan_Pos, EIC_Peak_Detect[[k]]$acute_mass, EIC_Peak_Detect[[k]]$Peak_Intensity,EIC_Peak_Detect[[k]]$Peak_SNR, EIC_Peak_Detect[[k]]$Left_Boundary, EIC_Peak_Detect[[k]]$Right_Boundary, EIC_Peak_Detect[[k]]$Peak_Area, EIC_Peak_Detect[[k]]$Peak_Significance), file = fid, sep=",")
					cat("\n",file = fid)
					row_index <-(row_index +1)
				}
			}
		}	
	}	
	
	
	
	close(fid)	
}






#Do Peak Detection process for all of the EIC Masstrace Array 
PeakDetect_Process_parallel<- function(EIC_Masstrace_Array, scan_num, Code_Path, sep_char)
{
	cat("In Parallel PeakDetect_Process \n")
	
	library(wmtsa)
	library(xcms)
	
	
#	Code_Path <- "/home/xdu/testLECO/LCMS"

	
	
	
	#Need to get the Extract_EIC function, Wavelet based peak picking function in each cluster node
	source(paste(Code_Path, "EIC_MassTrace.R", sep=sep_char))
	source(paste(Code_Path, "Misc_Process.R", sep=sep_char)) 
	source(paste(Code_Path, "WavCWT_PeakPicking.R",sep=sep_char))
	source(paste(Code_Path, "BaseLine_Remove.R",sep=sep_char)) 
	source(paste(Code_Path, "Peak_Annotation.R", sep=sep_char)) 
	
	
	peak_list <- data.frame(EIC_ID=as.integer(), Apex_Pos=as.integer(),acute_mass=as.numeric(), Peak_Intensity=as.integer(), Peak_SNR=as.numeric(), Left_Boundary=as.integer(), Right_Boundary=as.integer(), Peak_Area=as.numeric(), Peak_Significance=as.numeric())
	
	
	
	row_index <-1
	Sel_BaseLine_Remove <- 0

	for (m in 1:length(EIC_Masstrace_Array))
	{
		cat(paste("current EIC Index = ", m, "\n"))
		#Use Centwave method to do peak detection  
		
		
		if (length(EIC_Masstrace_Array[[m]]$intensity)>4) # Only consider the EIC/mass trace with more than 5 valid points.
		{
			EIC_Smooth <- EIC_FillMissing_Smoothing(EIC_Masstrace_Array[[m]]$scan_index, EIC_Masstrace_Array[[m]]$intensity,scan_num)
			EIC_Smooth <- Chromatogram_WindowFilter_Smooth(EIC_Smooth, 5)
			EIC_BaseLine_Remove <- BaseLine_Remove(EIC_Smooth, 100)
			Acute_MassList <- EIC_FillMissing_Smoothing(EIC_Masstrace_Array[[m]]$scan_index, EIC_Masstrace_Array[[m]]$mass,scan_num)
			
			if(Sel_BaseLine_Remove == 0)
     		   EIC_Peak_Detect<- PeakDetect_CentWave(EIC_Smooth, Acute_MassList)
			else
			   EIC_Peak_Detect<- PeakDetect_CentWave(EIC_BaseLine_Remove, Acute_MassList)   # Use the EIC after removing baseline to do centwave based peak detection
			
			if(length(EIC_Peak_Detect)>0)
			{
				# Save the corresponding peak information to csv file
				for (k in 1:length(EIC_Peak_Detect))
				{
					new_peak <- data.frame(EIC_ID=m, Apex_Pos=EIC_Peak_Detect[[k]]$Scan_Pos, acute_mass=EIC_Peak_Detect[[k]]$acute_mass, Peak_Intensity=EIC_Peak_Detect[[k]]$Peak_Intensity, Peak_SNR=EIC_Peak_Detect[[k]]$Peak_SNR, Left_Boundary=EIC_Peak_Detect[[k]]$Left_Boundary, Right_Boundary=EIC_Peak_Detect[[k]]$Right_Boundary, Peak_Area=EIC_Peak_Detect[[k]]$Peak_Area, Peak_Significance=EIC_Peak_Detect[[k]]$Peak_Significance)
					peak_list = rbind(peak_list, new_peak)
				}
			}
		}		
	}	
	peak_list

}

Detect_PlotPeaks_OneEIC<- function(EIC_Extracted,scan_number,scan_acquisition_time, EIC_Index)
{			
	
	
#	EIC_Extracted_mzdata,scan_number_mzdata,scan_acquisition_time_mzdata, EIC_Index
#	EIC_Extracted<- EIC_Extracted_mzdata
#	scan_number<- 
#	scan_acquisition_time
#	EIC_Index
	
	
	
#	browser()
	
	cat("Entering Detect_PlotPeaks_OneEIC \n")
	
	EIC_Smooth <- EIC_FillMissing_Smoothing(EIC_Extracted[[EIC_Index]]$scan_index, EIC_Extracted[[EIC_Index]]$intensity, scan_number)	
	#EIC_Smooth <- Chromatogram_WindowFilter_Smooth(EIC_Smooth, 5)
	
	Sel_BaseLine_Remove <- 0
	
	EIC_BaseLine_Remove <- BaseLine_Remove(EIC_Smooth, 100)  	
	Acute_MassList <- EIC_FillMissing_Smoothing(EIC_Extracted[[EIC_Index]]$scan_index, EIC_Extracted[[EIC_Index]]$mass,  scan_number)
	
#	browser()
	
	if (Sel_BaseLine_Remove==0) Peaklist <- PeakDetect_CentWave(EIC_Smooth, Acute_MassList) else Peaklist<- PeakDetect_CentWave(EIC_BaseLine_Remove, Acute_MassList)# Use the EIC after removing baseline to do centwave based peak detection
	
	cat(paste("Number of peaks for this EIC = ", length(PeakList), "\n"))
	
	if(length(Peaklist)>0)
	{
		peak_sim<-numeric(length(scan_number))
		for (i in 1:length(Peaklist))
		{
			for (j in (Peaklist[[i]]$Left_Boundary):(Peaklist[[i]]$Right_Boundary))
			{   
				peak_scale <- ((Peaklist[[i]]$Right_Boundary) - (Peaklist[[i]]$Left_Boundary))/2.0 
				# peak_sim[j]<- (Peaklist[[i]]$Peak_Intensity) *exp(-((j-(Peaklist[[i]]$Scan_Pos))*1.0/peak_scale)^2)    #Use exponent curve to simulate 
				peak_sim[j]<- (Peaklist[[i]]$Peak_Intensity) *(1-(((j-(Peaklist[[i]]$Scan_Pos))*1.0/peak_scale)^2))    #Use 2 order polynomonial curve to simulate
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
		
		if(Sel_BaseLine_Remove==1) EIC_Plot <- EIC_BaseLine_Remove else EIC_Plot <- EIC_Smooth
		
		playwith({
					lower_scan <- min(scan_number)
					upper_scan <- max(scan_number)
					lower_time <- scan_acquisition_time[lower_scan]
					upper_time <- scan_acquisition_time[upper_scan]				
					
				#	plot(EIC_Plot, type="l", xlab="time (min)", ylab="intensity")  #,
					plot(EIC_Plot, type="l", xlab="scans",lwd=2.5, ylab="intensity")  #,
					points(1:length(EIC_Plot), EIC_Plot, pch=1, col="black")
				#	axis(1, at=lower_scan:upper_scan-1, lab= scan_acquisition_time[lower_scan:upper_scan]/60.0)
			    #   axis(1, at=lower_scan:upper_scan-1, lab= lower_scan:upper_scan-1)
				#	legend(50,max(EIC_Extracted[[EIC_Index]]$intensity), paste("EIC_Extracted_", EIC_Index ," mz:", min(EIC_Extracted[[EIC_Index]]$mass),"~",max(EIC_Extracted[[EIC_Index]]$mass), sep=""))
					lines(bound_pos, bound_int, type="h", col ="green", lwd=2.5, cex =1)
				#	lines(peak_sim, type="l", col="blue")
					Point_x <- numeric()
					Point_y <- numeric()
					for (i in 1: length(Peaklist))
					{
						Point_x[i] <-Peaklist[[i]]$Scan_Pos
						Point_y[i] <-Peaklist[[i]]$Peak_Intensity
					}
					points(Point_x,Point_y,pch=16,col="red",cex=1.2)	
					
				}, new=TRUE)
	}
}

# Fill the missing point of a EIC Mass Trace and smooth the corresponding point with linear interpolation, or other polynorminal interpolation
EIC_FillMissing_Smoothing <- function(scan_array, intensity_array, scan_num)
{
   intensity_array_Filter <- intensity_array   
  #intensity_array_Filter <- Median_Filter_WidowLength_3(intensity_array) # Use a median filter to smooth the abrupt peak
	
	Scan_min <- min(scan_array)
	Scan_max <- max(scan_array)
	EIC_Filling_Smoothing <-numeric(length=max(scan_num)-min(scan_num)+1)
	
#	EIC_Filling_Smoothing <-c(rep(intensity_array[1],Scan_min-1),rep(0,Scan_max-Scan_min+1),rep(intensity_array[length(intensity_array)],max(scan_num)-min(scan_num)+1-Scan_max))
	
	for (i in 1:(length(scan_array)-1))
	{
		EIC_Filling_Smoothing[scan_array[i]]   <- intensity_array_Filter[i]  
		EIC_Filling_Smoothing[scan_array[i+1]] <- intensity_array_Filter[i+1]
		scan_gap <- scan_array[i+1]-scan_array[i]
		
		if(((scan_array[i+1])!=(scan_array[i]+1))&(scan_gap<5)) #Missing points in the scans, need to interpolate 
		{
			for (j in (scan_array[i]+1):(scan_array[i+1]-1))
			{
				lamda <- 1.0*(j-scan_array[i])/(scan_array[i+1]-scan_array[i])
				beta  <- 1.0*(scan_array[i+1]-j)/(scan_array[i+1]-scan_array[i])
				EIC_Filling_Smoothing[j] <-     beta*intensity_array_Filter[i] + lamda*intensity_array_Filter[i+1]
			}
		}	   
	}
	EIC_Filling_Smoothing	
}


# Using CentWave algorithm to do peak detection for one EIC Mass trace
PeakDetect_CentWave<-function(EIC_MassTrace, Acute_MassList)
{	
#	cat("In PeakDetect_CentWave \n")
	
#	browser()
	
#	6/1/2012: Xiuxia changed the upper scale from 32 to 64 to accommodate wide peaks
	EIC.cwt <-wavCWT_Modify(EIC_MassTrace, scale.range = deltat(EIC_MassTrace) * c(1, min(length(EIC_MassTrace), 64)))
#	EIC.cwt <-wavCWT_Modify(EIC_MassTrace, scale.range = deltat(EIC_MassTrace) * c(1, min(length(EIC_MassTrace), 32)))
	
#	browser()
	
	EIC_Peak_List = vector("list",0)
	if(length(which(is.na(EIC.cwt)==TRUE))==0)
	{
		EIC_Tree<-wavCWTTree_Modify(EIC.cwt,type="maxima")	    
		
		EIC_Peak<-wavCWTPeaks_Modify(EIC_Tree)		
		att<-attr(EIC_Peak, "peaks")		
		
		if (nrow(att)>0)  # Only when detect some valid peak, then save it
		{		
			#Need to convert the corresponding peak format into our know peak format to compatible with our exisiting peak detection module 
			peak_x  <- EIC_Peak$x
			peak_y  <- EIC_Peak$y
			peak_snr<- EIC_Peak$snr
			
			peak_scale<-att$scale
			
			for (i in 1:length(peak_x))
			{
				peakcent_pos  <- peak_x[i]
				Boundary_AfterExtending <- Boundary_Extending(peakcent_pos, peak_scale[i], EIC_MassTrace)
				Peak_Area_Tem    <-Peak_Area_Calculation(EIC_MassTrace, Boundary_AfterExtending[1], Boundary_AfterExtending[2])
				Peak_Sign_Tem    <-Peak_Significance_Calculation(EIC_MassTrace, peakcent_pos, Boundary_AfterExtending[1], Boundary_AfterExtending[2])
				#Peak_Area_Tem <- Peak_Area_Calculation(EIC_MassTrace, peakcent_pos-peak_scale[i], peakcent_pos+peak_scale[i])
				#EIC_Peak_Node<- data.frame(Scan_Pos=c(peakcent_pos), acute_mass=c(Acute_MassList[peakcent_pos]), Peak_Intensity=c(peak_y[i]), Peak_SNR=c(peak_snr[i]), Left_Boundary=c(max(1,peakcent_pos-peak_scale[i])), Right_Boundary=c(min(peakcent_pos+peak_scale[i], length(EIC_MassTrace))), Peak_Area=c(Peak_Area_Tem))
				
#				#########################################################################
#				Added on 6/1/2012 by Xiuxia: correct peak position and peak intensity
				correction_range <- 11 # scans to the right
				II <- which(EIC_MassTrace[peakcent_pos:(peakcent_pos+correction_range)] == max(EIC_MassTrace[peakcent_pos:(peakcent_pos+correction_range)])) 
				peakcent_pos_corrected <- peakcent_pos + II[1] - 1
				peak_intensity_corrected <- round(EIC_MassTrace[peakcent_pos_corrected])
#				#########################################################################
			
			
#				EIC_Peak_Node<- data.frame(Scan_Pos=c(peakcent_pos), acute_mass=c(Acute_MassList[peakcent_pos]), Peak_Intensity=c(peak_y[i]), Peak_SNR=c(peak_snr[i]), Left_Boundary=Boundary_AfterExtending[1], Right_Boundary=Boundary_AfterExtending[2], Peak_Area=c(Peak_Area_Tem), Peak_Significance=c(Peak_Sign_Tem))
#				EIC_Peak_List[[i]] <- EIC_Peak_Node
				EIC_Peak_List[[i]] <- data.frame(Scan_Pos=peakcent_pos_corrected, acute_mass=Acute_MassList[peakcent_pos_corrected], Peak_Intensity=peak_intensity_corrected, Peak_SNR=c(peak_snr[i]), Left_Boundary=Boundary_AfterExtending[1], Right_Boundary=Boundary_AfterExtending[2], Peak_Area=c(Peak_Area_Tem), Peak_Significance=c(Peak_Sign_Tem))
			}	
		}
	}		
	EIC_Peak_List	
}

#Boundary extending, Input the symmetry boundary but they are not equal in peak intensity. Through boundary extending and try to make the peak intensity at the left boundary and right boundary are equal. 
Boundary_Extending <-function(peakcent_pos, peak_scale, EIC_MassTrace)
{
	Boundary_AfterExtending <- c(max(1,peakcent_pos-peak_scale), min(peakcent_pos+peak_scale, length(EIC_MassTrace)))
	if((Boundary_AfterExtending[1]>1)&(EIC_MassTrace[Boundary_AfterExtending[1]]>EIC_MassTrace[Boundary_AfterExtending[2]]))  #Left_Boundary larger than Right Boundary, Extend Left boundary
	{
		i<- 1
		Boundary_Err <- EIC_MassTrace[Boundary_AfterExtending[1]]-EIC_MassTrace[Boundary_AfterExtending[2]]
		Extend_Err <-   EIC_MassTrace[Boundary_AfterExtending[1]]- EIC_MassTrace[Boundary_AfterExtending[1]-i]
		while ((i<20)&(Boundary_Err>40)&((Boundary_AfterExtending[1]-i)>1)&(Extend_Err>0)) # The maximum extending point number is 20 
		{
			Boundary_Err <- EIC_MassTrace[Boundary_AfterExtending[1]-i]-EIC_MassTrace[Boundary_AfterExtending[2]]
			Extend_Err   <- EIC_MassTrace[Boundary_AfterExtending[1]-i]-EIC_MassTrace[Boundary_AfterExtending[1]-i-1]
			i<- i+1
		}
		Boundary_AfterExtending <- c(Boundary_AfterExtending[1]-i, Boundary_AfterExtending[2])
	}else if((Boundary_AfterExtending[2]<length(EIC_MassTrace))&(EIC_MassTrace[Boundary_AfterExtending[1]]< EIC_MassTrace[Boundary_AfterExtending[2]]))  #Left_Boundary smaller than Right Boundary, Extend Right boundary 
	{
		i<- 1
		Boundary_Err <- EIC_MassTrace[Boundary_AfterExtending[2]]-EIC_MassTrace[Boundary_AfterExtending[1]]
		Extend_Err <-   EIC_MassTrace[Boundary_AfterExtending[2]]- EIC_MassTrace[Boundary_AfterExtending[2]+i]
		while ((i<20)&(Boundary_Err>40)&((Boundary_AfterExtending[2]+i)<length(EIC_MassTrace))&(Extend_Err>0)) # The maximum extending point number is 20 
		{
			Boundary_Err <- EIC_MassTrace[Boundary_AfterExtending[2]+i]-EIC_MassTrace[Boundary_AfterExtending[1]]
			Extend_Err   <- EIC_MassTrace[Boundary_AfterExtending[2]+i]-EIC_MassTrace[Boundary_AfterExtending[1]+i+1]
			i<- i+1
		}
		Boundary_AfterExtending <- c(Boundary_AfterExtending[1], Boundary_AfterExtending[2]+i)
	}
	Boundary_AfterExtending
}

Save_scan_data<- function(mass_values, intensity_values, scan_index, scannum_save, FilePath, data_name)
{
	Scan_File_Name <-file(paste(FilePath, "/",data_name,"_scandata_",scannum_save,".txt",sep=""), "w")
	
	ind_start <- scan_index[scannum_save]+1
	ind_end   <- scan_index[scannum_save+1]
	
	cat(c("mass_values","intensity"),file=Scan_File_Name, sep=",")
	cat("\n",file = Scan_File_Name)
	
	for (pos_ind in ind_start:ind_end)
	{
		cat(c(mass_values[pos_ind], intensity_values[pos_ind]), file=Scan_File_Name, sep=",")
		cat("\n", file = Scan_File_Name)
	}
	close(Scan_File_Name)
}
#When read into a data file, then generate and save the relation map of scan_acquisition_time(seconds) and scan number 
Save_RetentionTime_ScanNum_Map <- function(scan_acquisition_time, scan_number, File_Path, data_name)
{
	if(length(scan_acquisition_time)!= length(scan_number))
		stop("scan_acquisition_time and scan_number must have the same length!")

	RetentionTime_ScanNum_FileName <- file(paste(File_Path, "/", data_name, "_Retentiont_time_scannum_Map.csv",sep=""),"w")
	cat(c("scan_num","retention_time(m)"),file=RetentionTime_ScanNum_FileName, sep=",")
	cat("\n",file = RetentionTime_ScanNum_FileName)
	
	for (scan_index in min(scan_number):max(scan_number))
	{
		cat(c(scan_number[scan_index],scan_acquisition_time[scan_index]/60.0),file=RetentionTime_ScanNum_FileName, sep=",")
		cat("\n",file = RetentionTime_ScanNum_FileName)
	}
	close(RetentionTime_ScanNum_FileName)
}

#This function is used to filter the input Data by median filter and the filter length=3,  
Median_Filter_WidowLength_3 <- function(Data_Input)
{
	Data_Len <- length(Data_Input)
	Data_Output <- numeric(length=Data_Len)
	Data_Padding <- c(Data_Input[1],Data_Input,Data_Input[Data_Len])
	for (i in 1:Data_Len)
	{
		Data_Output[i] <- median(Data_Padding[i:(i+2)])
	}
	Data_Output
}

#This function is used to smooth the EIC Chromatogram Intensity by a weighting window filter [1,2,5,2,1]
Chromatogram_WindowFilter_Smooth <- function (EIC_Chromatogram, Smooth_Times)
{	
#	6/1/2012: Added by Xiuxia to prevent over-smoothing of "smooth" data
#	TODO: vec_coef should be data-dependent, so needs to add it as a parameter vector in PARAMS
#	vec_coef <- c(1, 2, 10, 2, 1)
	vec_coef <- c(0.0, 0.0, 1.0, 0.0, 0.0)
	
	if(Smooth_Times <1) Smooth_Times <-1
	if(Smooth_Times >10)Smooth_Times <-10
	
	EIC_Chromatogram_Smooth <- EIC_Chromatogram
	for (Smooth_Index in 1:Smooth_Times)
	{
		EIC_Chromatogram_Smooth_pre<- EIC_Chromatogram_Smooth # Update Chromatogram Every Smooth Filter Loop
		for (point_index in 3:(length(EIC_Chromatogram)-2))
		{
			EIC_Chromatogram_Smooth[point_index] <- (vec_coef[1]*EIC_Chromatogram_Smooth_pre[point_index-2]+ vec_coef[2]*EIC_Chromatogram_Smooth_pre[point_index-1] + vec_coef[3]*EIC_Chromatogram_Smooth_pre[point_index]+ vec_coef[4]*EIC_Chromatogram_Smooth_pre[point_index+1] + vec_coef[5]*EIC_Chromatogram_Smooth_pre[point_index+2])/sum(vec_coef)
#			EIC_Chromatogram_Smooth[point_index] <- (EIC_Chromatogram_Smooth_pre[point_index-2]+ 2.0*EIC_Chromatogram_Smooth_pre[point_index-1] + 10.0*EIC_Chromatogram_Smooth_pre[point_index]+ 2.0*EIC_Chromatogram_Smooth_pre[point_index+1] + EIC_Chromatogram_Smooth_pre[point_index+2])/16.0
		}
	}
	
	EIC_Chromatogram_Smooth	
}