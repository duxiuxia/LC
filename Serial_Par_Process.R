Par_Process_mzdata_File <- function(LC_MS_File, Code_Path, In_File_Path, Out_File_Path, PARAMS, sep_char)
{   
	library(wmtsa)
	library(xcms)

	#Need to get the Extract_EIC function, Wavelet based peak picking function in each cluster node
	source(paste(Code_Path, "EIC_MassTrace.R", sep=sep_char))
	source(paste(Code_Path, "Misc_Process.R", sep=sep_char)) 
	source(paste(Code_Path, "WavCWT_PeakPicking.R",sep=sep_char))
	source(paste(Code_Path, "BaseLine_Remove.R",sep=sep_char)) 
	source(paste(Code_Path, "Peak_Annotation.R", sep=sep_char))
	
	
	
	
	# Generate the full output file name including its path
	pos<-regexpr('.mzdata',LC_MS_File)
	DataSet_Name <-substring(LC_MS_File, first=1, last=pos-1)
	PeakFileName<-paste(DataSet_Name, '.csv', sep="")	
	Out_PeakFileName <- paste(Out_File_Path, sep_char, "Peak_List_", PeakFileName, sep= "")	
	fileName_mzdata<- paste(In_File_Path, LC_MS_File, sep=sep_char)  #" #blank1pos.mzdata.xml"  QC1pos.mzdata.xml" "2_7.mzdata.xml" "blank3.mzdata.xml"
	
##################################################################################################
# Code below was commented out by Xiuxia so raw data will be read from files within EIC_MassTrace.R and transfer of big data is avoided

#	xs<-xcmsRaw(filename=fileName_mzdata, includeMSn=FALSE)    # xs<-xcmsRaw(filename="/home/wenchao/Projects/Rproject/mix1_12.mzdata.xml", includeMSn=FALSE)
#	scan_acquisition_time_mzdata <- xs@scantime  # same to scantime_array<- slot(xs, "scantime")
#	
#	scan_index_mzdata <- xs@scanindex # same to scanindex_array<- slot(xs, "scanindex"), get the index for each scan
#	scan_number_mzdata <- c(1:length(scan_index_mzdata))
#	
#	mass_values_mzdata <- numeric() # initialize a null array for mass value all of scans 
#	intensity_values_mzdata <- numeric() # initialize a null array for intensity value all of scans 
#	
#	TIC_mzdata<- numeric()
#	for (i in scan_number_mzdata)
#	{
#		Mass_spectrum_scan_temp<- getScan(xs, scan=i)  # Read mass_spectrum data of the corresponding scan  Mass_spectrum_scan_temp<- getScan(xs, i, massrange=numeric())
#		
#		mass_temp<- Mass_spectrum_scan_temp[,1]
#		intensity_temp <- Mass_spectrum_scan_temp[,2]
#		
#		TIC_mzdata[i]<- sum(intensity_temp)
#		mass_values_mzdata[(length(mass_values_mzdata)+1):(length(mass_values_mzdata)+length(mass_temp))] <- mass_temp
#		intensity_values_mzdata[(length(intensity_values_mzdata)+1):(length(intensity_values_mzdata)+length(intensity_temp))] <- intensity_temp	
#	}	
##################################################################################################

#	current_variable_name <- paste(DataSet_Name, "_EIC_Extracted_mzdata", sep="")
#	assign(current_variable_name, Extract_EIC(fileName_mzdata, PARAMS))	
##	EIC_Extracted_mzdata <- Extract_EIC(mass_values_mzdata, intensity_values_mzdata, scan_index_mzdata, scan_number_mzdata, PARAMS$Scan_Start, PARAMS$Scan_End, PARAMS$ROI_ppm_Th, PARAMS$MinimumEIC_Count, PARAMS$Hit_Miss_Count, PARAMS$Cutoff_Intensity, PARAMS$Top_number, PARAMS$Quantile_Percentnumber, PARAMS$Cutoff_Mode)	
#	out_file_name <- paste(DataSet_Name, "_EIC_extraction.RData", sep="")
#	out_file_full_name <- paste(Out_File_Path, out_file_name, sep=sep_char)
#	save(list=paste(DataSet_Name, "_EIC_Extracted_mzdata", sep=""), file=out_file_full_name)

	
	
	
#	xs<-xcmsRaw(filename=fileName_mzdata, includeMSn=FALSE)    # 
#	scan_acquisition_time_mzdata <- xs@scantime  # same to scantime_array<- slot(xs, "scantime")	
#	scan_index_mzdata <- xs@scanindex # same to scanindex_array<- slot(xs, "scanindex"), get the index for each scan
#	scan_number_mzdata <- c(1:length(scan_index_mzdata))
#	rm(list="xs")
#	Save_RetentionTime_ScanNum_Map(scan_acquisition_time_mzdata, scan_number_mzdata, Out_File_Path, DataSet_Name)
	
	
	
	
	
	cat("I am entering peak detection!\n")
	in_file_name <- paste(DataSet_Name, "_EIC_extraction_filtered.RData", sep="")
	in_file_full_name <- paste(Out_File_Path, in_file_name, sep=sep_char)
	load(in_file_full_name)	
	PeakDetect_Process(get(paste(DataSet_Name, "_EIC_extraction_filtered", sep="")), scan_number_mzdata, Out_PeakFileName)
#	PeakDetect_Process(EIC_Extracted_mzdata, scan_number_mzdata, Out_PeakFileName)
	







#	Retention_Time_Scan_Map <- data.frame(scan_num=scan_number_mzdata,retention_time=scan_acquisition_time_mzdata/60.0)  # Generate a map for retention time from scan number to minutes  

	
	
	
	
	
	
	
#	Pipeline_Peak_Group_Annotate_Refinement(Out_PeakFileName, PARAMS$AdductFileName, PARAMS$Isotope_AdductFileName, PARAMS$hmdb_library_FileName, PARAMS$Max_Iso, PARAMS$Max_CS, PARAMS$mass_tol, EIC_Extracted_mzdata, Retention_Time_Scan_Map, PARAMS$Group_Scan_Shift, PARAMS$SNR_Th, PARAMS$PeakWidth_Low, PARAMS$PeakWidth_High, PARAMS$Group_Cluster_Angle_Th, PARAMS$Refinement_Cluster_Angle_Th)
}

Par_Process_cdf_File <- function(LC_MS_File, Code_Path, In_File_Path, Out_File_Path, PARAMS, sep_char)
{   
	library(wmtsa)
	library(xcms)
	library(ncdf)
	
	#Need to get the Extract_EIC function, Wavelet based peak picking function in each cluster node
	source(paste(Code_Path, "EIC_MassTrace.R", sep=sep_char))
	source(paste(Code_Path, "Misc_Process.R", sep=sep_char)) 
	source(paste(Code_Path, "WavCWT_PeakPicking.R",sep=sep_char))
	source(paste(Code_Path, "BaseLine_Remove.R",sep=sep_char)) 
	source(paste(Code_Path, "Peak_Annotation.R", sep=sep_char)) 
	# Generate the full output file name including its path
	pos<-regexpr('.cdf',LC_MS_File)
	DataSet_Name <-substring(LC_MS_File, first=1, last=pos-1)
	PeakFileName<-paste(DataSet_Name, '.csv', sep="")	
	Out_PeakFileName <- paste(Out_File_Path, sep_char,"Peak_List_", PeakFileName, sep= "")	
	
	fileName_cdf<- paste(In_File_Path, LC_MS_File, sep=sep_char)  #"/home/wenchao/Projects/Rproject/S1_1.cdf" #P_S_QC_0101.CDF
	ncid <- open.ncdf(fileName_cdf)
	print.ncdf(ncid)
	
	scan_number_cdf           <- get.var.ncdf(ncid, varid="scan_number")
	scan_acquisition_time_cdf <- get.var.ncdf(ncid, varid="scan_acquisition_time")
	total_intensity_cdf       <- get.var.ncdf(ncid, varid="total_intensity")
	scan_index_cdf            <- get.var.ncdf(ncid, varid="scan_index")
	point_count_cdf           <- get.var.ncdf(ncid, varid="point_count")
	mass_values_cdf           <- get.var.ncdf(ncid, varid="mass_values")
	intensity_values_cdf      <- get.var.ncdf(ncid, varid="intensity_values")
    	
	EIC_Extracted_cdf <- Extract_EIC(mass_values_cdf, intensity_values_cdf, scan_index_cdf, scan_number_cdf, PARAMS$Scan_Start, PARAMS$Scan_End, PARAMS$ROI_ppm_Th, PARAMS$MinimumEIC_Count, PARAMS$Hit_Miss_Count, PARAMS$Cutoff_Intensity, PARAMS$Top_number, PARAMS$Quantile_Percentnumber, PARAMS$Cutoff_Mode)	
	Save_RetentionTime_ScanNum_Map(scan_acquisition_time_cdf, scan_number_cdf, Out_File_Path, DataSet_Name)
	PeakDetect_Process(EIC_Extracted_cdf, scan_number_cdf, Out_PeakFileName)
	Retention_Time_Scan_Map <- data.frame(scan_num=scan_number_cdf,retention_time=scan_acquisition_time_cdf/60.0)  # Generate a map for retention time from scan number to minutes  
	Pipeline_Peak_Group_Annotate_Refinement(Out_PeakFileName, PARAMS$AdductFileName, PARAMS$Isotope_AdductFileName, PARAMS$hmdb_library_FileName, PARAMS$Max_Iso, PARAMS$Max_CS, PARAMS$mass_tol, EIC_Extracted_cdf, Retention_Time_Scan_Map, PARAMS$Group_Scan_Shift, PARAMS$SNR_Th, PARAMS$PeakWidth_Low, PARAMS$PeakWidth_High, PARAMS$Group_Cluster_Angle_Th, PARAMS$Refinement_Cluster_Angle_Th)
}










Serial_Process_mzdata_File <- function(LC_MS_File, In_File_Path, Out_File_Path, PARAMS, sep_char)
{   
	# Generate the full output file name including its path

	cat("Entering Serial_Process_mzdata_File \n")
#	browser()
	
	pos<-regexpr('.mzXML',LC_MS_File)
	DataSet_Name <-substring(LC_MS_File, first=1, last=pos-1)
	PeakFileName<-paste(DataSet_Name, '.csv', sep="")	
	Out_PeakFileName <- paste(Out_File_Path, sep_char,"Peak_List_", PeakFileName, sep= "")	

	fileName_mzdata<- paste(In_File_Path, LC_MS_File, sep=sep_char)  
	
##################################################################################################
	# Code below was commented out by Xiuxia so raw data will be read from files within EIC_MassTrace.R and transfer of big data is avoided
	
#	xs<-xcmsRaw(filename=fileName_mzdata, includeMSn=FALSE)    # 
#	scan_acquisition_time_mzdata <- xs@scantime  # same to scantime_array<- slot(xs, "scantime")
#	
#	scan_index_mzdata <- xs@scanindex # same to scanindex_array<- slot(xs, "scanindex"), get the index for each scan
#	scan_number_mzdata <- c(1:length(scan_index_mzdata))
#	
#	mass_values_mzdata <- numeric() # initialize a null array for mass value all of scans 
#	intensity_values_mzdata <- numeric() # initialize a null array for intensity value all of scans 
#	
#	TIC_mzdata<- numeric()
##	for (i in PARAMS$Scan_Start:PARAMS$Scan_End)
#	for (i in scan_number_mzdata)
#	{
#		Mass_spectrum_scan_temp<- getScan(xs, scan=i)  # Read mass_spectrum data of the corresponding scan  Mass_spectrum_scan_temp<- getScan(xs, i, massrange=numeric())
#		
#		mass_temp<- Mass_spectrum_scan_temp[,1]
#		intensity_temp <- Mass_spectrum_scan_temp[,2]
#		
#		TIC_mzdata[i]<- sum(intensity_temp)
#		mass_values_mzdata[(length(mass_values_mzdata)+1):(length(mass_values_mzdata)+length(mass_temp))] <- mass_temp
#		intensity_values_mzdata[(length(intensity_values_mzdata)+1):(length(intensity_values_mzdata)+length(intensity_temp))] <- intensity_temp	
#	}
##################################################################################################

	
	
#	current_variable_name <- paste(DataSet_Name, "_EIC_Extracted_mzdata", sep="")
#	assign(current_variable_name, Extract_EIC(fileName_mzdata, PARAMS))
#	
##	assign(current_variable_name, Extract_EIC(mass_values_mzdata, intensity_values_mzdata, scan_index_mzdata, scan_number_mzdata, PARAMS$Scan_Start, PARAMS$Scan_End, PARAMS$ROI_ppm_Th, PARAMS$MinimumEIC_Count, PARAMS$Hit_Miss_Count, PARAMS$Cutoff_Intensity, PARAMS$Top_number, PARAMS$Quantile_Percentnumber, PARAMS$Cutoff_Mode))
##	EIC_Extracted_mzdata <- Extract_EIC(mass_values_mzdata, intensity_values_mzdata, scan_index_mzdata, scan_number_mzdata, PARAMS$Scan_Start, PARAMS$Scan_End, PARAMS$ROI_ppm_Th, PARAMS$MinimumEIC_Count, PARAMS$Hit_Miss_Count, PARAMS$Cutoff_Intensity, PARAMS$Top_number, PARAMS$Quantile_Percentnumber, PARAMS$Cutoff_Mode)
#	
#	out_file_name <- paste(DataSet_Name, "_EIC_extraction.RData", sep="")
#	out_file_full_name <- paste(Out_File_Path, out_file_name, sep=sep_char)
#	save(list=paste(DataSet_Name, "_EIC_Extracted_mzdata", sep=""), file=out_file_full_name)
##	unlink(paste(DataSet_Name, "_EIC_extraction.RData", sep=""))
##	rm(paste(DataSet_Name, "_EIC_Extracted_mzdata", sep=""))
#	
#
#
#
#	xs<-xcmsRaw(filename=fileName_mzdata, includeMSn=FALSE)    # 
#	scan_acquisition_time_mzdata <- xs@scantime  # same to scantime_array<- slot(xs, "scantime")
#	scan_index_mzdata <- xs@scanindex # same to scanindex_array<- slot(xs, "scanindex"), get the index for each scan
#	scan_number_mzdata <- c(1:length(scan_index_mzdata))
#	rm(list="xs")
#	Save_RetentionTime_ScanNum_Map(scan_acquisition_time_mzdata, scan_number_mzdata, Out_File_Path, DataSet_Name)
	
	
	
	
	
	
	cat("I am entering peak detection!\n")
#	in_file_name <- paste(DataSet_Name, "_EIC_extraction_filtered.RData", sep="")
	in_file_name <- paste(DataSet_Name, "_EIC_extraction.RData", sep="")
	in_file_full_name <- paste(Out_File_Path, in_file_name, sep=sep_char)
	load(in_file_full_name)
	
	scan_number_mzdata <- c(1:PARAMS$Scan_End)
	
	
#	PeakDetect_Process(get(paste(DataSet_Name, "_EIC_extraction_filtered", sep="")), scan_number_mzdata, Out_PeakFileName)
	PeakDetect_Process(get(paste(DataSet_Name, "_EIC_Extracted_mzdata", sep="")), scan_number_mzdata, Out_PeakFileName)
#	PeakDetect_Process(current_variable_name, scan_number_mzdata, Out_PeakFileName)







#	Retention_Time_Scan_Map <- data.frame(scan_num=scan_number_mzdata,retention_time=scan_acquisition_time_mzdata/60.0)  # Generate a map for retention time from scan number to minutes  
#	
#	
#	
#	
#	
#	
## ToDo added by Xiuxia on 5/22/2012: pass the entire PARAMS to this function, rather than individual parameters.  	
#	cat("I am entering peak annotation!\n")
#	Pipeline_Peak_Group_Annotate_Refinement(Out_PeakFileName, PARAMS$AdductFileName, PARAMS$Isotope_AdductFileName, PARAMS$hmdb_library_FileName, PARAMS$Max_Iso, PARAMS$Max_CS, PARAMS$mass_tol, EIC_Extracted_mzdata, Retention_Time_Scan_Map, PARAMS$Group_Scan_Shift, PARAMS$SNR_Th, PARAMS$PeakWidth_Low, PARAMS$PeakWidth_High, PARAMS$Group_Cluster_Angle_Th, PARAMS$Refinement_Cluster_Angle_Th)
}









Serial_Process_cdf_File <- function(LC_MS_File, In_File_Path, Out_File_Path, PARAMS, sep_char)
{   
	# Generate the full output file name including its path
	pos<-regexpr('.cdf',LC_MS_File)
	DataSet_Name <-substring(LC_MS_File, first=1, last=pos-1)
	PeakFileName<-paste(DataSet_Name, '.csv', sep="")	
	Out_PeakFileName <- paste(Out_File_Path, sep_char,"Peak_List_", PeakFileName, sep= "")	
	
	fileName_mzdata<- paste(In_File_Path, LC_MS_File, sep=sep_char)  #
	
	pos<-regexpr('.cdf',LC_MS_File)
	DataSet_Name <-substring(LC_MS_File, first=1, last=pos-1)
	PeakFileName<-paste(DataSet_Name, '.csv', sep="")	
	Out_PeakFileName <- paste(Out_File_Path, sep_char,"Peak_List_", PeakFileName, sep= "")	
	
	fileName_cdf<- paste(In_File_Path, LC_MS_File, sep=sep_char)  #"/home/wenchao/Projects/Rproject/S1_1.cdf" #P_S_QC_0101.CDF
	ncid <- open.ncdf(fileName_cdf)
	print.ncdf(ncid)
	
	scan_number_cdf           <- get.var.ncdf(ncid, varid="scan_number")
	scan_acquisition_time_cdf <- get.var.ncdf(ncid, varid="scan_acquisition_time")
	total_intensity_cdf       <- get.var.ncdf(ncid, varid="total_intensity")
	scan_index_cdf            <- get.var.ncdf(ncid, varid="scan_index")
	point_count_cdf           <- get.var.ncdf(ncid, varid="point_count")
	mass_values_cdf           <- get.var.ncdf(ncid, varid="mass_values")
	intensity_values_cdf      <- get.var.ncdf(ncid, varid="intensity_values")
	EIC_Extracted_cdf <- Extract_EIC(mass_values_cdf, intensity_values_cdf, scan_index_cdf, scan_number_cdf, PARAMS$Scan_Start, PARAMS$Scan_End, PARAMS$ROI_ppm_Th, PARAMS$MinimumEIC_Count, PARAMS$Hit_Miss_Count, PARAMS$Cutoff_Intensity, PARAMS$Top_number, PARAMS$Quantile_Percentnumber, PARAMS$Cutoff_Mode)	
	Save_RetentionTime_ScanNum_Map(scan_acquisition_time_cdf, scan_number_cdf, Out_File_Path, DataSet_Name)
	PeakDetect_Process(EIC_Extracted_cdf, scan_number_cdf, Out_PeakFileName)
	Retention_Time_Scan_Map <- data.frame(scan_num=scan_number_cdf,retention_time=scan_acquisition_time_cdf/60.0)  # Generate a map for retention time from scan number to minutes  
	Pipeline_Peak_Group_Annotate_Refinement(Out_PeakFileName, PARAMS$AdductFileName, PARAMS$Isotope_AdductFileName, PARAMS$hmdb_library_FileName, PARAMS$Max_Iso, PARAMS$Max_CS, PARAMS$mass_tol, EIC_Extracted_cdf, Retention_Time_Scan_Map, PARAMS$Group_Scan_Shift, PARAMS$SNR_Th, PARAMS$PeakWidth_Low, PARAMS$PeakWidth_High, PARAMS$Group_Cluster_Angle_Th, PARAMS$Refinement_Cluster_Angle_Th)
}
