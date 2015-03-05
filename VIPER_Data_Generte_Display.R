# Generate the input file for VIPER program
Generate_VIPER_Import_Data<- function(mass_values, intensity_values, scan_index, scan_number)
{
	
	File_Viper_Im<-file("/home/wenchao/Projects/Rproject/w2_12_Temp_LCMS_Import_Threshold500.txt", "w")
	
	cat(c("scan_num	charge	abundance	mz	fit	average_mw	monoisotopic_mw	mostabundant_mw	fwhm	signal_noise	mono_abundance	mono_plus2_abundance	index"),file=File_Viper_Im)
	cat("\n",file = File_Viper_Im)
	
	unindex <- 1
	for(scan_in in 1:(length(scan_number)-1))
	{
		ind_start <- scan_index[scan_in]+1
		ind_end   <- scan_index[scan_in+1]
		
		# Get the top 200 intensity threshold
		inten_sort<-sort(intensity_values[ind_start:ind_end], decreasing = TRUE)
		if(length(inten_sort)>100)
		{
			inten_th <-500  #inten_sort[100]			
			for (pos_ind in ind_start:ind_end)
			{
				if(intensity_values[pos_ind]>=inten_th)
				{
					cat(c(scan_in, 1,intensity_values[pos_ind], mass_values[pos_ind], 0, mass_values[pos_ind], mass_values[pos_ind], mass_values[pos_ind], 0, 0, 0,0, unindex),file = File_Viper_Im, sep="\t")
					cat("\n", file = File_Viper_Im)
					unindex <- unindex+1
				}				
			}
		}
		else  # If the point number is less than 200, save all of the scan  
		{
			for (pos_ind in ind_start:ind_end)
			{
				cat(c(scan_in, 1,intensity_values[pos_ind], mass_values[pos_ind], 0, mass_values[pos_ind], mass_values[pos_ind], mass_values[pos_ind], 0, 0, 0,0, unindex),file = File_Viper_Im, sep="\t")
				cat("\n", file = File_Viper_Im)
				unindex <- unindex+1								
			}
		}
	}
	close(File_Viper_Im)	
}

# Read whole VIPER Inputting Data and outputting Feature data for whole display, 2D Display VIPER Input data and overlap the feature  
Read_VIPERData_Display <- function()
{	
	File_Import <- read.table("/home/wenchao/Projects/Rproject/w2_7_Temp_LCMS_Import_Threshold500.txt", header = TRUE, sep = "\t")
	File_Feature<- read.table("/home/wenchao/Projects/Rproject/w2_7_Temp_LCMS_Import_Threshold500_Features.txt", header = TRUE, sep = "\t")
	
	# Read the data and feature information
	data_scan_num <- File_Import$scan_num
	data_mz       <- File_Import$monoisotopic_mw   #mz  
	data_Intensity<- File_Import$abundance
	
	feature_scan_start<-File_Feature$Scan_Start
	feature_scan_end  <-File_Feature$Scan_End
	feature_mz        <-File_Feature$Monoisotopic_Mass  #Class_Rep_MZ  
	faeture_cs        <-File_Feature$Class_Rep_Charge
	
	# Get the data number and feature number 
	Data_Num   <-length(data_scan_num)
	Feature_Num<-length(feature_mz)		
	
	#Display the 2D point of the Import data file 
	playwith({
				
	            #pdf("/home/wenchao/Projects/Rproject/Mass_Scan_Map.pdf")
	            Start_Scan <- min(data_scan_num)
				End_Scan <-max(data_scan_num)
				
				ppm<-25
				
				Start_MZ<-min(data_mz)
				End_MZ  <-max(data_mz)
				plot(0, type="n", xlab="scan", ylab="mass", xlim= c(Start_Scan, End_Scan), ylim= c(Start_MZ, End_MZ))
				legend(Start_Scan, End_MZ, "VIPER_Data_Disply")
				# plot the 2D data points
				for (i in 1:Data_Num)	
				{
					points(data_scan_num[i], data_mz[i],pch=".",  col="blue")
				}
				
				#plot the feature rectangle zone				
				for (i in 1:Feature_Num)	
				{
					lines(c(feature_scan_start[i],feature_scan_end[i]),c(feature_mz[i]*(1-ppm/1000000),feature_mz[i]*(1-ppm/1000000)), type="l", col="red")
					lines(c(feature_scan_start[i],feature_scan_end[i]),c(feature_mz[i]*(1+ppm/1000000),feature_mz[i]*(1+ppm/1000000)), type="l", col="red")
					lines(c(feature_scan_start[i],feature_scan_start[i]),c(feature_mz[i]*(1-ppm/1000000),feature_mz[i]*(1+ppm/1000000)), type="l", col="red")
					lines(c(feature_scan_end[i],feature_scan_end[i]),c(feature_mz[i]*(1-ppm/1000000),feature_mz[i]*(1+ppm/1000000)), type="l", col="red")				   
				} 
				#dev.off()
				
			}, new=TRUE)
	
	
	
}

# Read whole VIPER Inputting Data and outputting Feature data for cropped display, 2D Display VIPER Input data and overlap the feature  
Read_VIPERData_Crop_Display <- function()
{		
	File_Import <- read.table("/home/wenchao/Projects/Rproject/w2_7_Temp_LCMS_Import_Threshold500.txt", header = TRUE, sep = "\t")
	File_Feature<- read.table("/home/wenchao/Projects/Rproject/w2_7_Temp_LCMS_Import_Threshold500_Features.txt", header = TRUE, sep = "\t")
	
	# Read the data and feature information
	data_scan_num <- File_Import$scan_num
	data_mz       <- File_Import$monoisotopic_mw   #mz  
	data_Intensity<- File_Import$abundance
	
	feature_scan_start<-File_Feature$Scan_Start
	feature_scan_end  <-File_Feature$Scan_End
	feature_mz        <-File_Feature$Monoisotopic_Mass  #Class_Rep_MZ  
	faeture_cs        <-File_Feature$Class_Rep_Charge
	
	# Get the data number and feature number 
	
	ppm<-25
	Scan_Crop_Left  <- 0
	Scan_Crop_Right <- 1100
	Mass_Crop_Top   <- 420.0
	Mass_Crop_Bottom<- 410.0
	
	#Data point crop filter 
	Ind_1<- which((data_scan_num>Scan_Crop_Left)&(data_scan_num<Scan_Crop_Right))
	
	data_scan_num_Crop_Scan <- data_scan_num[Ind_1]
	data_mz_Crop_Scan       <- data_mz[Ind_1]
	data_Intensity_Crop_Scan<- data_Intensity[Ind_1]
	
	Ind_2<- which((data_mz_Crop_Scan>Mass_Crop_Bottom)&(data_mz_Crop_Scan<Mass_Crop_Top))
	
	data_scan_num_Crop_Scan_Mz <- data_scan_num_Crop_Scan[Ind_2]
	data_mz_Crop_Scan_Mz       <- data_mz_Crop_Scan[Ind_2]
	data_Intensity_Crop_Scan_Mz<- data_Intensity_Crop_Scan[Ind_2]	
	
	#Feature point crop filter
	Ind_F1<- which((feature_scan_start>Scan_Crop_Left)|(feature_scan_end <Scan_Crop_Right))
	
	feature_scan_start_Crop_Scan <- feature_scan_start[Ind_F1]
	feature_scan_end_Crop_Scan   <- feature_scan_end[Ind_F1]
	feature_mz_Crop_Scan         <- feature_mz[Ind_F1]
	faeture_cs_Crop_Scan         <- faeture_cs[Ind_F1]
	
	Ind_F2<- which((feature_mz_Crop_Scan>Mass_Crop_Bottom)&(feature_mz_Crop_Scan <Mass_Crop_Top))
	
	feature_scan_start_Crop_Scan_Mz <- feature_scan_start_Crop_Scan[Ind_F2]
	feature_scan_end_Crop_Scan_Mz   <- feature_scan_end_Crop_Scan[Ind_F2]
	feature_mz_Crop_Scan_Mz         <- feature_mz_Crop_Scan[Ind_F2]
	faeture_cs_Crop_Scan_Mz         <- faeture_cs_Crop_Scan[Ind_F2]
	
	Data_Num   <-length(data_scan_num_Crop_Scan_Mz)
	Feature_Num<-length(feature_mz_Crop_Scan_Mz)	
	
	#Display the 2D point of the Import data file 
	playwith({
				Start_Scan <- min(data_scan_num_Crop_Scan_Mz)
				End_Scan <-max(data_scan_num_Crop_Scan_Mz)					
				
				Start_MZ<-min(data_mz_Crop_Scan_Mz)
				End_MZ  <-max(data_mz_Crop_Scan_Mz)
#				plot(0, type="n", xlab="scan", ylab="mass", xlim= c(Start_Scan, End_Scan), ylim= c(Start_MZ, End_MZ))
#			    legend(Start_Scan, End_MZ, "VIPER_Data_Disply_Crop")
				plot(0, type="n", xlab="scan", ylab="mass", xlim= c(Scan_Crop_Left, Scan_Crop_Right), ylim= c(Mass_Crop_Bottom, Mass_Crop_Top))
				#legend(Scan_Crop_Left, Mass_Crop_Top, "VIPER_Data_Disply_Crop")
				# plot the 2D data points
				for (i in 1:Data_Num)	
				{
					points(data_scan_num_Crop_Scan_Mz[i], data_mz_Crop_Scan_Mz[i],pch=".",  col="blue")
				}
				
				#plot the feature rectangle zone				
				for (i in 1:Feature_Num)	
				{
					lines(c(feature_scan_start_Crop_Scan_Mz[i],feature_scan_end_Crop_Scan_Mz[i]),c(feature_mz_Crop_Scan_Mz[i]*(1-ppm/1000000),feature_mz_Crop_Scan_Mz[i]*(1-ppm/1000000)), type="l", col="red")
					lines(c(feature_scan_start_Crop_Scan_Mz[i],feature_scan_end_Crop_Scan_Mz[i]),c(feature_mz_Crop_Scan_Mz[i]*(1+ppm/1000000),feature_mz_Crop_Scan_Mz[i]*(1+ppm/1000000)), type="l", col="red")
					lines(c(feature_scan_start_Crop_Scan_Mz[i],feature_scan_start_Crop_Scan_Mz[i]),c(feature_mz_Crop_Scan_Mz[i]*(1-ppm/1000000),feature_mz_Crop_Scan_Mz[i]*(1+ppm/1000000)), type="l", col="red")
					lines(c(feature_scan_end_Crop_Scan_Mz[i],feature_scan_end_Crop_Scan_Mz[i]),c(feature_mz_Crop_Scan_Mz[i]*(1-ppm/1000000),feature_mz_Crop_Scan_Mz[i]*(1+ppm/1000000)), type="l", col="red")				   
				} 
				
			}, new=TRUE)
	
}