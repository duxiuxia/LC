
Looking_for_MissingPeaks <- function(fileName_mzdata)
{
	cat("Entering Looking_for_MissingPeaks \n")
	
	
#	source("http://bioconductor.org/biocLite.R")
#	#biocLite("CAMERA")
#	biocLite("xcms")
	#
	library(xcms)
	
	
	

	
	xs<-xcmsRaw(filename=fileName_mzdata, includeMSn=FALSE)    # 
	scan_acquisition_time_mzdata <- xs@scantime  # same to scantime_array<- slot(xs, "scantime")
	
	scan_index_mzdata <- xs@scanindex # same to scanindex_array<- slot(xs, "scanindex"), get the index for each scan
	scan_number_mzdata <- c(1:length(scan_index_mzdata))
	
	#mass_values_mzdata <- numeric() # initialize a null array for mass value all of scans 
	#intensity_values_mzdata <- numeric() # initialize a null array for intensity value all of scans 
	
	lower_scan <- 200
	upper_scan <- 4000
	
#	lower_mass <- 782.0
#	upper_mass <- 783.0
	
#	lower_mass <- 782.55
#	upper_mass <- 782.6
	
#	lower_mass <- 782.7
#	upper_mass <- 782.8
	
	lower_mass <- 183.5
	upper_mass <- 184.5
	

	
	Mass_spectrum_scan_temp<- getScan(xs, scan=200)  # Read mass_spectrum data of the corresponding scan  Mass_spectrum_scan_temp<- getScan(xs, i, massrange=numeric())
	
	mass_temp<- Mass_spectrum_scan_temp[,1]
	intensity_temp <- Mass_spectrum_scan_temp[,2]
	
	browser()
	
	playwith({
				cat("Entering playwith \n")
				
				plot(0, type="n", xlab="scan", ylab="mass", xlim= c(lower_scan, upper_scan), ylim=c(lower_mass, upper_mass))
				
				for (i in lower_scan:upper_scan)
				{
					Mass_spectrum_scan_temp<- getScan(xs, scan=i)  # Read mass_spectrum data of the corresponding scan  Mass_spectrum_scan_temp<- getScan(xs, i, massrange=numeric())
					
					mass_temp<- Mass_spectrum_scan_temp[,1]
					intensity_temp <- Mass_spectrum_scan_temp[,2]
					
					scan_temp <- rep(i, length(mass_temp))
					#points(scan_temp, mass_temp, pch=".")
					points(scan_temp, mass_temp, type="p", cex=0.3)
				}
				
			}, new=TRUE)
	
	browser()
	
	playwith({
				cat("Entering playwith \n")
				
				EIC_mzdata<- numeric()
				j <- 1
				
				for (i in lower_scan:upper_scan)
				{
					Mass_spectrum_scan_temp<- getScan(xs, scan=i)  # Read mass_spectrum data of the corresponding scan  Mass_spectrum_scan_temp<- getScan(xs, i, massrange=numeric())
					
					mass_temp<- Mass_spectrum_scan_temp[,1]
					intensity_temp <- Mass_spectrum_scan_temp[,2]
					
#						lower_mass_forEIC <- 782.0
#						upper_mass_forEIC <- 783.0
					
#						lower_mass_forEIC <- 782.55
#						upper_mass_forEIC <- 782.6
#				
#						lower_mass_forEIC <- 782.7
#						upper_mass_forEIC <- 782.8
					
					lower_mass_forEIC <- 184.071
					upper_mass_forEIC <- 184.075
					
					temp_index <- which(mass_temp > lower_mass_forEIC & mass_temp < upper_mass_forEIC)
					EIC_mzdata[j]<- sum(intensity_temp[temp_index])
					j <- j + 1
					
				}
				
				plot(lower_scan:upper_scan, EIC_mzdata, type="p", cex=0.3)
				
				cat("Done with playwith \n")
				
			}, new=TRUE)
	
	
	cat("Out of playwith \n")
	
	browser()
}






Examine_EIC <- function(file_name, file_path){
	cat("In Examine_EIC \n")
	browser()
	
	file_full_name_EIC <- paste(file_path, file_name, sep="/")
	
	load(file_full_name_EIC)
	
	pos <- regexpr('extraction.RData',file_name)
	EIC_list_name_part1 <- substring(file_name, first=1, last=pos-1)
	EIC_list_name_part2 <- "Extracted_mzdata"
	EIC_list_name <- paste(EIC_list_name_part1, EIC_list_name_part2, sep="")
	
	EIC_Extracted_mzdata <- get(EIC_list_name)
	
	browser()
	
	i <- 498
	playwith({
				plot(EIC_Extracted_mzdata[[i]]$scan_index, EIC_Extracted_mzdata[[i]]$intensity, type="p", cex=.3)
			}, new=TRUE)
	
	browser()
	
	scan_number_mzdata <- 1:9600
	Out_PeakFileName <- "xx.csv"
	PeakDetect_Process(EIC_Extracted_mzdata, scan_number_mzdata, Out_PeakFileName)
	
	browser()
}
	





# main

library(playwith)

graphics.off()
#rm(list=ls())






browser()




###############################################
# Start -- Looking for missing peaks
##############################################

fileName_mzdata <- "/Users/xdu4/Downloads/testLECO/A_2.mzXML"
#fileName_mzdata <- "/home/xdu/testLECO/Input/A_2.mzXML"
Looking_for_MissingPeaks(fileName_mzdata)

browser()

in_file_name <- "A_2_EIC_extraction.RData"
in_file_path <- "/home/xdu/testLECO/Output"
Examine_EIC(in_file_name, in_file_path)


browser()


###############################################
# End -- Looking for missing peaks
##############################################


