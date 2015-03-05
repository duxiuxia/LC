# Made by Xiuxia Du on 5/19/2012


################################################################################
# Start -- merging EIC traces that are less than the mass measurement resolution
################################################################################
biocLite("xcms")

library(xcms)
library(playwith)

rm(list=ls())
graphics.off()




# get the EICs
In_File_Name <- "A_2.mzXML"

In_File_Path <- "/home/xdu/testLECO/Input"

In_File_Full_Name <- paste(In_File_Path, In_File_Name, sep="/")

xs<-xcmsRaw(filename=In_File_Full_Name, includeMSn=FALSE)
TIC <- xs@tic
scan_range <- 1:length(TIC)
playwith({
			plot(scan_range, TIC, type="l")
		}, new=TRUE)


current_scan <- getScan(xs, scan=358)
playwith({
			plot(current_scan[,1], current_scan[,2], type="h")
		}, new=TRUE)

playwith({
			plot(mass_values_tempscan, intensity_values_tempscan, type="h")
		}, new=TRUE)

