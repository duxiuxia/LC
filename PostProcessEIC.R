# Made by Xiuxia Du on 5/19/2012


################################################################################
# Start -- merging EIC traces that are less than the mass measurement resolution
################################################################################

rm(list=ls())

library(playwith)


graphics.off()




# get the EICs
in_file_name <- "A_2_EIC_extraction.RData"

in_file_path <- "/home/xdu/testLECO/Output"

in_file_full_name <- paste(in_file_path, in_file_name, sep="/")

load(in_file_full_name)

pos <- regexpr('extraction.RData',in_file_name)
EIC_list_name_part1 <- substring(in_file_name, first=1, last=pos-1)
EIC_list_name_part2 <- "Extracted_mzdata"
EIC_list_name <- paste(EIC_list_name_part1, EIC_list_name_part2, sep="")

EIC_Extracted_mzdata <- get(EIC_list_name)





####################################################################################
# Plot a particular EIC
#MM <- 3
#playwith({
#			plot(EIC_Extracted_mzdata[[MM]]$scan_index, EIC_Extracted_mzdata[[MM]]$intensity)
#		}, new=TRUE)
####################################################################################





####################################################################################
# Remove short EICs

vector_EIC_length <- NULL
for (i in 1:length(EIC_Extracted_mzdata)){
	current_mass_trace <- EIC_Extracted_mzdata[[i]]$mass
	vector_EIC_length <- c(vector_EIC_length, length(current_mass_trace))
}

II <- which(vector_EIC_length < PARAMS$MinimumEIC_Count)
ratio_short_EIC <- length(II) / length(vector_EIC_length) 

playwith({
			plot(sort(vector_EIC_length))
		}, new=TRUE)

EIC_Extracted_mzdata <- EIC_Extracted_mzdata[-II]
# ##################################################################################











# ##################################################################################
# Remove large-gap EICs and save to a file
vector_scan_gap <- NULL
for (i in 1:length(EIC_Extracted_mzdata))
{
	current_scan_index <- EIC_Extracted_mzdata[[i]]$scan_index
	
	vector_gaps_in_current_scan <- diff(current_scan_index, lag=1, difference=1)
	
	II <- which(vector_gaps_in_current_scan > PARAMS$Hit_Miss_Count)
	
	vector_scan_gap <- c(vector_scan_gap, length(II))
}

II <- which(vector_scan_gap > 0)
ratio_EICs_with_big_gap <- length(II) / length(EIC_Extracted_mzdata)

EIC_Extracted_mzdata <- EIC_Extracted_mzdata[-II]

ind <- regexpr(".RData", in_file_name)
dataset_name <- substring(in_file_name, first=1, last=ind-1)
out_file_name <- paste(dataset_name, "_filtered.RData", sep="")
out_file_full_name <- paste(in_file_path, out_file_name, sep="/")

EIC_list_filtered <- paste(dataset_name, "_filtered", sep="")
assign(EIC_list_filtered, EIC_Extracted_mzdata)
save(list=paste(dataset_name, "_filtered", sep=""), file=out_file_full_name)
# ##################################################################################











# ##################################################################################
# Identify EICs that are closer in mass than the mass measurement resolution allows
vector_mass_center <- NULL
for (i in 1:length(EIC_Extracted_mzdata)){
	current_mass_trace <- EIC_Extracted_mzdata[[i]]$mass
	
	if (length(current_mass_trace) > 0){
		current_mass_center <- median(current_mass_trace)
		
		vector_mass_center <- c(vector_mass_center, current_mass_center)
	}
}

vector_mass_center <- sort(vector_mass_center)
vector_mass_diff_in_ppm <- 1000000 * abs(diff(vector_mass_center, lag=1, difference=1)) / vector_mass_center[-length(vector_mass_center)]
JJ <- which(vector_mass_diff_in_ppm <= 2*PARAMS$ROI_ppm_Th)
ratio_of_EIC_to_merge <- length(JJ) / length(vector_mass_diff_in_ppm)
playwith({
			plot(sort(vector_mass_diff_in_ppm))
		}, new=TRUE)
# ##################################################################################















# ##################################################################################
# Merge EICs that are closer in mass than the mass measurement resolution allows
EIC_merged <- vector("list", 0)

for (j in 1:length(JJ))
{
	graphics.off()
	current_index <- JJ[j]
	
	cat(vector_mass_diff_in_ppm[[current_index]])
	
	current_EIC_to_merge_1 <- EIC_Extracted_mzdata[[current_index]]
	current_EIC_to_merge_2 <- EIC_Extracted_mzdata[[(current_index+1)]]
	
	playwith({
				plot(current_EIC_to_merge_1$scan_index, current_EIC_to_merge_1$mass)
			}, new=TRUE)
	playwith({
				plot(current_EIC_to_merge_1$scan_index, current_EIC_to_merge_1$intensity)
			}, new=TRUE)	
	
	playwith({
				plot(current_EIC_to_merge_2$scan_index, current_EIC_to_merge_2$mass, col="red")
			}, new=TRUE)
	playwith({
				plot(current_EIC_to_merge_2$scan_index, current_EIC_to_merge_2$intensity, col="red")
			}, new=TRUE)	
	
	EIC_merged[[j]]	<- rbind(current_EIC_to_merge_1, current_EIC_to_merge_2)
	
	playwith({
				plot(EIC_merged[[j]]$scan_index, EIC_merged[[j]]$mass, col="blue")
			}, new=TRUE)
	
	playwith({
				plot(EIC_merged[[j]]$scan_index, EIC_merged[[j]]$intensity, col="blue")
			}, new=TRUE)		
}





















