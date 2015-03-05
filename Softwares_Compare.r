# ====================================
# To compare the num of detected peaks
# from ADAP and XCMS
# May 4, 2012, by Yan
# ====================================


rm(list=ls())

browser()

# get peak detection results from seperate folders
filepath_xcms <- "/Users/xdu4/Downloads/testLECO"
#filepath_xcms <- "/home/xdu/testLECO/Compare_withXCMS_5-9-2012/3scanDifference/LECO_XCMSCentWave_Peaks"
#filepath_xcms <- current_file_path
filenames_xcms <- list.files(filepath_xcms,recursive=TRUE)
filenames_xcms <- filenames_xcms[grep("peaks",filenames_xcms)]

filepath_ADAP <- "/Users/xdu4/Downloads/testLECO"
#filepath_ADAP <- "/home/xdu/testLECO/Compare_withXCMS_5-29-2012"
filenames_ADAP <- list.files(filepath_ADAP,recursive=TRUE,full.names = TRUE)
filenames_ADAP <- filenames_ADAP[grep("Peak_List",filenames_ADAP)]

# initialize a data frame to save results
peakNum_output <- data.frame(Filename=as.character(),ADAP=as.integer(),cross=as.integer(),xcms=as.integer())
compare_mode <- 1
# 2 mode decides which pipeline produce less number of peaks first 
# 1 mode Always compare ADAP to XCMS, and output those peaks that ADAP missed

if (compare_mode == 2) {
	# loop on samples
	for (i in c(1:length(filenames_xcms))) {
		browser()
		
		this.Name <- strsplit(filenames_xcms[i],split="_peaks.csv")
		ADAP_pos <- grep(this.Name,filenames_ADAP)
		
		ADAP_thisSample <- read.csv(filenames_ADAP[ADAP_pos])[,c("Apex_Pos","acute_mass","Peak_Intensity")]
		xcms_thisSample <- read.csv(paste(filepath_xcms,filenames_xcms[i],sep="/"), sep=";")[,c("mz","scpos", "maxo")]
		
		min_mass <- min(ADAP_thisSample$acute_mass)
		max_mass <- max(ADAP_thisSample$acute_mass)
		xcms_thisSample <- subset(xcms_thisSample, (scpos>=200&scpos<=4000)&(mz>=min_mass&mz<=max_mass))
#		xcms_thisSample <- subset(read.csv(paste(filepath_xcms,filenames_xcms[i],sep="/")),scpos>=200&scpos<=4000)[,c("mz","scpos","maxo")]
		ADAP_thisSample_peakNum <- nrow(ADAP_thisSample)
		xcms_thisSample_peakNum <- nrow(xcms_thisSample)
		xcms_thisSample$index <- c(1:nrow(xcms_thisSample))
		ADAP_thisSample$index <- c(1:nrow(ADAP_thisSample))
		# by checking which pipeline produce less peaks
		#cross_num <- 0
		intersect_peaks <- data.frame()
		
		if (ADAP_thisSample_peakNum <= xcms_thisSample_peakNum) {
			# ADAP produce less peaks
			for (j in c(1:ADAP_thisSample_peakNum)) {
				
				browser()
				
				curApex <- ADAP_thisSample[j,]$"Apex_Pos"
				curMass <- ADAP_thisSample[j,]$"acute_mass"
				#filtered_peak <- subset(xcms_thisSample,mz<=(curMass+0.0002)&mz>=(curMass-0.0002)&scpos<=(curApex+10)&scpos>=(curApex-10))
				#filtered_peak <- subset(xcms_thisSample,scpos<=(curApex+3)&scpos>=(curApex-3))
				filtered_num1 <- which(abs((xcms_thisSample$mz-curMass)*1000000/curMass)<=10)
				if (length(filtered_num1)>0) {
					filtered_peak1 <- xcms_thisSample[filtered_num1,]
					filtered_num2 <- which(abs(filtered_peak1$scpos-curApex)<=10)
					if (length(filtered_num2)>0) {
						filtered_peak2 <- filtered_peak1[filtered_num2,]
						# attach the closest one
						intersect_peaks <- rbind(intersect_peaks, filtered_peak2[which.min(abs(filtered_peak2$mz-curMass)),])
						#cross_num <- cross_num +1
					}
				}
			}
			if (nrow(intersect_peaks) >=0) {
				missing_peaks <- xcms_thisSample[-unique(intersect_peaks$index),]
				out_file_name <- paste(as.character(this.Name),"_ADAP_missPeaks.csv", sep="")
				out_file_full_name <- paste(filepath_ADAP, out_file_name, sep="/")
				write.csv(missing_peaks, out_file_full_name)
			}else{
				print("no matching results")
			}
		}else{
			# XCMS produce less peaks
			for (j in c(1:xcms_thisSample_peakNum)) {
				curApex <- xcms_thisSample[j,]$"scpos"
				curMass <- xcms_thisSample[j,]$"mz"
				#filtered_peak <- subset(ADAP_thisSample,acute_mass<=(curMass+0.0001)&mz>=(acute_mass-0.0001)&Apex_Pos<=(curApex+5)&Apex_Pos>=(curApex-5))	
				#filtered_peak <- subset(ADAP_thisSample,Apex_Pos<=(curApex+5)&Apex_Pos>=(curApex-5))
				filtered_num1 <- which(abs((ADAP_thisSample$acute_mass-curMass)*1000000/curMass)<=10)
				
				if (length(filtered_num1)>0) {
					filtered_peak1 <- ADAP_thisSample[filtered_num1,]
					filtered_num2 <- which(abs(filtered_peak1$Apex_Pos-curApex)<=10)
					if (length(filtered_num2)>0) {
						filtered_peak2 <- filtered_peak1[filtered_num2,]
						# attach the closest one
						intersect_peaks <- rbind(intersect_peaks, filtered_peak2[which.min(abs(filtered_peak2$acute_mass-curMass)),])
						#cross_num <- cross_num +1
					}
				}
			}
			if (nrow(intersect_peaks) >0) {
				missing_peaks <- ADAP_thisSample[-unique(intersect_peaks$index),]
				out_file_name <- paste(as.character(this.Name),"_XCMS_missPeaks.csv",sep="")
				out_file_full_name <- paste(filepath_ADAP, out_file_name, sep="/")
				write.csv(missing_peaks, out_file_full_name)
			}else{
				print("no matching results")
			}
		}
		
		peakNum_output <- rbind(peakNum_output,data.frame(Filename=as.character(this.Name),
						ADAP=ADAP_thisSample_peakNum,cross=nrow(intersect_peaks),xcms=xcms_thisSample_peakNum))
		print(i)
	}
}

if (compare_mode == 1) {
	for (i in c(1:length(filenames_xcms))) {
		this.Name <- strsplit(filenames_xcms[i],split="_peaks.csv")
		ADAP_pos <- grep(this.Name,filenames_ADAP)
		
		ADAP_thisSample <- read.csv(filenames_ADAP[ADAP_pos])[,c("Apex_Pos","acute_mass","Peak_Intensity")]
		min_mass <- min(ADAP_thisSample$acute_mass)
		max_mass <- max(ADAP_thisSample$acute_mass)
		
		
		xcms_thisSample <- read.csv(paste(filepath_xcms,filenames_xcms[i],sep="/"), sep=";")
		xcms_thisSample <- xcms_thisSample[,c("mz", "scpos","maxo")]
		xcms_thisSample <- subset(xcms_thisSample, (scpos>=200&scpos<=4000)&(mz>=min_mass&mz<=max_mass))
#		xcms_thisSample <- subset(read.csv(paste(filepath_xcms,filenames_xcms[i],sep="/")),scpos>=200&scpos<=7000)[,c("mz","scpos","maxo")]


		ADAP_thisSample_peakNum <- nrow(ADAP_thisSample)
		xcms_thisSample_peakNum <- nrow(xcms_thisSample)
		xcms_thisSample$index <- c(1:nrow(xcms_thisSample))
		ADAP_thisSample$index <- c(1:nrow(ADAP_thisSample))
		intersect_peaks <- data.frame()
		
		# Always compare ADAP to XCMS, and output those peaks that ADAP missed
		for (j in c(1:ADAP_thisSample_peakNum)) {
			curApex <- ADAP_thisSample[j,]$"Apex_Pos"
			curMass <- ADAP_thisSample[j,]$"acute_mass"
			#filtered_peak <- subset(xcms_thisSample,mz<=(curMass+0.0002)&mz>=(curMass-0.0002)&scpos<=(curApex+10)&scpos>=(curApex-10))
			#filtered_peak <- subset(xcms_thisSample,scpos<=(curApex+3)&scpos>=(curApex-3))
			filtered_num1 <- which(abs((xcms_thisSample$mz-curMass)*1000000/curMass)<=10)
			if (length(filtered_num1)>0) {
				filtered_peak1 <- xcms_thisSample[filtered_num1,]
				filtered_num2 <- which(abs(filtered_peak1$scpos-curApex)<=10)
				if (length(filtered_num2)>0) {
					filtered_peak2 <- filtered_peak1[filtered_num2,]
					# attach the closest one
					intersect_peaks <- rbind(intersect_peaks, filtered_peak2[which.min(abs(filtered_peak2$mz-curMass)),])
					#cross_num <- cross_num +1
				}
			}
		}
		
		# Check
		if (nrow(intersect_peaks) >0) {
			missing_peaks <- xcms_thisSample[-unique(intersect_peaks$index),]
			out_file_name <- paste(as.character(this.Name),"_missPeaks.csv",sep="")
			write.csv(missing_peaks,paste(filepath_ADAP,out_file_name, sep="/"))
		}else{
			print("no matching results")
		}
		
		# output
		peakNum_output <- rbind(peakNum_output,data.frame(Filename=as.character(this.Name),
						ADAP=ADAP_thisSample_peakNum,cross=nrow(intersect_peaks),xcms=xcms_thisSample_peakNum))
		print(i)
	}
}


out_file_full_name <- paste(filepath_ADAP, "ADAP_XCMS_Comp_A_2.csv", sep="")
write.csv(peakNum_output, out_file_full_name)
