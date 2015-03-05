###In this file, from the detected peak list, do annotation based on the information of peak shape correlation and retention time
Display_LocalPeak <- function()
{
	
	Color_Array<- c("red","blue", "green","yellow", "brown", "black", "orange", "cyan")
	Peak_List  <- read.csv(file="/home/wenchao/Projects/Rproject/Peaklist_w2_07_new_20111011.csv",head=TRUE,sep=",")
	#Peak_List  <- read.csv(file="/home/wenchao/Projects/Rproject/Peaklist_w2_07_new_20111021.csv",head=TRUE,sep=",")
	Apex_pos   <- 225  #210, 215, 212, 240, 225, 333, 215, 640, 262, 206, 261, 486, 879, 204, 205
	Tol        <- 3
	Peak_Neighbor_Index <- which((Peak_List$Apex_Pos <= (Apex_pos+Tol))&(Peak_List$Apex_Pos >= (Apex_pos-Tol)))
	
	Peak_List_Sub<-Peak_List[Peak_Neighbor_Index,]
	Peak_List_Filter<-PeakQuality_filter(Peak_List_Sub, EIC_Extracted_mzdata)
	
	Peak_Data<-Generate_Peak_DataArray(Peak_List_Filter, EIC_Extracted_mzdata)
	Dist_Array<-as.dist(distCal(Peak_Data))
	hc <- hclust(Dist_Array, "ward")    #,"ward"
	memb <- cutree(hc, k=2)
	Peak_List_Filter_Group1<- Peak_List_Filter[which(memb==1),]
	Peak_List_Filter_Group2<- Peak_List_Filter[which(memb==2),]
	playwith({plot(hc)},new=TRUE)
	
	Peak_List_Filter_Group1$EIC_ID
	Peak_List_Filter_Group1$PeakID
	Peak_List_Filter_Group1$acute_mass
	Peak_List_Filter_Group1$Apex_Pos
	Peak_List_Filter_Group1$Peak_Intensity
	
	Peak_List_Filter_Group2$EIC_ID	
	Peak_List_Filter_Group2$PeakID
	Peak_List_Filter_Group2$acute_mass
	Peak_List_Filter_Group2$Apex_Pos
	Peak_List_Filter_Group2$Peak_Intensity
	
	Plot_PeakGroup(Peak_List_Filter_Group1, EIC_Extracted_mzdata, Color_Array,0,"Filtered Peak List Cluster Group1")
	Plot_PeakGroup(Peak_List_Filter_Group2, EIC_Extracted_mzdata, Color_Array,0,"Filtered Peak List Cluster Group2")
	
	Plot_PeakGroup(Peak_List_Filter, EIC_Extracted_mzdata, "Peak List Display --Before Cluster")
	
#	EIC_ID_Array         <- Peak_List_Sub$EIC_ID
#   acute_mass_Array     <- Peak_List_Sub$acute_mass
#	Peak_Intensity_Array <- Peak_List_Sub$Peak_Intensity
#	Peak_SNR_Array       <- Peak_List_Sub$Peak_SNR
#	Left_Boundary_Array  <- Peak_List_Sub$Left_Boundary
#	Right_Boundary_Array <- Peak_List_Sub$Right_Boundary	
	
#   The following codes are used to display the peak list after clustering, one group use one color and different color means different group. 
	EIC_ID_Array         <- Peak_List_Filter$EIC_ID
    acute_mass_Array     <- Peak_List_Filter$acute_mass
	Peak_Intensity_Array <- Peak_List_Filter$Peak_Intensity
	Peak_SNR_Array       <- Peak_List_Filter$Peak_SNR
	Left_Boundary_Array  <- Peak_List_Filter$Left_Boundary
	Right_Boundary_Array <- Peak_List_Filter$Right_Boundary		
	
	Left_Boundary_Median <- min(Left_Boundary_Array)
	Right_Boundary_Median<- max(Right_Boundary_Array)
	Intensity_Max        <- max(Peak_Intensity_Array)
	# plot a series of Peaks with the closer peak apex position from the corresponding EIC 
	playwith({
				sub_id<- 1
				EIC_Index<- EIC_ID_Array[sub_id]
				Group_Id <- memb[sub_id]
				Color<- Color_Array[Group_Id%%8] 
				EIC_Scan <- EIC_Extracted_mzdata[[EIC_Index]]$scan_index
				
				EIC_Intensity <- EIC_Extracted_mzdata[[EIC_Index]]$intensity
				Peak_Scan_Index <- which((EIC_Scan>=Left_Boundary_Median)&(EIC_Scan<=Right_Boundary_Median))
				# In order to compare the shape patterns of different peaks that come from the same group, the intensity normalization is required. 
				Normal_Scale<- Intensity_Max/max(EIC_Intensity[Peak_Scan_Index]) 
				
				plot(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale, xlim= c(Left_Boundary_Median, Right_Boundary_Median), ylim= c(0,1.2*Intensity_Max), type="l", col=Color)
				points(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale,pch=1, col=Color)	
				for (index in 2:length(EIC_ID_Array))
				{
					EIC_Index<- EIC_ID_Array[index]
					Group_Id <- memb[index]
					Color<- Color_Array[Group_Id%%8]
					EIC_Scan <- EIC_Extracted_mzdata[[EIC_Index]]$scan_index
					
					EIC_Intensity <- EIC_Extracted_mzdata[[EIC_Index]]$intensity
					Peak_Scan_Index <- which((EIC_Scan>=Left_Boundary_Median)&(EIC_Scan<=Right_Boundary_Median))
					Normal_Scale<- Intensity_Max/max(EIC_Intensity[Peak_Scan_Index]) 
					lines(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale, col=Color)				
					points(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale,pch=1,col=Color)	
				}	
				title(main="Peak List Display after Clustering, GroupID:1,2,3,4,5.., Color:Red,blue,green,yellow,brown...",col.main="blue", font.main=4, cex.lab=0.5)
				
			},new= TRUE)
	
}

# This function can be used to plot a compound's peaks with the specific compoundID  
Plot_CompoundPeak <- function(Peak_ListFileName, scan_acquisition_time, EIC_Extracted, CompoundID)	
{
	Peak_List  <- read.csv(file=Peak_ListFileName,head=TRUE, sep=",") 
	Peak_List_CompoundID <- Peak_List[which(Peak_List$CompoundID==CompoundID),]
	Color_Array<- c("red","blue", "green","yellow", "brown", "black", "orange", "cyan")
	Figure_Name<- paste("PeakList Figure of Compound", CompoundID, sep="") 
#	Plot_PeakGroup(Peak_List_CompoundID, EIC_Extracted, Color_Array,1,Figure_Name)  #Plot Compound's peaks in a figure, x axis is labled as scan 
	Plot_PeakGroup_Minutes(Peak_List_CompoundID, EIC_Extracted, scan_acquisition_time, Color_Array, 1, Figure_Name)  #Plot Compound's peaks in a figure, x axis is labled as minutes
}

Plot_PeakGroup_Minutes<- function(Peak_List, EIC_Extracted, scan_acquisition_time, Color_Array,Normalize,Figure_Name)
{	
	
	EIC_ID_Array         <- Peak_List$EIC_ID
	acute_mass_Array     <- Peak_List$acute_mass
	Peak_Intensity_Array <- Peak_List$Peak_Intensity
	Peak_SNR_Array       <- Peak_List$Peak_SNR
	Left_Boundary_Array  <- Peak_List$Left_Boundary
	Right_Boundary_Array <- Peak_List$Right_Boundary	
	
	Left_Boundary_Median <- min(Left_Boundary_Array)
	Right_Boundary_Median<- max(Right_Boundary_Array)
	Intensity_Max        <- max(Peak_Intensity_Array)
	
	Legend_Text_Array<- NULL
	Legend_Color_Array<- NULL
	Normal_Scale<- 1
	# plot a series of Peaks with the closer peak apex position from the corresponding EIC 
	#pdf(paste("/home/wenchao/Projects/Rproject/Output/",Figure_Name,".pdf"))
	playwith({ 
	
				sub_id<- 1
				EIC_Index<- EIC_ID_Array[sub_id]
				EIC_Scan <- EIC_Extracted[[EIC_Index]]$scan_index
				Color<- Color_Array[sub_id] 
				
				EIC_Intensity <- EIC_Extracted[[EIC_Index]]$intensity
				Peak_Scan_Index <- which((EIC_Scan>=Left_Boundary_Median)&(EIC_Scan<=Right_Boundary_Median))
				# In order to compare the shape patterns of different peaks that come from the same group, the intensity normalization is required. 
				if(Normalize==1) Normal_Scale<- Intensity_Max/max(EIC_Intensity[Peak_Scan_Index]) 
				
				Legend_Text_Array[sub_id] <- paste("mass:", acute_mass_Array[sub_id], "EIC_Index:",EIC_Index)
				Legend_Color_Array[sub_id]<- Color
				# If want to use the customed axis, you must use the parameter xaxt = "n"  
				plot(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale, lwd=2.5, xlim= c(Left_Boundary_Median, Right_Boundary_Median), xaxt = "n", ylim= c(0,1.2*Intensity_Max), xlab="time (minutes)", ylab="intesnity", col=Color, type="l")
				min_scan <- min(EIC_Scan[Peak_Scan_Index])
				max_scan <- max(EIC_Scan[Peak_Scan_Index])
				axis(1, at= EIC_Scan[Peak_Scan_Index], lab=scan_acquisition_time[EIC_Scan[Peak_Scan_Index]]/60.0)
				points(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale,pch=1)
				
				if(length(EIC_ID_Array)>=2)
				{
					for (index in 2:length(EIC_ID_Array))
					{
						EIC_Index<- EIC_ID_Array[index]
						EIC_Scan <- EIC_Extracted[[EIC_Index]]$scan_index
						Color<- Color_Array[(index-1)%%length(Color_Array)+1]
						EIC_Intensity <- EIC_Extracted[[EIC_Index]]$intensity
						Peak_Scan_Index <- which((EIC_Scan>=Left_Boundary_Median)&(EIC_Scan<=Right_Boundary_Median))
						if(Normalize==1) Normal_Scale<- Intensity_Max/max(EIC_Intensity[Peak_Scan_Index])
						lines(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale, lwd=2.5, col=Color)
						points(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale,pch=1)	
						Legend_Text_Array[index] <- paste("mass:", acute_mass_Array[index], "EIC_Index:",EIC_Index)
						Legend_Color_Array[index]<- Color
					}	
				}				
				title(main=Figure_Name, col.main="blue", font.main=4, cex.lab=0.5)	
				legend(min(Left_Boundary_Array), 1.3*max(Peak_Intensity_Array),Legend_Text_Array, lty=c(1,1), lwd=c(2.5,2.5), col=Legend_Color_Array)
			   #legend(880,1.28*Intensity_Max,Legend_Text_Array, lty=c(1,1), lwd=c(2.0,2.0), col=Legend_Color_Array)
	           #dev.off()
			},new= TRUE)
	
	#pdf(paste("/home/wenchao/Projects/Rproject/Output/PseudoSpectrum_",Figure_Name,".pdf"))		
	playwith({
			   plot(acute_mass_Array, Peak_Intensity_Array,lwd=2.5, type="h", xlab="m/z", ylab="intensity", col="black")	
	#		   dev.off()
           },new=TRUE)
}

# This function can be used to plot multi compound's peaks with the specific compoundIDs  
Plot_Multi_CompoundPeak <- function(Peak_ListFileName_1, scan_acquisition_time_1, EIC_Extracted_1, CompoundID_1, Peak_ListFileName_2, scan_acquisition_time_2, EIC_Extracted_2, CompoundID_2, Peak_ListFileName_3, scan_acquisition_time_3, EIC_Extracted_3, CompoundID_3)	
{
#	Peak_ListFileName_1<-compound_peak_list_file_1
#	scan_acquisition_time_1<- scan_acquisition_time_mzdata_1
#	Peak_ListFileName_2 <- compound_peak_list_file_2
#	scan_acquisition_time_2 <- scan_acquisition_time_mzdata_2
#	Peak_ListFileName_3 <- compound_peak_list_file_3
#	scan_acquisition_time_3 <- scan_acquisition_time_mzdata_3
	
	
	Peak_List  <- read.csv(file=Peak_ListFileName_1,head=TRUE, sep=",") 
	Peak_List_CompoundID_1 <- Peak_List[which(Peak_List$CompoundID==CompoundID_1),]
	
	Peak_List  <- read.csv(file=Peak_ListFileName_2,head=TRUE, sep=",") 
	Peak_List_CompoundID_2 <- Peak_List[which(Peak_List$CompoundID==CompoundID_2),]
	
	Peak_List  <- read.csv(file=Peak_ListFileName_3,head=TRUE, sep=",") 
	Peak_List_CompoundID_3 <- Peak_List[which(Peak_List$CompoundID==CompoundID_3),]
	
	Color_Array<- c("red","blue", "green","yellow", "brown", "black", "orange", "cyan")
	Figure_Name<- paste("PeakList Figure of", "Sample_1", "_", CompoundID_1," ","Sample_2", "_", CompoundID_2, "_","Sample_3", "_", CompoundID_3,sep="") 
	
	Plot_Multi_PeakGroup_Minutes(Peak_List_CompoundID_1, EIC_Extracted_1, scan_acquisition_time_1, Peak_List_CompoundID_2, EIC_Extracted_2, scan_acquisition_time_2, Peak_List_CompoundID_3, EIC_Extracted_3, scan_acquisition_time_3, Color_Array, Figure_Name)  #Plot Compound's peaks in a figure, x axis is labled as minutes
}

# The following function is used to plot multi peak groups
Plot_Multi_PeakGroup_Minutes<- function(Peak_List_1, EIC_Extracted_1, scan_acquisition_time_1, Peak_List_2, EIC_Extracted_2, scan_acquisition_time_2, Peak_List_3, EIC_Extracted_3, scan_acquisition_time_3, Color_Array,Figure_Name)
{	
	
	EIC_ID_Array_1         <- Peak_List_1$EIC_ID
	acute_mass_Array_1     <- Peak_List_1$acute_mass
	Peak_Intensity_Array_1 <- Peak_List_1$Peak_Intensity
	Peak_SNR_Array_1       <- Peak_List_1$Peak_SNR
	Left_Boundary_Array_1  <- Peak_List_1$Left_Boundary
	Right_Boundary_Array_1 <- Peak_List_1$Right_Boundary	
	
	Left_Boundary_Median_1 <- min(Left_Boundary_Array_1)
	Right_Boundary_Median_1<- max(Right_Boundary_Array_1)
	Intensity_Max_1        <- max(Peak_Intensity_Array_1)
	
	EIC_ID_Array_2         <- Peak_List_2$EIC_ID
	acute_mass_Array_2     <- Peak_List_2$acute_mass
	Peak_Intensity_Array_2 <- Peak_List_2$Peak_Intensity
	Peak_SNR_Array_2       <- Peak_List_2$Peak_SNR
	Left_Boundary_Array_2  <- Peak_List_2$Left_Boundary
	Right_Boundary_Array_2 <- Peak_List_2$Right_Boundary	
	
	Left_Boundary_Median_2 <- min(Left_Boundary_Array_2)
	Right_Boundary_Median_2<- max(Right_Boundary_Array_2)
	Intensity_Max_2        <- max(Peak_Intensity_Array_2)
	
	EIC_ID_Array_3         <- Peak_List_3$EIC_ID
	acute_mass_Array_3     <- Peak_List_3$acute_mass
	Peak_Intensity_Array_3 <- Peak_List_3$Peak_Intensity
	Peak_SNR_Array_3       <- Peak_List_3$Peak_SNR
	Left_Boundary_Array_3  <- Peak_List_3$Left_Boundary
	Right_Boundary_Array_3 <- Peak_List_3$Right_Boundary	
	
	Left_Boundary_Median_3 <- min(Left_Boundary_Array_3)
	Right_Boundary_Median_3<- max(Right_Boundary_Array_3)
	Intensity_Max_3        <- max(Peak_Intensity_Array_3)
	
		
	Legend_Text_Array<- NULL
	Legend_Color_Array<- NULL
	Normal_Scale<- 1.0
	Normalize=0
	# plot a series of Peaks with the closer peak apex position from the corresponding EIC 
	#pdf(paste("/home/wenchao/Projects/Rproject/Output/",Figure_Name,".pdf"))
	playwith({ 
				
				sub_id<- 1
				EIC_ID_Array <- EIC_ID_Array_1
				EIC_Extracted<- EIC_Extracted_1
				Intensity_Max<- Intensity_Max_1
				Left_Boundary_Median <- Left_Boundary_Median_1
				Right_Boundary_Median<- Right_Boundary_Median_1	
				scan_acquisition_time<- scan_acquisition_time_1
				
				EIC_Index<- EIC_ID_Array[sub_id]
				EIC_Scan <- EIC_Extracted[[EIC_Index]]$scan_index												
				EIC_Intensity <- EIC_Extracted[[EIC_Index]]$intensity	
				
				Peak_Scan_Index <- which((EIC_Scan>=min(Left_Boundary_Median_1, Left_Boundary_Median_2, Left_Boundary_Median_3))&(EIC_Scan<=max(Right_Boundary_Median_1, Right_Boundary_Median_2, Right_Boundary_Median_3)))				
				Color<- Color_Array[sub_id] 
				# In order to compare the shape patterns of different peaks that come from the same group, the intensity normalization is required. 
				if(Normalize==1) Normal_Scale<- Intensity_Max/max(EIC_Intensity[Peak_Scan_Index]) 
				
			#	Legend_Text_Array[sub_id] <- paste("mass:", acute_mass_Array[sub_id], "EIC_Index:",EIC_Index)
			#	Legend_Color_Array[sub_id]<- Color
				# If want to use the customed axis, you must use the parameter xaxt = "n"  
				plot(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale, lwd=2.5, xlim= c(min(Left_Boundary_Median_1, Left_Boundary_Median_2, Left_Boundary_Median_3),max(Right_Boundary_Median_1, Right_Boundary_Median_2, Right_Boundary_Median_3)), xaxt = "n", ylim= c(0,1.2*max(Intensity_Max_1,Intensity_Max_2, Intensity_Max_3)), xlab="time (minutes)", ylab="intesnity", col=Color, type="l")
				min_scan <- min(EIC_Scan[Peak_Scan_Index])
				max_scan <- max(EIC_Scan[Peak_Scan_Index])
				axis(1, at= EIC_Scan[Peak_Scan_Index], lab=scan_acquisition_time[EIC_Scan[Peak_Scan_Index]]/60.0)
				points(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale,pch=1)
				
				if(length(EIC_ID_Array)>=2)
				{
					for (index in 2:length(EIC_ID_Array))
					{
						EIC_Index<- EIC_ID_Array[index]
						EIC_Scan <- EIC_Extracted[[EIC_Index]]$scan_index
						Color<- Color_Array[(index-1)%%length(Color_Array)+1]
						EIC_Intensity <- EIC_Extracted[[EIC_Index]]$intensity
						Peak_Scan_Index <- which((EIC_Scan>=Left_Boundary_Median)&(EIC_Scan<=Right_Boundary_Median))
						if(Normalize==1) Normal_Scale<- Intensity_Max/max(EIC_Intensity[Peak_Scan_Index])
						lines(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale, lwd=2.5, col=Color)
						points(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale,pch=1)						
					}	
				}				
				
				#Plot compound_2                
                EIC_ID_Array <- EIC_ID_Array_2
                EIC_Extracted<- EIC_Extracted_2
                Intensity_Max<- Intensity_Max_2
                Left_Boundary_Median <- Left_Boundary_Median_2
                Right_Boundary_Median<- Right_Boundary_Median_2	
                scan_acquisition_time<- scan_acquisition_time_2

                EIC_Index<- EIC_ID_Array[sub_id]
                EIC_Scan <- EIC_Extracted[[EIC_Index]]$scan_index												
                EIC_Intensity <- EIC_Extracted[[EIC_Index]]$intensity	
                Peak_Scan_Index <- which((EIC_Scan>=Left_Boundary_Median)&(EIC_Scan<=Right_Boundary_Median))				
                
				if(length(EIC_ID_Array)>=1)
				{
					for (index in 1:length(EIC_ID_Array))
					{
						EIC_Index<- EIC_ID_Array[index]
						EIC_Scan <- EIC_Extracted[[EIC_Index]]$scan_index
						Color<- Color_Array[(index-1)%%length(Color_Array)+1]
						EIC_Intensity <- EIC_Extracted[[EIC_Index]]$intensity
						Peak_Scan_Index <- which((EIC_Scan>=Left_Boundary_Median)&(EIC_Scan<=Right_Boundary_Median))
						if(Normalize==1) Normal_Scale<- Intensity_Max/max(EIC_Intensity[Peak_Scan_Index])
						lines(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale, lwd=2.5, col=Color)
						points(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale,pch=1)						
					}	
				}	
				
                #Plot compound 3
                EIC_ID_Array <- EIC_ID_Array_3
                EIC_Extracted<- EIC_Extracted_3
                Intensity_Max<- Intensity_Max_3
                Left_Boundary_Median <- Left_Boundary_Median_3
                Right_Boundary_Median<- Right_Boundary_Median_3	
                scan_acquisition_time<- scan_acquisition_time_3

                EIC_Index<- EIC_ID_Array[sub_id]
                EIC_Scan <- EIC_Extracted[[EIC_Index]]$scan_index												
                EIC_Intensity <- EIC_Extracted[[EIC_Index]]$intensity	
                Peak_Scan_Index <- which((EIC_Scan>=Left_Boundary_Median)&(EIC_Scan<=Right_Boundary_Median))				

               if(length(EIC_ID_Array)>=1)
               {
	             for (index in 1:length(EIC_ID_Array))
	             {
		           EIC_Index<- EIC_ID_Array[index]
		           EIC_Scan <- EIC_Extracted[[EIC_Index]]$scan_index
		           Color<- Color_Array[(index-1)%%length(Color_Array)+1]
		           EIC_Intensity <- EIC_Extracted[[EIC_Index]]$intensity
		           Peak_Scan_Index <- which((EIC_Scan>=Left_Boundary_Median)&(EIC_Scan<=Right_Boundary_Median))
		           if(Normalize==1) Normal_Scale<- Intensity_Max/max(EIC_Intensity[Peak_Scan_Index])
		           lines(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale, lwd=2.5, col=Color)
		           points(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale,pch=1)						
	             }	
               }				
				
				title(main=Figure_Name, col.main="blue", font.main=4, cex.lab=0.5)	
				#legend(min(Left_Boundary_Array), 1.3*max(Peak_Intensity_Array),Legend_Text_Array, lty=c(1,1), lwd=c(2.5,2.5), col=Legend_Color_Array)
				#legend(880,1.28*Intensity_Max,Legend_Text_Array, lty=c(1,1), lwd=c(2.0,2.0), col=Legend_Color_Array)
				#dev.off()
			},new= TRUE)
	
	
}

# This function can plot a peak list group that defined in Peak_List, the peak color use the corresponding color defined in Color_Array
# Normalize, means whether normalize the peaks defined in Peak_List. Figure_Name Corresponds to the titile for the figure
Plot_PeakGroup<- function(Peak_List, EIC_Extracted,Color_Array,Normalize,Figure_Name)
{	
	
	EIC_ID_Array         <- Peak_List$EIC_ID
	acute_mass_Array     <- Peak_List$acute_mass
	Peak_Intensity_Array <- Peak_List$Peak_Intensity
	Peak_SNR_Array       <- Peak_List$Peak_SNR
	Left_Boundary_Array  <- Peak_List$Left_Boundary
	Right_Boundary_Array <- Peak_List$Right_Boundary	
	
	Left_Boundary_Median <- min(Left_Boundary_Array)
	Right_Boundary_Median<- max(Right_Boundary_Array)
	Intensity_Max        <- max(Peak_Intensity_Array)
	
	Legend_Text_Array<- NULL
	Legend_Color_Array<- NULL
	Normal_Scale<- 1
	# plot a series of Peaks with the closer peak apex position from the corresponding EIC 
	playwith({
				sub_id<- 1
				EIC_Index<- EIC_ID_Array[sub_id]
				EIC_Scan <- EIC_Extracted[[EIC_Index]]$scan_index
				Color<- Color_Array[sub_id] 
				
				EIC_Intensity <- EIC_Extracted[[EIC_Index]]$intensity
				Peak_Scan_Index <- which((EIC_Scan>=Left_Boundary_Median)&(EIC_Scan<=Right_Boundary_Median))
				# In order to compare the shape patterns of different peaks that come from the same group, the intensity normalization is required. 
				if(Normalize==1) Normal_Scale<- Intensity_Max/max(EIC_Intensity[Peak_Scan_Index]) 
				
				Legend_Text_Array[sub_id] <- paste("mass:", acute_mass_Array[sub_id], "EIC_Index:",EIC_Index)
				Legend_Color_Array[sub_id]<- Color
				
				plot(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale, xlim= c(Left_Boundary_Median, Right_Boundary_Median), ylim= c(0,1.2*Intensity_Max), col=Color, type="l")
				points(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale,pch=1)
				 
				if(length(EIC_ID_Array)>=2)
				{
					for (index in 2:length(EIC_ID_Array))
					{
						EIC_Index<- EIC_ID_Array[index]
						EIC_Scan <- EIC_Extracted[[EIC_Index]]$scan_index
						Color<- Color_Array[(index-1)%%length(Color_Array)+1]
						EIC_Intensity <- EIC_Extracted[[EIC_Index]]$intensity
						Peak_Scan_Index <- which((EIC_Scan>=Left_Boundary_Median)&(EIC_Scan<=Right_Boundary_Median))
						if(Normalize==1) Normal_Scale<- Intensity_Max/max(EIC_Intensity[Peak_Scan_Index])
						lines(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale, col=Color)
						points(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale,pch=1)	
						Legend_Text_Array[index] <- paste("mass:", acute_mass_Array[index], "EIC_Index:",EIC_Index)
						Legend_Color_Array[index]<- Color
					}	
				}				
				title(main=Figure_Name, col.main="blue", font.main=4, cex.lab=0.5)	
				legend(min(Left_Boundary_Array),max(Peak_Intensity_Array),Legend_Text_Array, lty=c(1,1), lwd=c(2.5,2.5), col=Legend_Color_Array)
			},new= TRUE)
}

#Display Isotopic Peak pair
Display_IsotopePeaks <- function()
{
	Color_Array<- c("red","blue", "green","yellow", "brown", "black", "orange", "cyan")
	Peak_List  <- read.csv(file="/home/wenchao/Projects/Rproject/Peaklist_w2_07_new_20111011.csv",head=TRUE,sep=",")
	PeakGroup_Index <- c(788, 813) #76.0399(),  (788, 813), (1288, 1322), 116.0706 (1698,1738), (1761, 1785), (1762, 1787), (2572,2628), 136.0432(), (4060,4063, 4065), (4387, 4391, 4394), (4516, 4522), (4896, 4907, 4925), (5033, 5062), (5295, 5314)
	
	EIC_ID_Array         <- Peak_List[PeakGroup_Index,]$EIC_ID
	acute_mass_Array     <- Peak_List[PeakGroup_Index,]$acute_mass
	Apex_pos_Array       <- Peak_List[PeakGroup_Index,]$Apex_Pos
	Peak_Intensity_Array <- Peak_List[PeakGroup_Index,]$Peak_Intensity
	Peak_SNR_Array       <- Peak_List[PeakGroup_Index,]$Peak_SNR
	Left_Boundary_Array  <- Peak_List[PeakGroup_Index,]$Left_Boundary
	Right_Boundary_Array <- Peak_List[PeakGroup_Index,]$Right_Boundary
	
	Left_Boundary_Median <- min(Left_Boundary_Array)-5
	Right_Boundary_Median<- max(Right_Boundary_Array)+5
	Intensity_Max        <- max(Peak_Intensity_Array)
	# plot a series of Peaks with the closer peak apex position from the corresponding EIC 
	playwith({
				EIC_Index<- EIC_ID_Array[1]
				Color<- Color_Array[1] 
				EIC_Scan <- EIC_Extracted_mzdata[[EIC_Index]]$scan_index
				
				EIC_Intensity <- EIC_Extracted_mzdata[[EIC_Index]]$intensity
				Peak_Scan_Index <- which((EIC_Scan>=Left_Boundary_Median)&(EIC_Scan<=Right_Boundary_Median))
				# In order to compare the shape patterns of different peaks that come from the same group, the intensity normalization is required. 
				Normal_Scale<- Intensity_Max/max(EIC_Intensity[Peak_Scan_Index]) 
				
				plot(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale, xlim= c(Left_Boundary_Median, Right_Boundary_Median), ylim= c(0,1.2*Intensity_Max), type="l", col=Color)
				points(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale,pch=1, col=Color)	
				for (index in 2:length(EIC_ID_Array))
				{
					EIC_Index<- EIC_ID_Array[index]
					Color<- Color_Array[index%%8]
					EIC_Scan <- EIC_Extracted_mzdata[[EIC_Index]]$scan_index
					
					EIC_Intensity <- EIC_Extracted_mzdata[[EIC_Index]]$intensity
					Peak_Scan_Index <- which((EIC_Scan>=Left_Boundary_Median)&(EIC_Scan<=Right_Boundary_Median))
					Normal_Scale<- Intensity_Max/max(EIC_Intensity[Peak_Scan_Index]) 
					lines(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale, col=Color)
					#	lines(EIC_Extracted_mzdata[[EIC_Index]]$scan_index, EIC_Extracted_mzdata[[EIC_Index]]$intensity, col=Color) 
					points(EIC_Scan[Peak_Scan_Index], EIC_Intensity[Peak_Scan_Index]*Normal_Scale,pch=1,col=Color)	
				}
				
				#legend(50,max(EIC_Extracted_mzdata[[EIC_Index]]$intensity), paste("EIC_Extracted_", EIC_Index ," mz:", min(EIC_Extracted_mzdata[[EIC_Index]]$mass),"~",max(EIC_Extracted_mzdata[[EIC_Index]]$mass), sep=""))
			},new= TRUE)
}

# Pipeline for Peak Group, Peak Annotate, GroupRefinement
Pipeline_Peak_Group_Annotate_Refinement <- function(PeakFileName, AdductFileName, Isotope_AdductFileName, hmdb_library_FileName, Max_Iso, Max_CS, mass_tol, EIC_Extracted, Retention_Time_Scan_Map, Group_Scan_Shift, SNR_Th, PeakWidth_Low, PeakWidth_High, Group_Cluster_Angle_Th, Refinement_Cluster_Angle_Th)
{   
	PeakFileName_Group<-Peak_Group(PeakFileName, EIC_Extracted, Group_Scan_Shift, SNR_Th, PeakWidth_Low, PeakWidth_High, Group_Cluster_Angle_Th)
	cat("Peak Group is done!\r")
	PeakFileName_Isotope_Annotate<- Isotope_Annotation(PeakFileName_Group, Max_Iso, Max_CS, mass_tol)
	cat("Isotope Annotate is done!\r")
	PeakFileName_Adduct_Annotate <- Adduct_Annotation(PeakFileName_Isotope_Annotate, AdductFileName, mass_tol)
	cat("Adduct Annotate is done!\r")
	PeakFileName_Isotope_Adduct_Annotate_Refinement<- Group_Refinement(PeakFileName_Adduct_Annotate, Isotope_AdductFileName, hmdb_library_FileName,mass_tol,EIC_Extracted, Refinement_Cluster_Angle_Th)
	cat("Group Refinement is done!\r")
	# The following codes are used to map the compound's retention time from scan number to minutes 
	Compound_List_In  <- read.csv(file=PeakFileName_Isotope_Adduct_Annotate_Refinement,head=TRUE,sep=",")
	Compound_List_Out <- Acquire_Compound_RetentionIime_Map(Compound_List_In, Retention_Time_Scan_Map)
	pos<-regexpr('.csv',PeakFileName_Isotope_Adduct_Annotate_Refinement)
	PeakFileName_Isotope_Adduct_Annotate_Refinement_RtMap<-paste(substring(PeakFileName_Isotope_Adduct_Annotate_Refinement, first=1, last=pos-1), '_RtMap', '.csv', sep="")
	write.csv(Compound_List_Out, file = PeakFileName_Isotope_Adduct_Annotate_Refinement_RtMap)
	cat("Compound Retention time Mapping is done!\r")
}

#Do peak grouping according to the inputting peaklist, one output group means a cluster of molecular compound peak and its corresponding isotoping peaks, adduct peaks,    
Peak_Group <- function (PeakFileName, EIC_Extracted, Group_Scan_Shift, SNR_Th, PeakWidth_Low, PeakWidth_High, Group_Cluster_Angle_Th)
{
	#PeakFileName <- "/home/wenchao/Projects/Rproject/Peaklist_w3_7_20111213.csv"   # Peaklist_w2_7_20111205.csv Peaklist_w2_07_new_20111021.csv  Peaklist_w2_07_new_20111021.csv
	#EIC_Extracted<- EIC_Extracted_mzdata
	
	Tol          <- Group_Scan_Shift
	Peak_List_Orignal  <- read.csv(file=PeakFileName,head=TRUE, sep=",") 
	Peak_List_Filter   <- PeakQuality_filter(Peak_List_Orignal, EIC_Extracted, SNR_Th, PeakWidth_Low, PeakWidth_High)	
		
	#write.csv(Peak_List_Filter, file = "/home/wenchao/Projects/Rproject/Peaklist_w3_7_20111213_filter.csv")  #Peaklist_w2_07_new_20111021_filter.csv
	pos<-regexpr('.csv',PeakFileName)
	PeakFileName_filter<-paste(substring(PeakFileName, first=1, last=pos-1), '_filter', '.csv', sep="")
	write.csv(Peak_List_Filter, file = PeakFileName_filter)
	
	Peak_List_Filter$GroupID <- rep(0, nrow(Peak_List_Filter))		
	GroupID <-1
	Peak_List_Filter_ForGroup<- Peak_List_Filter[which(Peak_List_Filter$GroupID==0),]	
	while(nrow(Peak_List_Filter_ForGroup)>0)
	{
		Peak_Intensity_Array     <- Peak_List_Filter_ForGroup$Peak_Intensity
		Peak_Max_Index           <- which.max(Peak_Intensity_Array)
		Max_Apex_Pos             <- Peak_List_Filter_ForGroup[Peak_Max_Index,]$Apex_Pos
		Max_PeakID               <- Peak_List_Filter_ForGroup[Peak_Max_Index,]$PeakID
		
		MaxPeak_Neighbor_Index <- which((Peak_List_Filter$GroupID==0)&(Peak_List_Filter$Apex_Pos <= (Max_Apex_Pos+Tol))&(Peak_List_Filter$Apex_Pos >= (Max_Apex_Pos-Tol)))
		MaxPeak_List           <- Peak_List_Filter[MaxPeak_Neighbor_Index, ]
		if(nrow(MaxPeak_List)>1){
			MaxPeak_Data           <- Generate_Peak_DataArray(MaxPeak_List, EIC_Extracted)
			Dist_Array<-as.dist(distCal(MaxPeak_Data))
			hc <- hclust(Dist_Array, "ward")
			memb <- cutree(hc, h=Group_Cluster_Angle_Th)
			Compound_GroupID <- memb[which(names(memb)==Max_PeakID)]
			Compound_Group   <- MaxPeak_List[which(memb==Compound_GroupID),]
		}else{
			Compound_Group   <- MaxPeak_List
		}
				
		#Record the Peak that has been grouped by cluster method as it own groupID 
		for(i in 1:nrow(Compound_Group))
		{
			Peak_List_Filter[which((Peak_List_Filter$PeakID)==(Compound_Group[i,]$PeakID)),]$GroupID <- GroupID
		}
		GroupID <- GroupID + 1
		
		#Update Peak_List_ForGroup
	    Peak_List_Filter_ForGroup<- Peak_List_Filter[which(Peak_List_Filter$GroupID==0),]
	}	
	
	#write.csv(Peak_List_Filter, file = "/home/wenchao/Projects/Rproject/Peaklist_w3_7_20111213_filter_Group.csv")    #Peaklist_w2_07_filter_Group.csv
    pos<-regexpr('.csv',PeakFileName_filter)
    PeakFileName_Group <-paste(substring(PeakFileName_filter, first=1, last=pos-1), '_Group', '.csv', sep="")
    write.csv(Peak_List_Filter, file = PeakFileName_Group)  #Peaklist_w2_07_filter_Group.csv
	PeakFileName_Group
}

#For the Peak List after groupping, they need to do Isotope annotation. So, the output Peak list will be added 1 column that named as IsotopeID.  
#For each group, do isotope annotation. Do loop from the smallest mass to the largest mass, and pairwisely judge the mass difference whether the peak pair meet the requirements of isotope pair or not.   
Isotope_Annotation <- function(PeakFileName, Max_Iso, Max_CS, mass_tol)
{
#   PeakFileName<- "/home/wenchao/Projects/Rproject/Peaklist_w3_7_20111213_filter_Group.csv"     #Peaklist_w2_07_filter_Group.csv
#   mass_tol    <- 0.011
#   
#   Max_Iso     <- 2   #Define the maximum Isotope Number that should be considered
#   Max_CS      <- 1
   
   Peak_List <- read.csv(file=PeakFileName,head=TRUE, sep=",")
     
   Peak_List$Isotope_Annotation <- rep(0, nrow(Peak_List))	
   GroupID_Array       <- Peak_List$GroupID
   Isotope_ID<-1
   for (group_Index in min(GroupID_Array):max(GroupID_Array))
   {
	   Peak_List_OneGroup  <- Peak_List[which(GroupID_Array==group_Index),]
	   mz_list             <- Peak_List_OneGroup$acute_mass 
	   PeakID_OneGroup     <- Peak_List_OneGroup$PeakID
	   
	   MonoIsotope_Record  <- rep(0, nrow(Peak_List_OneGroup))
	   Mass_Matrix_Estimate<- matrix(0, ncol=length(mz_list), nrow=(Max_Iso+1)*Max_CS)
	   # Generate the Estimating isotope Mass Matrix 
	   for (MZ_Index in 1:length(mz_list))
	   {
	      for(Iso_Index in 0:Max_Iso)
	      {
		     for (CS_Index in 1:Max_CS)  
		     {
			     Mass_Matrix_Estimate[Iso_Index*Max_CS+CS_Index, MZ_Index]<- mz_list[MZ_Index]*CS_Index- (CS_Index+Iso_Index)*1.007276
		     }
	      }
	   }
	   
	  #Dump the 2D Matrix into a vector and use hclust 
	  Mass_List <- NULL
	  for (Row_Index in 1:nrow(Mass_Matrix_Estimate))
	  {
		 Mass_List <- c(Mass_List, Mass_Matrix_Estimate[Row_Index, ]) 
	  }
	  
	  # Use hclust to cluster the Mass_List 
	  hc<- hclust(dist(Mass_List))   
	  memb<- cutree(hc, h=mass_tol)	   
	  # Try to find the optimal cluster group that have more than 1 nmuber    
	  count <- rep(0,max(memb))
	  for (i in 1:length(memb))
	  {		 
		 count[memb[i]] <- count[memb[i]]+1
	  }
	  for(count_index in 1:length(count))
	  {
		  if(count[count_index]>1) # Found some meaningful Isotope information, Do Isotope Annotation 
		  {			  
			  SameMass_Pos <- which(memb==count_index)
			  # From the index to deduce the detailed Row and Column 
			  SameMass_Row_Pos<-floor((SameMass_Pos-1)/ncol(Mass_Matrix_Estimate))+1
			  SameMass_Col_Pos<- (SameMass_Pos-1)%%ncol(Mass_Matrix_Estimate)+1
			  
			  # From the Col Pos to deduce the detailed CS information and Isotope Information
			  SameMass_Iso <-floor((SameMass_Row_Pos-1)/Max_CS)
			  SameMass_CS  <-(SameMass_Row_Pos-1)%%Max_CS+1
			  
			  #Try to avoid the case annotating some peak with monisotope peak, but this peak has been annotated before. I fdon't avoid this case, will overwrite the previous annotation reslut  
		      Mono_pos<-SameMass_Col_Pos[which.min(SameMass_Iso)] 	  
		
		      if((min(SameMass_Iso)==0)&&(length(which(SameMass_Iso==0))==1)&&(MonoIsotope_Record[Mono_pos]==0)) #Must ensure the Isotope Mass peak [M] Exist
			  {   
				  for(Annote_Index in 1:length(SameMass_Pos))
				  {
					  pos<- SameMass_Col_Pos[Annote_Index]
					  Annote_Pos<- which((Peak_List$PeakID)==PeakID_OneGroup[pos])				  
					  Peak_List[Annote_Pos,]$Isotope_Annotation <- paste("[",Isotope_ID,"]","[M+",SameMass_Iso[Annote_Index],"]",SameMass_CS[Annote_Index],"+")
					  # Record this Isotope Annotation Record, its corresponding mass peak can't be used as monoisotope peak. 
					  MonoIsotope_Record[pos]<-1
				  }				  
				  Isotope_ID<- Isotope_ID+1	 
			  }			  	  
		  }	 
	  }	  	   
   } 
   
   #write.csv(Peak_List, file = "/home/wenchao/Projects/Rproject/Peaklist_w3_7_20111213_filter_Group_Isotope_Annotate.csv") #Peaklist_w2_07_filter_Group_Isotope_Annotate.csv
   
   pos<-regexpr('.csv',PeakFileName)
   PeakFileName_Isotope_Annotate <-paste(substring(PeakFileName, first=1, last=pos-1), '_Annotate_Isotope', '.csv', sep="")
   write.csv(Peak_List, file = PeakFileName_Isotope_Annotate) 
   PeakFileName_Isotope_Annotate
}

Adduct_Annotation <- function(PeakFileName, AdductFileName, mass_tol)
{
#	PeakFileName      <- "/home/wenchao/Projects/Rproject/Peaklist_w3_7_20111213_filter_Group_Isotope_Annotate.csv"        #Peaklist_w2_07_filter_Group_Isotope_Annotate.csv
#	AdductFileName    <- "/home/wenchao/Projects/Rproject/Adduct_Infor_Config.csv"
#	mass_tol          <- 0.01
    
	Peak_List         <- read.csv(file=PeakFileName,head=TRUE, sep=",")
	Adduct_Infor_List <- read.csv(file=AdductFileName,head=TRUE, sep=",")
   # Adduct_Infor_List <- data.frame(Formula=c("[M+Na]+","[M+K]+","[M-C3H9N]+","[M+2Na-H]","[2M+Na]+"), Para_N=c(1,1,1,1,2), Para_CS=c(1,1,1,1,1), Para_Mass_Shift=c(22.98977, 38.963708, -59.073499,44.96563,22.98977))
   # write.csv(Adduct_Infor_List,file ="/home/wenchao/Projects/Rproject/Adduct_Infor_Config.csv")
	Peak_List$Adduct_Annotation <- rep("", nrow(Peak_List))	
	GroupID_Array               <- Peak_List$GroupID
	for (group_Index in min(GroupID_Array):max(GroupID_Array))
	{
		Peak_List_OneGroup  <- Peak_List[which(GroupID_Array==group_Index),]
		mz_list             <- Peak_List_OneGroup$acute_mass 
		PeakID_OneGroup     <- Peak_List_OneGroup$PeakID
		
		Mass_Matrix_Estimate<- matrix(0, ncol=length(mz_list), nrow=nrow(Adduct_Infor_List))
		#Generate the 2D Estimating Mass Matrix according to the current peak group and Adduct Parameter Information List 
		for (MZ_Index in 1:length(mz_list))
		{
			for(Adduct_Index in 1:nrow(Adduct_Infor_List))
			{				
				Para_N<- Adduct_Infor_List[Adduct_Index,]$Para_N
				Para_Z<- Adduct_Infor_List[Adduct_Index,]$Para_CS
				Para_a<- Adduct_Infor_List[Adduct_Index,]$Para_Mass_Shift
				Mass_Matrix_Estimate[Adduct_Index, MZ_Index]<- (mz_list[MZ_Index]*Para_Z- Para_a)/Para_N
			}
		}
		
		#Dump the 2D Matrix into a vector and use hclust 
		Mass_List <- NULL
		for (Row_Index in 1:nrow(Mass_Matrix_Estimate))
		{
			Mass_List <- c(Mass_List, Mass_Matrix_Estimate[Row_Index, ]) 
		}
		
		# Use hclust to cluster the Mass_List 
		hc<- hclust(dist(Mass_List))   
		memb<- cutree(hc, h=mass_tol)	   
		# Try to find the optimal cluster group that have more than 1 nmuber    
		count <- rep(0,max(memb))
		for (i in 1:length(memb))
		{		 
			count[memb[i]] <- count[memb[i]]+1
		}
		
		for(count_index in 1:length(count))
		{
			if(count[count_index]>1) # Found some meaningful Isotope information, Do Adduct Annotation 
			{			  
				SameMass_Pos <- which(memb==count_index)
				# From the index to deduce the detailed Row and Column 
				SameMass_Row_Pos   <-floor((SameMass_Pos-1)/ncol(Mass_Matrix_Estimate))+1
				SameMass_Col_Pos   <- (SameMass_Pos-1)%%ncol(Mass_Matrix_Estimate)+1				
				
				# Do Adduct Annotation for each peak according to the corresponding adduct information
		        for(Annote_Index in 1:length(SameMass_Pos))
		        {
			        pos_col<- SameMass_Col_Pos[Annote_Index]
					pos_row<- SameMass_Row_Pos[Annote_Index]
					Molecular_Mass <- Mass_Matrix_Estimate[pos_row,pos_col]
			        Annote_Pos<- which((Peak_List$PeakID)==PeakID_OneGroup[pos_col])	
					#Peak_List[Annote_Pos,]$Adduct_Annotation <- paste(Adduct_Infor_List[pos_row,]$Formula, " M:", Molecular_Mass)	
			        old_Annotation<-Peak_List[Annote_Pos,]$Adduct_Annotation
			        Peak_List[Annote_Pos,]$Adduct_Annotation <- paste(old_Annotation," ",Adduct_Infor_List[pos_row,]$Formula," M:", Molecular_Mass) # Add the new adduct annotation into the end of old annotation.  
		        }		
			 }
		 }		
	}
	
	#write.csv(Peak_List, file = "/home/wenchao/Projects/Rproject/Peaklist_w3_7_20111213_filter_Group_Isotope_Adduct_Annotate.csv")       #Peaklist_w2_07_filter_Group_Isotope_Adduct_Annotate.csv
    pos<-regexpr('.csv',PeakFileName)
    PeakFileName_Adduct_Annotate <-paste(substring(PeakFileName, first=1, last=pos-1), '_Annotate_Adduct', '.csv', sep="")
    write.csv(Peak_List, file = PeakFileName_Adduct_Annotate)
	PeakFileName_Adduct_Annotate
}

# A group can be acquired according to the closer retention time and closer peak shape(use dot product), but the peaks of a group maybe contain more than 1 pure compound. So, the group need another group refinement process
# Before do group refinement, we need do isotoping annotation and adduct annotation, if some peaks belong to a same compound group, they should have some relationship as isotoping or adduct. So, we can use the annotation information of 
# Isotoping Annotation and Adduction annotation to do group refinement. 
# After get the annotation information of Isotoping and Adduction for a group,  
Group_Refinement<- function(PeakFileName,Isotope_AdductFileName, hmdb_library_FileName,mass_tol,EIC_Extracted, Refinement_Cluster_Angle_Th)
{
#	EIC_Extracted             <- EIC_Extracted_mzdata
#	Color_Array               <- c("red","blue", "green","yellow", "brown", "black", "orange", "cyan", "purple", "orange", "magenta", "pink")
#	PeakFileName              <- "/home/wenchao/Projects/Rproject/Peaklist_w3_7_20111213_filter_Group_Isotope_Adduct_Annotate.csv"        #Peaklist_w2_07_filter_Group_Isotope_Adduct_Annotate.csv
#	Isotope_AdductFileName    <- "/home/wenchao/Projects/Rproject/Isotope_Adduct_Infor_Config.csv"
#   hmdb_library_FileName     <- "/home/wenchao/Projects/Rproject/hmdb.csv"
#	mass_tol                  <- 0.01
	
	Peak_List                 <- read.csv(file=PeakFileName,head=TRUE,sep=",")
	Isotope_Adduct_Infor_List <- read.csv(file=Isotope_AdductFileName,head=TRUE, sep=",")
	hmdb_library              <- read.csv(file=hmdb_library_FileName, head=TRUE, sep=",")

#	Peak_OneGroup <- Peak_List[which(Peak_List$GroupID==7), ]
#	Plot_PeakGroup(Peak_OneGroup, EIC_Extracted_mzdata, Color_Array,1,"Peak List Group_07")
			
	Peak_List_Trim            <- Peak_List[,4:16]	
	Peak_List_Trim$CompoundID <- rep(0, nrow(Peak_List_Trim))	# Add a new column for Compound ID Information 
	Peak_List_Trim$Compound_MolecularMass <- rep(0, nrow(Peak_List_Trim))	# Add a new column for Compound Molecular Mass Information 
	
	Compoud_ID                <- 1
	GroupID_Array             <- Peak_List_Trim$GroupID
	
	Peak_Compound_List        <- matrix(ncol=ncol(Peak_List_Trim))		
	Peak_Compound_List        <- data.frame(Peak_Compound_List)
	colnames(Peak_Compound_List)<- colnames(Peak_List_Trim)
	
	Peak_Compound_List_Index <- 1		
	for (group_index in min(GroupID_Array):max(GroupID_Array))
    {
#		cat(paste("group_index=",group_index))
		Peak_List_OneGroup  <- Peak_List_Trim[which(GroupID_Array==group_index),]
		mz_list             <- Peak_List_OneGroup$acute_mass 
		PeakID_OneGroup     <- Peak_List_OneGroup$PeakID
       	
		Mass_Matrix_Estimate<- matrix(0, ncol=length(mz_list), nrow=nrow(Isotope_Adduct_Infor_List))
		#Generate the 2D Estimating Mass Matrix according to the current peak group and Adduct Parameter Information List 
		for (MZ_Index in 1:length(mz_list))
		{
			for(Isotope_Adduct_Index in 1:nrow(Isotope_Adduct_Infor_List))
			{				
				Para_N<- Isotope_Adduct_Infor_List[Isotope_Adduct_Index,]$Para_N
				Para_Z<- Isotope_Adduct_Infor_List[Isotope_Adduct_Index,]$Para_CS
				Para_a<- Isotope_Adduct_Infor_List[Isotope_Adduct_Index,]$Para_Mass_Shift
				Mass_Matrix_Estimate[Isotope_Adduct_Index, MZ_Index]<- (mz_list[MZ_Index]*Para_Z- Para_a)/Para_N
			}
		}
		
		#Dump the 2D Matrix into a vector and use hclust 
		Mass_List <- NULL
		for (Row_Index in 1:nrow(Mass_Matrix_Estimate))
		{
			Mass_List <- c(Mass_List, Mass_Matrix_Estimate[Row_Index, ]) 
		}
		
		# Use hclust to cluster the Mass_List 
		hc<- hclust(dist(Mass_List))   
		memb<- cutree(hc, h=mass_tol)
		
		memb<- Group_contain_Process(memb, length(mz_list))
		# Try to find the optimal cluster group that have more than 1 nmuber    
		count <- rep(0,max(memb))
		for (i in 1:length(memb))
		{		 
			count[memb[i]] <- count[memb[i]]+1
		}
		
		# Try to distinguish the isolated peaks, they maybe a compound peak or a fragment of a compound group 
        IsolatePeak_Index_Array<- Find_IsolatedPeak(memb, count, length(mz_list), nrow(Isotope_Adduct_Infor_List))
		# Try to distinguish the compound group that conncected as isotope or adduct relationship
        Group_Compound_ID_Array<-numeric()
        group_array <- which(count>1)   # Each of the Group_array can be a valid compound group memb[which(count>1)]
		
		if(length(group_array)>0) # Only consider the group number is not empty
		{
			# Add some codes for get the isotope group information
			pure_isotope_group_array <-vector("list", 0)
			pure_isotope_group_array_index <- 1
			for (index in 1:length(group_array))
			{
				SameMass_Pos       <- which(memb==group_array[index])
				# From the index to deduce the detailed Row and Column 
				SameMass_Row_Pos   <- floor((SameMass_Pos-1)/ncol(Mass_Matrix_Estimate))+1
				SameMass_Col_Pos   <- (SameMass_Pos-1)%%ncol(Mass_Matrix_Estimate)+1
				
				if((min(SameMass_Row_Pos)==1)&((max(SameMass_Row_Pos)==3)|(max(SameMass_Row_Pos)==2))) 
				{
					pos_col   <- SameMass_Col_Pos[1]
					Mass_temp <- Mass_Matrix_Estimate[1,pos_col]
					if(length(SameMass_Col_Pos)==3) pure_isotope_group_array[[pure_isotope_group_array_index]]<- c(SameMass_Col_Pos, Mass_temp, 0)
					else pure_isotope_group_array[[pure_isotope_group_array_index]]<- c(SameMass_Col_Pos, 0, Mass_temp, 0)
					pure_isotope_group_array_index <- pure_isotope_group_array_index+1
				}
			 }
			 pure_isotope_group_matrix <- matrix(0, nrow=length(pure_isotope_group_array), ncol=5)
			 if(length(pure_isotope_group_array)>0)
			 {
				 for(index in 1:length(pure_isotope_group_array))
				 {
					 pure_isotope_group_matrix[index,]<- pure_isotope_group_array[[index]]						 
				 }	
			 }			 
			 
			 # Add some codes for get the isotope group information over
			
			# Try to get each compound group peaks 
			for (index in 1:length(group_array))
			{
				SameMass_Pos       <- which(memb==group_array[index])
				# From the index to deduce the detailed Row and Column 
				SameMass_Row_Pos   <- floor((SameMass_Pos-1)/ncol(Mass_Matrix_Estimate))+1
				SameMass_Col_Pos   <- (SameMass_Pos-1)%%ncol(Mass_Matrix_Estimate)+1	
				
				if(!((length(SameMass_Row_Pos)==2)&((min(SameMass_Row_Pos)==2)|(min(SameMass_Row_Pos)==3)))) # Must exclude the default isotop group case, in this case, only [M+1], [M+2] can be found, but the [M+0] can't be found.  
				{
					# Do Adduct Annotation for each peak according to the corresponding adduct information   #if(!((min(SameMass_Row_Pos)==1)&(max(SameMass_Row_Pos)==2)))
					if(!((min(SameMass_Row_Pos)==1)&((max(SameMass_Row_Pos)==3)|(max(SameMass_Row_Pos)==2)))) 	            
					{
						for(Annote_Index in 1:length(SameMass_Pos))
						{
							pos_col<- SameMass_Col_Pos[Annote_Index]
							pos_row<- SameMass_Row_Pos[Annote_Index]
							Molecular_Mass <- Mass_Matrix_Estimate[pos_row,pos_col]
							
							Peak_List_OneGroup[pos_col,]$CompoundID            <- Compoud_ID
							Peak_List_OneGroup[pos_col,]$Compound_MolecularMass<- Molecular_Mass 
							# Save this peak into Peak_Compound_List				
							Peak_Compound_List[Peak_Compound_List_Index, ]<-RowCopy_Data_Frame(Peak_Compound_List[Peak_Compound_List_Index, ], Peak_List_OneGroup[pos_col,])				
							Peak_Compound_List_Index <- Peak_Compound_List_Index+1	
							
							# Add codes for considering the isotoping peak for adduct peak.
							find_isotope<- which(pure_isotope_group_matrix[,1]==pos_col)
							if(length(find_isotope)==1) 
							{   
								isotope_Group_col <-  ncol(pure_isotope_group_matrix)
								for(isotope_Group_col_index in 2:(isotope_Group_col-2))
								{
									pos_col <- pure_isotope_group_matrix[find_isotope,isotope_Group_col_index]
									if(pos_col >0)
									{
										Peak_List_OneGroup[pos_col,]$CompoundID            <- Compoud_ID
										Peak_List_OneGroup[pos_col,]$Compound_MolecularMass<- Molecular_Mass 
										Peak_Compound_List[Peak_Compound_List_Index, ]<-RowCopy_Data_Frame(Peak_Compound_List[Peak_Compound_List_Index, ], Peak_List_OneGroup[pos_col,])				
										Peak_Compound_List_Index <- Peak_Compound_List_Index+1	
									}								
								}
								pure_isotope_group_matrix[find_isotope, isotope_Group_col]<- 1
							}
							# Add codes for considering the isotoping peak for adduct peak over					
						}
						Group_Compound_ID_Array<- c(Group_Compound_ID_Array, Compoud_ID)
						Compoud_ID             <- Compoud_ID+1	
					}
				}										
			 }
			 
			 # Add some codes for considering the isotope peaks that can be a compound group
             if(nrow(pure_isotope_group_matrix)>0)
			 {
				 for(index in 1:nrow(pure_isotope_group_matrix))
				 {
					 isotope_Group_col <-  ncol(pure_isotope_group_matrix)
					 if(pure_isotope_group_matrix[index,ncol(pure_isotope_group_matrix)]==0) # means this group may be a compound, have no relation with other compound's adduct peaks 
					 {
						 for(isotope_Group_col_index in 1:(isotope_Group_col-2))
						 {
							 pos_col <- pure_isotope_group_matrix[index,isotope_Group_col_index]
							 if(pos_col>0)
							 {
								 Peak_List_OneGroup[pos_col,]$CompoundID            <- Compoud_ID
								 Peak_List_OneGroup[pos_col,]$Compound_MolecularMass<- pure_isotope_group_matrix[index,isotope_Group_col-1] 
								 Peak_Compound_List[Peak_Compound_List_Index, ]<-RowCopy_Data_Frame(Peak_Compound_List[Peak_Compound_List_Index, ], Peak_List_OneGroup[pos_col,])				
								 Peak_Compound_List_Index <- Peak_Compound_List_Index+1	 
							 }							 
						 }
						 Group_Compound_ID_Array<- c(Group_Compound_ID_Array, Compoud_ID)
						 Compoud_ID        <- Compoud_ID+1
					 }
				 }	
			 }             
             # Add some codes for considering the isotope peaks that can be a compound group over			 
		 }
				  
		 #Decide each isolated peak is a valid compound or not by comparing its actuate mass with hmdb lib, if found some, mark it as a valid compound too.
         if(length(IsolatePeak_Index_Array)>0)
		 {
			 for(IsolatePeak_Index in 1:length(IsolatePeak_Index_Array))
			 {
				 temp<-IsolatePeak_Index_Array[IsolatePeak_Index]
				 IsolatePeak_actuateMass <- Peak_List_OneGroup[temp,]$acute_mass- 1.007276 # Suppose the isolate peak output MW [M+H]format
				 if(is.factor(hmdb_library$weight_mono)) lib_mass<- as.numeric(as.character(hmdb_library$weight_mono))
				 else lib_mass<- hmdb_library$weight_mono
					 
				 match_result<-which((lib_mass>(IsolatePeak_actuateMass-mass_tol))&(lib_mass<(IsolatePeak_actuateMass+mass_tol)))
				 if(length(match_result)>0) # Found some good match, mark the isolated peak as a vaild compound 
				 {
					 Peak_List_OneGroup[temp,]$CompoundID            <- Compoud_ID
					 Peak_List_OneGroup[temp,]$Compound_MolecularMass<- IsolatePeak_actuateMass
					 # Save this peak into Peak_Compound_List	
					 Peak_Compound_List[Peak_Compound_List_Index, ]<-RowCopy_Data_Frame(Peak_Compound_List[Peak_Compound_List_Index, ], Peak_List_OneGroup[temp,])				
					 Peak_Compound_List_Index <- Peak_Compound_List_Index+1
					 Compoud_ID        <- Compoud_ID+1
				 }else if(length(Group_Compound_ID_Array)==0){    # Can't find some good match,the compound group is empty, keep these peaks as unknown peaks  
					 Peak_Compound_List[Peak_Compound_List_Index, ]<-RowCopy_Data_Frame(Peak_Compound_List[Peak_Compound_List_Index, ], Peak_List_OneGroup[IsolatePeak_Index,])				
					 Peak_Compound_List_Index <- Peak_Compound_List_Index+1					 
				 }else{	 # Can't find some good match and the compound group is not empty, Check whether this peak can be clustered into the group_compound or not
					 for (Group_Compound_ID_Array_Index in 1:length(Group_Compound_ID_Array))
					 {
						 Peak_List_Temp_Group   <- Peak_Compound_List[which(Peak_Compound_List$CompoundID==Group_Compound_ID_Array[Group_Compound_ID_Array_Index]),]
						 Peak_List_Temp_Group[nrow(Peak_List_Temp_Group)+1,] <- Peak_List_OneGroup[temp,]
						 Peak_Data_Temp_Group   <- Generate_Peak_DataArray(Peak_List_Temp_Group, EIC_Extracted)
						 #Peak_Data_Temp_Isolate <- Generate_Peak_DataArray(Peak_List_OneGroup[IsolatePeak_Index,], EIC_Extracted)
						 angle_sum<- 0
						 for(Peak_List_Temp_Group_Index in 1:(nrow(Peak_List_Temp_Group)-1))
						 {
							 vec1 <- Peak_Data_Temp_Group[,Peak_List_Temp_Group_Index]
							 vec2 <- Peak_Data_Temp_Group[,nrow(Peak_List_Temp_Group)]
							 angle_sum <- angle_sum + Vectors_dotproc(vec1,vec2)
						 }
						 
						 if(angle_sum/(nrow(Peak_List_Temp_Group)-1)<Refinement_Cluster_Angle_Th)   # if(angle_sum/(nrow(Peak_List_Temp_Group)-1)<15)
						 {
							 Peak_List_OneGroup[temp,]$CompoundID            <- Group_Compound_ID_Array[Group_Compound_ID_Array_Index]
							 Peak_List_OneGroup[temp,]$Compound_MolecularMass<- 0
							 Peak_Compound_List[Peak_Compound_List_Index, ]<-RowCopy_Data_Frame(Peak_Compound_List[Peak_Compound_List_Index, ], Peak_List_OneGroup[temp,])				
							 Peak_Compound_List_Index <- Peak_Compound_List_Index+1
						 }					 
					 }
				 }
			 }
		 }         
		
    }# End for one group
	
	#write.csv(Peak_Compound_List, file = "/home/wenchao/Projects/Rproject/Peaklist_w3_7_20111213_filter_Group_Isotope_Adduct_Annotate_Refinement.csv")	#Peaklist_w2_07_filter_Group_Isotope_Adduct_Annotate_Refinement.csv
	pos<-regexpr('.csv',PeakFileName)
	PeakFileName_Isotope_Adduct_Annotate_Refinement <-paste(substring(PeakFileName, first=1, last=pos-1), '_Refinement', '.csv', sep="")
	write.csv(Peak_Compound_List, file = PeakFileName_Isotope_Adduct_Annotate_Refinement)
	PeakFileName_Isotope_Adduct_Annotate_Refinement
}

#This function can get the Compound's Retention time by LUT of the scan_retention_time map for compound's scan value. 
#The scan value of Compound_Apex, Compound_Left_Boundary and Compound_Right_Boundary can be got by the median scan value of each compound's annaotation peaks scan value.
# Here, we think the annotation peaks have more confidence than fragment peaks to a compound.
Acquire_Compound_RetentionIime_Map <- function(Compound_List_In, Retention_Time_Scan_Map)
{
	
	Compound_List_Out<- Compound_List_In
	CompoundID_Low <- max(1, min(Compound_List_Out$CompoundID))  # Must exclude the cases of CompoundID==0, becuase CompoundID==0 means the the peaks can't be grouped into any compound cluster, and also these peak can't be identified from the hmdb libary.  
	CompoundID_High<- max(Compound_List_Out$CompoundID)
	
	Compound_List_Out$Compound_Apex          <- rep(0, nrow(Compound_List_Out))
	Compound_List_Out$Compound_Left_Boundary <- rep(0, nrow(Compound_List_Out))
	Compound_List_Out$Compound_Right_Boundary<- rep(0, nrow(Compound_List_Out))
	Compound_List_Out$Compound_Intensity     <- rep(0, nrow(Compound_List_Out))
	
	for(Id_index in CompoundID_Low:CompoundID_High)
	{
		One_Compound <- Compound_List_Out[which(Compound_List_Out$CompoundID ==Id_index),]
		One_Compound_Valid <- One_Compound[which(One_Compound$Compound_MolecularMass>0),]
		Compound_Apex_Scan <- ceiling(median(One_Compound_Valid$Apex_Pos)) #Use ceiling function to avoid the fragment Apex_Scan case
		Compound_Left_Scan <- ceiling(median(One_Compound_Valid$Left_Boundary))
		Compound_Right_Scan<- ceiling(median(One_Compound_Valid$Right_Boundary))
		Compound_Intensity <- max(One_Compound_Valid$Peak_Intensity)
		
		Compound_Apex_minutes          <-Retention_Time_Scan_Map[which(Retention_Time_Scan_Map$scan_num == Compound_Apex_Scan),2]
		Compound_Left_Boundary_minutes <-Retention_Time_Scan_Map[which(Retention_Time_Scan_Map$scan_num == Compound_Left_Scan),2]
		Compound_Right_Boundary_minutes<-Retention_Time_Scan_Map[which(Retention_Time_Scan_Map$scan_num == Compound_Right_Scan),2]
		
		Compound_List_Out[which(Compound_List_Out$CompoundID ==Id_index),]$Compound_Apex          <-Compound_Apex_minutes
		Compound_List_Out[which(Compound_List_Out$CompoundID ==Id_index),]$Compound_Left_Boundary <-Compound_Left_Boundary_minutes
		Compound_List_Out[which(Compound_List_Out$CompoundID ==Id_index),]$Compound_Right_Boundary<-Compound_Right_Boundary_minutes
		Compound_List_Out[which(Compound_List_Out$CompoundID ==Id_index),]$Compound_Intensity     <-Compound_Intensity
	}		
	Compound_List_Out
	
}

# For some certain reasons, we can get some groups with different molecular mass, but maybe one small group can be contained wholely in the big group. So, we need to consider this case.
Group_contain_Process <- function (memb, Peak_Num)
{
	memb_new <- memb
	count <- rep(0,max(memb))
	for (i in 1:length(memb))
	{		 
		count[memb[i]] <- count[memb[i]]+1
	}
	
	group_array_tem <- which(count>1)
	if(length(group_array_tem)>1)
	{
		# Use a list to store the group that contain the peak information of each group 
		group_list      <- vector("list", 0)
		for(group_list_index in 1:length(group_array_tem))
		{
			group_list[[group_list_index]] <- (which(memb==group_array_tem[group_list_index])-1)%%Peak_Num+1
		}
		
		group_valid_flag<- rep(1,length(group_array_tem))
		# Use a 2-loop to determine the group that contained in a large group, and should be considered as a void group 
		for (group_list_index_1 in 1:(length(group_array_tem)-1))
		for (group_list_index_2 in (group_list_index_1+1):length(group_array_tem))	
		{
			Result <- Set_Relation_Decision(group_list[[group_list_index_1]], group_list[[group_list_index_2]])  
			if(Result==1) group_valid_flag[group_list_index_1]<- 0
			if(Result==2) group_valid_flag[group_list_index_2]<- 0
		}
	    
	    for(group_list_index in 1:length(group_array_tem))
		{
			if(group_valid_flag[group_list_index]==0) 
			{
				memb_index_arry <- which(memb==group_array_tem[group_list_index])
				for (memb_index_arry_index in 2:length(memb_index_arry))
				{
					memb_new_index <-memb_index_arry[memb_index_arry_index]
					memb_new[memb_new_index]<- max(memb_new)+1
				}				
			}					
		}	
	    
	}		
		
	memb_new
}

# Set relation decision, Result =0: Set1==Set2
#                        Result =1: Set2 Contain Set1   , Set1<Set2
#                        Result =2: Set1 Contain Set2   , Set1>Set2
#                        Result =3: Set1 InterCross Set2, Set1|=Set2 
Set_Relation_Decision<- function(Set1, Set2)
{	
	Result12<-1
	for(Set1_Index in 1:length(Set1))
	{
		if(length(which(Set2== Set1[Set1_Index]))==0) Result12 <- 0
	}
	
	Result21<-1
	for(Set2_Index in 1:length(Set2))
	{
		if(length(which(Set1== Set2[Set2_Index]))==0) Result21 <- 0
	}
	
	if((Result12==1)&(Result21==1)){ # Set1=Set2
	   Result<- 0
    }else if((Result12==1)&(Result21==0)){ # Set1<Set2	
	   Result<-1
    }else if((Result12==0)&(Result21==1)){ # Set1>Set2
	   Result<-2    
	}else{ # Set1|=Set2 
	   Result<-3
	}                                     
	
	Result
}

# According to the memb and count information to parse the isolated peaks 
Find_IsolatedPeak <- function(memb, count, Peak_Num, Adduct_Num)
{
		
	if(length(memb)!=Peak_Num*Adduct_Num) stop("memb size must equal to Peak_Num*Adduct_Num!")
	if(length(count)>Peak_Num*Adduct_Num) stop("count size must smaller than Peak_Num*Adduct_Num!")
	
	IsolatedPeaks<- numeric()
	
	for(Peak_Index in 1:Peak_Num)
	{
	   	Peak_Adduct_pos <- seq(Peak_Index, Peak_Num*Adduct_Num, by=Peak_Num)
		Peak_Adduct_gro <- memb[Peak_Adduct_pos]
		if(sum(count[Peak_Adduct_gro])==Adduct_Num) IsolatedPeaks<- c(IsolatedPeaks, Peak_Index)
	}	
	IsolatedPeaks
}

RowCopy_Data_Frame <- function(Data_Frame_Row_Dest, Data_Frame_Row_Sour)
{
	if(nrow(Data_Frame_Row_Dest)!=1) stop("Row copy for DataFrame must be 1 row!")
	if(nrow(Data_Frame_Row_Sour)!=1) stop("Row copy for DataFrame must be 1 row!")
	if(nrow(Data_Frame_Row_Dest)!=nrow(Data_Frame_Row_Sour)) stop("Row copy for DataFrame must be have row!")
	Out_Data_Frame_Row <- Data_Frame_Row_Dest
	for(i in 1:length(Data_Frame_Row_Sour))
	{
		if(is.factor(Data_Frame_Row_Sour[1,i]))
		{
			Out_Data_Frame_Row[1,i]<- as.character(Data_Frame_Row_Sour[1,i]) # Convert factor into character 			
		}else{
			Out_Data_Frame_Row[1,i]<- Data_Frame_Row_Sour[1,i]	
		}
	}
	Out_Data_Frame_Row	
}

#Peak filtering according to the peak quality such as peak boundary range, peak signifiance. Try to keep the good peaks
#PeakQuality_filter <- function(Peak_list_sub, EIC_Extracted)	
#{ 	
#	SNR_Th         <- 10.0
#	Boundary_Th_Min<- 7
#	Boundary_Th_Max<- 60
#	Significance_Th<- 1.1
#	j<- 1
#    Peak_list_filter <-data.frame(PeakID=0,EIC_ID=0,Apex_Pos=0,acute_mass=0,Peak_Intensity=0,Peak_SNR=0, Left_Boundary=0,Right_Boundary=0,Peak_Area=0,  Peak_Significance=0)	
#	
#	for (i in 1:nrow(Peak_list_sub))
#	{
#		SNR              <-Peak_list_sub[i,]$Peak_SNR
#		Left_Boundary    <-Peak_list_sub[i,]$Left_Boundary
#		Right_Boundary   <-Peak_list_sub[i,]$Right_Boundary
#		EIC_Index        <-Peak_list_sub[i,]$EIC_ID
#		Apex_pos         <-Peak_list_sub[i,]$Apex_Pos
#		Peak_Intensity   <-Peak_list_sub[i,]$Peak_Intensity		
#	#	Peak_Significance<-Peak_list_sub[i,]$Peak_Significance
#		
#		EIC_Scan         <-EIC_Extracted[[EIC_Index]]$scan_index
#		EIC_Intensity    <-EIC_Extracted[[EIC_Index]]$intensity
#				
#		Peak_Index <-which((EIC_Scan<Right_Boundary)&(EIC_Scan>Left_Boundary))
#    	Peak_Significance <- (length(Peak_Index)-1)*Peak_Intensity/(sum(EIC_Intensity[Peak_Index])- Peak_Intensity)
#				
#		Peak_Boundary  <-Right_Boundary-Left_Boundary		
#		
#		if((SNR>SNR_Th)&(Peak_Boundary>Boundary_Th_Min)&(Peak_Boundary<Boundary_Th_Max)&(Peak_Significance>Significance_Th))     #&(SNR<10000)
#		{
#			Peak_list_filter[j,] <- Peak_list_sub[i,]
#			j<-j+1
#		}
#	}
#	Peak_list_filter
#}

PeakQuality_filter <- function(Peak_list_sub, EIC_Extracted, SNR_Th, PeakWidth_Low, PeakWidth_High)	
{			
	#SNR_Th         <- 10.0
	Boundary_Th_Min<- PeakWidth_Low     #7
	Boundary_Th_Max<- PeakWidth_High    #60
	Significance_Th<- 1.5
	j<- 1
	Peak_list_filter <-data.frame(PeakID=0,EIC_ID=0,Apex_Pos=0,acute_mass=0,Peak_Intensity=0,Peak_SNR=0,Left_Boundary=0,Right_Boundary=0,Peak_Area=0,Peak_Significance=0)	
	
	for (i in 1:nrow(Peak_list_sub))
	{
		SNR              <-Peak_list_sub[i,]$Peak_SNR
		Left_Boundary    <-Peak_list_sub[i,]$Left_Boundary
		Right_Boundary   <-Peak_list_sub[i,]$Right_Boundary
		EIC_Index        <-Peak_list_sub[i,]$EIC_ID
		Apex_pos         <-Peak_list_sub[i,]$Apex_Pos
		Peak_Intensity   <-Peak_list_sub[i,]$Peak_Intensity		
		Peak_Significance<-Peak_list_sub[i,]$Peak_Significance
		Peak_Area        <-Peak_list_sub[i,]$Peak_Area
		
		EIC_Scan         <-EIC_Extracted[[EIC_Index]]$scan_index
		EIC_Intensity    <-EIC_Extracted[[EIC_Index]]$intensity

		Peak_Boundary         <-Right_Boundary-Left_Boundary+1		
		Area_Estimate_Triangle<-(Peak_Intensity*0.5*Peak_Boundary)		
		Area_Estimate_Bias    <- abs(Area_Estimate_Triangle-Peak_Area)/(Area_Estimate_Triangle+Peak_Area)
		
		if((SNR>SNR_Th)&(Peak_Boundary>Boundary_Th_Min)&(Peak_Boundary<Boundary_Th_Max)&(Peak_Significance>Significance_Th)&(Area_Estimate_Bias<0.25))     #&(SNR<100000)
		{
			Peak_list_filter[j,] <- Peak_list_sub[i,]
			j<-j+1
		}
	}
	Peak_list_filter
}

#calculate dotproductor of two vectors, require the two inputting vectors are of the same length and aligned in scan  
Vectors_dotproc <- function(vec1, vec2)
{
	norm_vec1<- vec1/sqrt(sum(vec1^2))
	norm_vec2<- vec2/sqrt(sum(vec2^2))
	dotproduct <- crossprod(norm_vec1,norm_vec2)
	angle<- acos(dotproduct)*180/pi
	angle
}

#Pairwisely Calculate the distance of the peak profile. The inputting x is a matrix in formats as following:
#[PeakData1, PeakData2, ....]. Each peak data have the same data point. If they are different, you should padding 
#zeros in all of the position that out off the range of [left_boundary, right_boundary]. 
#The return value is a distance matrix. 

distCal<-function (x) 
{
	ncy <- ncx <- ncol(x)
	if (ncx == 0) 
		stop("'x' is empty")
	r <- matrix(0, nrow = ncx, ncol = ncy)
	for (i in seq_len(ncx)) {
		for (j in seq_len(i)) {
			x2 <- x[, i]
			y2 <- x[, j]
			r[i, j] <- ifelse(i==j,0,Vectors_dotproc(x2, y2))
		}
	}
	r <- r + t(r) - diag(diag(r))
	rownames(r) <- colnames(x)
	colnames(r) <- colnames(x)
	r
}

# Generate a smooth peak data array according to the peak_list information and EIC_Extracted data. 
# The rule for smoothing: 1.)linear interpolation for the missing points that fall in the range [left_boundary, right_boundary] according its neighboring points.
#                         2.)Padding zero for the points that fall off the range [left_boundary, right_boundary], because the inputting peaks are usually not the same. 
# The outputting is a peak profile data matrix as [PeakData1, PeakData2, ....]. The col_num = Peak number and row_num= Max(Peak data length). 
Generate_Peak_DataArray<- function(Peak_List, EIC_Extracted)
{	
	Right<-max(Peak_List$Right_Boundary)
	Left <-min(Peak_List$Left_Boundary)
	Peak_Length   <- Right- Left+1  
	Peak_Num      <- nrow(Peak_List)
	Peak_DataArray<- matrix(0, nrow = Peak_Length, ncol = Peak_Num)
	
	for(i in 1:Peak_Num)
	{
	   	EIC_Index         <- Peak_List[i,]$EIC_ID
		Peak_LeftBoundary <- Peak_List[i,]$Left_Boundary
		Peak_RightBoundary<- Peak_List[i,]$Right_Boundary
		Peak_Data_Smooth  <- Peak_Smooth_Interpolate(EIC_Extracted[[EIC_Index]]$scan_index, EIC_Extracted[[EIC_Index]]$intensity, Peak_LeftBoundary, Peak_RightBoundary)
		Peak_DataArray[,i]<- c(rep(0, (Peak_LeftBoundary-Left)), Peak_Data_Smooth, rep(0, Right-Peak_RightBoundary))	
	}
	
	colnames(Peak_DataArray) <- Peak_List$PeakID
	Peak_DataArray
}
# Generate a linear interpolating peak for missing peaks that fall in the range of [Peak_LeftBoundary, Peak_RightBoundary]
Peak_Smooth_Interpolate <- function(scan_array, intensity_array, Peak_LeftBoundary, Peak_RightBoundary)
{		
	Peak_Interpolate_Smoothing <-numeric(length=(Peak_RightBoundary-Peak_LeftBoundary+1))
	
	for (i in Peak_LeftBoundary:Peak_RightBoundary)
	{
		ret<-which(scan_array==i)
		if(length(ret)==0) # can't find the point in scan_array. missing this point, use the linear interpolation of the neighboring point
		{   
			if((i>=min(scan_array))&(i<=max(scan_array)))
			{
			   scan_array_left <- scan_array[which(scan_array<i)]
			   scan_array_right<- scan_array[which(scan_array>i)]
			   ret_left        <- which.min(i-scan_array_left)
			   ret_right       <- which.min(scan_array_right-i)
			   scan_left       <- scan_array_left[ret_left]
			   scan_right      <- scan_array_right[ret_right]
			   index_left      <- which(scan_array==scan_left)
			   index_right     <- which(scan_array==scan_right)
			
			   lamda <- 1.0*(i-scan_left)/(scan_right-scan_left)
			   beta  <- 1.0*(scan_right-i)/(scan_right-scan_left)
						
			   Peak_Interpolate_Smoothing[i-Peak_LeftBoundary+1] <- beta*intensity_array[index_left]+ lamda*intensity_array[index_right]
		     }
			
		}else { # find the point in scan_array. use it 
			Peak_Interpolate_Smoothing[i-Peak_LeftBoundary+1] <- intensity_array[ret]
		}		
	}	
	Peak_Interpolate_Smoothing 
}
