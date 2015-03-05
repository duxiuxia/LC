#This files contain all of the function modules related to Compound Alignment Algorithm

AlignmentFile_Load <- function(FileFolderPath)
{
	Peak_FileName_List <-list.files(path=FileFolderPath, pattern ="Group_Annotate_Isotope_Annotate_Adduct_Refinement_RtMap.csv", full.names=TRUE)
	Peak_List <- read.csv(file=Peak_FileName_List[1],head=TRUE,sep=",")
	Compound_List_Temp <- Peak_List
	Compound_List_Temp$FileID <- rep(1, nrow(Peak_List))
	Compound_List <- Compound_List_Temp
	
	if(length(Peak_FileName_List)>1)
	{
		for (File_ID in 2:length(Peak_FileName_List))
		{
			Peak_List <- read.csv(file=Peak_FileName_List[File_ID],head=TRUE,sep=",")
			Compound_List_Temp <- Peak_List
			Compound_List_Temp$FileID <- rep(File_ID, nrow(Peak_List))
			Compound_List <- rbind(Compound_List, Compound_List_Temp)		
		}
	}	

     Compound_List
}

#This function can get the Compound's Retention time by LUT of the scan_retention_time map for compound's scan value. 
#The scan value of Compound_Apex, Compound_Left_Boundary and Compound_Right_Boundary can be got by the median scan value of each compound's annaotation peaks scan value.
# Here, we think the annotation peaks have more confidence than fragment peaks to a compound.
Acquire_Compound_RetentionIme_Map <- function(Compound_List_In, Retention_Time_Scan_Map)
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

# Alignment process is for aligning the retention time across samples. So, we need convert the scan number into rention time for each sample.  
Compound_Alignment <- function(Compound_List, Align_Window_Phase1, Align_Window_Phase2, Mass_Tol, Two_Phase, Aligned_FileName)
{	
	
	Compound_List$Align_ID <-rep(0, nrow(Compound_List))   # Align_ID is the final unique ID acorss samples. 
	Compound_List$Align_RT <-rep(0, nrow(Compound_List))   # Align_RT is the final aligned Retention time across samples.
	
	Compound_List_ForAlign <-Compound_List[which((Compound_List$Align_ID==0)&(Compound_List$CompoundID>0)),]
	Align_ID <- 1
	
	while(nrow(Compound_List_ForAlign)>0)     #&(Compound_List$Compound_MolecularMass>0)
	{
		cat(paste("Align_ID=",Align_ID, "\r", sep=""))
		Compound_Intensity_Array                <- Compound_List_ForAlign$Compound_Intensity 
		Compound_MolecularMass_Array            <- Compound_List_ForAlign$Compound_MolecularMass
		
		# Calculate the index of Align Target Peak  
		Index_Peak_MeetMass_Larger_Zero         <- which(Compound_MolecularMass_Array>0)
		Index_Peak_MaxIntensity_Mass_Larger_Zero<- Index_Peak_MeetMass_Larger_Zero[which.max(Compound_Intensity_Array[Index_Peak_MeetMass_Larger_Zero])]
		
		# Get the Compound_Mass and Retention time of the Align Target Compound
		FileID_Max_Intensity     <- Compound_List_ForAlign[Index_Peak_MaxIntensity_Mass_Larger_Zero,]$FileID
		CompoundID_Max_Intensity <- Compound_List_ForAlign[Index_Peak_MaxIntensity_Mass_Larger_Zero,]$CompoundID
		
		Target_Compound_index     <- which((Compound_List_ForAlign$FileID==FileID_Max_Intensity)&(Compound_List_ForAlign$CompoundID==CompoundID_Max_Intensity)&(Compound_List_ForAlign$Compound_MolecularMass>0))
		AlignTarget_Compound_Mass <- median(Compound_List_ForAlign[Target_Compound_index,]$Compound_MolecularMass)
		AlignTarget_RT            <- median(Compound_List_ForAlign[Target_Compound_index,]$Compound_Apex)
		
#		AlignTarget_Compound_Mass               <- Compound_List_ForAlign[Index_Peak_MaxIntensity_Mass_Larger_Zero,]$Compound_MolecularMass
#		AlignTarget_RT                          <- Compound_List_ForAlign[Index_Peak_MaxIntensity_Mass_Larger_Zero,]$Compound_Apex
		
		cat(paste("AlignTarget_Compound_Mass=",AlignTarget_Compound_Mass, "\r", sep=""))
		cat(paste("AlignTarget_RT=",AlignTarget_RT, "\r", sep=""))
		
		
		FileID_Min                              <- min(Compound_List_ForAlign$FileID)
		FileID_Max                              <- max(Compound_List_ForAlign$FileID)
		
		# Use Compound_Align_List to store the compounds that should be aligned to the group of AlignTarget  
		Compound_Align_List                     <- NULL
		Compound_Mass_Align_Array               <- NULL
		RT_Align_Array                          <- NULL
		
		for(index_FileID in FileID_Min:FileID_Max)
		{
			Compound_List_ForAlign_File <- Compound_List_ForAlign[which(Compound_List_ForAlign$FileID==index_FileID),]
			CompoundID_Min              <- min(Compound_List_ForAlign_File$CompoundID)
			CompoundID_Max              <- max(Compound_List_ForAlign_File$CompoundID)
			
			for(index_CompoundID in CompoundID_Min:CompoundID_Max)
			{
				Compound_List_ForAlign_File_Compound       <- Compound_List_ForAlign_File[which(Compound_List_ForAlign_File$CompoundID==index_CompoundID),]
				Compound_List_ForAlign_File_Compound_Mass  <- Compound_List_ForAlign_File_Compound[which(Compound_List_ForAlign_File_Compound$Compound_MolecularMass>0),]
				
				# If the compound has been aligned already, we need neglect it
		        if(nrow(Compound_List_ForAlign_File_Compound_Mass)>0)
				{
					Compound_MolecularMass_Temp                <- median(Compound_List_ForAlign_File_Compound_Mass$Compound_MolecularMass)
					Compound_RT_Temp                           <- median(Compound_List_ForAlign_File_Compound_Mass$Compound_Apex)
									
					#Compare the compound's retention time and molecular mass fall in the setting align window of the align target compound or not
					if((abs(Compound_RT_Temp-AlignTarget_RT)<(Align_Window_Phase1/2.0))&(abs(Compound_MolecularMass_Temp-AlignTarget_Compound_Mass)<Mass_Tol))
					{					
						Compound_Align_List       <- rbind(Compound_Align_List, Compound_List_ForAlign_File_Compound)
						Compound_Mass_Align_Array <- c(Compound_Mass_Align_Array, Compound_MolecularMass_Temp)
						RT_Align_Array            <- c(RT_Align_Array, Compound_RT_Temp)
					}
				}
		     }				
		 }
		
		 if(Two_Phase)
		 {
			# According to the Compound list that can be aligned to the target compound, calculate the statistical center of the cluster compounds,then align the compounds to the statistical center      
			AlignTarget_RT           <- median(RT_Align_Array)     #(Compound_Align_List$Compound_Apex)
			AlignTarget_Compound_Mass<- median(Compound_Mass_Align_Array) 
			Compound_Align_List      <- NULL
			Compound_Mass_Align_Array<- NULL
			RT_Align_Array           <- NULL
			
			for(index_FileID in FileID_Min:FileID_Max)
			{
				Compound_List_ForAlign_File <- Compound_List_ForAlign[which(Compound_List_ForAlign$FileID==index_FileID),]
				CompoundID_Min              <- min(Compound_List_ForAlign_File$CompoundID)
				CompoundID_Max              <- max(Compound_List_ForAlign_File$CompoundID)
				
				for(index_CompoundID in CompoundID_Min:CompoundID_Max)
				{
					Compound_List_ForAlign_File_Compound       <- Compound_List_ForAlign_File[which(Compound_List_ForAlign_File$CompoundID==index_CompoundID),]
					Compound_List_ForAlign_File_Compound_Mass  <- Compound_List_ForAlign_File_Compound[which(Compound_List_ForAlign_File_Compound$Compound_MolecularMass>0),]
					if(nrow(Compound_List_ForAlign_File_Compound_Mass)>0)
					{
						Compound_MolecularMass_Temp                <- median(Compound_List_ForAlign_File_Compound_Mass$Compound_MolecularMass)
						Compound_RT_Temp                           <- median(Compound_List_ForAlign_File_Compound_Mass$Compound_Apex)
						
						#Compare the compound's retention time and molecular mass fall in the setting align window of the align target compound or not
						if((abs(Compound_RT_Temp-AlignTarget_RT)<Align_Window_Phase2/2.0)&(abs(Compound_MolecularMass_Temp-AlignTarget_Compound_Mass)<Mass_Tol))
						{					
							Compound_Align_List       <- rbind(Compound_Align_List, Compound_List_ForAlign_File_Compound)	
							Compound_Mass_Align_Array <- c(Compound_Mass_Align_Array, Compound_MolecularMass_Temp)
							RT_Align_Array            <- c(RT_Align_Array, Compound_RT_Temp)
						}
					}					
				}
			 }
		 }		

		# Update the align result to Compound_List
		if(nrow(Compound_Align_List)>0)
		{
			for(index_Align in 1:nrow(Compound_Align_List))
			{
				FileID_Temp     <- Compound_Align_List[index_Align, ]$FileID
				PeakID_Temp     <- Compound_Align_List[index_Align, ]$PeakID
				CompoundID_Temp <- Compound_Align_List[index_Align, ]$CompoundID
				
				index_Peak  <- which((Compound_List$FileID == FileID_Temp)&(Compound_List$PeakID == PeakID_Temp)&(Compound_List$CompoundID ==CompoundID_Temp))
				Compound_List[index_Peak,]$Align_ID <-Align_ID
				Compound_List[index_Peak,]$Align_RT <-median(RT_Align_Array)
			}
			Align_ID <- Align_ID+1
		}		
		cat(paste("Num_Compound_List_ForAlign=",nrow(Compound_List_ForAlign), "\r", sep=""))
		Compound_List_ForAlign <-Compound_List[which((Compound_List$Align_ID==0)&(Compound_List$CompoundID>0)),]
	}	
	write.csv(Compound_List, file = Aligned_FileName)  #Peaklist_w2_07_filter_Group.csv
}

