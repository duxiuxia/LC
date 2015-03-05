# TODO: Add comment
# 
# Author: wenchao
# realize remove a node from list at poistion
remove_list_node<- function(existing_list, position)
#remove_list_node<- function(list, position)
{
	new_list <- NULL
	
	if (position == 1)
	{
		new_list <- existing_list[2:length(existing_list)]	
	}
	else if ((position > 1) & (position < length(existing_list)))
	{
		new_list <- c(existing_list[1:(position-1)], existing_list[(position+1):length(existing_list)])
	}
	else
	{
		new_list <- existing_list[1:(position-1)]
#		stop("Error: node to be removed is out of range \n")
	}
	
	return(new_list)
	
############################################################
# Code below was commented out by Xiuxia on 5/16/2012 to speed up the process
#	list_new <-vector("list", 0)
#	list_num=length(list)
#	
#	if (position ==1)
#		list_new= list[2:list_num]
#	else if((position>1)&(position<(list_num+1)))
#	{         
#		list_new= list[1:position-1]
#		for(i in position:list_num-1)
#		{
#			list_new[[i]] = list[[i+1]]
#		}
#	}
#	else
#		stop("POSITION FOR REMOVING NODE IS OUTRANGE!")
#	
#	list_new
############################################################	
}

# realize add a list_node(data_frame) at position 
add_list_node <- function(existing_list, position, new_node)
{
	new_list_node <- NULL
	new_list_node[[1]] <- new_node
	
	new_list <- NULL
	
	#browser()
	
	if (position == 1)
	{
		new_list <- c(new_list_node, existing_list)
	}
	else if ((position > 1) & (position < length(existing_list)))
	{
		new_list <- c(existing_list[1:(position-1)], new_list_node, existing_list[position:length(existing_list)])
	}
	else
	{
		new_list <- c(existing_list, new_list_node)
	}
	return(new_list)
	
##############################################################
# Code below was commented out by Xiuxia on 5/16/2012 to speed up the process
#	list_new <-vector("list", 0)
#	
#	list_num=length(list)
#	
#	if(position==1)
#	{
#		list_new[[1]] = list_node	  
#		for(i in 1:list_num)
#		{
#			list_new[[i+1]] = list[[i]]
#		}          
#	}
#	else if((position>1)&(position<=(list_num+1))) 
#	{   
#		list_new= list[1:position-1]
#		list_new[[position]] = list_node
#		if(position!=(list_num+1))
#		{ 
#			for(i in position:list_num)
#			{
#				list_new[[i+1]] = list[[i]]
#			}
#		}
#	}
#	else
#	{
#		cat("adding list node, postion is", position, "\r")
#		cat("list length", length(list), "\r")
#		stop("POSITION FOR ADDING NODE IS OUTRANGE!")
#	}	 
#	
#	list_new
##############################################################	
}

# Find the optimal matched EIC_Masstrace index for the new adding EICData_Node from existing EICInformation_list 
# the returned value (Masstrace index) 0 : Failed to find matched mass trace, need to insert EIC_Masstrace to list
#--------------------------------------------------------------------------------------------------------------
Find_optimal_matched_EIC_Masstrace <- function(EICData_Node, EICInformation_list, ppm_Th)
#Find_optimal_matched_EIC_Masstrace <- function(EICData_Node, EICInformation_list, delta_mass_Th)
{	
	current_node_mass      <- EICData_Node$mass
#	current_node_intensity <- EICData_Node$intensity
#	current_node_scanindex <- EICData_Node$scan_index
	
	vector_mass_cent <- NULL
	for (candidate_index in 1:length(EICInformation_list))
	{
		vector_mass_cent <- c(vector_mass_cent, EICInformation_list[[candidate_index]]$mass_cent)
	}
	
	
	vector_distance_in_ppm <- abs((current_node_mass - vector_mass_cent)) / vector_mass_cent * 1000000
	
	II <- which(vector_distance_in_ppm == min(vector_distance_in_ppm))
	
	if (vector_distance_in_ppm[II] <= ppm_Th)
	{
		return(II)
	}
	else
	{
		return(NULL)
	}
	
	
#############################################################
#	The code below looking for the nearest mass trace is commented out by Xiuxia on 5/15/2012 	
#	mass_lowerbound<- current_node_mass - delta_mass_Th
#	mass_upperbound<- current_node_mass + delta_mass_Th
#	dist<- 10000000000000.0 
#	position <- 0

#	for(candidate_index in 1:length(EICInformation_list))
#	{
#		temp_mass<-EICInformation_list[[candidate_index]]$mass_cent 
#		# first, use [mass_lowerbound, mass_upperbound] to select several mass trace candidates
#		if((temp_mass>=mass_lowerbound)&(temp_mass<=mass_upperbound))
#		{
#			newnode_mass      <- EICInformation_list[[candidate_index]]$newnode_mass
#			newnode_intensity <- EICInformation_list[[candidate_index]]$newnode_intensity
#			newnode_scanindex <- EICInformation_list[[candidate_index]]$newnode_scanindex
#			
#			# using weighing eludiance distance to select the best mass trace from the selceted candidiates 
#			# dist_temp<- 0.5*(current_node_mass- newnode_mass)^2/(current_node_mass+newnode_mass)^2 + 0.4*(current_node_intensity-newnode_intensity)^2/(current_node_intensity+newnode_intensity)^2+ 0.1*(current_node_scanindex-newnode_scanindex)^2/(current_node_scanindex+newnode_scanindex)^2
#			# dist_temp<- 0.8*(current_node_mass- newnode_mass)^2 + 0.1*(current_node_intensity-newnode_intensity)^2+ 0.1*(current_node_scanindex-newnode_scanindex)^2
#			Delta_Mass_PPM      <- 1000000*abs(current_node_mass- newnode_mass)/(newnode_mass)
#			Delta_Log_Intensity <- abs(log(current_node_intensity) -log(newnode_intensity))
#			Delta_Scan          <- abs(current_node_scanindex-newnode_scanindex)
#			
#			dist_temp<- 0.6*Delta_Mass_PPM + 0.3*Delta_Log_Intensity + 0.1*Delta_Scan			
#			
#			if(dist_temp<=dist)
#			{
#				position<- candidate_index
#				dist<- dist_temp				
#			}							
#		}
#	}
#	position		
	#############################################################	
}

# Find a position for the new adding EICDataNode, inset a new Mass_Trace, return the suitable position for new adding EICData_Node
Find_position_addingnode <- function(EICData_Node, EICInformation_list)
{
	current_mass <- EICData_Node$mass
	
	vector_mass_cent <- NULL
	
	for (i in 1:length(EICInformation_list))
	{
		vector_mass_cent <- c(vector_mass_cent, EICInformation_list[[i]]$mass_cent)	
	}
	
	II <- which(vector_mass_cent > current_mass)
	
	if (length(II) > 0)
	{
		return(II[1])
	}
	else
	{
		return((length(EICInformation_list)+1))
	}
	
	
#################################################################	
# Code below was commented out by Xiuxia on 5/16/2012 to speed up the process
#	CentMass_Row <-1
#	position <- 1
#	List_Num <- length(EICInformation_list)
#	for (position_index in 1: List_Num)
#	{
#		PreviousNode <- EICInformation_list[[position_index]]
#		if (PreviousNode[1, CentMass_Row]< EICData_Node[1,1])
#			position <- (position+1)
#	}
#	position
#################################################################	
}

# After adding a new EICDataNode into the mass trace, then need to update the corresponding Node's EICInformation 
Update_EICInformation_Node <- function(EICInformation_Node, EICData_Node)
{
	mass_cent_old        <- EICInformation_Node$mass_cent         #EICInformation_Node[1,1]
	mass_max_old         <- EICInformation_Node$mass_max          #EICInformation_Node[1,2]
	mass_min_old         <- EICInformation_Node$mass_min          #EICInformation_Node[1,3]
	point_num_old        <- EICInformation_Node$point_num         #EICInformation_Node[1,4]
	hit_miss_count_old   <- EICInformation_Node$hit_miss_count    #EICInformation_Node[1,5]
	newnode_mass_old     <- EICInformation_Node$newnode_mass      #EICInformation_Node[1,6]
	newnode_intensity_old<- EICInformation_Node$newnode_intensity #EICInformation_Node[1,7]
	newnode_scanindex_old<- EICInformation_Node$newnode_scanindex #EICInformation_Node[1,8]
	
	mass_new      <- EICData_Node$mass       #EICData_Node[1,1]
	intensity_new <- EICData_Node$intensity  #EICData_Node[1,2]
	scan_index_new<- EICData_Node$scan_index #EICData_Node[1,3]
	
	mass_cent_update <- (mass_cent_old*point_num_old+ mass_new)/(point_num_old+1)

	if (mass_max_old<mass_new)
	{
		mass_max_update<- mass_new
	}else
	{
		mass_max_update<- mass_max_old
	}
		
	
	if (mass_min_old>mass_new)
	{
		mass_min_update <- mass_new
	}else
	{
		mass_min_update <- mass_min_old
	}
		
	
	point_num_update <- point_num_old+1
	hit_miss_count_update<- -1 # after ++, can be reverted to 0. 
	
	newnode_mass_update <- mass_new
	newnode_intensity_update <-intensity_new
	newnode_scanindex_update <-scan_index_new
	
	Update_Node <- data.frame(mass_cent=mass_cent_update, mass_max=mass_max_update, mass_min=mass_min_update, point_num=point_num_update, hit_miss_count=hit_miss_count_update, newnode_mass=newnode_mass_update, newnode_intensity=newnode_intensity_update, newnode_scanindex=newnode_scanindex_update)
#	Update_Node <- data.frame(mass_cent=c(mass_cent_update), mass_max=c(mass_max_update), mass_min=c(mass_min_update), point_num=c(point_num_update), hit_miss_count=c(hit_miss_count_update), newnode_mass=c(newnode_mass_update), newnode_intensity=c(newnode_intensity_update), newnode_scanindex=c(newnode_scanindex_update))
	
	return(Update_Node)
}

#realize the EIC extraction based on ROI methods, mainly consider the all of mass trace mass centroid and mass dynamic range, the eludiance distance 
# among the latest adding node (mass, intensity, scan) and the considering node(mass, intensity, scan), If we can get the final deisotope data, each node 
# represents as a 4D vector ( mass, intensity, scan, desiotope fit). we only need a samll modification to suit with the ecliduance distance. 
#--------------------------------------------------------------------------------------------------------------------------------------------------------

#	EIC_Extracted_mzdata <- Extract_EIC(mass_values_mzdata, intensity_values_mzdata, scan_index_mzdata, scan_number_mzdata, PARAMS$Scan_Start, PARAMS$Scan_End, PARAMS$ROI_ppm_Th, PARAMS$MinimumEIC_Count, PARAMS$Hit_Miss_Count, PARAMS$Cutoff_Intensity, PARAMS$Top_number, PARAMS$Quantile_Percentnumber, PARAMS$Cutoff_Mode)
Extract_EIC <- function(fileName_mzdata, PARAMS)
#Extract_EIC<-function(mass_values, intensity_values, scan_index, scan_number, Scan_start, Scan_end, ROI_ppm_Th, point_num_Th, hit_miss_count_Th, cutoff_intensity, Top_number, Quantile_Percentnumber, cutoff_mode)
{
	cat("In Extract_EIC \n")
	
	xs<-xcmsRaw(filename=fileName_mzdata, includeMSn=FALSE)    # 

#	Start_Scan <- max(Scan_start, min(scan_number))
#	End_Scan <-   min(Scan_end,(max(scan_number)-1))
	
	#library(doBy) # Ensure we can use meas2<- orderBy(~-wt, data=meas2)
	EICData_list = vector("list",0)	#create a null list for storing each EICData
	EICInformation_list =vector("list", 0) #create a null list for saving each EICInformation
	
#	Start_MassTrace <- 0	# Commented out by Xiuxia on 5/16/2012 since this appears to be unnecessary
	
	for (i in PARAMS$Scan_Start:PARAMS$Scan_End)
#	for(i in Start_Scan:End_Scan)
	{
		# get the mass-spectrum data for the current scan 
		# Use a database to store the raw spectral data later on to speed up the query
		
		if (i%%1 == 0) 
		{
			cat("current scan is", i, "\r")
		}
		
		browser()
		
		current_scan <- getScan(xs, scan=i)
		mass_values_tempscan1 <- current_scan[,1]
		intensity_values_tempscan1 <- current_scan[,2]
		
		
#		ind_start<- scan_index[i]+1
#		ind_end <- scan_index[i+1]
		# mass_values_tempscan <- mass_values[ind_start:ind_end]
		# intensity_values_tempscan <- intensity_values[ind_start:ind_end]	 
		
#		mass_values_tempscan1 <- mass_values[ind_start:ind_end]
#		intensity_values_tempscan1 <- intensity_values[ind_start:ind_end]	  
		
		# Add some codes for cutoff threshold determination 
#		cutoff_th<- cutoff_intensity
#		if(cutoff_mode==2)
#		{
#			intensity_values_tempscan1_sort<- sort(intensity_values_tempscan1, decreasing=TRUE)
#			cutoff_th <- intensity_values_tempscan1_sort[length(intensity_values_tempscan1_sort)]
#			if(length(intensity_values_tempscan1_sort)>Top_number) cutoff_th<- intensity_values_tempscan1_sort[Top_number]
#		}
#		
#		if(cutoff_mode==3)
#		{
#			if((Quantile_Percentnumber>=0)&(Quantile_Percentnumber<=1)) cutoff_th <- quantile(intensity_values_tempscan1, prob=Quantile_Percentnumber)
#			if(Quantile_Percentnumber<0) cutoff_th <- quantile(intensity_values_tempscan1, prob=0)
#			if(Quantile_Percentnumber>1) cutoff_th <- quantile(intensity_values_tempscan1, prob=1)						
#		} 
		# Add over for the cutoff threshold calculation 
		
#		ind <- which(intensity_values_tempscan1 >= cutoff_th) 
		ind <- which(intensity_values_tempscan1 >= PARAMS$Cutoff_Intensity)
		
		mass_values_tempscan <- mass_values_tempscan1[ind]
		intensity_values_tempscan <- intensity_values_tempscan1[ind]	 	  
		
		# need to sort the current scan according to mass asscending sequence 
        if(length(mass_values_tempscan)>0) # Only sort when the data point number is more than 1
		{
			II <- order(mass_values_tempscan)
			mass_values_tempscan <- mass_values_tempscan[II]
			intensity_values_tempscan <- intensity_values_tempscan[II]
			
			if (i==PARAMS$Scan_Start)
#			if((i==PARAMS$Scan_Start)|(Start_MassTrace==0)) #the first scan, Create EICData_List and EICInformation_list 
			{			 

#				Start_MassTrace <-1    # You can start a valid mass trace
				
				for (index in 1:length(mass_values_tempscan))
				{
					EICData_list[[index]] <- data.frame(mass=mass_values_tempscan[index], intensity=intensity_values_tempscan[index], scan_index=i)
#					EICData_Node<- data.frame(mass=c( mass_values_tempscan[index]), intensity=c(intensity_values_tempscan[index]), scan_index=c(i))
#					EICData_list[[index]] <- EICData_Node
					
					EICInformation_list[[index]] <- data.frame(mass_cent=mass_values_tempscan[index], mass_max=mass_values_tempscan[index], mass_min=mass_values_tempscan[index], point_num=1, hit_miss_count=-1, newnode_mass=mass_values_tempscan[index], newnode_intensity=intensity_values_tempscan[index], newnode_scanindex=i)
#					EICInformation_Node<- data.frame(mass_cent=c( mass_values_tempscan[index]), mass_max=c( mass_values_tempscan[index]), mass_min=c( mass_values_tempscan[index]), point_num=c(1), hit_miss_count=c(0), newnode_mass=c( mass_values_tempscan[index]), newnode_intensity=c(intensity_values_tempscan[index]), newnode_scanindex=c(i))
#					EICInformation_list[[index]] <-EICInformation_Node
			
# Xiuxia Note: Some of the initial masses could be less than the specified mass tolerance			
				}
							
			}#over initialization of ROI(EICData_list and EICInformation_list)  
			else  #Append feature point into EICData_list and EICInformation_list 
			{
				
				for (index in 1:length(mass_values_tempscan))
				{							
					current_mass <- mass_values_tempscan[index]
					current_intensity <- intensity_values_tempscan[index]
					
					EICData_Node<- data.frame(mass=current_mass, intensity=current_intensity, scan_index=i)
#					EICData_Node<- data.frame(mass=c(mass_values_tempscan[index]), intensity=c(intensity_values_tempscan[index]), scan_index=c(i))
					
#					temp_mass <- mass_values_tempscan[index]   # Commented out by Xiuxia
#					delta_mass_Th <- ((ROI_ppm_Th*temp_mass)/1000000)	# Commented out by Xiuxia
					
					#trying to find out the matched EIC_MassTrace to select the corresponding EICData_Node 			 
					MatchedTrace_Index <- Find_optimal_matched_EIC_Masstrace(EICData_Node, EICInformation_list, PARAMS$ROI_ppm_Th) 	
#					MatchedTrace_Index <- Find_optimal_matched_EIC_Masstrace(EICData_Node, EICInformation_list, delta_mass_Th)	# Commented out by Xiuxia on 5/15/2012		 														 			 			 			
					
					if (is.null(MatchedTrace_Index)) # Need to insert a new EIC trace
#					if (MatchedTrace_Index <=0) # Commented out by Xiuxia on 5/15/2012 # Can't find a matched mass trace, need to insert a new mass trace in the exisiting mass list at its correponding mass position 
					{
						position <- Find_position_addingnode(EICData_Node, EICInformation_list) # find the suitable position in the mass trace list for the new mass trace, mainly according to the mass value order  
						EICData_list        <- add_list_node(EICData_list, position, EICData_Node)   
						
						EICInformation_Node <- data.frame(mass_cent=current_mass, mass_max=current_mass, mass_min=current_mass, point_num=1, hit_miss_count=0, newnode_mass=current_mass, newnode_intensity=current_intensity, newnode_scanindex=i) 
#						EICInformation_Node <- data.frame(mass_cent=c( mass_values_tempscan[index]), mass_max=c( mass_values_tempscan[index]), mass_min=c( mass_values_tempscan[index]), point_num=c(1), hit_miss_count=c(0), newnode_mass=c( mass_values_tempscan[index]), newnode_intensity=c(intensity_values_tempscan[index]), newnode_scanindex=c(i)) 
						
						EICInformation_list <- add_list_node(EICInformation_list, position, EICInformation_Node)
						
						##################################################################################
# 						Code below was commented out by Xiuxia on 5/16/2012 since this appears to be unnecessary
#						if(length(EICData_list)!= length(EICInformation_list))
#							stop("EICDatalist and EICInformation_list are not equal in length after adding new mass trace!")				
						##################################################################################
						
					} 
					else #find a matched mass trace, then add the feature point into the corresponding EICData Mass trace and update the corresponding EICInfomration Node
					{
						EICData_list[[MatchedTrace_Index]] <- rbind(EICData_list[[MatchedTrace_Index]], EICData_Node)
						# by row binding, add a EICData_Node into the corresponding mass trace.   
						EICInformation_list[[MatchedTrace_Index]]<- Update_EICInformation_Node(EICInformation_list[[MatchedTrace_Index]], EICData_Node)
						# add EICData node into EICMass trace, need to update the EICInformation Node,inclduing recalculating mass parameters for the corresponding EIC Mass trace and add 1 to mass trace count and clear hit_miss_count.                                    	         
					}                           
					
				} # end of for
									 
			}	# end of else # over appending
			

			
#################################################################			
#	Commented out the sorting loop below by Xiuxia to speed up the sorting			
#			for(j in 1:(length(mass_values_tempscan)-1))
#			{
#				for(k in (j+1):length(mass_values_tempscan))
#				{
#					if (mass_values_tempscan[j]> mass_values_tempscan[k])
#					{
#						# swap mz 
#						temp_mass <- mass_values_tempscan[j]
#						mass_values_tempscan[j] <- mass_values_tempscan[k]
#						mass_values_tempscan[k] <- temp_mass
#						# sawp intensity
#						temp_intensity <- intensity_values_tempscan[j]
#						intensity_values_tempscan[j] <- intensity_values_tempscan[k]
#						intensity_values_tempscan[k] <- temp_intensity
#					}	  					  
#				}		  
#			}#over sorting as ascending	 
#################################################################	
		}	# end of one scan
			
		cat("After one scan \n")
		browser()
		
		vector_node_to_remove <- NULL
		
		if (length(EICInformation_list) > 0)
		{
			temp <- 1:length(EICInformation_list)
			
			#for(k in 1:length(EICInformation_list)){
			for(k in temp){			
				EICInformation_list[[k]]$hit_miss_count <- EICInformation_list[[k]]$hit_miss_count + 1
				
				if ((EICInformation_list[[k]]$point_num < PARAMS$MinimumEIC_Count) & (EICInformation_list[[k]]$hit_miss_count > PARAMS$Hit_Miss_Count))
				{
					vector_node_to_remove <- c(vector_node_to_remove, k)
#					EICInformation_list <- remove_list_node(EICInformation_list, k)
#					EICData_list <- remove_list_node(EICData_list, k)
				}
			}
		}	
		
		if (!is.null(vector_node_to_remove))
		{
			EICData_list <- EICData_list[-vector_node_to_remove]
			EICInformation_list <- EICInformation_list[-vector_node_to_remove]	
		}
	
		
##################################################################################				
# This section below was commented out by Xiuxia

		# After a loop for a whole scan, need to remove some mass traces from the EIC mass trace list if their length is smaller than a fix threshold, also the length keep unchanged for several scans, 
#		temp_EICData_list <- vector("list",0)	     #create a null list for saving the 
#		temp_EICInformation_list <- vector("list", 0) #create a null list for saving each EICInformation
#		temp_position <- 1	
#		
#		if(length(EICInformation_list)>0)  #Ensure the 
#		{
#			for (k in 1:length(EICInformation_list)) 
#			{
#				
#			
#				temp_EICInformation<- EICInformation_list[[k]]
#				temp_EICInformation$hit_miss_count<- temp_EICInformation$hit_miss_count+1 # temp_EICInformation[1,5], miss count +1 for each scan loop, if the current mass trace find a node, it must be revert 0 by add 1
#		
##				if(((temp_EICInformation$point_num) <point_num_Th) && ((temp_EICInformation$hit_miss_count)>hit_miss_count_Th))
##				{
##					  EICInformation_list<- remove_list_node(EICInformation_list, k) 
##					  EICData_list       <- remove_list_node(EICData_list, k)			  
#					
##					if(length(EICData_list)!= length(EICInformation_list))
##						stop("EICDatalist and EICInformation_list are not equal in length after removing random mass trace!")
##				} 	
##				else
##				{
#					# EICInformation_list[[k]]$hit_miss_count<- temp_EICInformation$hit_miss_count
#			
#			
#					temp_EICData  <- EICData_list[[k]]
#					temp_EICInformation_list[[temp_position]]<- temp_EICInformation		       
#					temp_EICData_list[[temp_position]]       <- temp_EICData
#					temp_position <- temp_position + 1
#					
#				
##					if(length(temp_EICData_list)!= length(temp_EICInformation_list))
##						stop("temp_EICData_list and temp_EICInformation_list are not equal in length after removing random mass trace!")
#			
##				}			  
#			}
#		}
		
#		EICInformation_list<- temp_EICInformation_list
#		EICData_list       <- temp_EICData_list
####################################################################################	
		
	} # end of all scans # over scan loop	
	
	#print out the extracted mass trace (EICs)
	EICData_list
	
}# over Extract_EIC function
# 

# Summing all of EICs(Mass trace) to get the new TIC
TIC_Summing_EIC<-function(EICData_list, scan_max)
{
	TIC_bysummingEIC <- numeric(length=scan_max)
	for (i in 1:length(EICData_list))
	{
		inten_arr<- EICData_list[[i]]$intensity
		scan_arr<-  EICData_list[[i]]$scan_index
		for(j in 1:length(inten_arr))
		{
			scan_ind <- scan_arr[j] 
			TIC_bysummingEIC[scan_ind] <- TIC_bysummingEIC[scan_ind] + inten_arr[j]
		}
	}
	TIC_bysummingEIC
}

