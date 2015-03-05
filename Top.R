#Top File for LC-MS project, Written by Wenchao Zhang, 02/14/2012
#From Top File surface, you need to configure:
# Code_Path          : The R code directory path, usually the same to Top file  
# Input_Data_Path    : The LC-MS Data file directory path. You need to copy the files in this path  
# Output_Data_Path   : The Path that storing the output peaklist files from each processing step, also include the final compound_align_list file. 
# Configure_FilePath : The path that need to include 3 files, Adduct_Infor_Config.csv, Isotope_Adduct_Infor_Config.csv, hmdb.csv,   
# PARAMS             : The PARAMS is a data.frame, which need to be configured parameters such as the start_scan, end_scan, mass trace mass tolerance, minimum 
#                       point number of a valid mass trace, miss hit count, intensity_threshold, Maximum conisdering Isotope number, charge state, mass tolerance 
#                       for the mass bias, align window for phase1, align window for phase2, whether adopt two phase alignment. Of course, you can use the default
#                       parameters. These default paramters are used for Xie standard sample data.  
############################################################################################################################################################################

source("http://bioconductor.org/biocLite.R")
#biocLite("CAMERA")
biocLite("xcms")

library(xcms)
#library(CAMERA)

window_linux <- 0 #1: windows system, 0: linux system.
if(window_linux) sep_char <- "\\"  else  sep_char <- "/"


#Define a variable for the Code folder       
if (window_linux)
{
	Code_Path          <- "C:\\Duxiuxia\\LCMS\\LC_Code_testLeco"             # "C:\\lc_R_test\\test" 
	Input_Data_Path    <- "C:\\Duxiuxia\\LCMS\\Input"      # "C:\\lc_R_test\\lc_ms_file"
	Output_Data_Path   <- "C:\\Duxiuxia\\LCMS\\Output"           # "C:\\lc_R_test\\lc_ms_output"
	Configure_FilePath <- "C:\\Duxiuxia\\LCMS\\Configure"  # "C:\\lc_R_test\\lc_ms_configure"
}else{
	Code_Path          <- "/home/xdu/testLECO/LCMS"             # "C:\\lc_R_test\\test" 
	Input_Data_Path    <- "/home/xdu/testLECO/Input"      # "C:\\lc_R_test\\lc_ms_file"
	Output_Data_Path   <- "/home/xdu/testLECO/Output"           # "C:\\lc_R_test\\lc_ms_output"
	Configure_FilePath <- "/home/xdu/testLECO/Configure"  # "C:\\lc_R_test\\lc_ms_configure"
}


FileFormat_mzdata_cdf <- 1  # 1 Mean the sample data are as format of.mzdata.xml file, 0 means the sample data are as format of .CDF file 
if(FileFormat_mzdata_cdf) file_extend <- ".mzXML" else file_extend<- ".CDF"

source(paste(Code_Path,"Misc_Process.R", sep=sep_char)) 
source(paste(Code_Path,"EIC_MassTrace.R", sep=sep_char)) 
source(paste(Code_Path,"EIC_Binning.R", sep=sep_char)) 
source(paste(Code_Path,"WavCWT_PeakPicking.R", sep=sep_char)) 
source(paste(Code_Path,"BaseLine_Remove.R", sep=sep_char))
source(paste(Code_Path,"Peak_Annotation.R", sep=sep_char)) 
source(paste(Code_Path,"Serial_Par_Process.R", sep=sep_char))  
source(paste(Code_Path,"Compound_Alignment.R", sep=sep_char)) 

# The following code segments are related to global processing 
graphics.off()

if (window_linux==0) {
	library(playwith)
	library(cairoDevice)
	library(RGtk2)
}

#browser()

# wavelet
library(waveslim)
library(wmtsa)

library(ncdf)

#Define some Folder for inputting LC-MS files and the internal files that from peak detecting modules, peak group modules, peak annotation modules, peak refinement modules, and the final alignment modules

AdductFileName        <- paste(Configure_FilePath, "Adduct_Infor_Config.csv", sep=sep_char)
Isotope_AdductFileName<- paste(Configure_FilePath, "Isotope_Adduct_Infor_Config.csv", sep=sep_char) 
hmdb_library_FileName <- paste(Configure_FilePath, "hmdb.csv", sep=sep_char)

# Define the default parameters for LC Project
PARAMS <- data.frame(Scan_Start=50, Scan_End=1200, ROI_ppm_Th=40, MinimumEIC_Count=10, Hit_Miss_Count=3, Cutoff_Intensity=500, Top_number=500, Quantile_Percentnumber=0.1, Cutoff_Mode=1, AdductFileName="", Isotope_AdductFileName="", hmdb_library_FileName="",Max_Iso=0.011, Max_CS=2, mass_tol=0.1, Group_Scan_Shift=3, SNR_Th=10.0, PeakWidth_Low =7, PeakWidth_High =60, Group_Cluster_Angle_Th=30, Refinement_Cluster_Angle_Th=15, Align_Window_Phase1= 0.2, Align_Window_Phase2= 0.2,Two_Phase=1, Cluster_num=2)  

# The following parameters are related to the module of Mass trace based EIC_Extraction 
PARAMS$Scan_Start<-200                        # unit(scan)
PARAMS$Scan_End  <-4000                      # unit(scan)
PARAMS$ROI_ppm_Th<-20
PARAMS$MinimumEIC_Count<-10                   # unit(scan)
PARAMS$Hit_Miss_Count<-3                      # unit(scan)
PARAMS$Cutoff_Intensity<-50                  # Cutoff Intensity, Only the data point with intensity>500 of a scan are considered 
PARAMS$Top_number<- 200                       # Top 500 data points of a scan are considered 
PARAMS$Quantile_Percentnumber<-0.1            # Quantile point, 10 percent of the lowest data point of each scan are not considered.
PARAMS$Cutoff_Mode <-1                        # Cutoff Mode, 1: CutOff directly by specific value Cutoff_Intensity, 2: Only the Top_Number data point for each spectrum are kept; 3: The data point of the specific percenttage of a spectrum are not considered.   

# The following parameters are related to the module of Peak filtering
PARAMS$SNR_Th                       <- 10.0
PARAMS$PeakWidth_Low                <- 7      # unit(scan)
PARAMS$PeakWidth_High               <- 160     # unit(scan)

# The following parameters are related to the modules of Peak group, Annotation, and refinement 
PARAMS$Group_Scan_Shift             <- 3      # unit(scan)                
PARAMS$Group_Cluster_Angle_Th       <- 30     # unit(degree)   

# The following parameters are related to the modules of Peak Annotation, and Refinement 
PARAMS$AdductFileName <- AdductFileName
PARAMS$Isotope_AdductFileName <- Isotope_AdductFileName
PARAMS$hmdb_library_FileName <- hmdb_library_FileName
PARAMS$Max_Iso<-2
PARAMS$Max_CS <-1
PARAMS$mass_tol<-0.011  # Da
PARAMS$Refinement_Cluster_Angle_Th  <- 15  # unit(degree) 

# The following parameters are related to the module of Compound Alignment
PARAMS$Align_Window_Phase1<- 0.3   # unit(minutes)
PARAMS$Align_Window_Phase2<-0.2    # unit(minutes)
PARAMS$Two_Phase <- 1              # 1-> Two phase alignment, 0-> One phase alignment


PARAMS$Cluster_num<- 1

# Load two packages for MPI, Only work in Linux 
if(window_linux == 0 & PARAMS$Cluster_num > 1)
{
	library(Rmpi)
	library(snow)	
}

# ======================================================
# read the LC-MS data files as .cdf file or mzdata.xml file
# ======================================================
FileName_List <-list.files(path=Input_Data_Path, pattern = file_extend)  #,full.names=TRUE   ,

cat(paste("FileName_List is", FileName_List, sep=""))
 
if((window_linux==1)||(PARAMS$Cluster_num<2))  #Windows system or linux system but Cluster node Number <2, use serial process to run the LC-MS data 
{
	cat("Running in Serial model", "\r")
	cat("Before Run Serial Process, The current time is:")
	Sys.time()
	T1<- Sys.time()

	cat(FileFormat_mzdata_cdf)

	for(File_Index in 1:length(FileName_List))
	{
		if(FileFormat_mzdata_cdf) Serial_Process_mzdata_File(FileName_List[File_Index], Input_Data_Path, Output_Data_Path, PARAMS, sep_char)
		else Serial_Process_cdf_File(FileName_List[File_Index], Input_Data_Path, Output_Data_Path, PARAMS, sep_char)
	}
	cat("After Run Serial Process, The current time is:")
	Sys.time()
	T2<-Sys.time()
	TimeDif <- T2 -T1
	TimeDif
}else{
	#cat(paste("Before Run MPI, The current time is:", format(Sys.time(), "%H:%M:%S> "),sep=""))
	cat("Before Run Par Process, The current time is:")
	Sys.time()
	T1<- Sys.time()
	cl<-makeCluster(PARAMS$Cluster_num,type="MPI")
	if(FileFormat_mzdata_cdf)
	{
		clusterExport(cl,"Par_Process_mzdata_File")
		clusterApply(cl,FileName_List,Par_Process_mzdata_File, Code_Path, Input_Data_Path, Output_Data_Path, PARAMS, sep_char)#parallel version
	}else{
		clusterExport(cl,"Par_Process_cdf_File")
		clusterApply(cl,FileName_List,Par_Process_cdf_File, Code_Path, Input_Data_Path, Output_Data_Path, PARAMS, sep_char)#parallel version
	}	
	stopCluster(cl)
    #cat(paste("After Run MPI, The current time is:", format(Sys.time(), "%H:%M:%S> "),sep=""))
	cat("After Run Par Process, The current time is:")
	Sys.time()
	T2<-Sys.time()
	TimeDif <- T2 -T1
	TimeDif
}

#After finished the Parallel processing for EIC_Extraction, Peak Picking, and Peak filter, Peak Group, Peak Annotation, Group Refinement(By MPI). Go to Compound Alignment processing.
cat("Before Run Alignment, The current time is:")
Sys.time()
T3<-Sys.time()
Compound_List <-AlignmentFile_Load(Output_Data_Path)
Compound_Align_FileName <- paste(Output_Data_Path, "Compound_List_Align.csv", sep=sep_char)
Compound_Alignment(Compound_List, PARAMS$Align_Window_Phase1, PARAMS$Align_Window_Phase2, PARAMS$mass_tol, PARAMS$Two_Phase, Compound_Align_FileName)
cat("After run Alignment, The current time is:")
Sys.time()
T4<-Sys.time()
TimeDif <- T4 -T3
TimeDif
