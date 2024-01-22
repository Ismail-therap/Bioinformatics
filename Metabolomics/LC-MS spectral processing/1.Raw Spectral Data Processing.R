# Ref: https://www.metaboanalyst.ca/resources/vignettes/LCMSMS_Raw_Spectral_Processing.html


## LC-MS/MS data processing workflow in MetaboAnalystR involves several steps:

# 1. Raw spectral data import, 
# 2. MS data processing (auto-optimized peak picking, alignment, gap filling and annotation)
# 3. Data-dependent acquisition (DDA)/SWATH-data-independent acquisition (DIA) data deconvolution, 
# 4. Spectrum consensus from replicates, 
# 5. MS/MS reference library searching, results export, and
# integration into functional prediction.

## Loading the packages 

# Load MetaboAnalystR
library(MetaboAnalystR)
# Load OptiLCMS
library(OptiLCMS)


#############################################################
####### 1. Example raw spectral data processing (MS1) #######
#############################################################

# 1.0  Data download and extraction:

# Run this part to download and unzip the data (run only once)
#download.file("https://www.xialab.ca/api/download/metaboanalyst/malaria_r_example.zip",
#              destfile = "malaria_raw.zip",
#              method = "curl")
#unzip("malaria_raw.zip", exdir = "upload")


# 1.1 Region of Interest (ROI) extraction

# Here, we extract ROIs from 3 QC samples.
DataFiles <- list.files("upload/QC/", full.names = TRUE)
mSet <- PerformROIExtraction(datapath = DataFiles, rt.idx = 0.9, rmConts = TRUE) # rt.idx should be > 0.5 and rmCounts = TRUE to remove instrumental noise 



### Printed messages:

# Step 0/6: Scanning ROIs for parameters optimization...
# QC_001.mzML
# QC_003.mzML
# QC_005.mzML
# Data Loading...
# Raw file import begin...
# Import data: QC_001.mzML finished!                                         
#   Import data: QC_003.mzML finished!                                         
#   Import data: QC_005.mzML finished!                                         
#   Data Loaded !
#   Empty Scan detecting...
# No Empty scan found !
#   Identifying regions of interest (ROI)...
# Identifying regions of potential contaminants........Done!
#   70 potential contaminamts will not be used for parameters optimization !
#   Going to the next step...
# MS data Preparing...
# MS Data ready !
#   Identifying ROIs in m/z dimensions...
# Identifying ROIs in m/z dimensions Done !                                       
#   Identifying ROIs in RT dimensions...
# Chromatogram Plotting Begin...
# Identification on ROIs Finished!
#   Warning message:
#   In MulticoreParam(4L) :
#   MulticoreParam() not supported on Windows, use SnowParam()






# 1.2 Auto-optimization of parameters

# Here we use PerformParamsOptimization to optimize parameters based on 
# the extracted ROI (stored in 'mSet') before process the entire dataset


# Recommendations for parallel core setting, 64GB RAM ~ 6 cores; 32GB RAM ~ 4 cores; 16GB RAM ~ 2 cores. Minimum 16GB RAM is required for raw spectral processing.
best_params <- PerformParamsOptimization(mSet, param = SetPeakParam(platform = "UPLC-Q/E"), ncore = 6)

### Printed messages:

# Step 1/6: Start to optimize parameters! 
#   
#   This step may take a long time...
# DoE Optimization Starting Now... (2023-11-02 18:00:07.205273)
# Evaluating Noise level...
# ppm is estimated as : 4.3                                                  
# Done!
#   Preparing Parameters for optimization finished !
#   Data Spliting Finished !
#   Peak Preparing Begin...
# Peak Preparing Done !
#   Round:1
# DoE Running Begin...
# Your OS is Windows, there might be unexpected errors.
# 
# If there is some unexpected bugs, please reduce the 'core' as 1.
# 
# 100% |                                                                     
#   Round 1 Finished !
#   Model Parsing...
# Gaussian peak ratio (%): 39.3
# Model Parsing Done !
#   
#   Round:2
# DoE Running Begin...
# Your OS is Windows, there might be unexpected errors.
# 
# If there is some unexpected bugs, please reduce the 'core' as 1.
# 
# 100% |                                                                     
#   Round 2 Finished !
#   Model Parsing...
# Gaussian peak ratio (%): 42.6
# Model Parsing Done !
#   
#   Round:3
# DoE Running Begin...
# Your OS is Windows, there might be unexpected errors.
# 
# If there is some unexpected bugs, please reduce the 'core' as 1.
# 
# 100% |                                                                     
#   Round 3 Finished !
#   Model Parsing...
# Gaussian peak ratio (%): 34.3
# Model Parsing Done !
#   
#   No Increase Stopping !
#   best parameter settings:
#   min_peakwidth:  5
# max_peakwidth:  15.75
# mzdiff:  0.036
# snthresh:  7.1
# bw:  2
# Peak_method:  centWave
# RT_method:  loess
# ppm:  4.3
# noise:  6772
# prefilter:  2
# value_of_prefilter:  11556.95
# minFraction:  0.8
# minSamples:  1
# maxFeatures:  100
# fitgauss:  FALSE
# mzCenterFun:  wMean
# integrate:  1
# extra:  1
# span:  0.4
# smooth:  loess
# family:  gaussian
# verbose.columns:  FALSE
# polarity:  negative
# perc_fwhm:  0.6
# mz_abs_iso:  0.005
# max_charge:  2
# max_iso:  2
# corr_eic_th:  0.85
# mz_abs_add:  0.001
# rmConts:  TRUE
# BlankSub:  TRUE
# Step 1/6: Parameters Optimization Finished ! (2023-11-02 18:35:47.809023)
# Time Spent In Total:35.7mins




# 1.3 Importing example data

# "path" is used to specify the path to the folder containing the raw MS spectra to be processed.
# BPI and TIC plotting can be enabled with parameter, 
# "plotSettings = SetPlotParam(Plot = T)", or disabled by changing "T" into "F";

mSet <- ImportRawMSData(path = c("upload"), plotSettings = SetPlotParam(Plot = T))

# 1.4 Raw spectral data processing

# "mSet" include complete raw MS spectra to be processed.
# "params" is using the "best_params" generated above
# Plotting functions can be enabled with parameter, 
# "plotSettings = SetPlotParam(Plot = T)", or disabled by changing "T" into "F";
mSet <- PerformPeakProfiling(mSet, Params = best_params, plotSettings = SetPlotParam(Plot=TRUE))

# 1.5 Feature annotation

# We firstly define the parameters for feature annotation

# 'polarity' is required, can be either 'negative' or 'positive';
# 'perc_fwhm' is used to set the percentage of the width of the FWHM for peak grouping. 
#              Default is set to 0.6;
# 'mz_abs_iso' is used to set the allowed variance for the search (for isotope annotation). 
#              The default is set to 0.005;
# 'max_charge' is set the maximum number of the isotope charge. 
#              For example, the default is 2, therefore the max isotope charge is 2+/-;
# 'max_iso' is used to set the maximum number of isotope peaks.
#              For example, the default is 2, therefore the max number of isotopes per peaks is 2;
# 'corr_eic_th' is used to set the threshold for intensity correlations across samples. 
#              Default is set to 0.85.
# 'mz_abs_add' is used to set the allowed variance for the search (for adduct annotation). 
#              Default is set to 0.001.
# 'adducts' is used to specify the adducts based on your instrument settings.

annParams <- SetAnnotationParam(polarity = 'positive', mz_abs_add = 0.015);

# "mSet" include processed raw MS spectra to be processed.
# "annParams" is the parameters used for annotation

mSet <- PerformPeakAnnotation(mSet, annParams)


# 1.6 Feature table generation

# Here we format and filter the peak list for following analysis with MetaboAnalystR

# Parameters are explained as below,
# annParams, is the object created using the SetAnnotationParam function above;
# filtIso, is used to decide to filter out all isotopes (TRUE) or not (FALSE);
# filtAdducts, is used to decide to filter out all adducts (TRUE) or not (FALSE);
# missPercent, specify the threshold to remove features missing in a certain percentage
#              of all samples in a group.

mSet <- FormatPeakList(mSet, annParams, filtIso = FALSE, filtAdducts = FALSE, missPercent = 1)

# Export annotation results, the annotation will be save as "annotated_peaklist.csv";
Export.Annotation(mSet)

# Export complete feature table. It will be saved as "metaboanalyst_input.csv";
# This table can be used for statistic analysis, functional analysis, biomarker analysis module directly.
Export.PeakTable(mSet)

# Export a summary table (peak_result_summary.txt) to summarize the information of all peaks in
# different samples. The columns are sample names, groups, retention time range, m/z range of all peaks,
# number of all peaks and percentage of missing features.
Export.PeakSummary(mSet)



#############################################################
####### 2.  Example of DDA MS/MS data processing ############
#############################################################

# 3.0 Download example DDA MS/MS spectra


# There are 6 raw spectral files (.mzML) and 1 MS feature table included in the downloaded folder
# The MS feature table was generated with the MS1 data processing pipeline (as described in section 2)

download.file("https://www.xialab.ca/api/download/metaboanalyst/ms2_dda_example.zip",
              destfile = "ms2_dda_example.zip",
              method = "curl")
unzip("ms2_dda_example.zip", exdir = "ms2_dda")

# 3.1 Download MS/MS spectra reference database

## Here we are downloading biology MS/MS reference database.

download.file("https://www.xialab.ca/api/download/metaboanalyst/MS2ID_Bio_v09102023.zip",
              destfile = "MS2ID_Bio.zip",
              method = "curl")
unzip("MS2ID_Bio.zip", exdir = "ms2_db")

# 3.2 Importing DDA MS/MS spectra

## We load MetaboAnalystR and OptiLCMS
library(MetaboAnalystR)
library(OptiLCMS)

## Clean the environment first
rm(list = ls())

## Read the MS1 feature table as the target features for MS/MS data processing
## This table include four columns (mzmin, mzmax, rtmin, rtmax)
## mzmin and mzmax is the minimum and the maximum value of m/z for the feature;
## rtmin and rtmax is the minimum and the maximum value of retention time for the feature;
ft_dt <- qs::qread("ms2_dda/ms2_dda_example/ft_dt.qs")

## Here we use function PerformMSnImport to read MS/MS data
## This step may take seconds to minutes depending on the size of your dataset
mSet <- PerformMSnImport(filesPath = c(list.files("ms2_dda/ms2_dda_example/",
                                                  pattern = ".mzML",
                                                  full.names = T, recursive = T)),
                         targetFeatures = ft_dt,
                         acquisitionMode = "DDA")

# 3.3 DDA MS/MS spectra Deconvolution

# In this step, we are processing the DDA spectra deconvolution;

# The deconvolution is based on the Biology MS/MS spectral database;
# Parallel computing is supported, users are encouraged to use multiple cores to speed up;
# This step may take minutes to hours to finish, depending on the size of datasets

system.time(mSet <- PerformDDADeconvolution(mSet,
                                            ppm1 = 5,
                                            ppm2 = 10,
                                            sn = 12,
                                            filtering = 0,
                                            window_size = 1.5,
                                            intensity_thresh = 1.6e5,
                                            database_path = "ms2_db/MS2ID_Bio_v09102023.sqlite",
                                            ncores = 6L))

# 3.4 Spectrum consensus of replicates

mSet <- PerformSpectrumConsenus (mSet,
                                 ppm2 = 15,
                                 concensus_fraction = 0.2,
                                 database_path = "",
                                 use_rt = FALSE,
                                 user_dbCorrection = FALSE)

# 3.5 Reference database searching

# PerformDBSearchingBatch is used to seatching MS/MS reference library
# Results will be scored based on the similarity rules above
# Parallel computing is allowed. CPU cores used are controlled by argument "ncores".

mSet <- PerformDBSearchingBatch (mSet,
                                 ppm1 = 10,
                                 ppm2 = 25,
                                 rt_tol = 5,
                                 database_path = "ms2_db/MS2ID_Bio_v09102023.sqlite",
                                 use_rt = FALSE,
                                 enableNL = FALSE,
                                 ncores = 6L)

# This code is not giving us any ouput...
