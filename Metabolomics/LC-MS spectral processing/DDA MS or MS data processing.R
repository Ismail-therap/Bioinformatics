
#############################################################
####### 2.  Example of DDA MS/MS data processing ############
#############################################################

# In this case study, we are implementing a small DDA dataset from a COVID-19 metabolomics study.
# DDA acquires MS/MS spectra by fragmentation of precursor ions selected using a relatively narrow MS/MS isolation window (e.g., 1 m/z).



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
# Stopped here because of error!
# https://omicsforum.ca/t/error-during-reference-database-searching-dda-ms-ms-data-processing/2967
