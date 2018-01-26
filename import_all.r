# This script runs everything.
setwd('/home/dan/projects/imagen/')

source('./analysis/import_routines.r')
source('./analysis/quality_control.r')
source('./analysis/impute_data.r')
source('./plotting/plot_wrappers.r')
source('./plotting/plot_all.r')
source('./misc/misc_tools.r')

#---------------- Do Everything For Age 14 Data ---------------#
# Import the basic data
data <- import_non_fmri('local')

# Import the standard fmri statistics
fmri_std <- readRDS('data/std_stats.rds')
data$raw <- merge(data$raw, fmri_std, on = 'Subject', all = TRUE)

# Import the fmri ROI covariance stats
fmri_cov <- readRDS('data/covariance_stats.rds')
data$raw <- merge(data$raw, fmri_cov, on = 'Subject', all = TRUE)

# Import Rob's data
rob      <- import_rob_data('./data/rob_data.csv')
data$raw <- merge(data$raw, rob, on = 'Subject')

# Quality control and z-score stuff
data <- gen_addnl_flds(data)
data <- quality_control(data)

#stop("You don't want to process the 18 yo data now, do you?")
#---------------- Do Everything For Age 18 Data ---------------#
# Import the basic data
data_18 <- import_non_fmri('local', time = 'FU2')

# Import the standard fmri statistics
fmri_std_18 <- readRDS('data/fmri_stats_SST_FU2.rds')
data_18$raw <- merge(data_18$raw, fmri_std_18, on = 'Subject', all = TRUE)

# Quality control and z-score stuff
data_18 <- gen_addnl_flds(data_18, skip_tci = TRUE)
data_18 <- quality_control(data_18)

#---------------- Merge Age 14 Data w/ Age 18 Data ---------------#
# Save the un-suffixed names
names   <- data$names

# Suffix data sets properly
data    <- append_to_data_names(data   , '_14')
data_18 <- append_to_data_names(data_18, '_18')

# Combine the data sets
data$raw <- merge(data$raw, data_18$raw, on = 'Subject')

# Replace the 'names' fields to have _14, _18, and un-suffixed
data$names_18 <- data_18$names
data$names_14 <- data$names
data$names    <- names

# Get that memory back!
rm(data_18)

#--------------------------------------------------------------#
#data <- import_fmri('local', data, 1, 396, process_fmri = FALSE)
#ag_results <- ag_analysis()      # Perform ancestral graph analysis
#data <- impute_data(data)
#plots(data)

# The end.
