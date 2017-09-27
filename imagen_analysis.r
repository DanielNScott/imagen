# This script runs everything.

source('./analysis/import_routines.r')
source('./analysis/quality_control.r')
source('./analysis/impute_data.r')
source('./fmri/fmri_routines.r')
source('./plotting/plot_wrappers.r')
source('./plotting/plot_all.r')

# Import the basic data
data <- import_non_fmri('local')

# Import the standard fmri statistics
fmri_std <- readRDS('data/std_stats.rds')
data$raw <- merge(data$raw, fmri_std, on = 'Subject', all = TRUE)

# Import the fmri ROI covariance stats
fmri_cov <- readRDS('data/covariance_stats.rds')
data$raw <- merge(data$raw, fmri_cov, on = 'Subject', all = TRUE)

# Import Rob's data
rob  <- import_rob_data('./data/rob_data.csv')
data$raw <- merge(data$raw, rob, on = 'Subject')

# Quality control and z-score stuff
data <- gen_addnl_flds(data)
data <- quality_control(data)

#data <- import_fmri('local', data, 1, 396, process_fmri = FALSE)
#ag_results <- ag_analysis()      # Perform ancestral graph analysis
#data <- impute_data(data)
#plots(data)

# The end.
