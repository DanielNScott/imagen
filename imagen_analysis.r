# This script runs everything.

source('./analysis/import_all.r')
source('./analysis/quality_control.r')
source('./analysis/fmri_routines.r')
source('./analysis/impute_data.r')
source('./plotting/plot_all.r')
#source('./plotting/plot_results.r')
data <- import_all('local')
data <- gen_addnl_flds(data)
#ag_results <- ag_analysis()      # Perform ancestral graph analysis
data <- attach_fmri_results(data) #
data <- quality_control(data)
#data <- impute_data(data)
#analysis <- analyze_data(data)
#plot_results(data, analysis)

# The end.
