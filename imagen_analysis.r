# This script runs everything.

source('./analysis/import_all.r')
source('./analysis/quality_control.r')
source('./analysis/fmri_routines.r')
source('./analysis/impute_data.r')
source('./analysis/cca_wrapper.r')
source('./plots.r')

data <- import_all('local')
data <- gen_addnl_flds(data)
#ag_results <- ag_analysis()      # Perform ancestral graph analysis
data <- quality_control(data)
#data <- impute_data(data)
#plots(data)

# The end.
