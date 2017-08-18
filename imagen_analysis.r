# This script runs everything.

source('./preproc/import_all.r')
source('./preproc/quality_control.r')
source('./preproc/impute_data.r')
source('./analysis/analyze_data.r')
source('./plotting/plot_results.r')

# Import everything. This calls some preprocessing steps as well such as the fmri analysis
data <- import_all()

# Construct additional features from the raw set
data <- gen_addnl_flds(data)

# Perform ancestral graph analysis
ag_results <- ag_analysis()

#
data <- attach_fmri_results(data)

#
data <- quality_control(data)

#
data <- impute_data(data)

#
analysis <- analyze_data(data)

plot_results(data, analysis)

# The end.
