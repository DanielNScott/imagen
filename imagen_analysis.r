# This script runs everything.

source('./preproc/import_all.r')
source('./preproc/quality_control.r')
source('./preproc/impute_data.r')
source('./analysis/analyze_data.r')
source('./plotting/plot_results.r')

data <- import_all()
data <- quality_control(data)
data <- impute_data(data)

analysis <- analyze_data(data)

plot_results(data, analysis)

# The end.
