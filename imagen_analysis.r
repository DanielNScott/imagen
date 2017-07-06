# This script runs everything.

source('./preproc/import_all.r')
data <- import_all()

source('./preproc/quality_control.r')
data <- quality_control(data)

source('./preproc/impute_data.r')
data <- impute_data(data)

source('./analysis/analyze_data.r')
analysis <- analyze_data(data)

source('./plotting/plot_results.r')
plot_results(data, analysis)

# The end.