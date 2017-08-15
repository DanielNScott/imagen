# This script runs everything.

source('./preproc/import_all.r')
source('./preproc/quality_control.r')
source('./preproc/impute_data.r')
source('./analysis/analyze_data.r')
source('./plotting/plot_results.r')

data <- import_all()
data <- quality_control(data)
data <- impute_data(data)

ag_results <- ag_analysis()
ag_subjs   <- readRDS('data/ag_subjs.rds')

data$raw$Subject <- read.csv('data/sandbox_subject_list.csv', sep = '')
stop_betas <- data.frame(cbind(ag_subjs, ag_results$betas$connectivity_ST))
colnames(stop_betas)[1] <- 'Subject'
imputed    <- merge(imputed, stop_betas, by = 'Subject', all = FALSE)

analysis <- analyze_data(data)

plot_results(data, analysis)

# The end.
