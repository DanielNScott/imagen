
# Plot feature densities
options(repr.plot.width=9, repr.plot.height=4)
ggplot( melt(data$raw[data$sst_names], id.vars=NULL), aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram(bins=30, na.rm=TRUE) + ggtitle("Recovered SST Parameter Densities")
ggplot( melt(data$raw[data$mid_names], id.vars=NULL), aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram(bins=30, na.rm=TRUE) + ggtitle("Recovered MID Parameter Densities")
ggplot( melt(data$raw[data$cgt_names], id.vars=NULL), aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram(bins=30, na.rm=TRUE) + ggtitle("CGT Data")
ggplot( melt(data$raw[data$agn_names], id.vars=NULL), aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram(bins=30, na.rm=TRUE) + ggtitle("Affective Go/No-Go Data")
ggplot( melt(data$raw[data$msc_names], id.vars=NULL), aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram(bins=30, na.rm=TRUE) + ggtitle("Other Task Data")

# Plot survey densities
options(repr.plot.width=9, repr.plot.height=5)
ggplot( melt(data$raw[data$survey_names_14]), aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram(bins=30, na.rm=TRUE) + ggtitle("Survey Data")

#

# Display strip plots of multiply imputed data (imp)
options(repr.plot.width=15, repr.plot.height=15)
stripplot(imp, pch = 20, cex = 0.75, subset = 1:(198*3))

###
# Display histograms of 'complete' data
###
library(reshape2)
library(ggplot2)

options(repr.plot.width=9, repr.plot.height=5)

#d <- melt(data_14_task.imputed)
#ggplot(d, aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram()

d <- melt(data_14_task.raw)
ggplot(d, aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram()

#d <- melt(data_14_survey.raw)
#ggplot(d, aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram()

d <- melt(data_18_survey.raw)
ggplot(d, aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram()

