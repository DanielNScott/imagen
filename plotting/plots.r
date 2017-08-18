
# Plot feature densities
#options(repr.plot.width=9, repr.plot.height=4)
density_plot <- function(data, names, title) {
  ggplot( melt(data.frame(data[names]), id.vars = NULL), aes(x = value)) +
  facet_wrap(~variable, scales = "free") + geom_histogram(bins = 30, na.rm = TRUE) +
  ggtitle(title)
}

density_plot(data$raw, data$names$sst   , 'SST Parameter Densities')
density_plot(data$raw, data$names$mid   , 'MID Parameter Densities')
density_plot(data$raw, data$names$CANTAB, 'CANTAB Parameter Densities')

#
d <- melt(data.frame(data$raw[data$names$sst], set = as.factor(rep(c('test','train'), each = 198)) ), id.vars = 'set')
ggplot(d, aes(x = value, fill = set)) +
  facet_wrap(~variable, scales = 'free') +
  geom_histogram(alpha = 0.2, position = 'identity') +
  ggtitle('SST Params by Train / Test Set')

# Plot survey densities
#options(repr.plot.width=9, repr.plot.height=5)
#ggplot( melt(data$raw[data$survey_names_14]), aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram(bins=30, na.rm=TRUE) + ggtitle("Survey Data")

# Display strip plots of multiply imputed data (imp)
#options(repr.plot.width=15, repr.plot.height=15)
#stripplot(imp, pch = 20, cex = 0.75, subset = 1:(198*3))

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


plotdata <- data.frame(x = imputed$ifg..stn, y = imputed$presma..stn, z = imputed$stn..gpi)
ggplot(plotdata, aes(x, y, colour = z)) + scale_color_viridis(option = 'plasma') +
  geom_point() + theme_bw() +
  labs(x = 'IFG to STN', y = 'preSMA to STN', colour = 'STN to GPi', title = 'STN Connectivity') +
  theme(legend.position = c(0.9, 0.75)) +
  geom_abline(slope = -0.399, intercept = -0.385, lty = 2) +
  annotate("text", x = 0.1, y = -0.6, label = "Slope: -0.4 \n R: 0.32")

z <- imputed$stn..gpi
npts <- 2
colors <- viridis(npts)
zcolor <- colors[(z - min(z))/diff(range(z))*npts + 1]
plotdata <- data.frame(x = imputed$ifg..stn, y = imputed$presma..stn)
ggplot(plotdata, aes(x, y, colour = zcolor)) + geom_point() + theme_bw() +
  labs(x = 'IFG to STN', y = 'preSMA to STN', colour = 'STN to GPi', title = 'STN Connectivity') +     theme(legend.position = c(0.85, 0.8)) +
  scale_color_identity('STN to GPi', labels = c('-0.64 < val < 0.03', ' 0.03 < val < 0.58'), breaks = colors, guide = "legend")



