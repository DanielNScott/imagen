
# This is ok for now, but some information about the standard deviatins should also be encoded.
# I'm uncertain what at this point, but distributional overlap seems important.
imp_ind <- (data$train['mu_go_14'] - data$train['mu_stop_14']) / data$train['sigma_go_14']
colnames(imp_ind) <- c('imp_ind')

options(repr.plot.width=4, repr.plot.height=2)
dmelt  <- melt(imp_ind, id.vars = NULL)

probs  <- c(0.25, 0.5, 0.75)
quants <- quantile(imp_ind, prob=probs, na.rm=TRUE)

ggplot(dmelt, aes(x = value)) +
   facet_wrap(~variable, scales = "free") +
   geom_line(aes(y = ..density.., colour = 'Empirical'), stat = 'density') +  
   geom_histogram(aes(y = ..density..), alpha = 0.4) +                        
   scale_colour_manual(name = 'Density', values = c('red', 'blue')) + 
   theme(legend.position = c(0.85, 0.85)) +
   geom_vline(data=dmelt, aes(xintercept=quants[1],), linetype="dashed", size=0.5) +
   geom_vline(data=dmelt, aes(xintercept=quants[2],), linetype="dashed", size=0.5) +
   geom_vline(data=dmelt, aes(xintercept=quants[3],), linetype="dashed", size=0.5)

dmelt  <- melt(data$raw['stn_std'], id.vars = NULL)

probs  <- c(0.25, 0.5, 0.75)
stn_quants <- quantile(data$raw['stn_std'], prob=probs, na.rm=TRUE)

options(repr.plot.width=4, repr.plot.height=2)
ggplot(dmelt, aes(x = value)) +
   facet_wrap(~variable, scales = "free") +
   geom_line(aes(y = ..density.., colour = 'Empirical'), stat = 'density', na.rm=TRUE) +  
   geom_histogram(aes(y = ..density..), alpha = 0.4, bins=30, na.rm=TRUE) +                        
   scale_colour_manual(name = 'Density', values = c('red', 'blue')) + 
   theme(legend.position = c(0.85, 0.85)) +
   geom_vline(data=dmelt, aes(xintercept=stn_quants[1],), linetype="dashed", size=0.5) +
   geom_vline(data=dmelt, aes(xintercept=stn_quants[2],), linetype="dashed", size=0.5) +
   geom_vline(data=dmelt, aes(xintercept=stn_quants[3],), linetype="dashed", size=0.5) +
   ggtitle("STN Activation StdDevs & Quartiles")