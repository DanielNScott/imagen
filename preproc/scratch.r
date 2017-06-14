options(repr.plot.width=8, repr.plot.height=3)

nplts <- 3
for (subj_num in 1:nplts) {
    dmelt <- melt(cbind(lm_list[[subj_num]]$conv_regs[1:40,c(1,3,8)], stn_series[1:40,subj_num]/12.5), id.vars = 'TR', variable.name = 'Condition')
    plt <- ggplot(dmelt, aes(TR, value)) + geom_line(aes(colour=Condition))
    print(plt)
}

##### DESIGN MATRIX w/ DRIFTS
#dmatdf <- data.frame(design_mat)
#colnames(dmatdf) <- c(labels[setdiff(1:7, bad_cond)], 'Intercept', 'Linear Drift', 'Sq. Drift')
#dmatdf['TR'] <- 1:444

#dmatdf.m <- melt(dmatdf, id.vars = 'TR')
#dmatdf.m <- ddply(dmatdf.m, .(variable), transform, rescale = scale(value))

#options(repr.plot.width=6, repr.plot.height=5)

#p <- ggplot(dmatdf.m, aes(variable, TR)) +
#      geom_tile(aes(fill = value)) + 
#      scale_fill_gradient(low = "white", high = "steelblue") + 
#      labs(x = "") + scale_x_discrete(expand = c(0, 0)) + #scale_y_discrete(expand = c(0, 0)) +
#      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#              panel.background = element_blank(),
#              legend.position = "none",
#              axis.ticks = element_blank(),
#              axis.text.x = element_text(size = base_size * 0.8, angle = 330, hjust = 0, colour = "grey50"))

#print(p)

options(repr.plot.width=6, repr.plot.height=3)

dmelt <- melt(betas[c('STOP_SUCCESS', 'STOP_FAILURE')], id.vars = 'STOP_SUCCESS', na.rm=TRUE)
ggplot(dmelt, aes(STOP_SUCCESS, value)) + geom_point(colour='dark cyan') + scale_y_continuous('Stop Failure Beta') + scale_x_continuous('Stop Success Beta') + theme(aspect.ratio = 1)

dmelt <- melt(betas[c('GO_SUCCESS', 'STOP_SUCCESS')], id.vars = 'GO_SUCCESS', na.rm=TRUE)
ggplot(dmelt, aes(GO_SUCCESS, value)) + geom_point(colour='dark cyan') + scale_x_continuous('GO Success Beta') + scale_y_continuous('Stop Success Beta') + theme(aspect.ratio = 1)

dmelt <- melt(betas[c('GO_SUCCESS', 'STOP_FAILURE')], id.vars = 'GO_SUCCESS', na.rm=TRUE)
ggplot(dmelt, aes(GO_SUCCESS, value)) + geom_point(colour='dark cyan') + scale_x_continuous('Go Success Beta') + scale_y_continuous('Stop Failure Beta') + theme(aspect.ratio = 1)

#ggplot(melt(data.frame(contrast=contrasts), id.vars=NULL), aes(x=value)) + geom_histogram(bins=50)

if (FU) {
    betas_BLFU <- array(dim=c(198,8))
    contrasts_BLFU <- array(dim=c(198,2))

    cols <- c('design_matGO_SUCCESS','design_matSTOP_SUCCESS','design_matSTOP_FAILURE','design_matGO_TOO_LATE')
    fus_lost <- 0
    for (subj_num in 1:50) {
        if (subj_num %in% subs_lost) {
            fus_lost <- fus_lost + 1
            next
        }
        betas_BLFU[subj_num,1:4] <- lm_list[[subj_num]]$lm$coefficients[cols]
        betas_BLFU[subj_num,5:8] <- lm_list_fu2[[subj_num-fus_lost]]$lm$coefficients[cols]

        # Contrast go conds with stop conds ignoring 'failures'
        #contr <- as.vector(c(-1, 0.5, 0.5))
        #contrasts[subj_num,1] <- contr %*% as.vector(betas[subj_num,1:3])
    }
    betas_BLFU <- data.frame(betas_BLFU)
    colnames(betas_BLFU) <- c('GO_SUCCESS','STOP_SUCCESS','STOP_FAILURE','GO_TOO_LATE','GO_SUCCESS_FU2','STOP_SUCCESS_FU2','STOP_FAILURE_FU2','GO_TOO_LATE_FU2')
}

### 

dmelt <- melt(betas_BLFU[c('STOP_SUCCESS', 'STOP_SUCCESS_FU2')], id.vars = 'STOP_SUCCESS', na.rm=TRUE)
ggplot(dmelt, aes(STOP_SUCCESS, value)) + geom_point(colour='dark cyan') + scale_y_continuous('Stop Success Beta FU2') + scale_x_continuous('Stop Success Beta') + theme(aspect.ratio = 1)

dmelt <- melt(betas_BLFU[c('GO_SUCCESS', 'GO_SUCCESS_FU2')], id.vars = 'GO_SUCCESS', na.rm=TRUE)
ggplot(dmelt, aes(GO_SUCCESS, value)) + geom_point(colour='dark cyan') + scale_x_continuous('GO Success Beta') + scale_y_continuous('Go Success Beta FU2') + theme(aspect.ratio = 1)

dmelt <- melt(betas_BLFU[c('STOP_FAILURE', 'STOP_FAILURE_FU2')], id.vars = 'STOP_FAILURE', na.rm=TRUE)
ggplot(dmelt, aes(STOP_FAILURE, value)) + geom_point(colour='dark cyan') + scale_x_continuous('Stop Failure Beta') + scale_y_continuous('Stop Failure Beta FU2') + theme(aspect.ratio = 1)








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