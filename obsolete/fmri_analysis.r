shuffle <- FALSE
save    <- FALSE
FU      <- FALSE

use_cached <- FALSE

source('import_fmri.r')

roi    <- 'Left_preSMA'
fuc    <- 'BL'
subset <- 'train'

data   <- import_fmri(subset,fuc,roi)

stn_series  <- data$stn_series
sst_times   <- data$sst_times
sst_outcome <- data$sst_outcome

#data <- import_stn_series('FU2')

#stn_series_fu2  <- data$stn_series
#sst_times_fu2   <- data$sst_times
#sst_outcome_fu2 <- data$sst_outcome

#subs_lost <- data$subs_lost

if (shuffle){
    stn_series[,1] <- sample(stn_series[1:dim(stn_series)[1],1], replace=FALSE)
}



library(ggplot2)
library(reshape2)

options(repr.plot.width=8, repr.plot.height=2)
dmelt <- melt(stn_series[c('TR', 'Subj_1', 'Subj_198')], id.vars = 'TR', variable.name = 'series')
ggplot(dmelt, aes(TR, value)) + geom_line(aes(colour=series))

options(repr.plot.width=8, repr.plot.height=2)
dmelt <- melt(stn_series[c('TR', 'Subj_2', 'Subj_102')], id.vars = 'TR', variable.name = 'series')
ggplot(dmelt, aes(TR, value)) + geom_line(aes(colour=series))

options(repr.plot.width=8, repr.plot.height=2)
dmelt <- melt(stn_series[c('TR', 'Subj_3', 'Subj_58')], id.vars = 'TR', variable.name = 'series')
ggplot(dmelt, aes(TR, value)) + geom_line(aes(colour=series))

options(repr.plot.width=8, repr.plot.height=2)
dmelt <- melt(stn_series[c('TR', 'Subj_4', 'Subj_138')], id.vars = 'TR', variable.name = 'series')
ggplot(dmelt, aes(TR, value)) + geom_line(aes(colour=series))


labels <- c("GO_SUCCESS", "GO_FAILURE", 'STOP_SUCCESS', 'STOP_FAILURE', 'GO_TOO_LATE', 'GO_WRONG_KEY_RESPONSE', 'STOP_TOO_EARLY_RESPONSE')

dmelt <- melt(data.frame('Outcome' = sst_outcome[,1], 'Time' = sst_times[,1]/1000/2.2), id.vars = 'Time', variable.name = 'series', na.rm=TRUE)
ggplot(dmelt, aes(Time, value)) + geom_point(colour='red') + scale_y_continuous('outcome', breaks=seq(1,7,1), labels=labels, c(1,7))

dmelt <- melt(data.frame('Outcome' = sst_outcome[,194], 'Time' = sst_times[,194]/1000/2.2), id.vars = 'Time', variable.name = 'series', na.rm=TRUE)
ggplot(dmelt, aes(Time, value)) + geom_point(colour='dark cyan') + scale_y_continuous('outcome', breaks=seq(1,7,1), labels=labels, c(1,7))





## Fit GLM to STN Activity ##
#
#For each participant:
#- Remove non-existent conditions from outcome set
#- Determine onset times (in TRs)
#- Set event durations (event-based --> 0s)
#- Create expected activation timeseries via HRF - regressor convolution
#- Create design matrix with second order polynomial drift
#- Fit the model to the data

source('get_stn_betas.r')

# --------------------------------------------#
#           THIS CAN TAKE A WHILE!!           #
# --------------------------------------------#
n_fit <- 198
if (use_cached) {
    lm_list <- readRDS(paste('lm_list_198_',fuc,'_',roi,'_',subset,'.rds',sep=''))
} else {
    #clumping <- list(c(1,2,5,6), c(3,4))
    clumping <- list(c(1), c(2), c(3), c(4), c(5), c(6), c(7))
    lm_list  <- get_stn_betas(stn_series, sst_times, sst_outcome, n_fit, clumping, labels)
}

if (save) {
    saveRDS(lm_list,paste('lm_list_198_',fuc,'_',roi,'_',subset,'.rds',sep=''))
}

###
#
#

if (use_cached) {
    contrasts <- readRDS(contrasts, paste('contrasts_198_',fuc,'_',roi,'_',subset,'.rds',sep=''))
    betas <- readRDS(betas, paste('betas_198_',fuc,'_',roi,'_',subset,'.rds',sep=''))  
} else {
    n_cols <- 4
    betas     <- array(dim=c(n_fit,n_cols))
    contrasts <- array(dim=c(n_fit,1))

    cols <- c('design_matGO_SUCCESS','design_matSTOP_SUCCESS','design_matSTOP_FAILURE','design_matGO_TOO_LATE')#,'design_matGO_WRONG_KEY_RESPONSE','design_matSTOP_TOO_EARLY_RESPONSE')
    cols <- cols[1:n_cols]
    print(cols)
    for (subj_num in 1:n_fit) {
        betas[subj_num,1:n_cols] <- lm_list[[subj_num]]$lm$coefficients[cols]

        # Contrast go conds with stop conds ignoring 'failures'
        contr <- as.vector(c(-1, 1, 0, 0))
        contrasts[subj_num,1] <- contr %*% as.vector(betas[subj_num,c(1:n_cols)])
    }
    betas <- data.frame(betas)
    colnames(betas) <- c('GO_SUCCESS','STOP_SUCCESS','STOP_FAILURE','GO_TOO_LATE') #,'GO_WRONG_KEY_RESPONSE','STOP_TOO_EARLY_RESPONSE')
}

if (save) {
    saveRDS(contrasts, paste('contrasts_198_',fuc,'_',roi,'_',subset,'.rds',sep=''))
    saveRDS(betas, paste('betas_198_',fuc,'_',roi,'_',subset,'.rds',sep=''))
}

options(repr.plot.width=8, repr.plot.height=2)
ggplot(melt(data.frame(contrasts=contrasts), id.vars=NULL), aes(x=value)) + geom_histogram(bins=50)

options(repr.plot.width=8, repr.plot.height=3)
#ggplot(melt(data.frame(betas_go=betas[,1], betas_stop=betas[,2], betas_stop_fail=betas[,3], betas_go_late=betas[,4]), id.vars=NULL), aes(x=value)) + facet_wrap(~variable, scales = "free") + geom_histogram(bins=50)
ggplot(melt(data.frame(betas_go=betas[,1], betas_stop=betas[,2]), id.vars=NULL), aes(x=value)) + facet_wrap(~variable, scales = "free") + geom_histogram(bins=50)

#t.test(x = betas[,2], y = betas[,1], alternative='greater')
t.test(x = contrasts, alternative='greater')


options(repr.plot.width=6, repr.plot.height=3)

dmelt <- melt(betas[c('STOP_SUCCESS', 'STOP_FAILURE')], id.vars = 'STOP_SUCCESS', na.rm=TRUE)
ggplot(dmelt, aes(STOP_SUCCESS, value)) + geom_point(colour='dark cyan') + scale_y_continuous('Stop Failure Beta') + scale_x_continuous('Stop Success Beta') + theme(aspect.ratio = 0.7)

dmelt <- melt(betas[c('GO_SUCCESS', 'STOP_SUCCESS')], id.vars = 'GO_SUCCESS', na.rm=TRUE)
ggplot(dmelt, aes(GO_SUCCESS, value)) + geom_point(colour='dark cyan') + scale_x_continuous('GO Success Beta') + scale_y_continuous('Stop Success Beta') + theme(aspect.ratio = 0.7)

dmelt <- melt(betas[c('GO_SUCCESS', 'STOP_FAILURE')], id.vars = 'GO_SUCCESS', na.rm=TRUE)
ggplot(dmelt, aes(GO_SUCCESS, value)) + geom_point(colour='dark cyan') + scale_x_continuous('Go Success Beta') + scale_y_continuous('Stop Failure Beta') + theme(aspect.ratio = 0.7)


###

conv_reg.m <- melt(lm_list[[1]]$conv_regs, id.vars = 'TR')

options(repr.plot.width=6, repr.plot.height=5)

p <- ggplot(conv_reg.m, aes(variable, TR)) +
      geom_tile(aes(fill = value)) + 
      scale_fill_gradient(low = "white", high = "steelblue") + 
      labs(x = "") + scale_x_discrete(expand = c(0, 0)) + #scale_y_discrete(expand = c(0, 0)) +
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              legend.position = "none",
              axis.ticks = element_blank(),
              axis.text.x = element_text(angle = 330, hjust = 0, colour = "grey50"))

print(p)

##### DESIGN MATRIX w/ DRIFTS
dmatdf <- data.frame(design_mat)
colnames(dmatdf) <- c(labels[setdiff(1:7, bad_cond)], 'Intercept', 'Linear Drift', 'Sq. Drift')
dmatdf['TR'] <- 1:444

dmatdf.m <- melt(dmatdf, id.vars = 'TR')
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





###
###
###

source('fir.r') 
library('corpcor')
#lm_list[[1]]$onsets

onsets    <- lm_list[[1]]$onsets[[1]]
durations <- integer(length(onsets))

activations <- stn_series[,1]
n_scans <- 444
n_conds <- 1
n_voxels <- 1
TR <- 1
    
# Scale parameter - sets resolution
scale <- 10

# Rescale onsets
onsets    <- onsets    * scale
durations <- durations * scale

n_scans <- n_scans * scale
res     <- TR / scale
stims   <- rep(0, ceiling(n_scans))

## ESTIMATE HRF USING FIR BASIS SET

# CREATE FIR DESIGN MATRIX
# WE ASSUME HRF IS 16 TRS LONG
hrf_len <- 8
 
# BASIS SET FOR EACH CONDITOIN IS A TRAIN OF INPULSES
fir_bases <- matrix(0, n_scans, hrf_len*n_conds)
 
for (cond in 1:n_conds) {
    col_subset <- ((cond - 1)* hrf_len + 1):(cond*hrf_len)
    
    #for (onset in 1:length(onsets[[cond]]) ) {
    for (onset in 1:length(onsets) ) {
        #impulse_times <- onsets(onset):onsets(onset) + hrf_len - 1;
        impulse_times <- seq(onsets[onset], onsets[onset] + hrf_len*scale - 1, scale)
        
        for (impulse in 1:length(impulse_times)) {
            fir_bases[impulse_times[impulse], col_subset[impulse]] <- 1;
        }
    }
}

options(repr.plot.width=6, repr.plot.height=2)

acts <- matrix(0, length(activations)*scale, 1)
acts[1,1] <- activations[1]
for (i in 1:(length(activations)-1)) {
    acts[ ((i-1)*(scale)+1):(i*(scale)), 1] <- as.matrix(approx(c(1,11), c(activations[i], activations[i+1]), n=11)$y[1:10])
}

# ESTIMATE HRF FOR EACH CONDITION AND VOXEL
fir_hrf_est <- pseudoinverse( t(fir_bases) %*% fir_bases) %*% (t(fir_bases) %*% acts)

# RESHAPE HRFS
#hHatFIR <- reshape(fir_hrf_est, hrf_len, n_conds, n_voxels)

ggplot(melt(data.frame(TR=1:hrf_len, hrf=fir_hrf_est), id.vars='TR'), aes(TR,value)) + geom_point()



###
###
###
plot_avg_hrf <- function(sst_times, subj_num, cond, roi_series, look_ahead, color) {

    # Masks for getting trial times
    big_cond_msk <- apply(sst_outcome, c(1,2), function(x) x %in% cond)
    cond_msk     <- big_cond_msk[,subj_num]

    # If this condition doesn't occur, skip it.
    if (sum(cond_msk) == 0) {
        bad_cond <- cbind(bad_cond,cond)
        next
    }

    # For some reason fmri.stimulus can't handle the zero.
    # Only one of thousands of pts, so just 'fix' it.
    if (sst_times[cond_msk,subj_num][1] == 0) {
        sst_times[cond_msk,subj_num][1] <- 500
    }

    # Trial onset times and their durations - 0s for event design
    onsets    <- sst_times[cond_msk, subj_num]
    onsets    <- as.vector(na.omit(onsets/1000/2.2))
    durations <- double(length(onsets))

    #print(onsets)
        
    # Get convolved regressor for condition
    #conv_regs[,cond] <- fmri.stimulus(scans = n_scans, onsets = onsets, duration = durations)

    n_ons <- length(onsets)
    blah  <- array(dim=c(n_ons,look_ahead+1))
    for (i in 1:n_ons) {
        blah[i,1] <- ceiling(onsets[i]) - onsets[i]
        for (j in 1:look_ahead) {
            blah[i,j + 1] <- roi_series[ceiling(onsets[i]) + j - 1]
        }
    }

    blah_2 <- array(dim=c(n_ons*look_ahead,2))
    j <- 0
    for (i in 1:n_ons){
        blah_2[j + 1:look_ahead, 1] <- blah[i, 1] + 0:(look_ahead-1)
        blah_2[j + 1:look_ahead, 2] <- blah[i, 2:(look_ahead+1)]
        j <- j + look_ahead
    }

    blah_3 <- data.frame('Times' = blah_2[,1], 'Activation' = blah_2[,2])
    plt <- ggplot(melt(blah_3, id.vars = 'Times'), aes(Times,value)) + geom_point(colour=color)
    print(plt)
    
    blah_3_sorted <- blah_3[order(blah_3['Times']),]

    library(signal)
    bf <- butter(2, 1/800, type="low")
    b1 <- filtfilt(bf, blah_3_sorted[['Activation']])

    df <- data.frame('Times' = blah_3_sorted[,1], 'Activation' = b1)
    plt <- ggplot(melt(df, id.vars = 'Times'), aes(Times,value)) + geom_point(colour=color)
    print(plt)
    
}

subj <- 10
cond <- 1
look_ahead <- 9
                          
options(repr.plot.width=7, repr.plot.height=2)

plot_avg_hrf(sst_times, subj, cond, lm_list[[subj]]$conv_regs[,cond],look_ahead, 'red')
plot_avg_hrf(sst_times, subj,cond, stn_series[,subj],look_ahead, 'dark cyan')

len <- length(stn_series[,subj])
rand <- sample(stn_series[,subj], size=len, replace=FALSE)
plot_avg_hrf(sst_times, subj,cond, rand,look_ahead, 'dark green')

