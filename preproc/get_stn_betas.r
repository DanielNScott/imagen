get_stn_betas <- function(stn_series, sst_times, sst_outcome, n_to_fit, contrast, cnames) {

# fMRI analysis functions
library(fmri)
#library(MASS)

# Proportion of subjects w/ significant STN betas for STOP_SUCCESS
frac_go_sig   <- 0
frac_stop_sig <- 0

# Output variable in which everything is saved
lm_list  <- vector("list", n_to_fit) 

# Loop through subjects
for (subj_num in 1:n_to_fit) {

    # Number of scans in this series -- varies by individual.
    n_scans <- length(stn_series[!is.na(stn_series[,1]),1])

    # Array for events convolved w/ HRF
    conv_regs <- array(dim = c(n_scans, length(contrast)))
    
    # Some conditions do not occur - they will be removed
    bad_cond  <- c()

    # List of onset times by condition
    onset_list <- list(c(), c(), c(), c(), c(), c(), c())
    
    # Fill the conv_regressors array for each condition
    for (cond in 1:length(contrast)) {
        
        # Masks for getting trial times
        big_cond_msk <- apply(sst_outcome, c(1,2), function(x) x %in% contrast[[cond]])
        cond_msk     <- big_cond_msk[,subj_num]
        #print(sum(cond_msk) == 0)

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
        conv_regs[,cond] <- fmri.stimulus(scans = n_scans, onsets = onsets, duration = durations)

        #print(conv_regs)
        onset_list[[cond]] <- onsets
    }
    #print(bad_cond)
    #print(conv_regs)
    
    # Create design matrix w/ 2nd deg drift, remove any bad conditions
    design_mat <- fmri.design(conv_regs[, setdiff(1:7, bad_cond)], order = 2)
    colnames(design_mat) <- c(cnames[setdiff(1:7, bad_cond)], 'Intercept', 'Linear Drift', 'Sq. Drift')

    #design_mat <- fmri.design(conv_regs, order = 2)
    #colnames(design_mat) <- c(cnames, 'Intercept', 'Linear Drift', 'Sq. Drift')

    #print(design_mat)
    #design_mat <- design_mat[,setdiff(1:7, bad_cond)]
    #print(design_mat)

    # Clean up conv_regs a bit, append 'TR' column, label things
    conv_regs[,bad_cond] <- integer(n_scans)
    conv_regs  <- data.frame(cbind(conv_regs, 1:n_scans))
    colnames(conv_regs) <- c(cnames,'TR')

    # Fit the linear model
    linear_model <- lm(stn_series[1:n_scans,1] ~ design_mat)

    # Save the linear model, list of conds removed, 
    lm_list[[subj_num]] <- list(lm=linear_model, bad_cond=bad_cond, conv_regs=conv_regs, onsets = onset_list)

    # ------------------------- #
    # Actual processing is done #
    # ------------------------- #
    # Now stuff is just getting printed
    lm_desc <- summary(linear_model)
    coeff   <- lm_desc$coefficients
    
    # STOP_SUCCESS is coeff. 4 originally (cond 3)
    stop_ind <- 4
    if (2 %in% bad_cond) {
        stop_ind <- stop_ind - 1
    }

    # Information
    print(paste('Subject ', subj_num,
                ' GO_SUCCESS beta: ', round(coeff[2,'Estimate'], digits=2),
                '. Sig: ', coeff[2,'Pr(>|t|)'] < 0.05,
                ' STOP_SUCCESS beta: ', round(coeff[stop_ind,'Estimate'], digits=2),
                '. Sig: ', coeff[stop_ind,'Pr(>|t|)'] < 0.05, sep=''))
    
    if (coeff[2,'Pr(>|t|)'] < 0.05) {
        frac_go_sig <- frac_go_sig + 1
    }
    if (coeff[stop_ind,'Pr(>|t|)'] < 0.05) {
        frac_stop_sig <- frac_stop_sig + 1
    }
}

# Report on betas:
frac_go_sig   <- frac_go_sig  /n_to_fit
frac_stop_sig <- frac_stop_sig/n_to_fit

print('')
print(paste('Fraction of GO_SUCCESS   betas which are significant:', frac_go_sig))
print(paste('Fraction of STOP_SUCCESS betas which are significant:', frac_stop_sig))

return (lm_list)
}