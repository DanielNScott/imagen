# ------------------------------------------------------------------------------ #
#                      Subroutine for reading the fmri data                      #
# ------------------------------------------------------------------------------ #
read_fmri_data <- function(subj_id, timeseries_dir, taskdata_dir, movement_dir, data_dir) {

  # Set filename prefixes and suffixes
  task_pfx <- 'ss_'
  fmri_sfx <- paste('MaskDumpOutput',sep = '')
  move_pfx <- 'EPI_stop_'

  # Specify ordinal coding for trial outcome
  task_key  <- c("GO_SUCCESS" = "1", "GO_FAILURE" = '2', 'STOP_SUCCESS' = '3', 'STOP_FAILURE' = '4',
                 'GO_TOO_LATE' = '5','GO_WRONG_KEY_RESPONSE' = '6', 'STOP_TOO_EARLY_RESPONSE' = '7')

  # Convert subject id to padded string & read
  subj_id_str <- formatC(subj_id, width = 12, format = 'd', flag = '0')

  # Filename beurocracy...
  fmri_fname <- paste(subj_id_str, '_', fmri_sfx, sep = '')
  task_fname <- paste(task_pfx, subj_id_str, '.csv', sep = '')
  move_fname <- paste(move_pfx, subj_id_str, '.txt', sep = '')

  fmri_full_name <- paste(timeseries_dir, fmri_fname, sep = '')
  task_full_name <- paste(taskdata_dir  , task_fname, sep = '')
  move_full_name <- paste(movement_dir  , move_fname, sep = '')

  # Read the fmri activation series
  if (file.exists(fmri_full_name)) {
    activations <- read.csv(file = fmri_full_name, header = FALSE, sep = '')
    activations <- t(activations[,4:447]) # 1st 3 cols spurious, 444 TRs
    rownames(activations) <- NULL
    activations <- as.matrix(activations)

    print(paste('Successful read of ', fmri_fname))
  } else {
    print(paste('Failure to read ', fmri_fname))
    return
  }

  # Read the SST series
  if (file.exists(task_full_name)) {
    # Only want to read 2 columns, time and outcome
    cols <- c('NULL','NULL','numeric','NULL','NULL',NA    ,'NULL','NULL',
              'NULL','NULL','NULL'   , NA   ,'NULL','NULL','NULL','numeric')
    task_data <- read.csv(file = task_full_name, header = TRUE, sep = '\t', colClasses = cols, skip = 1)
    print(paste('Successful read of ', task_fname))
  } else {
    print(paste('Failure to read ', task_fname))
    return
  }

  # Reading the fmri series
  if (file.exists(move_full_name)) {
    movement <- read.csv(file = move_full_name, header = FALSE, sep = '')
    print(paste('Successful read of ', move_fname))
  } else {
    print(paste('Failure to read ', move_fname))
    return
  }

  # Req. library for revalue()
  library(plyr)

  # Extract the trial times and outcomes
  task_times   <- task_data['Trial.Start.Time..Onset.'][,1]
  task_outcome <- as.character(task_data['Response.Outcome'][,1])

  # Old, numeric coding.
  #task_data['Response.Outcome'][,1] <- revalue(task_data['Response.Outcome'][,1], task_key, warn_missing = FALSE)
  #task_outcome <- as.integer(as.character(task_data['Response.Outcome'][,1]))

  # Get stimulus identities (e.g. 'LeftArrow') and get rid of 'Arrow'
  stim <- as.character(task_data['Stimulus.Presented'][,])
  stim <- sapply(stim, function(x){x <- substr(x, 1, nchar(x) - 5) }, USE.NAMES = FALSE)
  stim <- toupper(stim)

  # If we're looking at just the rIFG (a QC check)...
  if (FALSE) {
    rIFG_voxels <- read.csv(file = '../data/rIFGClusterRows.csv', header = FALSE)
    activations <- activations[, as.logical(rIFG_voxels)]
    activations <- cbind(activations, rowMeans(activations, na.rm = TRUE))
  }

  roi_voxels  <- read.csv(file = paste(data_dir, '/roi_voxels.csv', sep = ''), header = TRUE)
  rois <- colnames(roi_voxels)

  n_voxels_by_roi <- colSums(roi_voxels)
  n_voxels_total  <- sum(n_voxels_by_roi)

  # Whole ROI Averaging
  # n_voxels_total  <- length(rois)
  # colnames(tmp) <- rois


  #tmp <- matrix(NA, dim(activations)[1], n_voxels_total)
  tmp <- matrix(NA, dim(activations)[1], 0)
  for (roi in rois) {
    # Whole roi average
    #tmp[,roi]  <- rowMeans(activations[, as.logical(roi_voxels[[roi]]) ], na.rm = TRUE)

    # Lazy alternative to mask FIX THIS LATER
    tmp <- cbind(tmp, activations[, as.logical(roi_voxels[[roi]]) ])
  }

  # Keeping each roi: This doesn't work because masks aren't mutually exclusive...
  # activations <- activations[, rowSums( roi_voxels)]

  # Indiidual STN voxels...
  # stn <- activations[, as.logical(roi_voxels[['rSTN']]) ]
  # tmp <- cbind(tmp, stn)
  # rois <- append(rois, paste('STN', 1:sum(roi_voxels[['rSTN']]), sep = ''))

  activations <- tmp

  # For plotting out ROI voxels...
  #d <- data.frame(activations[, as.logical(roi_voxels[['rSTN']]) ], TR = 1:444)
  #ggplot(melt(d, id.vars = 'TR'), aes(TR, value)) + geom_point() +
  #  facet_wrap(~ variable) + ggtitle('Raw STN Timeseries')

  # Save the fmri_series data (for raw <-> processed comparison)
  #saveRDS(activations, 'activations.rds')

  # Z-Score the data for each participant
  mean <- colMeans(activations, na.rm = TRUE)
  std  <- apply(activations, 2, sd, na.rm = TRUE)

  activations <- sweep(activations, 2, mean, FUN = "-")
  activations <- sweep(activations, 2, std , FUN = "/")

  # Initial spike removal.
  # - Zero is preferred to NA so that the number of non NA entries is retained.
  # - These will be set to 0 again after being perturbed by corrections.
  spike_inds <- which(abs(activations) >= 3)
  activations[spike_inds] <- 0

  # Drift correction
  # In principle filter this is redundant w/ filter...
  #trs         <- as.vector(seq(1,444)) # used right here...
  #polyfit     <- lm( as.matrix(activations) ~ poly(trs, 2))
  #activations <- activations - polyfit$fitted.values

  # High-Pass Filter
  bf <- signal::butter(2, 1/128, type = "high")
  filter <- function(x) {signal::filtfilt(bf, x)}
  activations <- apply(activations, 2, filter)

  # Re Z-Score and secondary spike removal
  mean <- colMeans(activations, na.rm = TRUE)
  std  <- apply(activations, 2, sd, na.rm = TRUE)

  activations <- sweep(activations, 2, mean, FUN = "-")
  activations <- sweep(activations, 2, std , FUN = "/")

  spike_inds <- which(abs(activations) >= 3)
  activations[spike_inds] <- 0 # See note on spike removal above.

  whole_roi = FALSE
  if (whole_roi) {
    tmp <- matrix(NA, dim(activations)[1], length(rois))
    prev_ind <- 1
    for (roi in 1:length(rois)) {
      # Whole roi average
      cur_ind    <- cumsum(n_voxels_by_roi)[[roi]]
      tmp[,roi]  <- rowMeans(activations[, prev_ind:cur_ind], na.rm = TRUE)
      prev_ind   <- cur_ind + 1
    }
    activations <- tmp

    # Z-Score the data for each participant
    mean <- colMeans(activations, na.rm = TRUE)
    std  <- apply(activations, 2, sd, na.rm = TRUE)

    activations <- sweep(activations, 2, mean, FUN = "-")
    activations <- sweep(activations, 2, std , FUN = "/")
  }

  # Initial spike removal.
  # - Zero is preferred to NA so that the number of non NA entries is retained.
  # - These will be set to 0 again after being perturbed by corrections.

  # Return data
  data <- list(acts = activations, spikes = spike_inds, task_outcome = task_outcome,
               task_times = task_times, task_key = task_key, stim = stim, movement = movement,
               n_voxels_by_roi = n_voxels_by_roi, rois = rois)
  return(data)

  # For looking at ROIs after preprocessing...
  #d <- data.frame(activations[, 8:13], TR = 1:444)
  #ggplot(melt(d, id.vars = 'TR'), aes(TR, value)) + geom_point() +
  #  facet_wrap(~ variable) + ggtitle('Processed STN Voxel Timeseries')
} # EOF
# ------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------ #
#            Fits a robust general linear model to raw fMRI activations          #
# ------------------------------------------------------------------------------ #
fit_fmri_glm <- function(fmri_data, seperate, sig_voxels = NULL, core_par = TRUE) {

  # Required libraries
  library(parallel)
  #library(robust)
  library(robustarima)
  library(fmri)

  # Number of scans (varies by individual) is num of non-na rows in activation matrix.
  n_scans  <- length(fmri_data$acts[!is.na(fmri_data$acts[,1]), 1])
  n_voxels <- dim(fmri_data$acts)[2]

  # In case of ancestral graph analysis, need book keeping vars for
  # additional by-trial regressors and to overwrite some things.
  if (seperate) {

    # Need to modify the task outcome fields to reflect handedness
    hands      <- c('LEFT', 'RIGHT')
    cond_names <- names(fmri_data$task_key)
    cond_names <- apply(expand.grid(cond_names, hands),1,function(x) paste(x,collapse = "_"))
    n_conds    <- length(cond_names)
    cond_nums  <- 1:n_conds

    outcomes   <- paste(fmri_data$task_outcome, fmri_data$stim, sep = '_')

    # Book-keeping
    total_trials  <- length(fmri_data$task_outcome[!(fmri_data$task_outcome == '')])
    conv_regs     <- array(NA, dim = c(n_scans, total_trials))
    n_trials_prev <- 0
    offset        <- 0
    index_list    <- list()
    trial_names   <- c()

  } else {
    # Extract numbers and names for task conditions, and their outcomes
    cond_nums  <- as.integer( unique( fmri_data$task_key) )
    cond_names <- names(fmri_data$task_key)
    n_conds    <- length(cond_nums)
    outcomes   <- fmri_data$task_outcome

    # Initialize the array for events convolved w/ HRF
    conv_regs <- array(dim = c(n_scans, n_conds))

    # FIR variables (but currently the FIR model doesn't work)
    #fir_res <- 1
    #fir_dmat  <- matrix(nrow = n_scans*(1/fir_res), ncol = 0)
  }

  # Some conditions do not occur and will be removed
  bad_cond  <- c()

  # Convolve HRF w/ task indicators to create regressors for each condition
  for (cond in cond_nums) {

    # Want to extract only trial in current condition
    cond_msk <- outcomes %in% cond_names[cond]

    # If this condition doesn't occur, skip it.
    n_trials <- sum(cond_msk)
    if (n_trials == 0) {
        bad_cond <- cbind(bad_cond, cond)
        next
    }

    # For some reason fmri.stimulus can't handle an initial-zero onset time.
    # Only one of thousands of pts, so just 'fix' it.
    if (fmri_data$task_times[cond_msk][1] == 0) {
        fmri_data$task_times[cond_msk][1] <- 500
    }

    # Trial onset times and their durations - 0s for event design
    # Onsets are converted from [ms] to units of [TR]
    onsets    <- fmri_data$task_times[cond_msk]
    onsets    <- as.vector(na.omit(onsets/1000/2.2))
    durations <- double(length(onsets))

    # If we're doing the AG-analysis we need to consider each trial as it's own condition.
    if (seperate) {
      # Function to apply fmri.stimulus to each condition onset
      conv_hrf <- function(x) { fmri.stimulus(scans = n_scans, onsets = x, duration = 0) }

      # In ancestral graph case, each trial is given a seperate regressor
      offset  <- offset + n_trials_prev
      ind_beg <- offset + 1
      ind_end <- offset + n_trials

      conv_regs[,ind_beg:ind_end] <- simplify2array(lapply(onsets, conv_hrf))
      trial_names[ind_beg:ind_end] <- paste(cond_names[cond], 1:n_trials, sep = ':')

      n_trials_prev <- n_trials
      index_list[[cond]]    <- c(ind_beg, ind_end)

    } else {

       # Standard fmri analysis - convolve hrf w/ event times:
       conv_regs[,cond] <- fmri.stimulus(scans = n_scans, onsets = onsets, duration = durations)

       # Extract FIR HRFs
       #fir_dmat <- cbind(fir_dmat, fir_design(onsets = onsets, n_scans = n_scans, res = fir_res))
    }
  }

  # Assign appropriate names to the columns in the conv_regs matrix.
  regressor_names <- ifelse(seperate, list(trial_names), list(cond_names))[[1]]
  colnames(conv_regs) <- regressor_names

  # FIR stuff - Can't seem to get this working reasonably
  #approx_fun <- function(y) {approx(1:444, y, xout = seq(fir_res, 444, fir_res), rule = 2)$y }
  #    For testing:
  #interp <- rowSums(apply(conv_regs[,c(1,3:7)], 2, approx_fun ))
  #fir_hrf_est <- corpcor::pseudoinverse( t(fir_dmat) %*% fir_dmat) %*% t(fir_dmat) %*% interp
  #    For use:
  #fir_hrf_est <- corpcor::pseudoinverse( t(fir_dmat) %*% fir_dmat) %*% t(fir_dmat) %*% fmri_data$acts

  # Throw away any bad conditions
  cond_nums    <- setdiff(cond_nums, bad_cond)
  cond_names   <- cond_names[cond_nums]
  n_good_conds <- length(cond_nums)

  # Build the design matrices for use with the linear models.
  if (seperate) {
    # De-index bad conditions
    index_list <- index_list[cond_nums]

    # Get design matrix but dith the linear drift term
    design_mat <- fmri.design(conv_regs, order = 0)
    design_mat <- design_mat[,1:(dim(design_mat)[2] - 1)]

  } else {
    # Create design matrix from conv_regs, less any bad conditions, with 2nd order drift regressor.
    #design_mat <- fmri.design(conv_regs[, cond_nums], order = 2)
    design_mat <- conv_regs[, cond_nums]
    regressor_names <- regressor_names[cond_nums]

    # Remove the unnecessary intercept column (it's automatically fit in lm() )
    #n_good_conds <- dim(design_mat)[2] - 3
    #design_mat  <- design_mat[, c(1:n_good_conds, (n_good_conds + 2):(n_good_conds + 3)) ]
    #colnames(design_mat) <- c(regressor_names, 'Linear Drift', 'Sq. Drift')
    #regressor_names <- colnames(design_mat)
    #n_regs <- n_good_conds + 2
  }

  # Add the motion parameters to the set of regressors
  design_mat <- cbind(design_mat, as.matrix(fmri_data$movement))
  colnames(design_mat) <- c(regressor_names, 'Motion_x', 'Motion_y', 'Motion_z', 'Motion_pitch', 'Motion_yaw', 'Motion_roll')
  regressor_names <- colnames(design_mat)
  n_regs <- dim(design_mat)[2]

  ctrst_names  <- c('Go-Stop', 'Stop-StopFail')
  cois <- c('GO_SUCCESS', 'STOP_SUCCESS', 'STOP_FAILURE', 'Go-Stop', 'Stop-StopFail')

  ###
  # Find pool for core-level SIMD parallelism.
  # - Using makeClust() here will break this on some clusters w/o socket based core communication
  # - mclapply defaults to fork-based dispatching.
  ###
  cores <- detectCores()
  cores <- ifelse(core_par, min(detectCores(), n_voxels), 4)
  flog.info('Cores found: %d', cores[1])
  flog.info('Using %d of them (for %d rois/voxels).', cores[1], n_voxels)

  ## If we want to test with just a few of the ROIs for time-sake:
  #fmri_data$rois <- c('rPreSMA','rIFG','rSTN')
  #fmri_data$acts <- fmri_data$acts[,c(1:76,730:1231,1231:1237)]
  #fmri_data$rois <- c('rPreSMA', 'rSTN')
  #fmri_data$acts <- fmri_data$acts[,c(1:76,1231:1237)]
  #fmri_data$n_voxels_by_roi <- fmri_data$n_voxels_by_roi[fmri_data$rois]

  design_mat <<- design_mat

  # Definitions for the regressions...
  if (seperate) {
    # For AG, need user defined function,
    coefficients <- matrix(NA, total_trials, n_voxels )
    colnames(coefficients) <- colnames(fmri_data$acts)

    # The function to fit an activation time-series
    fit_acts <- function(x) {
      movement <- as.matrix(fmri_data$movement)
      return( fit_ag_lm(conv_regs, movement, x, index_list) )
    }

  } else {
    # for Standard analysis need robust linear model
    coefficients <- matrix(NA, n_regs, n_voxels )

    # The function to fit an activation time-series
    fit_acts <- function(x) {
      # Standard linear regression, but with AR(1) errors
      # model <- lm(x ~ design_mat)
      # model <- orcutt::cochrane.orcutt(model)

      # Robust linear regression
      # model  <- lmRob(x ~ design_mat)

      # Extract coefficients and p-values
      #coef   <- model$coefficients
      #p_vals <- summary(model)$coefficients[,4]

      ###
      # Robust linear regression with ARIMA errors
      # arima.rob has some bug requiring global assignment of it's args...
      #
      response   <<- as.matrix(x)
      model      <- arima.rob(response ~ design_mat, p = 1)

      # Extract coefficients and p-values
      coef   <- model$regcoef
      p_vals <- summary(model)$reg.coef[,4]
      ###

      # Compute contrast vector and it's p-value, assuming white residuals
      # (which is not quite right, but empirically it is conservative enough to multiply by 2)

      # model$df.residual should be n_obs - rank(design_mat)
      residuals <- cbind(rep.int(1, 444), design_mat) %*% as.matrix(model$regcoef)
      residual_df  <- 430
      residual_var <- t(residuals) %*% residuals / residual_df

      n_reg     <- dim(design_mat)[2]
      ctrst_vec <- double(n_reg)
      ctrst_vec[1:2] <- c(1,-1)
      ctrst_var <- 1/sum(abs(ctrst_vec)) * t(ctrst_vec) %*% solve(t(design_mat) %*% design_mat) %*% ctrst_vec * residual_var
      ctrst     <- ctrst_vec %*% model$regcoef[2:(n_reg+1)]
      ctrst_t   <- ctrst / sqrt(ctrst_var)
      ctrst_p   <- 2*( 2*pt(-abs(ctrst_t), df = residual_df) )
      coef      <- c(coef, ctrst)
      p_vals    <- c(p_vals, ctrst_p)

      n_reg     <- dim(design_mat)[2]
      ctrst_vec <- double(n_reg)
      ctrst_vec[2:3] <- c(1,-1)
      ctrst_var <- 1/sum(abs(ctrst_vec)) * t(ctrst_vec) %*% solve(t(design_mat) %*% design_mat) %*% ctrst_vec * residual_var
      ctrst     <- ctrst_vec %*% model$regcoef[2:(n_reg+1)]
      ctrst_t   <- ctrst / sqrt(ctrst_var)
      ctrst_p   <- 2*( 2*pt(-abs(ctrst_t), df = residual_df) )
      coef      <- c(coef, ctrst)
      p_vals    <- c(p_vals, ctrst_p)

      names(coef)   <- c('Int.', colnames(design_mat), ctrst_names)
      names(p_vals) <- c('Int.', colnames(design_mat), ctrst_names)

      return(list(coef, p_vals))
    }
  }

  browser()
  # Perform regression and time it
  time <- system.time(
    coef_and_ps <- mclapply(data.frame(fmri_data$acts), fit_acts, mc.cores = cores[1], mc.silent = TRUE)
    #coef_and_ps <- lapply(data.frame(fmri_data$acts), fit_acts)
  )

  # Recovery test:
  # design_mat <- design_mat[,1:8]
  # n_fake_vox <- 1
  # noise      <- matrix(rnorm(444* n_fake_vox, mean = 0, sd = 0.1), 444)
  # betas      <- rnorm(8, mean = 2, sd = 2)
  # fake_vox   <- t(apply(design_mat, 1, function(x) {betas * x} ))
  # fake_acts  <- rowSums(fake_vox) + noise
  # model      <- summary(lmRob(fake_acts ~ design_mat))
  # rbind(model$coefficients[2:9,1], betas, model$coefficients[2:9, 4])

  # Report time required for regression
  flog.info('Computation time for voxel betas:')
  print(time)

  # Checkpoint for debugging with long computations
  #saveRDS(coefficients, 'betas2.rds')
  #coef_and_ps <- readRDS('coef_and_ps.rds')

  ### Return coefficients to desired format, using only significant voxels ###
  # Extraction from coef_and_ps, thresholding, and some set-up
  coefficients <- data.frame(lapply(coef_and_ps, function(x){x[[1]]} ))
  p_values     <- data.frame(lapply(coef_and_ps, function(x){x[[2]]} ))

  # if (seperate) { do something else here }
  rownames(coefficients) <- c('Int.', colnames(design_mat), ctrst_names)
  rownames(p_values)     <- c('Int.', colnames(design_mat), ctrst_names)

  n_rois   <- length(fmri_data$rois)
  n_ctrsts <- length(ctrst_names)
  n_cois   <- length(cois)

  roi_beg_end  <- c(0, cumsum(fmri_data$n_voxels_by_roi))
  sig_beg_end  <- matrix(0, n_cois, n_rois + 1)
  mean_betas   <- matrix(0, n_cois, n_rois    )
  beta_dists   <- list()
  pval_dists   <- list()
  relax_thresh <- FALSE

  #save(list = ls(), file = 'fmri_checkpoint.rds')
  #browser()
  whole_rois = FALSE
  re_average = TRUE
  if (!whole_rois) {
  for (j in 1:n_cois) {

    # Get roi-corrected significant voxels
    # Index the set of significant voxels by ROI
    for (i in 1:n_rois) {
      vox_first <- roi_beg_end[[i]] + 1
      vox_last  <- roi_beg_end[[i + 1]]

      sig_voxels <- 0
      factor     <- 0
      cur_n_vox  <- fmri_data$n_voxels_by_roi[[i]]

      # Accept increasing amounts of error
      #while ((sum(sig_voxels) < max(cur_n_vox/100, 5)) & (factor <= 0.05)) {
      #  sig_voxels <- p_values[cois,vox_first:vox_last][j,] < (0.0001/cur_n_vox + factor)
      #  factor <- factor + 0.0005
      #}

      # Std. conservative p-vals
      std_thresh <- 0.05/cur_n_vox
      #max_thresh <-  0.05/cur_n_vox
      sig_voxels <- p_values[cois,vox_first:vox_last][j,] < std_thresh
      #n_sig_vox  <- sum(sig_voxels)
      #min_acc    <- ceiling(cur_n_vox / 50)

      # If too few, relax threshold up to 0.05/cur_n_vox
      #while ((n_sig_vox < min_acc ) & (factor <= max_thresh) & relax_thresh) {
      #  sig_voxels <- p_values[cois,vox_first:vox_last][j,] < (0.0001/cur_n_vox + factor)
      #  factor <- factor + 0.0005/cur_n_vox
      #}

      # Take avg and save
      sig_betas  <- coefficients[cois, colnames(sig_voxels)][j, sig_voxels]
      pval_betas <- p_values[cois, colnames(sig_voxels)][j, sig_voxels]
      mean_betas[j,i] <- mean( as.matrix(sig_betas) )

      beta_dists[[paste(fmri_data$rois[i], cois[j], sep = '-')]] <- as.double(sig_betas)
      pval_dists[[paste(fmri_data$rois[i], cois[j], sep = '-')]] <- as.double(pval_betas)

      if (re_average) {
        # Get average time-series of significant voxels only, then fit.
        sig_voxel_avg <- data.frame(rowMeans(as.matrix(fmri_data$acts[, sig_voxels][,vox_first:vox_last])))
        new_res  <- mclapply(sig_voxel_avg, fit_acts, mc.cores = cores[1], mc.silent = TRUE)
        new_coef <- as.matrix(new_res[[1]][[1]])
        rownames(new_coef) <- rownames(coefficients)

        mean_betas[j,i] <- new_coef[cois,][j]
        sig_betas <- coefficients[cois,][j, sig_voxels][vox_first:vox_last]
      }
    }
  }
  }  else {
    mean_betas <- coefficients[cois,]
  }
  #browser()

  colnames(mean_betas) <- fmri_data$rois
  mean_betas <- data.frame(mean_betas)
  rownames(mean_betas) <- cois

  # Return names...
  #if (!seperate) {
  #}

  # Some diagnostic plots...
  #ggplot(melt(rob_model$fitted.values), aes(1:444,value)) + geom_point(col = 'red')
  #+ geom_point(aes(1:444,melt(fmri_data$acts[,503])), col = 'blue')

  # Save the linear model, list of conds removed,
  lm_list <- list(coef = mean_betas, beta_dists = beta_dists, pval_dists = pval_dists,
                  bad_cond = bad_cond, design = design_mat,
                  sig_voxels = sig_voxels, sig_beg_end = sig_beg_end)
  return(lm_list)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#         This performs the linear model fitting for the ancestral graphs        #
# ------------------------------------------------------------------------------ #
fit_ag_lm <- function(regressors, movement, activations, index_list) {

  # Initialize the output variable:
  # Linear model coefficients for each trial
  trial_betas <- rep(NA, dim(regressors)[2])
  names(trial_betas) <- colnames(regressors)

  # Index into trial_betas for use in loops below
  global_trial_num <- 0

  # Construct the "whole condition" regressors. Also see NOTE below.
  whole_cond_regs <- sapply(index_list, function(x) {apply(regressors[, x[1]:x[2], drop = FALSE],1,sum)})

  conds <- 1:length(index_list)
  for (cond in conds) {

    # If not, get the set of trials in this condition
    beg <- index_list[[cond]][1]
    end <- index_list[[cond]][2]
    trials <- beg:end

    # Columns for the other whole condition regressors
    other_conds <- setdiff(conds, cond)

    # Construct a regressor for each trial and fit lm with everything else
    # as nnuisance vars
    for (trial in trials) {
      global_trial_num <- global_trial_num + 1

      # NOTE: take care when partial_cond_inds is not a sequence...
      # ... drop = FALSE handles case where it is only 1 number
      # ... apply returns a list of zeros when partial_cond_inds is empty
      partial_cond_inds <- setdiff(trials, trial)
      partial_cond_reg  <- apply(regressors[,partial_cond_inds, drop = FALSE],1,sum)

      design <- cbind(regressors[,trial], partial_cond_reg, whole_cond_regs[,other_conds], movement)
      fit    <- lm(activations ~ design)

      # Save the coefficients
      trial_betas[global_trial_num] <- fit$coefficients[2]
    }
  }
  return(trial_betas)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#                          A finite impulse response model                       #
# ------------------------------------------------------------------------------ #
fir <- function(onsets, durations, activations, n_scans, n_conds, n_voxels, TR = 1, ...) {

  # Scale parameter - sets resolution
  scale <- 1

  # Rescale onsets
  onsets    <- onsets    * scale
  durations <- durations * scale

  n_scans <- n_scans * scale
  res     <- TR / scale

  #if (type == "user")
  #   shrf <- sum(hrf(0:(ceiling(scans) - 1)/scale))
  #   no  <- length(onsets)
  #if (length(durations) == 1) {
  #  durations <- rep(durations, no)
  #}
  #else if (length(durations) != no) {
  #  stop("Length of duration vector does not match the number of onsets!")
  #}
  stims <- rep(0, ceiling(n_scans))

  ## ESTIMATE HRF USING FIR BASIS SET

  # CREATE FIR DESIGN MATRIX
  # WE ASSUME HRF IS 16 TRS LONG
  hrf_len <- 16

  # BASIS SET FOR EACH CONDITOIN IS A TRAIN OF INPULSES
  fir_bases <- matrix(0, n_scans, hrf_len*n_conds)

  for (cond in 1:n_conds) {
    #
    fir_bases[round(onsets) , 1]  <- 1

    for (i in 2:hrf_len) {
      # i is the number of the 'time-point' in the hrf.
      fir_bases[, i] <- c(rep(0, i - 1), fir_bases[1:(444 - (i - 1)), 1])
    }
  }

  # ESTIMATE HRF FOR EACH CONDITION AND VOXEL
  fir_hrf_est <- corpcor::pseudoinverse( t(fir_bases) %*% fir_bases) %*% t(fir_bases) %*% activations

  # RESHAPE HRFS
  hHatFIR <- reshape(fir_hrf_est, hrf_len, n_conds, n_voxels)

  return(hHatFIR)

  # d <- melt(data.frame(fir_hrf_est[,1:7], TR = 1:16), id.vars = 'TR')
  # ggplot(d, aes(TR, value)) + facet_wrap(~variable) + geom_point()

}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#                     Creates a design matrix for an FIR analysis                #
# ------------------------------------------------------------------------------ #
fir_design <- function(onsets, n_scans, res = 0.1, hrf_len = 16) {

  scale      <- 1/res
  fir_design <- matrix(0, n_scans*scale, hrf_len)

  fir_design[round(onsets*scale) , 1]  <- 1
  for (i in 2:hrf_len) {
    # i is the number of the 'time-point' in the hrf.
    zero_pad <- rep(0, (i - 1)*scale)
    fir_design[, i] <- c(zero_pad, fir_design[1:(n_scans*scale - (i - 1)*scale), 1])
  }

  return(fir_design)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#                     Old means of getting ROI beta value                        #
# ------------------------------------------------------------------------------ #
function(){
  for (j in 1:n_cois) {
    sig_voxels <- p_values[cois, ][j, ] < (0.05/fmri_data$n_voxels_by_roi[j])

    # Index the set of significant voxels by ROI
    for (i in 1:n_rois) {
      vox_first <- roi_beg_end[i] + 1
      vox_last  <- roi_beg_end[i+1]
      sig_beg_end[j, i + 1] <- sum(sig_voxels[vox_first:vox_last])
    }

    # Go from # of significant per roi to a set of indexes
    sig_beg_end[j, ] <- cumsum(sig_beg_end[j,])

    # Take means over significant voxels to use as ROI activation series
    for (i in 1:n_rois) {
      if (sig_beg_end[j,i] != sig_beg_end[j,i + 1]) {
         vox_first <- sig_beg_end[j,i] + 1
         vox_last  <- sig_beg_end[j,i+1]

         sig_voxel_avg <- data.frame(rowMeans(as.matrix(fmri_data$acts[, sig_voxels][,vox_first:vox_last])))
         new_res  <- mclapply(sig_voxel_avg, fit_acts, mc.cores = cores[1], mc.silent = TRUE)
         new_coef <- as.matrix(new_res[[1]][[1]])
         rownames(new_coef) <- rownames(coefficients)

         mean_betas[j,i] <- new_coef[cois,][j]
         sig_betas <- coefficients[cois,][j, sig_voxels][vox_first:vox_last]

         dists[[paste(fmri_data$rois[i], cois[j], sep='-')]] <- as.double(sig_betas)
      }
    }
  }
}
# ------------------------------------------------------------------------------ #
