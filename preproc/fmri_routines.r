# ------------------------------------------------------------------------------ #
#                      Subroutine for reading the fmri data                      #
# ------------------------------------------------------------------------------ #
import_fmri <- function(subj_id, timeseries_dir, taskdata_dir, movement_dir, data_dir) {

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
    # tmp[,roi]    <- rowMeans(activations[, as.logical(roi_voxels[[roi]]) ], na.rm = TRUE)

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
fit_fmri_glm <- function(fmri_data, seperate) {

  # Required libraries
  library(parallel)
  library(robust)
  library(fmri)

  # Number of scans in this series -- varies by individual.
  n_scans <- length(fmri_data$acts[!is.na(fmri_data$acts[,1]),1])

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
    # Task conditions
    cond_nums  <- as.integer(unique(fmri_data$task_key))
    cond_names <- names(fmri_data$task_key)
    n_conds    <- length(cond_nums)
    outcomes   <- fmri_data$task_outcome

    # Array for events convolved w/ HRF
    conv_regs <- array(dim = c(n_scans, n_conds))

    # FIR
    #fir_res <- 1
    #fir_dmat  <- matrix(nrow = n_scans*(1/fir_res), ncol = 0)
  }

  # Some conditions do not occur - they will be removed
  bad_cond  <- c()

  # Fill the conv_regressors array for each condition
  for (cond in cond_nums) {

    # Masks for getting trial times
    cond_msk <- outcomes %in% cond_names[cond]

    # If this condition doesn't occur, skip it.
    n_trials <- sum(cond_msk)
    if (n_trials == 0) {
        bad_cond <- cbind(bad_cond, cond)
        next
    }

    # For some reason fmri.stimulus can't handle the zero.
    # Only one of thousands of pts, so just 'fix' it.
    if (fmri_data$task_times[cond_msk][1] == 0) {
        fmri_data$task_times[cond_msk][1] <- 500
    }

    # Trial onset times and their durations - 0s for event design
    # Onsets are converted from [ms] to units of [TR]
    onsets    <- fmri_data$task_times[cond_msk]
    onsets    <- as.vector(na.omit(onsets/1000/2.2))
    durations <- double(length(onsets))

    # Get convolved regressor for condition if 'seperate' is flagged
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

  # Assign appropriate names: Need to list-ify and un-listify for if-else.
  regressor_names <- ifelse(seperate, list(trial_names), list(cond_names))[[1]]
  colnames(conv_regs) <- regressor_names

  # FIR stuff - Can't seem to get this working reasonably
  #approx_fun <- function(y) {approx(1:444, y, xout = seq(fir_res, 444, fir_res), rule = 2)$y }
  #    For testing:
  #interp <- rowSums(apply(conv_regs[,c(1,3:7)], 2, approx_fun ))
  #fir_hrf_est <- corpcor::pseudoinverse( t(fir_dmat) %*% fir_dmat) %*% t(fir_dmat) %*% interp
  #    For use:
  #fir_hrf_est <- corpcor::pseudoinverse( t(fir_dmat) %*% fir_dmat) %*% t(fir_dmat) %*% fmri_data$acts

  # Throw away bad conditions
  cond_nums  <- setdiff(cond_nums, bad_cond)
  cond_names <- cond_names[cond_nums]
  if (seperate) {index_list <- index_list[cond_nums]}

  if (seperate) {
    # Get design matrix but dith the linear drift term
    design_mat <- fmri.design(conv_regs, order = 0)
    design_mat <- design_mat[,1:(dim(design_mat)[2] - 1)]

  } else {
    # Create design matrix by adding 2nd deg drift and remove any bad conditions from conv_regs
    design_mat <- fmri.design(conv_regs[, cond_nums], order = 2)
    regressor_names <- regressor_names[cond_nums]

    # Remove the unnecessary intercept column which is redundant intercept est. in lm()
    n_good_conds <- dim(design_mat)[2] - 3
    design_mat  <- design_mat[, c(1:n_good_conds, (n_good_conds + 2):(n_good_conds + 3)) ]
    colnames(design_mat) <- c(regressor_names, 'Linear Drift', 'Sq. Drift')
    regressor_names <- colnames(design_mat)
  }

  # Add the motion parameters to the set of regressors
  design_mat <- cbind(design_mat, as.matrix(fmri_data$movement))
  colnames(design_mat) <- c(regressor_names, 'Motion_x', 'Motion_y', 'Motion_z', 'Motion_pitch', 'Motion_yaw', 'Motion_roll')
  regressor_names <- colnames(design_mat)

  # Setup for regressions...
  n_voxels     <- dim(fmri_data$acts)[2]
  n_regressors <- dim(design_mat)[2]

  # Indicate pool for core-level SIMD parallelism:
  # Note: Using makeClust() here will break this on some clusters which do not let R
  #       use socket based core communication. mclapply defaults to fork-based dispatching.
  cores <- detectCores()
  flog.info('Cores found: %d', cores[1])
  cores <- min(detectCores(), n_voxels)
  flog.info('Using %d of them (for %d rois/voxels).', cores[1], n_voxels)

  # For AG, need user defined function,
  # for std analysis need robust linear model
  if (seperate) {
    coefficients <- matrix(NA, total_trials, n_voxels )
    colnames(coefficients) <- colnames(fmri_data$acts)

    fit_cols <- function(x) {
      return(fit_ag_lm(conv_regs, as.matrix(fmri_data$movement), x, index_list))
    }

  } else {
    coefficients <- matrix(NA, n_regressors, n_voxels )

    fit_cols <- function(x) {
      model  <- lmRob(x ~ design_mat)
      coef   <- model$coefficients[2:(n_regressors + 1)]
      p_vals <- summary(model)$coefficients[,4]

      return(list(coef, p_vals))
    }
  }

  # Perform regression and time it
  time <- system.time(
    coef_and_ps <- mclapply(data.frame(fmri_data$acts), fit_cols, mc.cores = cores[1], mc.silent = TRUE)
    #coefficients <- mclapply(data.frame(fmri_data$acts), fit_cols, mc.cores = cores[1], mc.silent = TRUE)
    #coefficients <- lapply(data.frame(fmri_data$acts), fit_cols)
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
  #saveRDS(coefficients, 'betas2.rds')
  #coef_and_ps <- readRDS('coef_and_ps.rds')

  # Return coefficients to desired format, using only significant voxels
  coefficients <- data.frame(lapply(coef_and_ps, function(x){x[[2]]} ))
  p_values     <- data.frame(lapply(coef_and_ps, function(x){x[[2]]} ))
  sig_voxels   <- colSums(p_values[2:4,]) < 0.05
  roi_beg_end  <- c(1, cumsum(fmri_data$n_voxels_by_roi))
  sig_beg_end  <- integer(8)
  for (i in 2:8) { sig_beg_end[i] <- sum(sig_voxels[roi_beg_end[i - 1]:roi_beg_end[i]]) }
  sig_beg_end <- cumsum(sig_beg_end)
  sig_beg_end[1] <- 1
  y <- matrix(0, 15, 7)
  for (i in 2:8) {y[,i - 1] <- rowSums(coefficients[ 2:15, sig_beg_end[i - 1]:sig_beg_end[i]]) }
  coefficients <- y
  rownames(coefficients) <- colnames(design_mat)
  colnames(coefficients) <- fmri_data$rois
  coefficients <- data.frame(coefficients)

  # Return names...
  if (!seperate) { rownames(coefficients) <- colnames(design_mat)}

  # Some diagnostic plots...
  #ggplot(melt(rob_model$fitted.values), aes(1:444,value)) + geom_point(col = 'red')
  #+ geom_point(aes(1:444,melt(fmri_data$acts[,503])), col = 'blue')

  # Save the linear model, list of conds removed,
  lm_list <- list(coef = coefficients, bad_cond = bad_cond, design = design_mat)
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
#               A finite impulse response model, currently unused.               #
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
# ------------------------------------------------------------------------------ #
attach_fmri_results <- function(data) {
  # Collect the data
  ag_subjs     <- readRDS('/home/dan/projects/imagen/data/ag_subjs.rds')
  ag_results   <- readRDS('/home/dan/projects/imagen/data/ag_results.rds')

  fmri_subjs   <- readRDS('/home/dan/projects/imagen/data/fmri_subjs.rds')
  fmri_results <- readRDS('/home/dan/projects/imagen/data/fmri_betas.rds')

  # Name the connection strengths according to condition
  cnames  <- colnames(ag_results$betas$connectivity_ST)

  conn_st <- ag_results$betas$connectivity_ST
  colnames(conn_st) <- paste(colnames(conn_st), 'st', sep = '_')

  conn_sr <- ag_results$betas$connectivity_SR
  colnames(conn_sr) <- paste(colnames(conn_sr), 'sr', sep = '_')

  conn_ctrst <- conn_st - conn_sr
  colnames(conn_ctrst) <- paste(cnames, 'st_sr', sep = '_')

  # Attach AG results to the data frame
  stop_betas <- data.frame('Subject' = ag_subjs, conn_st, conn_sr, conn_ctrst)
  data$raw   <- merge(data$raw, stop_betas, by = 'Subject', all = TRUE)

  # Compute Stop > Stop respond contrasts and attach them
  ctrst <- t( sapply(fmri_results, function(x){ unlist(x[2,] - x[3,]) }))
  colnames(ctrst) <- paste(colnames(ctrst), 'st_sr', sep = '_')
  ctrst_frame <- data.frame('Subject' = fmri_subjs, data.matrix(ctrst))
  data$raw   <- merge(data$raw, ctrst_frame, by = 'Subject', all = TRUE)

  # Compute Stop > Go contrasts and attach them
  ctrst <- t( sapply(fmri_results, function(x){ unlist(x[2,] - x[1,]) }))
  colnames(ctrst) <- paste(colnames(ctrst), 'st_go', sep = '_')
  ctrst_frame <- data.frame('Subject' = fmri_subjs, data.matrix(ctrst))
  data$raw   <- merge(data$raw, ctrst_frame, by = 'Subject', all = TRUE)

  # Attach stop, stop respond, and go betas
  ctrst <- t( sapply(fmri_results, function(x){ unlist(x[2,]) }))
  colnames(ctrst) <- paste(colnames(ctrst), 'st', sep = '_')
  ctrst_frame <- data.frame('Subject' = fmri_subjs, data.matrix(ctrst))
  data$raw   <- merge(data$raw, ctrst_frame, by = 'Subject', all = TRUE)

  ctrst <- t( sapply(fmri_results, function(x){ unlist(x[3,]) }))
  colnames(ctrst) <- paste(colnames(ctrst), 'sr', sep = '_')
  ctrst_frame <- data.frame('Subject' = fmri_subjs, data.matrix(ctrst))
  data$raw   <- merge(data$raw, ctrst_frame, by = 'Subject', all = TRUE)

  ctrst <- t( sapply(fmri_results, function(x){ unlist(x[1,]) }))
  colnames(ctrst) <- paste(colnames(ctrst), 'go', sep = '_')
  ctrst_frame <- data.frame('Subject' = fmri_subjs, data.matrix(ctrst))
  data$raw   <- merge(data$raw, ctrst_frame, by = 'Subject', all = TRUE)

  return(data)
}
# ------------------------------------------------------------------------------ #
