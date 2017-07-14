# ------------------------------------------------------------------------------ #
#                      Subroutine for reading the fmri data                      #
# ------------------------------------------------------------------------------ #
import_fmri <- function(subj_id, timeseries_dir, taskdata_dir, movement_dir) {

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
    cols <- c('NULL','NULL','numeric','NULL','NULL','NULL','NULL','NULL',
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
  task_times <- task_data['Trial.Start.Time..Onset.'][,1]
  task_data['Response.Outcome'][,1] <- revalue(task_data['Response.Outcome'][,1], task_key, warn_missing = FALSE)
  task_outcome <- as.numeric(as.character(task_data['Response.Outcome'][,1]))

  # If we're looking at just the rIFG (a QC check)...
  if (TRUE) {
    rIFG_voxels <- read.csv(file = '../data/rIFGClusterRows.csv', header = FALSE)
    activations <- activations[, as.logical(rIFG_voxels)]
    activations <- cbind(activations, rowMeans(activations, na.rm = TRUE))
  }

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
               task_times = task_times, task_key = task_key, movement = movement)
  return(data)

} # EOF
# ------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------ #
#            Fits a robust general linear model to raw fMRI activations          #
# ------------------------------------------------------------------------------ #
fit_fmri_glm <- function(fmri_data) {

  # Contrast groupings
  conditions <- as.integer(unique(fmri_data$task_key))
  cnames     <- names(fmri_data$task_key)
  n_conds    <- length(conditions)

  # fMRI analysis functions
  library(fmri)

  # Proportion of subjects w/ significant STN betas for STOP_SUCCESS
  frac_go_sig   <- 0
  frac_stop_sig <- 0

  # Number of scans in this series -- varies by individual.
  n_scans <- length(fmri_data$acts[!is.na(fmri_data$acts[,1]),1])

  # Array for events convolved w/ HRF
  conv_regs <- array(dim = c(n_scans, n_conds))

  # Some conditions do not occur - they will be removed
  bad_cond  <- c()

  # List of onset times by condition
  onset_list <- list(c(), c(), c(), c(), c(), c(), c())

  # Fill the conv_regressors array for each condition
  for (cond in conditions) {

    # Masks for getting trial times
    cond_msk <- fmri_data$task_outcome %in% cond

    # If this condition doesn't occur, skip it.
    if (sum(cond_msk) == 0) {
        bad_cond <- cbind(bad_cond,cond)
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
    #print(onsets)

    # Get convolved regressor for condition
    conv_regs[,cond] <- fmri.stimulus(scans = n_scans, onsets = onsets, duration = durations)

    #print(conv_regs)
    onset_list[[cond]] <- onsets
  }

  # Clean up conv_regs a bit, append 'TR' column, label things
  #conv_regs[,bad_cond] <- integer(n_scans)
  #conv_regs <- data.frame(cbind(conv_regs, 1:n_scans))
  #colnames(conv_regs) <- c(cnames,'TR')

  # Create design matrix w/ 2nd deg drift, remove any bad conditions
  conds <- setdiff(1:7, bad_cond)
  design_mat <- fmri.design(conv_regs[, conds], order = 2)

  # Remove the unnecessary intercept column (which will otherwise be identical to col. 1)
  nconds     <- length(conds)
  design_mat <- design_mat[, c(1:nconds, (nconds + 2):(nconds + 3)) ]
  colnames(design_mat) <- c(cnames[conds], 'Linear Drift', 'Sq. Drift')
  cnames <- colnames(design_mat)

  # Add the motion parameters to the set of regressors
  design_mat <- cbind(design_mat, as.matrix(fmri_data$movement))
  colnames(design_mat) <- c(cnames, 'Motion_x', 'Motion_y', 'Motion_z', 'Motion_pitch', 'Motion_yaw', 'Motion_roll')
  cnames <- colnames(design_mat)

  # Indicate pool for core-level SIMD parallelism:
  # Note: Using makeClust() here will break this on some clusters which do not allow
  #       socket assignment and the like. mclapply defaults to fork-based dispatching.
  library(parallel)
  cores <- detectCores()
  flog.info('Using %d cores', cores[1])

  # Specifics for fitting the linear model with robust regression
  library(robust)
  n_voxels <- dim(fmri_data$acts)[2]
  n_regressors <- dim(design_mat)[2]
  coefficients <- matrix(NA, n_regressors, n_voxels)

  # Function to apply lmRob (named y here) to each column of the activation data
  fit_cols <- function(x) {capture.output(
    lmRob(x ~ design_mat)$coefficients[2:(n_regressors + 1)]
  )}

  # Actually do it...
  time <- system.time(
    coefficients <- mclapply(data.frame(fmri_data$acts), fit_cols)
  )
  flog.info('Computation time for voxel betas:')
  print(time)

  # Return coefficients to desired format
  coefficients <- data.frame(coefficients)
  rownames(coefficients) <- cnames

  # Some diagnostic plots...
  #ggplot(melt(rob_model$fitted.values), aes(1:444,value)) + geom_point(col = 'red')
  #+ geom_point(aes(1:444,melt(fmri_data$acts[,503])), col = 'blue')

  # Save the linear model, list of conds removed,
  lm_list <- list(coef = coefficients, bad_cond = bad_cond, design = design_mat, onsets = onset_list)

  return(lm_list)

  ##################################################
  # The rest of this function needs dealing with!
  ##################################################

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

  # Report on betas:
  #frac_go_sig   <- frac_go_sig  /n_to_fit
  #frac_stop_sig <- frac_stop_sig/n_to_fit

  print('')
  print(paste('Fraction of GO_SUCCESS   betas which are significant:', frac_go_sig))
  print(paste('Fraction of STOP_SUCCESS betas which are significant:', frac_stop_sig))

  return (lm_list)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#               A finite impulse response model, currently unused.               #
# ------------------------------------------------------------------------------ #
fir <- function(onsets, durations, activations, n_scans, n_conds, n_voxels) {
  # Scale parameter - sets resolution
  scale <- 10

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
  stims <- rep(0, ceiling(scans))

  ## ESTIMATE HRF USING FIR BASIS SET

  # CREATE FIR DESIGN MATRIX
  # WE ASSUME HRF IS 16 TRS LONG
  hrf_len <- 16

  # BASIS SET FOR EACH CONDITOIN IS A TRAIN OF INPULSES
  fir_bases <- zeros(n_scans, hrf_len*n_conds)

  for (cond in 1:n_conds) {
      col_subset <- ((cond - 1)* hrf_len + 1):(cond*hrf_len)

      for (onset in 1:numel(onsets[,cond]) ) {
          #impulse_times <- onsets(onset):onsets(onset) + hrf_len - 1;
          impulse_times <- seq(onsets(onset), onsets(onset) + hrf_len*scale - 1, scale)

          for (impulse in 1:numel(impulse_times)) {
              fir_bases(impulse_times(impulse), col_subset(impulse)) <- 1;
          }
      }
  }

  # ESTIMATE HRF FOR EACH CONDITION AND VOXEL
  fir_hrf_est <- pseudoinverse(fir_bases %*% fir_bases) * fir_bases %*% activations

  # RESHAPE HRFS
  hHatFIR <- reshape(fir_hrf_est, hrf_len, n_conds, n_voxels)
}
# ------------------------------------------------------------------------------ #
