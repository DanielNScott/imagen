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

  tmp <- matrix(NA, dim(activations)[1], length(rois))
  colnames(tmp) <- rois
  for (roi in rois) {
    #activations[roi] <- rowMeans(activations[, as.logical(roi_voxels)], na.rm = TRUE)
    tmp[,roi]    <- rowMeans(activations[, as.logical(roi_voxels[[roi]]) ], na.rm = TRUE)
  }
  activations <- tmp

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
               task_times = task_times, task_key = task_key, stim = stim, movement = movement)
  return(data)

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
     }
  }

  # Assign appropriate names: Need to list-ify and un-listify for if-else.
  regressor_names <- ifelse(seperate, list(trial_names), list(cond_names))[[1]]
  colnames(conv_regs) <- regressor_names

  # Throw away bad conditions
  cond_nums  <- setdiff(cond_nums, bad_cond)
  cond_names <- cond_names[cond_nums]
  index_list <- index_list[cond_nums]

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


  # Indicate pool for core-level SIMD parallelism:
  # Note: Using makeClust() here will break this on some clusters which do not let R
  #       use socket based core communication. mclapply defaults to fork-based dispatching.
  cores <- detectCores()
  flog.info('Using %d cores', cores[1])


  # Setup for regressions...
  n_voxels     <- dim(fmri_data$acts)[2]
  n_regressors <- dim(design_mat)[2]

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
      lmRob(x ~ design_mat)$coefficients[2:(n_regressors + 1)]
    }
  }

  # Perform regression and time it
  time <- system.time(
    #coefficients <- mclapply(data.frame(fmri_data$acts), fit_cols, mc.cores = cores[1], mc.silent = TRUE)
    coefficients <- lapply(data.frame(fmri_data$acts), fit_cols)
  )

  # Report time required for regression
  flog.info('Computation time for voxel betas:')
  print(time)
  #saveRDS(coefficients, 'betas2.rds')

  # Return coefficients to desired format
  coefficients <- data.frame(coefficients)
  #if (!seperate) {
  #  rownames(coefficients) <- cond_names
  #}

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
