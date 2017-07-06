fit_fmri_glm <- function(fmri_data) {

  # Contrast groupings
  conditions <- as.integer(unique(fmri_data$task_key))
  cnames     <- names(fmri_data$task_key)
  n_conds    <- length(conditions)

  # fMRI analysis functions
  library(fmri)
  #library(MASS)

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
    onsets    <- fmri_data$task_times[cond_msk]
    onsets    <- as.vector(na.omit(onsets/1000/2.2))
    durations <- double(length(onsets))
    #print(onsets)

    # Get convolved regressor for condition
    conv_regs[,cond] <- fmri.stimulus(scans = n_scans, onsets = onsets, duration = durations)

    #print(conv_regs)
    onset_list[[cond]] <- onsets
  }

  # Create design matrix w/ 2nd deg drift, remove any bad conditions
  design_mat <- fmri.design(conv_regs[, setdiff(1:7, bad_cond)], order = 2)
  colnames(design_mat) <- c(cnames[setdiff(1:7, bad_cond)], 'Intercept', 'Linear Drift', 'Sq. Drift')

  # Remove the unnecessary intercept column (which will otherwise be identical to col. 1)
  design_mat <- design_mat[, 1:(dim(design_mat)[2] - 2) ]

  # Clean up conv_regs a bit, append 'TR' column, label things
  conv_regs[,bad_cond] <- integer(n_scans)
  conv_regs  <- data.frame(cbind(conv_regs, 1:n_scans))
  #colnames(conv_regs) <- c(cnames,'TR')

  # Fit the linear model
  linear_model <- lm(fmri_data$acts ~ design_mat)

  # Save the linear model, list of conds removed,
  lm_list <- list(lm = linear_model, bad_cond = bad_cond, conv_regs = conv_regs, onsets = onset_list)

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
