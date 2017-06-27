# ------------------------------------------------------------------------------ #
#                Subroutine for reading the SST parameters
# ------------------------------------------------------------------------------ #
read_sst_param_file <- function(filename){

  # Read the parameter file
  flog.info('Reading SST params from file: %s', filename)
  raw <- read.csv(file = filename, header = TRUE, sep = ';');

  # Determine params being read and subject numbers
  regex <- paste('*_subj.*', sep = '')
  param_names <- grep(glob2rx(regex), names(raw), value = TRUE)
  subject_ids <- sort(strtoi(unique( gsub('.+_subj.', '', param_names) )))

  # Strip subj nums
  param_names <- unique(gsub('subj.+', 'subj', param_names))
  flog.info('SST params identified: %s', toString(param_names))

  # Set up the dataframe to accept mean parameter values.
  n_subj   <- max(subject_ids)
  n_params <- length(param_names)
  data     <- data.frame(matrix(ncol = n_params, nrow = n_subj))
  colnames(data) <- param_names

  # Get mean chain values for each param
  trace_means <- colMeans(raw[2000:3000,])

  # Save the parameters:
  for (subj_num in subject_ids) {
    #flog.info('Reading params for subject num: %d', subj_num)
    for (param in param_names) {
      long_name <- paste(param, '.', subj_num, sep = '')
      data[subj_num, param] <- trace_means[long_name]
    }
  }

  # Remove the 'subj' part of param names
  names(data) <- gsub('_subj', '', names(data))

  return(data)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#                Subroutine for reading the MID parameters
# ------------------------------------------------------------------------------ #
read_mid_param_file <- function(filename){

  mid_params_raw  <- read.csv(file = filename, header = TRUE, sep = ',');

  mid_params       <- data.frame(matrix(ncol = 12, nrow = dim(unique(mid_params_raw['Subject']))[1]))
  mid_params_names <- c('Target Duration Coeff', 'Target is Left Coeff', 'Cue is Rewarded Coeff',
                        'Cue is High Reward Coeff', 'StDev of 1/RT', 'Intercept')
  mid_params_type  <- c('Expectation', 'Variance')

  colnames(mid_params) <- c( outer( mid_params_names, c('Expectation', 'Variance'), FUN=paste ,sep=" "))

  mid_params['Subject'] <- unique(mid_params_raw['Subject'])
  for (row in  1:dim(mid_params_raw)[1]) {
    subject <- mid_params_raw[row, 'Subject']
    param   <- mid_params_raw[row, 'Parameter']

    subj_msk <- mid_params['Subject'] == subject

    mid_params[subj_msk, paste(param,'Expectation',sep=' ')] <- mid_params_raw[row,'Expectation']
    mid_params[subj_msk, paste(param,'Variance'   ,sep=' ')] <- mid_params_raw[row,'Variance'   ]
  }
  return(mid_params)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#                      Subroutine for reading the fmri data                      #
# ------------------------------------------------------------------------------ #
import_fmri <- function(subj_id, timeseries_dir, taskdata_dir) {

  # Set filename prefixes and suffixes
  task_pfx <- 'ss_'
  fmri_sfx <- paste('MaskDumpOutput',sep='')

  # Specify ordinal coding for trial outcome
  task_key  <- c("GO_SUCCESS"="1", "GO_FAILURE"='2', 'STOP_SUCCESS'='3', 'STOP_FAILURE'='4',
                   'GO_TOO_LATE'='5','GO_WRONG_KEY_RESPONSE'='6', 'STOP_TOO_EARLY_RESPONSE'='7')

  #failed_fmri_reads <- 0
  #failed_task_reads <- 0

  #failed_fmri_cols  <- c()
  #failed_task_cols  <- c()

  # Convert subject id to padded string & read
  subj_id_str <- formatC(subj_id, width = 12, format = 'd', flag = '0')

  fmri_fname <- paste(subj_id_str, '_', fmri_sfx, sep='')
  task_fname <- paste(task_pfx, subj_id_str, '.csv', sep='')

  fmri_full_name <- paste(timeseries_dir, fmri_fname, sep='')
  task_full_name <- paste(taskdata_dir  , task_fname, sep='')

  # Reading the fmri Series
  if (file.exists(fmri_full_name)) {
    activations <- read.csv(file=fmri_full_name, header=FALSE, sep=' ')
    activations <- t(activations[,4:447])
    rownames(activations) <- NULL
    print(paste('Successful read of ', fmri_fname))
  } else {

    print(paste('Failure to read ', fmri_fname))
    return
    #failed_fmri_reads <- failed_fmri_reads + 1
    #failed_fmri_cols  <- c(failed_fmri_cols, subj_num)
  }

  # Reading the SST Series
  if (file.exists(task_full_name)) {
    # Only want to read 2 columns, time and outcome
    cols <- c('NULL','NULL','numeric','NULL','NULL','NULL','NULL','NULL',
              'NULL','NULL','NULL'   , NA   ,'NULL','NULL','NULL','numeric')
    task_data <- read.csv(file = task_full_name, header = TRUE, sep = '\t', colClasses = cols, skip = 1)
    print(paste('Successful read of ', task_fname))
  } else {
    print(paste('Failure to read ', task_fname))
    return
    #failed_task_reads <- failed_task_reads + 1
    #failed_task_cols  <- c(failed_task_cols, subj_num)
  }


  # Req. library for revalue()
  library(plyr)

  # Extract the trial times and outcomes
  task_times <- task_data['Trial.Start.Time..Onset.'][,1]
  task_data['Response.Outcome'][,1] <- revalue(task_data['Response.Outcome'][,1], task_key, warn_missing=FALSE)
  task_outcome <- as.numeric(as.character(task_data['Response.Outcome'][,1]))


  #print(paste('Completed ', n_train - failed_fmri_reads, 'out of', n_train, 'fmri reads.', sep=' '))
  #print(paste('Completed ', n_train - failed_task_reads, 'out of', n_train, 'SST reads.', sep=' '))

  #lost_cols  <- c(failed_fmri_cols, failed_task_cols)
  #valid_cols <- setdiff(1:subj_num, lost_cols)

  # Format the fmri data as a dataframe and add TRs
  activations <- as.matrix(activations)
  trs <- seq(1,444)

  # If we're looking at just the rIFG (a QC check)...
  if (TRUE) {
    rIFG_voxels <- read.csv(file = '../data/rIFGClusterRows.csv', header = FALSE)
    activations <- activations[, as.integer(rIFG_voxels)]
    activations <- cbind(activations, rowMeans(activations))
  }

  # Save the fmri_series data (for raw <-> processed comparison)
  #saveRDS(activations, 'activations.rds')

  # activations <- activations[,1:3]
  # Z-Score the data for each participant
  mean <- colMeans(activations, na.rm=TRUE)
  std  <- apply(activations, 2, sd, na.rm = TRUE)

  activations <- sweep(activations, 2, mean, FUN = "-")
  activations <- sweep(activations, 2, std , FUN = "/")

  # Spike removal
  spike_inds <- which(abs(activations) >= 3)
  activations[spike_inds] <- 0

  # Drift correction
  polyfit     <- lm( as.matrix(activations) ~ poly( as.vector(trs), 2))
  activations <- activations - polyfit$fitted.values

  # High-Pass Filter
  bf <- signal::butter(2, 1/128, type="high")
  filter <- function(x) {signal::filtfilt(bf, x)}
  activations <- apply(activations, 2, filter)

  # Should drift correct and Z-Score again after filtering

  data <- list(acts = activations, task_outcome = task_outcome, task_times = task_times, task_key = task_key)
  return (data)

} # EOF
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#                Master routine for creating 'data' variable
# ------------------------------------------------------------------------------ #
import_all <- function(){
  library('futile.logger')

  # Where is data located, and who are the subjects?
  base_dir      <- '/gpfs/home/dscott3/projects/imagen/data/'
  sbx_subj_list <- paste(base_dir, 'sandbox_subject_list.csv', sep = '')

  # Read subject list & demarcate test/training split
  data <- c()
  subj_list <- read.csv(sbx_subj_list, header = TRUE, sep = ',')
  data$train_inds <- 1:198
  data$test_inds  <- 199:396

  ### ---------------------- ###
  ### SST AND MID PARAMETERS ###
  ### ---------------------- ###
  # Files to read
  # NOTE: This is just reading the test data twice for debugging purposes!!
  sst_params_files <- paste(base_dir, c('sst_data_test/parameters1.csv', 'sst_data_test/parameters1.csv'), sep = '')
  mid_params_file  <- paste(base_dir, 'MIDT Subject Fits/MIDT_SubjectFits_BL_All.csv', sep = '')

  # Read the SST parameters files
  subj_ids_1   <- as.matrix(subj_list[data$train_inds, 'Subject'])
  sst_params_1 <- read_sst_param_file(sst_params_files[1])
  colnames(subj_ids_1) <- 'Subject'

  subj_ids_2   <- as.matrix(subj_list[data$train_inds, 'Subject'])
  sst_params_2 <- read_sst_param_file(sst_params_files[2])
  colnames(subj_ids_2) <- 'Subject'

  # Package into an data matrix
  data$raw <- rbind( cbind(subj_ids_1, sst_params_1), cbind(subj_ids_2, sst_params_2) )

  # Read the MID fit file in:
  flog.info('Reading MID parameters...')
  mid_params <- read_mid_param_file(mid_params_file)

  # Match rows by 'Subject' column
  tmp <- merge(sst_params, mid_params, by = "Subject")


  ### ---------------------- ###
  ###       FMRI Stuff       ###
  ### ---------------------- ###
  # Get SST FMRI beta values
  #
  # Reading in more than one subjects complete FMRI data at a time would be
  # impossibly memory intensive
  visit <- 'BL'
  task  <- 'SST'
  taskdata_dir   <- paste(base_dir, task, '_', visit, '_behavioraldata/', sep='')
  timeseries_dir <- paste(base_dir, visit, '_VoxelLevel_', task, '/', sep='')

  source('fit_fmri_glm.r')
  for (id in subj_list[1:5]) {
    fmri_data <- import_fmri(id, timeseries_dir, taskdata_dir)
    sst_betas <- fit_fmri_glm(fmri_data)
  }

  #jpeg('rplot.jpg')
  #plot(x,y)
  #dev.off()

  #### DOES NOT EXIST YET ###
  # Get MID FMRI beta values
  #mid_fmri_data <- import_fmri()
  #mid_betas <- fit_fmri_glm(mid_fmri_data)


  ### DOES NOT EXIST YET
  ### ---------------------- ###
  ###   Clinical & Genetic   ###
  ### ---------------------- ###
  #clinical_data <- import_clinical()
  #genetic_data  <- import_genetic()


  ### NEEDS UPDATE
  # Convert annoying names into better ones
  #data <- replace_bad_names(sst_params, mid_params, raw_df)

  return(data)
}
# ------------------------------------------------------------------------------ #
