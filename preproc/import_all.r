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
#                Master routine for creating 'data' variable
# ------------------------------------------------------------------------------ #
import_all <- function(){
  library('futile.logger')

  # Where is data located, and who are the subjects?
  base_dir      <- '/home/dan/projects/imagen/data/'
  sbx_subj_list <- paste(base_dir, 'sandbox_subject_list.csv', sep = '')

  # Read subject list & demarcate test/training split
  data <- data.frame()
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
  taskdata_ dir  <- paste('/home/dan/documents/lncc/From Catherine/Round 2/SST/SST_', visit,'_behavioraldata/',sep='')
  timeseries_dir <- paste('/home/dan/documents/lncc/From Catherine/', visit, '/',roi,'_',visit,'_SSTtimeseries/', sep='')

  for (id in subj_list) {
    fmri_data <- import_fmri(id, 'SST', timeseries_dir, taskdata_dir)
    sst_betas <- get_fmri_betas(fmri_data)
  }

  # Get MID FMRI beta values #### DOES NOT EXIST YET ###
  #mid_fmri_data <- import_fmri()
  #mid_betas <- get_fmri_betas(mid_fmri_data)


  ### ---------------------- ###
  ###   Clinical & Genetic   ### DOES NOT EXIST YET
  ### ---------------------- ###
  #clinical_data <- import_clinical()
  #genetic_data  <- import_genetic()


  # Convert annoying names into better ones
  data <- replace_bad_names(sst_params, mid_params, raw_df)

  return(data)
}
# ------------------------------------------------------------------------------ #