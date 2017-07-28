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

    mid_params[subj_msk, paste(param,'Expectation',sep = ' ')] <- mid_params_raw[row,'Expectation']
    mid_params[subj_msk, paste(param,'Variance'   ,sep = ' ')] <- mid_params_raw[row,'Variance'   ]
  }
  return(mid_params)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#                Subroutine for reading the genetic data
# ------------------------------------------------------------------------------ #
import_genes <- function(dose_file, rs_file, subj_ids){

  # Read the dosages
  flog.info('Reading gene dosages from file: %s', dose_file)
  dose_data <- data.frame( read.csv(file = dose_file, header = TRUE) );

  # Read the file indicating which genes to keep
  flog.info('Reading list of genes to keep from file: %s', rs_file)
  rs_list <- data.frame( read.csv(file = rs_file, header = TRUE) )

  # Subset the gene data for dopamine genes and the appropriate subjects
  dopamine_genes <- intersect(rs_list$gene, names(dose_data))
  dose_data <- dose_data[, c('Subject', dopamine_genes) ]
  dose_data <- dose_data[dose_data$Subject %in% subj_ids$Subject, ]

  return(dose_data)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#                Subroutine for reading generic data
# ------------------------------------------------------------------------------ #
import_generic <- function(data_file, subj_ids){

  data <- c()

  # Read the dosages
  flog.info('Reading generic file: %s', data_file)
  raw <- data.frame( read.csv(file = data_file, header = TRUE) );

  # Subset the gene data for dopamine genes and the appropriate subjects
  raw <- raw[raw$Subject %in% subj_ids$Subject, ]

  # Save names to an index. Have to remove path and suffix
  names <- list()
  subset_name <- gsub( '.+/', '', gsub('_.+', '', data_file) )
  names[subset_name] <- list(setdiff(colnames(raw), c('Subject', 'Gender')))

  data$raw   <- raw
  data$names <- names
  return(data)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#                Master routine for creating 'data' variable
# ------------------------------------------------------------------------------ #
import_all <- function(loc){

  library('futile.logger')
  source('fmri_routines.r')

  # Where is data located, and who are the subjects?
  if (loc == 'local') {
    base_dir      <- '/home/dan/projects/imagen/data/'
  } else {
    base_dir      <- '/users/dscott3/projects/imagen/data/'
  }
  sbx_subj_list <- paste(base_dir, 'sandbox_subject_list.csv', sep = '')

  # Read subject list & demarcate test/training split
  data <- c()
  data$subj_list  <- read.csv(sbx_subj_list, header = TRUE, sep = ',')
  data$train_inds <- 1:198
  data$test_inds  <- 199:396

  if (FALSE) {
  ### ---------------------- ###
  ### SST AND MID PARAMETERS ###
  ### ---------------------- ###
  # Files to read
  sst_params_files <- paste(base_dir, c('test_parameters.csv', 'train_parameters.csv'), sep = '')

  # Read the SST parameters files
  subj_ids_1   <- as.matrix(data$subj_list[data$test_inds, 'Subject'])
  sst_params_1 <- read_sst_param_file(sst_params_files[1])
  colnames(subj_ids_1) <- 'Subject'

  subj_ids_2   <- as.matrix(data$subj_list[data$train_inds, 'Subject'])
  sst_params_2 <- read_sst_param_file(sst_params_files[2])
  colnames(subj_ids_2) <- 'Subject'

  # Package into an data matrix
  data$raw <- rbind( cbind(subj_ids_1, sst_params_1), cbind(subj_ids_2, sst_params_2) )
  data$names$sst <- setdiff(colnames(sst_params_1), c('Subject'))

  # Read the MID fit file in:
  mid_params_file  <- paste(base_dir, 'MIDT_SubjectFits_BL_All.csv', sep = '')
  flog.info('Reading MID parameters...')
  mid_params <- read_mid_param_file(mid_params_file)

  # Match rows by 'Subject' column
  data$raw <- merge(data$raw, mid_params, by = "Subject", all = TRUE)
  data$names$mid <- setdiff(colnames(mid_params), c('Subject'))

  }
  ### ---------------------- ###
  ###       FMRI Stuff       ###
  ### ---------------------- ###
  # Get SST FMRI beta values
  #
  # Reading in more than one subjects complete FMRI data at a time would be
  # impossibly memory intensive
  visit <- 'BL'
  task  <- 'SST'

  if (loc == 'local') {
    # For cluster use:
    taskdata_dir   <- paste(base_dir, '/fmri/', sep = '')
    timeseries_dir <- paste(base_dir, '/fmri/', sep = '')
    movement_dir   <- paste(base_dir, '/fmri/', sep = '')
  } else {
    # For local testing:
    taskdata_dir   <- paste(base_dir, '/BL_SST_task/', sep = '')
    timeseries_dir <- paste(base_dir, '/BL_SST_fMRI/', sep = '')
    movement_dir   <- paste(base_dir, '/BL_SST_move/', sep = '')
  }

  ag_num  <- 1
  ag_data <- list()
  ag_subj <- c()

  for (id in data$subj_list[,]) {
    tryCatch({

      # Retrieve and preprocess rois and trial data
      fmri_data <- import_fmri(id, timeseries_dir, taskdata_dir, movement_dir, base_dir)

      # Fit linear models to roi activations for std and AG analysis
      #sst_roi_coef  <- fit_fmri_glm(fmri_data, seperate = FALSE)
      ag_trial_coef <- fit_fmri_glm(fmri_data, seperate = TRUE)

      # Set up data structure for use in AG analysis
      subj_fldname <- paste('subj-', sprintf('%04i', ag_num), sep = '')
      ag_input[[subj_fldname]] <- ag_trial_coef$coef

      # Book keeping
      ag_subj <- c(ag_subj, id)
      ag_num  <- ag_num + 1

      print(paste('AG preprocessing success for subj', id))

    }, error = function(e) {

      # If something went wrong...
      print(paste('AG preprocessing failure for subj', id))
      print(e)
    })

    # Sometimes the error message 'e' is not very helpful so check warnings...
    if (exists('w')) print(w)
  }
  #saveRDS(ag_input, paste(base_dir, '/data/', 'ag_input.rds', sep = ''))
  #saveRDS(ag_subj , paste(base_dir, '/data/', 'ag_subjs.rds', sep = ''))

  #### DOES NOT EXIST YET ###
  # Get MID FMRI beta values
  #mid_fmri_data <- import_fmri()
  #mid_betas <- fit_fmri_glm(mid_fmri_data)

  ### ---------------------- ###
  ###     Genetic Data       ###
  ### ---------------------- ###
  gene_file <- paste(base_dir, 'snps_dosage.csv', sep = '')
  gene_list <- paste(base_dir, 'snps_dopamine_list.csv', sep = '')
  gene_data <- import_genes(gene_file, gene_list, data$subj_list)

  # Match rows by 'Subject' column
  data$raw <- merge(data$raw, gene_data, by = "Subject", all = TRUE)
  data$names$genes <- setdiff(colnames(gene_data), c('Subject'))


  ### ---------------------- ###
  ###   Other Instruments    ###
  ### ---------------------- ###
  file_list <- c('DelayDiscounting_K_14.csv', 'ESPAD_Life_14.csv', 'CANTAB_14.csv', 'Age_IQ_Etc_14.csv',
                 'SURPS_14.csv', 'DAWBA_SDQ_14.csv', 'TCI_14.csv', 'Misc_14.csv')
  for (file in file_list) {
    file_to_import <- paste(base_dir, file, sep = '')
    file_data <- import_generic(file_to_import, data$subj_list)

    # Match rows by 'Subject' column
    data$raw <- merge(data$raw, file_data$raw, by = "Subject", all = TRUE)

    # Add names to groupings
    data$names[names(file_data$names)] <- file_data$names
  }

  ### NEEDS UPDATE
  # Convert annoying names into better ones
  #data <- replace_bad_names(sst_params, mid_params, raw_df)


  ### Generate additional summary fields and the like
  source('gen_addnl_flds.r')
  data <- gen_addnl_flds(data)

  return(data)
}
# ------------------------------------------------------------------------------ #
