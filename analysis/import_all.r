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

  colnames(mid_params) <- c( outer( mid_params_names, c('Expectation', 'Variance'), FUN = paste ,sep = " "))

  mid_params['Subject'] <- unique(mid_params_raw['Subject'])
  for (row in  1:dim(mid_params_raw)[1]) {
    subject <- mid_params_raw[row, 'Subject']
    param   <- mid_params_raw[row, 'Parameter']

    subj_msk <- mid_params['Subject'] == subject

    mid_params[subj_msk, paste(param,'Expectation',sep = ' ')] <- mid_params_raw[row,'Expectation']
    mid_params[subj_msk, paste(param,'Variance'   ,sep = ' ')] <- mid_params_raw[row,'Variance'   ]
  }

  new_param_names      <- c('mid_dur', 'mid_loc', 'mid_rew', 'mid_high_rew', 'mid_tau_stdev', 'mid_int')
  new_expanded_names   <- c( outer( new_param_names, c('', '_var'), FUN = paste ,sep = ""))
  colnames(mid_params) <- c(new_expanded_names, 'Subject')

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
#  Reads and processes fMRI data into a format useful for storing in primary DF  #
# ------------------------------------------------------------------------------ #
import_fmri <- function(data, base_dir, subj_beg, subj_end, loc = 'local', process_fmri = FALSE, seperate = FALSE, sig_voxels = NULL) {

  if (process_fmri) {
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

    fmri_betas <- list()
    fmri_sig   <- list()
    subj       <- c()
    subj_num   <- 1

    for (id in data$subj_list[subj_beg:subj_end,]) {
      tryCatch({

        # Retrieve and preprocess rois and trial data
        fmri_data <- read_fmri_data(id, timeseries_dir, taskdata_dir, movement_dir, base_dir)

        # Fit linear models to roi activations for std and AG analysis
        fmri_lm_fit  <- fit_fmri_glm(fmri_data, seperate = seperate, sig_voxels)

        # Set up data structure for use in AG analysis
        subj_fldname <- paste('subj-', sprintf('%04i', subj_num), sep = '')
        fmri_betas[[subj_fldname]] <- fmri_lm_fit$coef
        fmri_sig[[subj_fldname]] <- list(sig_voxels = fmri_lm_fit$sig_voxels, sig_beg_end = fmri_lm_fit$sig_beg_end)

        # Book keeping
        subj <- c(subj, id)
        subj_num  <- subj_num + 1

        print(paste('fMRI preprocessing success for subj', id))

      }, error = function(e) {

        # If something went wrong...
        print(paste('fMRI preprocessing failure for subj', id))
        print(e)
      })

      # Sometimes the error message 'e' is not very helpful so check warnings...
      if (exists('w')) print(w)
    }

    # Prefix to create different filenames if recovering AG betas or not.
    prfx <- ifelse(seperate, 'std', 'ag')

    saveRDS(fmri_betas, paste(base_dir, prfx, '_fmri_betas_', subj_beg, '_', subj_end, '.rds', sep = ''))
    saveRDS(subj      , paste(base_dir, prfx, '_fmri_subjs_', subj_beg, '_', subj_end, '.rds', sep = ''))
    saveRDS(fmri_sig  , paste(base_dir, prfx, '_fmri_sig_'  , subj_beg, '_', subj_end, '.rds', sep = ''))

    # Load and attach already-processed fmri data
    data <- attach_fmri_results(data, base_dir)
    return(data)

  } else{

    # Just load and attach already-processed fmri data
    data <- attach_fmri_results(data, base_dir)
    return(data)
  }
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#                Master routine for creating 'data' variable
# ------------------------------------------------------------------------------ #
import_all <- function(loc, subj_beg = 1, subj_end = 10, fmri_only = FALSE){

  library('futile.logger')
  #source('fmri_routines.r')

  # Different base directories for cluster or local use
  if (loc == 'local') { base_dir <- '/home/dan/projects/imagen/data/'
  } else {              base_dir <- '/users/dscott3/projects/imagen/data/'}

  # Sandbox' subjects are training set once analyses pipeline is done
  sbx_subj_list <- paste(base_dir, 'sandbox_subject_list.csv', sep = '')

  # Read subject list & demarcate test/training split
  data <- c()
  data$subj_list  <- read.csv(sbx_subj_list, header = TRUE, sep = ',')
  data$train_inds <- 1:198
  data$test_inds  <- 199:396

  # Initialize principle data frame (data$raw) with subjects and their partition.
  data$raw <- data.frame(Subject = data$subj_list, set = rep( as.factor(c('train', 'test')), each = 198))

  if (!fmri_only) {
  ### ---------------------- ###
  ### SST AND MID PARAMETERS ###
  ### ---------------------- ###
  # Files to read
  sst_params_files <- paste(base_dir, c('test_parameters.csv', 'train_parameters.csv'), sep = '')

  # Read the SST parameters files
  # Test partition
  subj_ids_1   <- as.matrix(data$subj_list[data$test_inds, 'Subject'])
  sst_params_1 <- read_sst_param_file(sst_params_files[1])
  colnames(subj_ids_1) <- 'Subject'

  # Training partition
  subj_ids_2   <- as.matrix(data$subj_list[data$train_inds, 'Subject'])
  sst_params_2 <- read_sst_param_file(sst_params_files[2])
  colnames(subj_ids_2) <- 'Subject'

  # Package into an data matrix
  # Columns are in different orders, but row bind matches on them.
  sst_data <- rbind( cbind(subj_ids_1, sst_params_1), cbind(subj_ids_2, sst_params_2))
  data$raw <- merge(data$raw, sst_data, by = "Subject", all = TRUE)
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
  data <- import_fmri(data, base_dir, subj_beg, subj_end, loc, seperate = FALSE, process_fmri = TRUE)
  #data <- import_fmri(data, base_dir, subj_beg, subj_end, loc, seperate = TRUE , process_fmri = FALSE, sig_voxels = sig_vox)

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

  return(data)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#             Generate some additional fields (i.e. feature construction)        #
# ------------------------------------------------------------------------------ #
gen_addnl_flds <- function(data) {

  # TCI impulsivity really looks like it falls into categories 'general' and 'financial' to me...
  # Also, it looks like taking rowSums is fine -- the summary vars get NAs if there are NAs in row.
  tci_gen_imp <- c('tci010', 'tci014', 'tci047', 'tci071', 'tci102', 'tci123', 'tci193', 'tci210')
  tci_fin_imp <- c('tci024', 'tci059', 'tci105', 'tci215')

  data$raw['tci_gen_imp'] <- rowSums(data$raw[, tci_gen_imp])
  data$raw['tci_fin_imp'] <- rowSums(data$raw[, tci_fin_imp])
  data$names$TCI <- c(data$names$TCI, 'tci_gen_imp', 'tci_fin_imp')

  #SURPS impulsivity just needs to be summed.
  # Question 22 about being manipulative is not included here.
  surps_imp   <- c('surps5', 'surps2', 'surps11', 'surps15')

  data$raw['surps_imp']   <- rowSums(data$raw[, surps_imp])
  data$names$SURPS <- c(data$names$SURPS, 'surps_imp')

  # Need some ESPAD summaries and other stuff...
  # Individual fields:
  # C.18i    - likelihood of regret
  # C.prev31 - scaled num drinks on typical day in which drinking
  # C.21     - scaled num drinks needed to get drunk
  # C.19b    - scaled num times drunk in last year
  # C.6      - scaled num times smoked in life -- should match 'espad_6_life_nic'
  #
  #data$raw['binge'] <- data$raw[,'C.prev31'] / data$raw[,'bmi']
  #data$names$ESPAD  <- c(data$names$ESPAD, 'binge')

  # Reconstruct stop and go distributions and their variances
  data$raw$SSRTGo    <- data$raw$mu_go   + data$raw$tau_go
  data$raw$SSRTStop  <- data$raw$mu_stop + data$raw$tau_stop

  data$raw$SSRTVarGo   <- data$raw$sigma_go^2   + data$raw$tau_go^2
  data$raw$SSRTVarStop <- data$raw$sigma_stop^2 + data$raw$tau_stop^2

  data$names$sst <- c(data$names$sst, 'SSRTGo', 'SSRTStop', 'SSRTVarGo', 'SSRTVarStop')

  return(data)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #
attach_fmri_results <- function(data, base_dir) {
  # Collect the data
  ag_subjs     <- readRDS(paste(base_dir, '/ag_subjs.rds'  , sep = ''))
  ag_results   <- readRDS(paste(base_dir, '/ag_results.rds', sep = ''))

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

  # Agglomerate output from node level parallelism
  subj_files <- list(c(1,25), c(26,50), c(51,75), c(76,100), c(101,125),
                     c(126,150), c(151,175), c(176,200), c(201,250), c(251,300),
                     c(301,350), c(351,396))

  attach <- function(data, subjs, results, suffix, fn, merge_flag) {
    ctrst <- t( sapply(results, fn))
    colnames(ctrst) <- paste(colnames(ctrst), suffix, sep = '_')
    if (merge_flag) {
      ctrst_frame <- data.frame('Subject' = subjs, data.matrix(ctrst))
      data$raw   <- merge(data$raw, ctrst_frame, by = 'Subject', all = TRUE)
    } else {
      row_indices <- match(subjs, data$raw[['Subject']])
      data$raw[colnames(ctrst)][row_indices, ] <- ctrst
    }
    return(data$raw)
  }

  for (i in 1:length(subj_files)) {
    beg <- subj_files[[i]][1]
    end <- subj_files[[i]][2]
    tmp_results <- readRDS(paste(base_dir, '/fmri/fmri_betas_', beg, '_', end, '.rds', sep = ''))
    tmp_subjs   <- readRDS(paste(base_dir, '/fmri/fmri_subjs_', beg, '_', end, '.rds', sep = ''))

    data$raw <- attach(data, tmp_subjs, tmp_results, 'st_sr', function(x){ unlist(x[2,] - x[3,]) }, i == 1)
    data$raw <- attach(data, tmp_subjs, tmp_results, 'st_go', function(x){ unlist(x[2,] - x[1,]) }, i == 1)
    data$raw <- attach(data, tmp_subjs, tmp_results, 'st'   , function(x){ unlist(x[2,]        ) }, i == 1)
    data$raw <- attach(data, tmp_subjs, tmp_results, 'sr'   , function(x){ unlist(x[3,]        ) }, i == 1)
    data$raw <- attach(data, tmp_subjs, tmp_results, 'go'   , function(x){ unlist(x[1,]        ) }, i == 1)
  }

  sfx <- c('_go', '_st', '_sr', '_st_go', '_st_sr')
  rois <- c('rPreSMA', 'rIFG', 'rCaudate', 'rSTN', 'rGPe', 'rGPi', 'rThalamus')
  titles <- c('SST Go Weights', 'SST Stop Weights', 'SST Stop Respond Weights',
              'SST Stop > Go Contrasts', 'SST Stop > Stop Respond Contrasts')
  std_fmri_feats <- as.vector(outer(rois, sfx, function(x,y) {paste(x,y, sep = '')}))

  # Convert zeros to NAs
  data$raw[ , std_fmri_feats][data$raw[,std_fmri_feats] == 0] <- NA
  return(data)
}
# ------------------------------------------------------------------------------ #