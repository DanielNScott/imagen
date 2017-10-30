# ------------------------------------------------------------------------------ #
#         This function performs misc. manually specified transforms             #
# ------------------------------------------------------------------------------ #
misc_transforms <- function(data) {

  # Age is in stupid units
  # data$raw['Age'] <- data$raw['Age']*365.24

  # Gender and Sex are factors
  data$raw$Gender.x <- as.integer(data$raw$Gender.x)
  data$raw$Sex      <- as.integer(data$raw$Sex)

  # -1 is a no-data code... I think?
  data$raw[data$raw == -1] <- NA

  return(data)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#                                   Z-Scoring...                                 #
# ------------------------------------------------------------------------------ #
z_score_data <- function(data, ignore) {
  # NOTICE: Here the input variable 'data' corresponds with 'data$raw' generally.
  # This is because the output is saved as 'data$scored' rather than 'data$raw'

  feature_list <- colnames(data)
  for (feature in setdiff(feature_list, ignore)) {
    flog.info('Z-scoring %s', feature)
    data[feature] <- data[feature] - colMeans(data[feature], na.rm = TRUE)
    data[feature] <- data[feature] / apply(data[feature], 2, sd, na.rm = TRUE)
  }

  return(data)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#               Subroutine for getting rid of useless stuff                      #
# ------------------------------------------------------------------------------ #
drop_useless_flds <- function(data){

  # Identify columns with nothing but zeros
  colsum_is_zero <- colSums(data$raw != 0, na.rm = TRUE) == 0
  colsum_is_zero <- names(colsum_is_zero[colsum_is_zero == TRUE])
  flog.info('These all zero cols will be dropped: %s', toString(colsum_is_zero))

  # Drop some other features that are suspect
  #drop_list <- c('sig_int_14', 'sig_std_inv_rt_14', 'mu_std_inv_rt_14')
  useless   <- c('QR_flag', 'Gender.y')
  drop_list <- c(useless, colsum_is_zero)
  data$raw  <- data$raw[, !names(data$raw) %in% drop_list ]

  # Remove names from indices
  for (fld in names(data$names) ) {
    data$names[[fld]] <- setdiff(data$names[[fld]], drop_list)
  }

  #data$task_names_14 <- setdiff(c(data$task_names_14, 'stn_std'), drop_list)
  #n_feat_task_14     <- length(data$task_names_14)

  # Indicate which features are still a part of the data:
  #flog('The following features were dropped:\n')
  #setdiff(raw_cols, colnames(data$raw))

  return(data)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#                       Subroutine for removing outliers                         #
# ------------------------------------------------------------------------------ #
remove_outliers <- function(data, ignore, thresh) {

  for (feature in names(data$raw)) {
    if (feature %in% ignore) {
      next
    }
    flog.info('Checking for outliers: %s', feature)

    var   <- data$raw[feature]
    stdev <- apply(var, 2, sd    , na.rm = TRUE)
    dmed  <- apply(var, 2, median, na.rm = TRUE)
    #dmean<- apply(var, 2, mean  , na.rm = TRUE)

    outmsk <- apply(var, 2, (function(x) ifelse(abs(dmed - x) < thresh*stdev | is.na(x) ,FALSE ,TRUE)))

    if (sum(outmsk) > 0) {
      print(paste('Removing ', toString(sum(outmsk)), ' points from ', feature))
      print(data$raw[feature][outmsk])
      cat('\n')

      data$raw[feature][outmsk] <- NA
    }
  }
  return(data)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# This function controls the
# ------------------------------------------------------------------------------ #
prepare_data <- function(re_read = TRUE, shuffle = FALSE, n_imputed = 30, max_iters = 1000) {

    #
    if (re_read) {
        # Actually read all the files
        # (Time consuming!)
        source('import_all.r')
        data <- import_all()
        saveRDS(data,'../data/raw_data_from_import_all.rds')
    } else {
        # Read saved R data file:
        data  <- readRDS('../data/raw_data_from_import_all.rds')
        nsubj <- nrow(data$raw)
    }

   # Might want to keep track separate copies as well as merged which is why I'm doing it this way...
   cleanup <- function(dframe){
       # These fields don't need to go into analysis:
       drop_flds <- c('Site_14'   , 'Group_14'      , 'Label_14_x', 'Label_14_y', 'Gender_14_x', 'Gender_14_y',
                      'QR_flag_14', 'Reliability_16', '', 'espad_6.life.nic_14', 'espad_8a.alc.life_14',
                      'agn_mean_correct_latency_neutral_14', 'agn_total_omissions_neutral_14')

       # Drop them:
       dframe <- dframe[,!(names(dframe) %in% drop_flds)]

       # Convert Age_14 data to days
       # Don't know exact conversion since don't know years, but this is fine
       #dframe['Age_14'] = dframe['Age_14']*365.24
   }



    # Quality control...
    source('remove_outliers.r')
    stdev_threshold <- 3
    exceptions <- c('nic_14', 'prm_perc_cor_14'     , 'agn_total_omissions_negative_18',
                    'agn_total_omissions_neutral_18', 'agn_total_omissions_positive_18')
    check_flds <- setdiff(c(data$cgt_names, data$agn_names, data$mid_names, data$sst_names, data$msc_names),exceptions)
    data$raw   <- remove_outliers(data$raw, check_flds, stdev_threshold)


    # Shuffle the data
    if (shuffle){
        # NOTE: I DON'T KNOW WHAT I WAS TALKING ABOUT HERE, SO THIS BLOCK MAY NOT BE RIGHT
        # Need to exclude the age 18 survey features so their absence remains co-occurant.
        # Not doing so will "actually" break the data for the 14-18 CCA.
        # Shuffling for that test will need to be done at that time.
        for (name in names(data$raw)) {
            col  <- cnames[i]
            perm <- sample(1:length(data$raw)[1], replace=TRUE)
            data$raw[col] <- data$raw[perm,col]
        }
    }

    #ncols_1 <- dim(data_14_task.imputed)[2]
    #ncols_2 <- dim(data_14_survey.imputed)[2]

    #nrows   <- dim(data_14_task.imputed)[1]

    #imputed_df <- data.frame(matrix(ncol=(ncols_1 + ncols_2), nrow=nrows))

    #imputed_df[,1:ncols_1]                     <- data_14_task.imputed
    #imputed_df[,(ncols_1+1):(ncols_1+ncols_2)] <- data_14_survey.imputed

    #colnames(imputed_df) <- c(colnames(data_14_task.imputed), colnames(data_14_survey.imputed))

    #dim(imputed_df)
    #dim(data_18_survey.raw)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#
# ------------------------------------------------------------------------------ #
replace_bad_names <- function(data){

  # Dictionary of old names and better names
  dict <- list(
    c('agn_mean_correct_latency_negative', 'agn_cor_lat_neg'),
    c('agn_mean_correct_latency_positive', 'agn_cor_lat_pos'),
    c('agn_total_omissions_negative', 'agn_tot_miss_neg'),
    c('agn_total_omissions_positive', 'agn_tot_miss_pos'),
    c('cgt_delay_aversion',  'cgt_delay_avers'),
    c('cgt_deliberation_time', 'cgt_delib_time'),
    c('cgt_quality_of_decision_making', 'cgt_decision_qual'),
    c('cgt_overall_proportion_bet', 'cgt_prop_bet'),
    c('cgt_risk_adjustment', 'cgt_risk_adjust'),
    c('cgt_risk_taking', 'cgt_risk_success'),
    c('prm_percent_correct', 'prm_perc_corr'),
    c( 'swm_between_errors',  'swm_errors'),
    c('log10.k.', 'discount'),
    c('sj1a', 'adhd_teacher'),
    c('sj1b', 'adhd_parent'),
    c('sj1c', 'adhd_child'),
    c('C.18i', 'alc_regret'),
    c('espad_6_life_nic', 'nic_use'),
    c('espad_8a_alc_life', 'alc_use')
    )

  to_replace <- sapply(dict, function(x){x[1]} )

  # Clunky nested loops but whatever...
  for (set in names(data$names) ) {
    for (fld in data$names[[set]]) {
      if (fld %in% to_replace) {
        old_names <- colnames(data$raw)
        new_name  <- dict[to_replace == fld][[1]][2]

        # Modify indexing
        data$names[[set]][data$names[[set]] == fld] <- new_name

        # Modify actual data frame
        colnames(data$raw) <- replace(old_names, colnames(data$raw) == fld, new_name)
      }
    }
  }

  return(data)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#         This function manages the quality control processes                    #
# ------------------------------------------------------------------------------ #
quality_control <- function(data) {
  library(futile.logger)

  # Order here matters...

  flog.info('Applying misc. transforms...')
  data <- misc_transforms(data)

  flog.info('Dropping useless fields...')
  data <- drop_useless_flds(data)

  flog.info('Replacing bad names...')
  data <- replace_bad_names(data)

  ignore    <- c(data$names$ESPAD, data$names$genes, data$names$Misc, 'Subject', 'set')
  threshold <- 4.5
  flog.info('Removing %s std outliars but ignoring %s ...', toString(threshold), toString(ignore))
  data   <- remove_outliers(data = data, ignore = ignore, thresh = threshold)

  flog.info('Z-scoring features...')
  scores <- z_score_data(data = data$raw, ignore = ignore)
  data$scores <- scores

  return(data)
}
# ------------------------------------------------------------------------------ #
