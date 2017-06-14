
import_all <- function(){

   # Set up all the files that need to be read in one spot...
   sst_files <- '/home/dan/documents/lncc/synthetic data/sst_round_2/safe_test_data.csv'
               '/home/dan/documents/lncc/synthetic data/sst_round_2/safe_train_data.csv'

   # Read the training data file
   fname    <- '/home/dan/documents/lncc/synthetic data/sst_round_2/safe_test_data.csv'
   df_train <- read.csv(file=fname, header=TRUE, sep=',');
    
   print('df_train dims & num(NA) ')
   print( dim(df_train) )
   print( t(sapply(df_train, function(x) sum(is.na(x)))) )

   # Read the training data file
   fname   <- '/home/dan/documents/lncc/synthetic data/sst_round_2/safe_train_data.csv'
   df_test <- read.csv(file=fname, header=TRUE, sep=',');

   print('df_test dims & num(NA) ')
   print( dim(df_test) )
   print( t(sapply(df_test, function(x) sum(is.na(x)))) )

   # Measures that exist for age 16: Currently not being used.
   #features_16 <- c('Age.for.timestamp_16', 'X6.life.nic_16', 'X8a.life.alc_16', 'life_amphet_16', 'life_anabolic_16',
   #    'life_coke_16', 'life_crack_16', 'life_ghb_16', 'life_glue_16', 'life_hash_16', 'life_heroin_16',
   #    'life_ketamine_16', 'life_lsd_16', 'life_mdma_16', 'life_mushrooms_16', 'life_tranq_16', 'log10.k._16')

   sst_params <- read_sst_param_files(df_train['Subject'], df_test['Subject'])

   print('sst_params dims & num(NA) ')
   print( dim(sst_params) )
   print( t(sapply(sst_params, function(x) sum(is.na(x)))) )

   # Read the MID fit file in:
   mid_params_file <- '/home/dan/documents/lncc/From Nick Jan 20/MIDT_SubjectFits_BL_All.csv'
   mid_params_raw  <- read.csv(file=mid_params_file, header=TRUE, sep=',');

   #
   mid_params       <- data.frame(matrix(ncol=12, nrow=dim(unique(mid_params_raw['Subject']))[1]))
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

   print('mid_params dims & num(NA) ')
   print( dim(mid_params) )
   print( t(sapply(mid_params, function(x) sum(is.na(x)))) )

   tmp <- merge(sst_params, mid_params, by="Subject")

   ### Read the STN timeseries data
   directory   <- '/home/dan/documents/lncc/From Catherine/BL_STN_SSTtimeseries/'
   file_suffix <- 'STNboth_BL_SSTtimeseries'

   n_train    <- dim(df_train['Subject'])[1]
   stn_series <- array(dim=c(444,n_train))

   for (subj_num in 1:n_train) {
      subj_id_num <- df_train['Subject'][subj_num,]
      subj_id_str <- formatC(subj_id_num, width=12, format='d', flag='0')
     
      fname   <- paste(directory, subj_id_str, '_', file_suffix, sep='')
     
      cur_series = tryCatch({
         read.csv(file=fname, header=TRUE, sep='\t', colClasses=c("NULL", "NULL", NA))
         }, warning = function(w) { print(paste('Failure to read ', fname)) }
      )
     
      if(is.character(cur_series)){
         next
      }
     
      if (dim(cur_series)[1] == 444){
         stn_series[,subj_num] <- cur_series$Mean_1
      }
   }
   stn_stats_train <- data.frame('Subject' = df_train['Subject'],
                           'stn_std' = apply(stn_series, 2, sd, na.rm = TRUE))

   n_test    <- dim(df_test['Subject'])[1]
   for (subj_num in 1:n_test) {
      subj_id_num <- df_test['Subject'][subj_num,]
      subj_id_str <- formatC(subj_id_num, width=12, format='d', flag='0')
     
      fname   <- paste(directory, subj_id_str, '_', file_suffix, sep='')
     
      cur_series = tryCatch({
         read.csv(file=fname, header=TRUE, sep='\t', colClasses=c("NULL", "NULL", NA))
         }, warning = function(w) { print(paste('Failure to read ', fname)) }
      )
     
      if(is.character(cur_series)){
         next
      }
     
      if (dim(cur_series)[1] == 444){
         stn_series[,subj_num] <- cur_series$Mean_1
      }
   }
   stn_stats_test <- data.frame('Subject' = df_test['Subject'],
                           'stn_std' = apply(stn_series, 2, sd, na.rm = TRUE))

   stn_stats <- rbind(stn_stats_train,stn_stats_test)

   print('stn_stats dims & num(NA) ')
   print( dim(stn_stats) )
   print( t(sapply(stn_stats, function(x) sum(is.na(x)))) )

   tmp    <- merge(tmp, stn_stats, by="Subject", all=TRUE)
   raw_df <- merge(raw_df, tmp, by="Subject", all=TRUE)

   train_inds <- 1:dim(df_train)[1]
   test_inds  <- (dim(df_train)[1]+1):nrow(raw_df)

   data <- replace_bad_names(sst_params, mid_params, raw_df)
   data$train_inds <- train_inds
   data$test_inds  <- test_inds

   print('data$raw dims & num(NA) ')
   print( dim(data$raw) )
   print( t(sapply(data$raw, function(x) sum(is.na(x)))) )

   return(data)
}

read_sst_param_files <- function(train_subjects, test_subjects){
    # Read the SST parameters' file:
    sst_params_file  <- '/home/dan/documents/lncc/sst_data_train/parameters1.csv'

    sst_params_train <- read_sst_params(sst_params_file, train_subjects, TRUE)

    # Read the SST parameters' file:
    sst_params_file  <- '/home/dan/documents/lncc/sst_data_test_100_200/parameters1.csv'

    sst_params_test <- read_sst_params(sst_params_file, test_subjects, FALSE)

    # Read the SST parameters' file:
    sst_params_file  <- '/home/dan/documents/lncc/sst_data_test_0_50/parameters1.csv'
    sst_params_test2 <- read_sst_params(sst_params_file, test_subjects, FALSE)

    sst_params_test <- rbind(sst_params_test2[1:50,], sst_params_test[51:198,])
    sst_params      <- rbind(sst_params_train, sst_params_test)
    #dim(sst_params)
    #dim(raw_df)

    print('sst_params dims & num(NA)')
    dim(sst_params)
    t(sapply(sst_params, function(x) sum(is.na(x))))

    return(sst_params)

}

read_sst_params <- function(filename, subjects, skip){
    sst_params_raw   <- read.csv(file=filename, header=TRUE, sep=';');

    # What parameters are being read?
    sst_params_names <- c('mu_go_subj' ,'mu_stop_subj' ,'sigma_go_subj', 'sigma_stop_subj',
                          'tau_go_subj','tau_stop_subj','p_tf_subjpt')

    # Set up the dataframe to accept mean parameter values.
    sst_params           <- data.frame(matrix(ncol=7, nrow=198))
    colnames(sst_params) <- sst_params_names 
    trace_means          <- colMeans(sst_params_raw[2000:3000,])
    sst_params['Subject']<- subjects

    # Save the parameters:
    for (subj_num in 1:198) {
        for (param_name in sst_params_names) {

            param_name_long <- paste(param_name,'.',subj_num, sep='')

            # Currently missing some data...
            if (any(subj_num == 166:185) & skip){            
            } else{
                sst_params[subj_num, param_name] <- trace_means[param_name_long]
            }
        }
    }

    return(sst_params)
}

replace_bad_names <- function(sst_params, mid_params, raw_df, new_features_14_task_names){
   # Task data at age 14:
   features_14_task_raw <- c('IQ_PR_14', 'IQ_VC_14', 'GB_SSRT_14',
        'agn_mean_correct_latency_negative_14',
        'agn_mean_correct_latency_positive_14',
        'agn_total_omissions_negative_14',
        'agn_total_omissions_positive_14',
        'cgt_delay_aversion_14',
        'cgt_deliberation_time_14',
        'cgt_quality_of_decision_making_14',
        'cgt_overall_proportion_bet_14',
        'cgt_risk_adjustment_14',
        'cgt_risk_taking_14',
        'prm_percent_correct_14',
        'rvp_a_14',
        'swm_between_errors_14',
        'swm_strategy_14',
        'log10.k._14', setdiff(colnames(sst_params),'Subject'), setdiff(colnames(mid_params),'Subject'))


   new_sst_names <- c('mu_go_14', 'mu_stop_14', 'sigma_go_14','sigma_stop_14','tau_go_14','tau_stop_14','p_tf_14')
   new_mid_names <- c('mu_targ_dur_co_14', 'mu_targ_left_co_14',  'mu_rewarded_co_14', 
                      'mu_high_rewarded_co_14', 'mu_std_inv_rt_14', 'mu_int_14',
                      'sig_targ_dur_co_14','sig_targ_left_co_14','sig_rewarded_co_14',
                      'sig_high_rewarded_co_14','sig_std_inv_rt_14','sig_int_14')           

   new_agn_names <- c('agn_mean_corr_lat_neg_14', 'agn_mean_corr_lat_pos_14',
                      'agn_num_omis_neg_14', 'agn_num_omis_pos_14')

   new_cgt_names <- c('cgt_delay_avers_14','cgt_delib_14','cgt_quality_14','cgt_prop_bet_14',
                      'cgt_risk_adjust_14','cgt_risk_taking_14')

   new_espad_names <- c('alc_14', 'nic_14', 'amphet_14', 'coke_14', 'crack_14', 'ghb_14',
                        'glue_14', 'hash_14', 'ketamine_14', 'lsd_14', 'mushrooms_14', 'narc_14')

   features_14_task <- c('IQ_PR_14', 'IQ_VC_14', 'GB_SSRT_14', new_agn_names, new_cgt_names,
                         'prm_perc_corr_14', 'rvp_a_14', 'swm_btwn_errs_14', 'swm_strategy_14',
                         'log10.k._14', new_sst_names, new_mid_names)

   misc_task_names <- c('IQ_PR_14', 'IQ_VC_14', 'GB_SSRT_14', 'prm_perc_corr_14', 'rvp_a_14',
                        'swm_btwn_errs_14', 'swm_strategy_14', 'log10.k._14')

    features_14_survey_raw <- c('Sex_best_M0_14', 'PDS_14','All_Alc_14', 'All_Nic_14',
           'espad_life_amphet_14',
           'espad_life_coke_14',
           'espad_life_crack_14',
           'espad_life_ghb_14',
           'espad_life_glue_14',
           'espad_life_hash_14',
           'espad_life_ketamine_14',
           'espad_life_lsd_14',
           'espad_life_mushrooms_14',
           'espad_life_narcotic_14')

   # Features left out:
   #       'espad_life_anabolic_14',
   #       'espad_life_heroin_14',
   #       'espad_life_mdma_14',
   #       'espad_life_tranq_14'

   features_14_survey <- c('bio_sex_14', 'pds_14', new_espad_names)

   features_18_task <- c('log10.k._18', 'agn_mean_correct_latency_negative_18', 'agn_mean_correct_latency_neutral_18',
     'agn_mean_correct_latency_positive_18', 'agn_total_omissions_negative_18', 'agn_total_omissions_neutral_18',
     'agn_total_omissions_positive_18', 'cgt_delay_aversion_18', 'cgt_deliberation_time_18', 
     'cgt_overall_proportion_bet_18', 'cgt_quality_of_decision_making_18', 'cgt_risk_adjustment_18',
     'cgt_risk_taking_18', 'prm_percent_correct_18', 'rvp_a_18', 'swm_between_errors_18', 'swm_strategy_18')
        
    #features_18_survey <- c('age_18', 'X6.life.nic_18', 'X8a.life.alc_18', 'Life.amph_18', 'Life.anab_18',
    #    'Life.coke_18', 'Life.crack_18', 'Life.hash.thc_18', 'Life.heroin_18', 'Life.GHB_18', 'Life.glue_18', 
    #    'Life.ketamine_18', 'Life.lsd_18', 'Life.MDMA_18', 'Life.mushrooms_18', 'Lif.narc_18', 'Life.tranq_18')

    features_18_survey <- c('X6.life.nic_18', 'X8a.life.alc_18', 'Life.amph_18', 'Life.anab_18',
        'Life.coke_18', 'Life.crack_18', 'Life.hash.thc_18', 'Life.heroin_18',
        'Life.ketamine_18', 'Life.lsd_18', 'Life.MDMA_18', 'Life.mushrooms_18', 'Lif.narc_18')


    n_feat_task   <- length(features_14_task)
    n_feat_survey <- length(features_14_survey)

    for (i in 1:n_feat_task){
        old_name <- features_14_task_raw[[i]]
        new_name <- features_14_task[[i]]
        colnames(raw_df)[colnames(raw_df) == old_name] <- new_name
    }
    for (i in 1:n_feat_survey){
        old_name <- features_14_survey_raw[[i]]
        new_name <- features_14_survey[[i]]
        names(raw_df)[names(raw_df) == old_name] <- new_name
    }


    output <- list('raw' = raw_df, 'survey_names_14' = features_14_survey,
                   'task_names_14' = features_14_task, 'sst_names' = new_sst_names,
                   'mid_names' = new_mid_names, 'espad_names' = new_espad_names,
                   'cgt_names' = new_cgt_names, 'agn_names' = new_agn_names,
                   'msc_names' = misc_task_names,
                   'task_names_18' = features_18_task,
                   'survey_names_18' = features_18_survey
                   )
    return(output)
}
