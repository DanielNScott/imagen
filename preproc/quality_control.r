remove_outliers <- function(dframe, feat_list, thresh) {
    
    for (i in 1:length(feat_list)){
        feature <- feat_list[[i]]
        
        var   <- dframe[feature]
        stdev <- apply(var, 2, sd    , na.rm=TRUE)
        dmean <- apply(var, 2, mean  , na.rm=TRUE)
        dmed  <- apply(var, 2, median, na.rm=TRUE)

        outmsk <- apply(var, 2, (function (x) ifelse(abs(dmed-x) < thresh*stdev | is.na(x) ,FALSE ,TRUE)))        

        if (sum(outmsk) > 0) {
            print(paste('Removing ', toString(sum(outmsk)), ' points from ', feature))
            print(dframe[feature][outmsk])
            cat('\n')

            dframe[feature][outmsk] <- NA
        }
    }
    return(dframe)
}



# Remove the wierd negative one values from the continuous task variables.
data$raw[data$task_names_14][data$raw[data$task_names_14] == -1] <- NA
data$raw[data$task_names_18][data$raw[data$task_names_18] == -1] <- NA

# Drop columns with all zeros, since they're uninformative.
raw_cols <- colnames(data$raw)
data$raw <- data$raw[, colSums(data$raw != 0, na.rm=TRUE) > 0]

# Drop some other features that are suspect
drop_list <- c('sig_int_14', 'sig_std_inv_rt_14', 'mu_std_inv_rt_14', 'p_tf_14', 'tau_stop_14')
data$raw  <- data$raw[, !names(data$raw) %in% drop_list ]

data$task_names_14 <- setdiff(c(data$task_names_14, 'stn_std'), drop_list)
n_feat_task_14     <- length(data$task_names_14)

# Indicate which features are still a part of the data:
flog('The following features were dropped:\n')
setdiff(raw_cols, colnames(data$raw))

cat('Hence, the following features are retained:\n')
print(colnames(data$raw))




# ------------------------------------------------------------------------------ #
#                Who knows what this is for...
# ------------------------------------------------------------------------------ #
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

# This function controls the 
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

    # Change annoying properties of all kinds

    # Clean up data frames:
    df_train <- cleanup(df_train)
    df_train['Age_14'] = df_train['Age_14']*365.24

    df_test  <- cleanup(df_test)
    df_test['Age_14'] = df_test['Age_14']*365.24

    # Stack them to create a big copy
    raw_df <- rbind(df_train, df_test)
    raw_df['Age_14'] = raw_df['Age_14']*365.24

    print('raw_df dims & num(NA) ')
    print( dim(raw_df) )
    print( t(sapply(raw_df, function(x) sum(is.na(x)))) )


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
    stdev_threshold <- 4
    exceptions <- c('nic_14', 'prm_perc_cor_14'     , 'agn_total_omissions_negative_18',
                    'agn_total_omissions_neutral_18', 'agn_total_omissions_positive_18')
    check_flds <- setdiff(c(data$cgt_names, data$agn_names, data$mid_names, data$sst_names, data$msc_names),exceptions)
    data$raw   <- remove_outliers(data$raw, check_flds, stdev_threshold)

    # Quick access to train and test data...
    data$train <- data$raw[data$train_inds,]
    data$test  <- data$raw[data$test_inds ,]

    # Z-score the features
    for (i in 1:n_feat_task_14){
        feature <- data$task_names_14[[i]]
        data$raw[feature] <- data$raw[feature] - colMeans(data$raw[feature], na.rm=TRUE)
        data$raw[feature] <- data$raw[feature] / apply(data$raw[feature], 2, sd, na.rm=TRUE)
    }

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

    # Multiple Imputation by Chained Equations (MICE)
    library(Rcpp)
    library(mice)
    imp <- mice(data$raw[data$task_names_14][data$train_inds,], print = FALSE)
    imp <- mice(data$raw[data$task_names_14][data$train_inds,], pred = imp$predictorMatrix, seed=23109, m=n_imputed, maxit=max_iters, print = FALSE)

    data$imp <- imp
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