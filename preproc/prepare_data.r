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