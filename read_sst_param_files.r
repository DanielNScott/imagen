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
            if (any(subj_num == 157:185) & skip){            
            } else{
                sst_params[subj_num, param_name] <- trace_means[param_name_long]
            }
        }
    }

    return(sst_params)
}

read_sst_param_files <- function(data){
    # Read the SST parameters' file:
    sst_params_file  <- '/home/dan/documents/lncc/sst_data_train/parameters1.csv'
    subjects <- data$raw['Subject'][1:data$train_indices,]
    #subjects <- df_train['Subject']

    sst_params_train <- read_sst_params(sst_params_file, subjects, TRUE)

    # Read the SST parameters' file:
    sst_params_file  <- '/home/dan/documents/lncc/sst_data_test_100_200/parameters1.csv'
    subjects <- data$raw['Subject'][data$train_indices+1:nrow(data$raw),]
    #subjects <- df_test['Subject']

    sst_params_test <- read_sst_params(sst_params_file, subjects, FALSE)

    # Read the SST parameters' file:
    sst_params_file  <- '/home/dan/documents/lncc/sst_data_test_0_50/parameters1.csv'
    sst_params_test2 <- read_sst_params(sst_params_file, subjects, FALSE)

    sst_params_test <- rbind(sst_params_test2[1:50,], sst_params_test[51:198,])
    sst_params      <- rbind(sst_params_train, sst_params_test)
    #dim(sst_params)
    #dim(raw_df)

    train_indices <- dim(sst_params_train)[1]
    data <- list("raw" = sst_params, "train_indices" = train_indices)

    return(data)
}