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
