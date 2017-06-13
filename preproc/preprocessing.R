#show_views <- FALSE  # Toggle slice views throughout
#shuffle    <- FALSE  # Shuffle data for shuffle testing
#re_read    <- FALSE  # Re-read the data
#re_imp     <- FALSE  # Re-impute the data

#library(ggplot2)
#library(reshape2)
#library(mice)

# Actually read all the files: (Time consuming!)
if (re_read) {
    source('import_data.r')
    data <- import_data()
    saveRDS(data,'raw_data.rds')
} else {
    # Read saved R data file:
    data  <- readRDS('raw_data.rds')
    nsubj <- nrow(data$raw)
}

# View into data
if(show_views){
    view_rows <- c(1:5, (nsubj-5):nsubj)
    data$raw[view_rows,c('mu_targ_dur_co_14', 'mu_targ_left_co_14',  'mu_rewarded_co_14', 
           'mu_high_rewarded_co_14', 'mu_int_14',
           'sig_targ_dur_co_14','sig_targ_left_co_14','sig_rewarded_co_14',
           'sig_high_rewarded_co_14')]
}

options(repr.plot.width=6, repr.plot.height=2)
ggplot(melt(data$raw['stn_std'], id.vars = NULL), aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram(bins=100, na.rm=TRUE)

options(repr.plot.width=9, repr.plot.height=4)
ggplot( melt(data$raw[data$sst_names], id.vars=NULL), aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram(bins=30, na.rm=TRUE) + ggtitle("Recovered SST Parameter Densities")
ggplot( melt(data$raw[data$mid_names], id.vars=NULL), aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram(bins=30, na.rm=TRUE) + ggtitle("Recovered MID Parameter Densities")
ggplot( melt(data$raw[data$cgt_names], id.vars=NULL), aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram(bins=30, na.rm=TRUE) + ggtitle("CGT Data")
ggplot( melt(data$raw[data$agn_names], id.vars=NULL), aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram(bins=30, na.rm=TRUE) + ggtitle("Affective Go/No-Go Data")
ggplot( melt(data$raw[data$msc_names], id.vars=NULL), aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram(bins=30, na.rm=TRUE) + ggtitle("Other Task Data")

options(repr.plot.width=9, repr.plot.height=5)
ggplot( melt(data$raw[data$survey_names_14]), aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram(bins=30, na.rm=TRUE) + ggtitle("Survey Data")

# Bad STN data:
bad_stn_msk <- !is.na(data$raw['stn_std']) & data$raw['stn_std'] > 50

print(c('Num bad stn data being dropped:', sum(bad_stn_msk)))
#data$raw[c('Subject','stn_std')][bad_stn_msk,]
data$raw['stn_std'][!is.na(data$raw['stn_std']) & data$raw['stn_std'] > 100,] <- NA

source('remove_outliers.r')
stdev_threshold <- 4

exceptions <- c('nic_14', 'prm_perc_cor_14'     , 'agn_total_omissions_negative_18',
                'agn_total_omissions_neutral_18', 'agn_total_omissions_positive_18')
check_flds <- setdiff(c(data$cgt_names, data$agn_names, data$mid_names, data$sst_names, data$msc_names, 'stn_std'),exceptions)

data$raw <- remove_outliers(data$raw, check_flds, stdev_threshold)



# Quick access to train and test data...
data$train <- data$raw[data$train_inds,]
data$test  <- data$raw[data$test_inds ,]

# Just the continuous features though:
for (i in 1:n_feat_task_14){
    feature <- data$task_names_14[[i]]
    data$raw[feature] <- data$raw[feature] - colMeans(data$raw[feature], na.rm=TRUE)
    data$raw[feature] <- data$raw[feature] / apply(data$raw[feature], 2, sd, na.rm=TRUE)
}


if (shuffle){
    dim(raw_df)
    # Need to exclude the age 18 survey features so their absence remains co-occurant.
    # Not doing so will "actually" break the data for the 14-18 CCA.
    # Shuffling for that test will need to be done at that time.
    cnames <- setdiff(colnames(raw_df), features_18_survey)
    for (i in 1:length(cnames)) {
        col  <- cnames[i]
        perm <- sample(1:dim(raw_df)[1], replace=TRUE)
        raw_df[col] <- raw_df[perm,col]
    }
}


library(Rcpp)
library(mice)

# Display the number of missing values for each feature in train & test sets:
sapply(data$raw[data$task_names_14][data$train_inds,], function(x) sum(is.na(x)))
sapply(data$raw[data$task_names_14][data$test_inds,], function(x) sum(is.na(x)))

#---------------------------------------------------------------------#
#---                     This may take hours!                    -----#
#---------------------------------------------------------------------#
# Impute the features:
if (re_imp) {
    imp <- mice(data$raw[data$task_names_14][data$train_inds,], print = FALSE)
    imp <- mice(data$raw[data$task_names_14][data$train_inds,], pred = imp$predictorMatrix, seed=23109, m=30, maxit=1000, print = FALSE)
    saveRDS(imp,'imputed.rds')
} else{
    imp <- readRDS('imputed.rds')
}

#print(imp)
options(repr.plot.width=15, repr.plot.height=15)
stripplot(imp, pch = 20, cex = 0.75, subset = 1:(198*3))



library(reshape2)
library(ggplot2)

options(repr.plot.width=9, repr.plot.height=5)

#d <- melt(data_14_task.imputed)
#ggplot(d, aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram()

d <- melt(data_14_task.raw)
ggplot(d, aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram()

#d <- melt(data_14_survey.raw)
#ggplot(d, aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram()

d <- melt(data_18_survey.raw)
ggplot(d, aes(x = value)) + facet_wrap(~variable, scales = "free") + geom_histogram()



# Creating imputed_df
# I must be missing something... R only hase merge for horiz. concat of data frames?
# And it requires the to-be-merged data frames to have an identical column on which to merge?
# Can there really be no horzcat()?

ncols_1 <- dim(data_14_task.imputed)[2]
ncols_2 <- dim(data_14_survey.imputed)[2]

nrows   <- dim(data_14_task.imputed)[1]

imputed_df <- data.frame(matrix(ncol=(ncols_1 + ncols_2), nrow=nrows))

imputed_df[,1:ncols_1]                     <- data_14_task.imputed
imputed_df[,(ncols_1+1):(ncols_1+ncols_2)] <- data_14_survey.imputed

colnames(imputed_df) <- c(colnames(data_14_task.imputed), colnames(data_14_survey.imputed))

dim(imputed_df)
dim(data_18_survey.raw)
#missing_imputed_df <- get_missing_fracs(imputed_df)