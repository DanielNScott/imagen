import_fmri <- function(subset, visit, roi) {

# Req. library for revalue()
library(plyr)

# Get training data subject IDs
sst_train_fname <- paste('/home/dan/documents/lncc/synthetic data/sst_round_2/safe_',subset,'_data.csv', sep='')

# Set directories for fMRI and stop signal data
stn_dir <- paste('/home/dan/documents/lncc/From Catherine/', visit, '/',roi,'_',visit,'_SSTtimeseries/', sep='')
sst_dir <- paste('/home/dan/documents/lncc/From Catherine/Round 2/SST/SST_', visit,'_behavioraldata/',sep='')

# Set filename prefixes and suffixes
sst_pfx <- 'ss_'
stn_sfx <- paste(roi,'_',visit,'_SSTtimeseries',sep='')

# Read the training data
sst_train_data  <- read.csv(file=sst_train_fname, header=TRUE, sep=',');

# Specify ordinal coding for trial outcome
factor_key  <- c("GO_SUCCESS"="1", "GO_FAILURE"='2', 'STOP_SUCCESS'='3', 'STOP_FAILURE'='4',
                 'GO_TOO_LATE'='5','GO_WRONG_KEY_RESPONSE'='6', 'STOP_TOO_EARLY_RESPONSE'='7')

# Num. of subjects in training data
n_train     <- dim(sst_train_data['Subject'])[1]

# Init. arrays w/ conservative size estimates
stn_series  <- array(dim=c(444,n_train))
sst_times   <- array(dim=c(800,n_train))
sst_outcome <- array(dim=c(800,n_train))

failed_stn_reads <- 0
failed_sst_reads <- 0

failed_stn_cols  <- c()
failed_sst_cols  <- c()

# Loop through subjects and collect their data
for (subj_num in 1:n_train) {
    
    # File reading setup
    subj_id_num <- sst_train_data['Subject'][subj_num,]
    subj_id_str <- formatC(subj_id_num, width=12, format='d', flag='0')

    stn_fname <- paste(subj_id_str, '_', stn_sfx, sep='')
    sst_fname <- paste(sst_pfx, subj_id_str, '.csv', sep='')

    stn_full_name <- paste(stn_dir, stn_fname, sep='')
    sst_full_name <- paste(sst_dir, sst_fname, sep='')

    # Reading the STN Series
    if (file.exists(stn_full_name)) {
        # Only want to read one column, activations
        cur_stn_series <- read.csv(file=stn_full_name, header=TRUE, sep='\t', colClasses=c("NULL", "NULL", NA))
        print(paste('Successful read of ', stn_fname))
    } else {
        print(paste('Failure to read ', stn_fname))
        failed_stn_reads <- failed_stn_reads + 1
        failed_stn_cols  <- c(failed_stn_cols, subj_num)
    }

    # Reading the SST Series
    if (file.exists(sst_full_name)) {
        #print(sst_full_name)
        #print(file.exists(sst_full_name))
        # Only want to read 2 columns, time and outcome
        cols <- c('NULL','NULL','numeric','NULL','NULL','NULL','NULL','NULL',
                  'NULL','NULL','NULL'   , NA   ,'NULL','NULL','NULL','numeric')
        cur_sst_series <- read.csv(file=sst_full_name, header=TRUE, sep='\t', colClasses=cols, skip=1)
        print(paste('Successful read of ', sst_fname))
    } else {
        print(paste('Failure to read ', sst_fname))
        failed_sst_reads <- failed_sst_reads + 1
        failed_sst_cols  <- c(failed_sst_cols, subj_num)
        next
    }
    
    # Extract the activations
    len <- length(cur_stn_series$Mean_1)
    stn_series[1:len,subj_num] <- cur_stn_series$Mean_1
    
    # Extract the trial times and outcomes
    len <- dim(cur_sst_series['Trial.Start.Time..Onset.'])[1]
    sst_times[1:len, subj_num] <- cur_sst_series['Trial.Start.Time..Onset.'][,1]
    
    # Code outcome values
    cur_sst_series['Response.Outcome'][,1] <- revalue(cur_sst_series['Response.Outcome'][,1], factor_key, warn_missing=FALSE)
    sst_outcome[1:len, subj_num] <- as.numeric(as.character(cur_sst_series['Response.Outcome'][,1]))
}
print(paste('Completed ', n_train - failed_stn_reads, 'out of', n_train, 'STN reads.', sep=' '))
print(paste('Completed ', n_train - failed_sst_reads, 'out of', n_train, 'SST reads.', sep=' '))

lost_cols  <- c(failed_stn_cols, failed_sst_cols)
#print(lost_cols)
valid_cols <- setdiff(1:subj_num, lost_cols)
#print(valid_cols)

# Standard deviatons - Was using this as a proxy for a little while. 
# stn_train <- data.frame('Subject' = sst_train_data['Subject'], 'stn_std' = apply(stn_series, 2, sd, na.rm = TRUE))

# Format the STN data as a dataframe and add TRs
stn_series <- data.frame(stn_series[,valid_cols])
colnames(stn_series) <- paste('Subj_', valid_cols, sep='');
stn_series['TR'] <- 1:444

# Need to save the polynomial function since 'signal' has 'poly' as well.
polynomial <- poly
library(signal)

# Clean activation data
for (feature in setdiff(colnames(stn_series), 'TR') ) {

    # Z-Score the data for each participant
    mean <- colMeans(stn_series[feature], na.rm=TRUE)
    std  <- apply(stn_series[feature], 2, sd, na.rm = TRUE)
    stn_series[feature] <- (stn_series[feature] - mean) / std

    # Spike removal
    cur_series <- stn_series[feature]
    cur_series[abs(cur_series) >= 3,] <- 0

    # Drift correction
    polyfit <- lm(cur_series[,1] ~ polynomial(stn_series[['TR']],2))
    stn_series[feature] <- cur_series - polyfit$fitted.values

    # High-Pass Filter
    bf <- butter(2, 1/128, type="high")
    stn_series[feature] <- filtfilt(bf, as.vector(stn_series[[feature]]))

    # Drift correction
    polyfit <- lm(cur_series[,1] ~ polynomial(stn_series[['TR']],2))
    stn_series[feature] <- cur_series - polyfit$fitted.values

    # Re-Z-Score
    mean <- colMeans(stn_series[feature], na.rm=TRUE)
    std  <- apply(stn_series[feature], 2, sd, na.rm = TRUE)
    stn_series[feature] <- (stn_series[feature] - mean) / std
}

detach("package:signal", unload=TRUE)

data <- list(stn_series  = stn_series, 
             sst_outcome = sst_outcome[,valid_cols],
             sst_times   = sst_times[,valid_cols],
             subs_lost   = lost_cols)
return (data)

} # EOF