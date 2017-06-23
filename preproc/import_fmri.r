import_fmri <- function(id, task, timeseries_dir, taskdata_dir) {

# subset
# visit
# roi

# Req. library for revalue()
library(plyr)

# Get training data subject IDs
sst_train_fname <- paste('/home/dan/documents/lncc/synthetic data/sst_round_2/safe_',subset,'_data.csv', sep='')

# Read the training data
sst_train_data  <- read.csv(file=sst_train_fname, header=TRUE, sep=',');

# Num. of subjects in training data
n_train     <- dim(sst_train_data['Subject'])[1]

# Set filename prefixes and suffixes
sst_pfx <- 'ss_'
fmri_sfx <- paste(roi,'_',visit,'_SSTtimeseries',sep='')

# Specify ordinal coding for trial outcome
factor_key  <- c("GO_SUCCESS"="1", "GO_FAILURE"='2', 'STOP_SUCCESS'='3', 'STOP_FAILURE'='4',
                 'GO_TOO_LATE'='5','GO_WRONG_KEY_RESPONSE'='6', 'STOP_TOO_EARLY_RESPONSE'='7')

# Init. arrays w/ conservative size estimates
fmri_series  <- array(dim=c(444,n_train))
sst_times   <- array(dim=c(800,n_train))
sst_outcome <- array(dim=c(800,n_train))

failed_fmri_reads <- 0
failed_sst_reads <- 0

failed_fmri_cols  <- c()
failed_sst_cols  <- c()

# Loop through subjects and collect their data
for (subj_num in 1:n_train) {
    
    # File reading setup
    subj_id_num <- sst_train_data['Subject'][subj_num,]
    subj_id_str <- formatC(subj_id_num, width=12, format='d', flag='0')

    fmri_fname <- paste(subj_id_str, '_', fmri_sfx, sep='')
    sst_fname <- paste(sst_pfx, subj_id_str, '.csv', sep='')

    fmri_full_name <- paste(timeseries_dir, fmri_fname, sep='')
    sst_full_name <- paste(taskdata_dir  , sst_fname, sep='')

    # Reading the fmri Series
    if (file.exists(fmri_full_name)) {
        # Only want to read one column, activations
        cur_fmri_series <- read.csv(file=fmri_full_name, header=TRUE, sep='\t', colClasses=c("NULL", "NULL", NA))
        print(paste('Successful read of ', fmri_fname))
    } else {
        print(paste('Failure to read ', fmri_fname))
        failed_fmri_reads <- failed_fmri_reads + 1
        failed_fmri_cols  <- c(failed_fmri_cols, subj_num)
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
    len <- length(cur_fmri_series$Mean_1)
    fmri_series[1:len,subj_num] <- cur_fmri_series$Mean_1
    
    # Extract the trial times and outcomes
    len <- dim(cur_sst_series['Trial.Start.Time..Onset.'])[1]
    sst_times[1:len, subj_num] <- cur_sst_series['Trial.Start.Time..Onset.'][,1]
    
    # Code outcome values
    cur_sst_series['Response.Outcome'][,1] <- revalue(cur_sst_series['Response.Outcome'][,1], factor_key, warn_missing=FALSE)
    sst_outcome[1:len, subj_num] <- as.numeric(as.character(cur_sst_series['Response.Outcome'][,1]))
}
print(paste('Completed ', n_train - failed_fmri_reads, 'out of', n_train, 'fmri reads.', sep=' '))
print(paste('Completed ', n_train - failed_sst_reads, 'out of', n_train, 'SST reads.', sep=' '))

lost_cols  <- c(failed_fmri_cols, failed_sst_cols)
#print(lost_cols)
valid_cols <- setdiff(1:subj_num, lost_cols)
#print(valid_cols)

# Standard deviatons - Was using this as a proxy for a little while. 
# fmri_train <- data.frame('Subject' = sst_train_data['Subject'], 'fmri_std' = apply(fmri_series, 2, sd, na.rm = TRUE))

# Format the fmri data as a dataframe and add TRs
fmri_series <- data.frame(fmri_series[,valid_cols])
colnames(fmri_series) <- paste('Subj_', valid_cols, sep='');
fmri_series['TR'] <- 1:444

# Need to save the polynomial function since 'signal' has 'poly' as well.
polynomial <- poly
library(signal)

# Clean activation data
for (feature in setdiff(colnames(fmri_series), 'TR') ) {

    # Z-Score the data for each participant
    mean <- colMeans(fmri_series[feature], na.rm=TRUE)
    std  <- apply(fmri_series[feature], 2, sd, na.rm = TRUE)
    fmri_series[feature] <- (fmri_series[feature] - mean) / std

    # Spike removal
    cur_series <- fmri_series[feature]
    cur_series[abs(cur_series) >= 3,] <- 0

    # Drift correction
    polyfit <- lm(cur_series[,1] ~ polynomial(fmri_series[['TR']],2))
    fmri_series[feature] <- cur_series - polyfit$fitted.values

    # High-Pass Filter
    bf <- butter(2, 1/128, type="high")
    fmri_series[feature] <- filtfilt(bf, as.vector(fmri_series[[feature]]))

    # Drift correction
    polyfit <- lm(cur_series[,1] ~ polynomial(fmri_series[['TR']],2))
    fmri_series[feature] <- cur_series - polyfit$fitted.values

    # Re-Z-Score
    mean <- colMeans(fmri_series[feature], na.rm=TRUE)
    std  <- apply(fmri_series[feature], 2, sd, na.rm = TRUE)
    fmri_series[feature] <- (fmri_series[feature] - mean) / std
}

detach("package:signal", unload=TRUE)

data <- list(fmri_series  = fmri_series, 
             sst_outcome = sst_outcome[,valid_cols],
             sst_times   = sst_times[,valid_cols],
             subs_lost   = lost_cols)
return (data)

} # EOF