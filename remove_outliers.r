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