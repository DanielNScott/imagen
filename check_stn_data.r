directory  <- '/home/dan/documents/lncc/From Catherine/BL_STN_SSTtimeseries/'
file.names <- dir(directory)
out.file   <- ""

stn_series <- array(dim=c(444,2076))

for(i in 1:length(file.names)){   
   cur_series = tryCatch({
      read.csv(file=paste(directory,file.names[i],sep=''), header=TRUE, sep='\t', colClasses=c("NULL", "NULL", NA))
      }, warning = function(w) { print(paste('Failure to read ', file.names[i])) }
   )
  
   if(is.character(cur_series)){
      next
   }
  
   if (dim(cur_series)[1] == 444){
      stn_series[,i] <- cur_series$Mean_1
   }
}
stn_data <- data.frame('filename' = file.names,
                       'stn_std'  = apply(stn_series, 2, sd, na.rm = TRUE))

bad_stn_msk <- !is.na(stn_data['stn_std']) & stn_data['stn_std'] > 50

print(c('Num bad stn data:', sum(bad_stn_msk)))
print(stn_data[c('filename','stn_std')][bad_stn_msk,])