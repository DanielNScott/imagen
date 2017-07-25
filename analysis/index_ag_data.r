index_ag_data <- function(data) {
 # Initialize lists of trial numbers for diff conds.
  Go    <- list()
  ST    <- list()
  SR    <- list()
  miss  <- list()
  error <- list();

  # Populate lists
  for (subj in 1:length(data)) {
  	conditions <- rownames(data[[subj]])

  	Go[[subj]]    <- grep('GO_SUCCESS' , conditions)
  	ST[[subj]]    <- grep('STOP_SUCCESS', conditions)
  	SR[[subj]]    <- grep('STOP_FAILURE', conditions)
  	#miss[[subj]]  <- grep('', yDat[[c]]$type)
  	#error[[subj]] <- grep('', yDat[[c]]$type)
  }

  # Copy subject names
  names(Go)    <- names(data)
  names(ST)    <- names(data)
  names(SR)    <- names(data)
  #names(error) <- names(data)
  #names(miss)  <- names(data)

  return(list(Go = Go, ST = ST, SR = SR))
}
